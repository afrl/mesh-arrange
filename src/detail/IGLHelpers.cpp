#include "detail/IGLHelpers.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <deque>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace mesh_arrange::detail
{

namespace
{

struct SortedTriKeyHash
{
    std::size_t operator()(const std::array<Index, 3>& k) const noexcept
    {
        std::size_t seed = 0;
        for (int m = 0; m < 3; ++m)
        {
            seed ^= std::hash<Index>{}(k[m]) + 0x9e3779b9ULL + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

Eigen::Vector3d closest_point_on_segment_vec(const Eigen::Vector3d& a,
                                             const Eigen::Vector3d& b,
                                             const Eigen::Vector3d& p) noexcept
{
    const Eigen::Vector3d ab = b - a;
    const double len2        = ab.squaredNorm();
    if (len2 <= 1e-30)
    {
        return a;
    }
    const double t = std::clamp((p - a).dot(ab) / len2, 0.0, 1.0);
    return a + t * ab;
}

Eigen::Vector3d closest_point_on_triangle(const Eigen::Vector3d& a,
                                          const Eigen::Vector3d& b,
                                          const Eigen::Vector3d& c,
                                          const Eigen::Vector3d& p) noexcept
{
    // Ericson, Real-Time Collision Detection, Section 5.1.5 (Voronoi region method)
    const Eigen::Vector3d ab = b - a;
    const Eigen::Vector3d ac = c - a;
    const Eigen::Vector3d ap = p - a;

    const double d1 = ab.dot(ap);
    const double d2 = ac.dot(ap);
    if (d1 <= 0.0 && d2 <= 0.0)
    {
        return a;
    }

    const Eigen::Vector3d bp = p - b;
    const double d3          = ab.dot(bp);
    const double d4          = ac.dot(bp);
    if (d3 >= 0.0 && d4 <= d3)
    {
        return b;
    }

    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
    {
        const double denom = d1 - d3;
        return (denom > 1e-30) ? (a + (d1 / denom) * ab) : a;
    }

    const Eigen::Vector3d cp = p - c;
    const double d5          = ab.dot(cp);
    const double d6          = ac.dot(cp);
    if (d6 >= 0.0 && d5 <= d6)
    {
        return c;
    }

    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
    {
        const double denom = d2 - d6;
        return (denom > 1e-30) ? (a + (d2 / denom) * ac) : a;
    }

    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
    {
        const double denom = (d4 - d3) + (d5 - d6);
        return (denom > 1e-30) ? (b + ((d4 - d3) / denom) * (c - b)) : b;
    }

    const double denom = va + vb + vc;
    if (denom > 1e-30)
    {
        const double v = vb / denom;
        const double w = vc / denom;
        return a + ab * v + ac * w;
    }

    // Degenerate (collinear) triangle fallback: minimum over the 3 edges
    const Eigen::Vector3d q0 = closest_point_on_segment_vec(a, b, p);
    const Eigen::Vector3d q1 = closest_point_on_segment_vec(b, c, p);
    const Eigen::Vector3d q2 = closest_point_on_segment_vec(c, a, p);
    double best_d2           = (p - q0).squaredNorm();
    Eigen::Vector3d best_q   = q0;
    const double d2_1        = (p - q1).squaredNorm();
    if (d2_1 < best_d2)
    {
        best_d2 = d2_1;
        best_q  = q1;
    }
    const double d2_2 = (p - q2).squaredNorm();
    if (d2_2 < best_d2)
    {
        best_q = q2;
    }
    return best_q;
}

double aabb_point_sq_dist(const Eigen::Vector3d& bmin,
                          const Eigen::Vector3d& bmax,
                          const Eigen::Vector3d& p) noexcept
{
    double d2 = 0.0;
    for (int k = 0; k < 3; ++k)
    {
        if (p[k] < bmin[k])
        {
            const double d = bmin[k] - p[k];
            d2 += d * d;
        }
        else if (p[k] > bmax[k])
        {
            const double d = p[k] - bmax[k];
            d2 += d * d;
        }
    }
    return d2;
}

void bfs_orient_manifold_patches(const Eigen::MatrixXi& F,
                                 Eigen::MatrixXi& FF,
                                 Eigen::VectorXi& C)
{
    const int m = static_cast<int>(F.rows());
    FF          = F;
    C.setConstant(m, -1);
    if (m == 0)
    {
        return;
    }

    // Map each undirected edge {min_v, max_v} to incident {face_idx, dir (+1 if u<v else -1)}
    std::unordered_map<uint64_t, std::vector<std::pair<int, int>>> edge_to_faces;
    edge_to_faces.reserve(static_cast<std::size_t>(m) * 2);

    for (int f = 0; f < m; ++f)
    {
        for (int e = 0; e < 3; ++e)
        {
            const int u       = F(f, e);
            const int v       = F(f, (e + 1) % 3);
            const uint32_t lo = static_cast<uint32_t>(std::min(u, v));
            const uint32_t hi = static_cast<uint32_t>(std::max(u, v));
            const uint64_t k  = static_cast<uint64_t>(lo) | (static_cast<uint64_t>(hi) << 32);
            edge_to_faces[k].emplace_back(f, (u < v) ? +1 : -1);
        }
    }

    // Build 2-manifold face adjacency (do not cross non-manifold edges with degree != 2)
    std::vector<std::vector<int>> adj(static_cast<std::size_t>(m));
    for (const auto& kv : edge_to_faces)
    {
        if (kv.second.size() == 2)
        {
            const int f0 = kv.second[0].first;
            const int f1 = kv.second[1].first;
            if (f0 != f1)
            {
                adj[static_cast<std::size_t>(f0)].push_back(f1);
                adj[static_cast<std::size_t>(f1)].push_back(f0);
            }
        }
    }

    int comp_id = 0;
    std::vector<bool> flipped(static_cast<std::size_t>(m), false);
    for (int seed = 0; seed < m; ++seed)
    {
        if (C(seed) != -1)
        {
            continue;
        }
        C(seed) = comp_id;
        std::deque<int> q = { seed };
        while (!q.empty())
        {
            const int cur = q.front();
            q.pop_front();

            for (const int nb : adj[static_cast<std::size_t>(cur)])
            {
                if (C(nb) != -1)
                {
                    continue;
                }
                C(nb) = comp_id;

                // Check if `cur` and `nb` share a half-edge in the SAME direction (requiring flip)
                bool same_he_dir = false;
                for (int ec = 0; ec < 3 && !same_he_dir; ++ec)
                {
                    const int cu = FF(cur, ec);
                    const int cv = FF(cur, (ec + 1) % 3);
                    for (int en = 0; en < 3; ++en)
                    {
                        if (FF(nb, en) == cu && FF(nb, (en + 1) % 3) == cv)
                        {
                            same_he_dir = true;
                            break;
                        }
                    }
                }
                if (same_he_dir)
                {
                    std::swap(FF(nb, 1), FF(nb, 2));
                    flipped[static_cast<std::size_t>(nb)] = true;
                }
                q.push_back(nb);
            }
        }
        ++comp_id;
    }
}

} // namespace

void TriangleAABB::init(const Eigen::MatrixX3d& V, const Eigen::MatrixX3i& F)
{
    nodes_.clear();
    const int m = static_cast<int>(F.rows());
    tri_order_.resize(static_cast<std::size_t>(m));
    if (m == 0)
    {
        return;
    }

    std::iota(tri_order_.begin(), tri_order_.end(), 0);
    std::vector<Eigen::Vector3d> centroids(static_cast<std::size_t>(m));
    for (int i = 0; i < m; ++i)
    {
        centroids[static_cast<std::size_t>(i)] =
          (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))).transpose() / 3.0;
    }

    nodes_.reserve(static_cast<std::size_t>(2 * m));
    build_recursive(V, F, centroids, 0, m);
}

int TriangleAABB::build_recursive(const Eigen::MatrixX3d& V,
                                  const Eigen::MatrixX3i& F,
                                  const std::vector<Eigen::Vector3d>& centroids,
                                  int start,
                                  int count)
{
    const int node_idx = static_cast<int>(nodes_.size());
    nodes_.emplace_back();

    Eigen::Vector3d bmin = Eigen::Vector3d::Constant(std::numeric_limits<double>::infinity());
    Eigen::Vector3d bmax = Eigen::Vector3d::Constant(-std::numeric_limits<double>::infinity());
    Eigen::Vector3d cmin = bmin;
    Eigen::Vector3d cmax = bmax;

    for (int i = start; i < start + count; ++i)
    {
        const int tid = tri_order_[static_cast<std::size_t>(i)];
        for (int c = 0; c < 3; ++c)
        {
            const Eigen::Vector3d pv = V.row(F(tid, c)).transpose();
            bmin                     = bmin.cwiseMin(pv);
            bmax                     = bmax.cwiseMax(pv);
        }
        const Eigen::Vector3d& ctr = centroids[static_cast<std::size_t>(tid)];
        cmin                       = cmin.cwiseMin(ctr);
        cmax                       = cmax.cwiseMax(ctr);
    }

    nodes_[static_cast<std::size_t>(node_idx)].bmin  = bmin;
    nodes_[static_cast<std::size_t>(node_idx)].bmax  = bmax;
    nodes_[static_cast<std::size_t>(node_idx)].start = start;
    nodes_[static_cast<std::size_t>(node_idx)].count = count;

    if (count <= 4)
    {
        return node_idx;
    }

    int split_axis            = 0;
    const Eigen::Vector3d ext = cmax - cmin;
    if (ext[1] > ext[0] && ext[1] >= ext[2])
    {
        split_axis = 1;
    }
    else if (ext[2] > ext[0] && ext[2] > ext[1])
    {
        split_axis = 2;
    }

    const int mid = start + count / 2;
    std::nth_element(tri_order_.begin() + start,
                     tri_order_.begin() + mid,
                     tri_order_.begin() + start + count,
                     [&](int a, int b)
                     {
                         return centroids[static_cast<std::size_t>(a)][split_axis] <
                                centroids[static_cast<std::size_t>(b)][split_axis];
                     });

    const int left_idx  = build_recursive(V, F, centroids, start, mid - start);
    const int right_idx = build_recursive(V, F, centroids, mid, start + count - mid);

    nodes_[static_cast<std::size_t>(node_idx)].left  = left_idx;
    nodes_[static_cast<std::size_t>(node_idx)].right = right_idx;
    return node_idx;
}

Eigen::Vector3d TriangleAABB::closest_point(const Eigen::MatrixX3d& V,
                                            const Eigen::MatrixX3i& F,
                                            const Eigen::Vector3d& p) const
{
    if (nodes_.empty())
    {
        throw std::runtime_error("TriangleAABB::closest_point called on an empty tree.");
    }

    double best_d2         = std::numeric_limits<double>::infinity();
    Eigen::Vector3d best_q = p;

    int stack[128];
    int sp    = 0;
    stack[sp++] = 0;

    while (sp > 0)
    {
        const int n_idx   = stack[--sp];
        const Node& node  = nodes_[static_cast<std::size_t>(n_idx)];
        const double box2 = aabb_point_sq_dist(node.bmin, node.bmax, p);
        if (box2 >= best_d2)
        {
            continue;
        }

        if (node.left < 0)
        {
            for (int i = node.start; i < node.start + node.count; ++i)
            {
                const int tid          = tri_order_[static_cast<std::size_t>(i)];
                const Eigen::Vector3d q =
                  closest_point_on_triangle(V.row(F(tid, 0)).transpose(),
                                            V.row(F(tid, 1)).transpose(),
                                            V.row(F(tid, 2)).transpose(),
                                            p);
                const double d2 = (p - q).squaredNorm();
                if (d2 < best_d2)
                {
                    best_d2 = d2;
                    best_q  = q;
                }
            }
        }
        else
        {
            const Node& ln   = nodes_[static_cast<std::size_t>(node.left)];
            const Node& rn   = nodes_[static_cast<std::size_t>(node.right)];
            const double dl2 = aabb_point_sq_dist(ln.bmin, ln.bmax, p);
            const double dr2 = aabb_point_sq_dist(rn.bmin, rn.bmax, p);

            if (dl2 < dr2)
            {
                if (dr2 < best_d2 && sp < 127)
                {
                    stack[sp++] = node.right;
                }
                if (dl2 < best_d2 && sp < 127)
                {
                    stack[sp++] = node.left;
                }
            }
            else
            {
                if (dl2 < best_d2 && sp < 127)
                {
                    stack[sp++] = node.left;
                }
                if (dr2 < best_d2 && sp < 127)
                {
                    stack[sp++] = node.right;
                }
            }
        }
    }

    return best_q;
}

IGLTriSurf triangulation_to_igl(const Triangulation& t, const Marker& marker_to_extract)
{
    Index curr_index              = 0;
    Index n_elements_in_component = 0;
    std::vector<Index> m_to_c_node_map(t.vertices.size(), NOT_SET_Index);
    std::vector<Index> c_to_m_node_map;
    c_to_m_node_map.reserve(t.vertices.size());

    for (std::size_t i = 0; i < t.triangles.size(); ++i)
    {
        if (t.tri_to_marker_map[i] == marker_to_extract)
        {
            n_elements_in_component++;
            const auto& conn = t.triangles[i];
            for (int c = 0; c < 3; ++c)
            {
                if (m_to_c_node_map[conn[c]] == NOT_SET_Index)
                {
                    m_to_c_node_map[conn[c]] = curr_index++;
                    c_to_m_node_map.push_back(conn[c]);
                }
            }
        }
    }

    IGLTriSurf surf;
    surf.vertices.resize(static_cast<Eigen::Index>(curr_index), 3);
    for (std::size_t i = 0; i < c_to_m_node_map.size(); ++i)
    {
        const auto& p = t.vertices[c_to_m_node_map[i]];
        for (int j = 0; j < 3; ++j)
        {
            surf.vertices(static_cast<Eigen::Index>(i), j) = p[j];
        }
    }

    Index out_row = 0;
    surf.triangles.resize(static_cast<Eigen::Index>(n_elements_in_component), 3);
    for (std::size_t i = 0; i < t.triangles.size(); ++i)
    {
        if (t.tri_to_marker_map[i] == marker_to_extract)
        {
            const auto& conn = t.triangles[i];
            for (int c = 0; c < 3; ++c)
            {
                surf.triangles(static_cast<Eigen::Index>(out_row), c) =
                  static_cast<int>(m_to_c_node_map[conn[c]]);
            }
            out_row++;
        }
    }

    if (surf.triangles.rows() > 0)
    {
        surf.tree.init(surf.vertices, surf.triangles);
    }

    return surf;
}

std::vector<IGLTriSurf> decompose_into_igl_facets(const Triangulation& combined_surface_mesh)
{
    std::vector<IGLTriSurf> facets;
    facets.reserve(combined_surface_mesh.n_markers);
    for (Marker fid = 0; fid < combined_surface_mesh.n_markers; ++fid)
    {
        facets.push_back(triangulation_to_igl(combined_surface_mesh, fid));
    }
    return facets;
}

IGLToTriangulationResult igl_to_triangulation(const Eigen::MatrixXd& vertices,
                                              const Eigen::MatrixXi& triangles)
{
    IGLToTriangulationResult r;
    r.t.vertices.resize(static_cast<std::size_t>(vertices.rows()));
    for (Eigen::Index i = 0; i < vertices.rows(); ++i)
    {
        r.t.vertices[static_cast<std::size_t>(i)] = { vertices(i, 0), vertices(i, 1), vertices(i, 2) };
    }

    r.element_to_origin_tris_map.reserve(static_cast<std::size_t>(triangles.rows()));
    std::unordered_map<std::array<Index, 3>, Index, SortedTriKeyHash> inserted_triangles;
    r.t.triangles.reserve(static_cast<std::size_t>(triangles.rows()));

    for (Eigen::Index i = 0; i < triangles.rows(); ++i)
    {
        std::array<Index, 3> conn  = { static_cast<Index>(triangles(i, 0)),
                                       static_cast<Index>(triangles(i, 1)),
                                       static_cast<Index>(triangles(i, 2)) };
        std::array<Index, 3> s_key = conn;
        std::sort(s_key.begin(), s_key.end());
        auto it = inserted_triangles.find(s_key);
        if (it == inserted_triangles.end())
        {
            inserted_triangles[s_key] = static_cast<Index>(r.t.triangles.size());
            r.t.triangles.push_back(conn);
            r.element_to_origin_tris_map.push_back({ static_cast<Index>(i) });
        }
        else
        {
            r.element_to_origin_tris_map[it->second].push_back(static_cast<Index>(i));
        }
    }

    r.t.tri_to_marker_map.assign(r.t.triangles.size(), 0);
    update_edges_and_maps(r.t);
    return r;
}

Eigen::Vector3d igl_closest_point(const IGLTriSurf& igl_surf, const Eigen::Vector3d& p)
{
    if (igl_surf.triangles.rows() == 0 || igl_surf.tree.empty())
    {
        throw std::runtime_error("Cannot compute closest point on an empty IGLTriSurf.");
    }
    return igl_surf.tree.closest_point(igl_surf.vertices, igl_surf.triangles, p);
}

double signed_volume(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
    double volume_total             = 0.0;
    const Eigen::RowVector3d origin = Eigen::RowVector3d::Zero();
    for (int f = 0; f < F.rows(); ++f)
    {
        const Eigen::RowVector3d t0 = V.row(F(f, 0));
        const Eigen::RowVector3d t1 = V.row(F(f, 1));
        const Eigen::RowVector3d t2 = V.row(F(f, 2));
        volume_total += (t1 - t0).dot((t2 - t0).cross(origin - t0)) / 6.0;
    }
    return -volume_total;
}

void orient_closed_surface(const Eigen::MatrixXd& V, Eigen::MatrixXi& F, bool positive)
{
    if ((positive ? 1.0 : -1.0) * signed_volume(V, F) < 0.0)
    {
        for (int f = 0; f < F.rows(); ++f)
        {
            F.row(f) = F.row(f).reverse().eval();
        }
    }
}

void orient_volume(const Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    if (F.rows() == 0)
    {
        return;
    }

    Eigen::MatrixXi F_tmp;
    Eigen::VectorXi C;
    bfs_orient_manifold_patches(F, F_tmp, C);

    const int n_comp = (C.size() > 0) ? (C.maxCoeff() + 1) : 0;
    if (n_comp <= 1)
    {
        orient_closed_surface(V, F_tmp, true);
        F = F_tmp;
        return;
    }

    std::vector<std::vector<int>> comp_faces(static_cast<std::size_t>(n_comp));
    for (int f = 0; f < F_tmp.rows(); ++f)
    {
        comp_faces[static_cast<std::size_t>(C(f))].push_back(f);
    }

    for (int c = 0; c < n_comp; ++c)
    {
        double vol_c = 0.0;
        for (int f : comp_faces[static_cast<std::size_t>(c)])
        {
            const Eigen::RowVector3d p0 = V.row(F_tmp(f, 0));
            const Eigen::RowVector3d p1 = V.row(F_tmp(f, 1));
            const Eigen::RowVector3d p2 = V.row(F_tmp(f, 2));
            vol_c += p0.dot(p1.cross(p2)) / 6.0;
        }
        if (vol_c < 0.0)
        {
            for (int f : comp_faces[static_cast<std::size_t>(c)])
            {
                std::swap(F_tmp(f, 1), F_tmp(f, 2));
            }
        }
    }

    // Exact Van Oosterom-Strackee solid-angle shell nesting parity:
    // even depth -> exterior boundary (+volume), odd depth -> interior cavity (-volume).
    constexpr double kFourPi = 12.56637061435917295385;
    for (int c = 0; c < n_comp; ++c)
    {
        const auto& cf = comp_faces[static_cast<std::size_t>(c)];
        if (cf.empty())
        {
            continue;
        }
        const int f0 = cf[0];
        const Eigen::Vector3d query =
          (V.row(F_tmp(f0, 0)) + V.row(F_tmp(f0, 1)) + V.row(F_tmp(f0, 2))).transpose() / 3.0;

        int nesting_depth = 0;
        for (int other = 0; other < n_comp; ++other)
        {
            const auto& of = comp_faces[static_cast<std::size_t>(other)];
            if (other == c || of.empty())
            {
                continue;
            }
            double solid_angle = 0.0;
            for (int f : of)
            {
                const Eigen::Vector3d a = V.row(F_tmp(f, 0)).transpose() - query;
                const Eigen::Vector3d b = V.row(F_tmp(f, 1)).transpose() - query;
                const Eigen::Vector3d d = V.row(F_tmp(f, 2)).transpose() - query;
                const double la         = a.norm();
                const double lb         = b.norm();
                const double ld         = d.norm();
                if (la < 1e-14 || lb < 1e-14 || ld < 1e-14)
                {
                    continue;
                }
                const double num = a.dot(b.cross(d));
                const double den = la * lb * ld + a.dot(b) * ld + b.dot(d) * la + d.dot(a) * lb;
                solid_angle += 2.0 * std::atan2(num, den);
            }
            if (solid_angle > 0.5 * kFourPi)
            {
                ++nesting_depth;
            }
        }

        if ((nesting_depth % 2) == 1)
        {
            for (int f : cf)
            {
                std::swap(F_tmp(f, 1), F_tmp(f, 2));
            }
        }
    }

    F = F_tmp;
}

} // namespace mesh_arrange::detail
