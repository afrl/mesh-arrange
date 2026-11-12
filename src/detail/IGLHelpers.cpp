#include "detail/IGLHelpers.h"

#include "igl/bfs_orient.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_map>

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
} // namespace

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
    if (igl_surf.triangles.rows() == 0)
    {
        throw std::runtime_error("Cannot compute closest point on an empty IGLTriSurf.");
    }
    Eigen::MatrixX3d point(1, 3);
    point.row(0) = p.transpose();

    Eigen::VectorXd squared_distance;
    Eigen::VectorXi closest_triangle;
    Eigen::MatrixX3d closest_point;
    igl_surf.tree.squared_distance(igl_surf.vertices,
                                   igl_surf.triangles,
                                   point,
                                   squared_distance,
                                   closest_triangle,
                                   closest_point);
    return closest_point.row(0).transpose();
}

double signed_volume(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
    double volume_total = 0.0;
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

    Eigen::MatrixXi F_tmp = F;
    Eigen::VectorXi C;
    igl::bfs_orient(F, F_tmp, C);

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
