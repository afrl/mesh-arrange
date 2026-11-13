#pragma once

#include "mesh_arrange/MeshArrange.h"
#include "mesh_arrange/PlcBackend.h"

#include <gmsh.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace mesh_arrange
{
namespace gmsh_sdk
{

namespace detail
{

struct ScopedGmshPlcState
{
    std::string prev_model;
    std::string scratch_model;
    double abort_err = 2.0;
    double alg3d     = 1.0;
    double opt       = 1.0;
    double optnet    = 0.0;
    double maxsz     = 1e22;
    double minsz     = 0.0;
    double extbnd    = 1.0;
    double renum     = 1.0;
    bool active      = false;

    explicit ScopedGmshPlcState(const std::string& scratch)
        : scratch_model(scratch)
    {
        gmsh::model::getCurrent(prev_model);
        gmsh::option::getNumber("General.AbortOnError", abort_err);
        gmsh::option::getNumber("Mesh.Algorithm3D", alg3d);
        gmsh::option::getNumber("Mesh.Optimize", opt);
        gmsh::option::getNumber("Mesh.OptimizeNetgen", optnet);
        gmsh::option::getNumber("Mesh.MeshSizeMax", maxsz);
        gmsh::option::getNumber("Mesh.MeshSizeMin", minsz);
        gmsh::option::getNumber("Mesh.MeshSizeExtendFromBoundary", extbnd);
        gmsh::option::getNumber("Mesh.Renumber", renum);

        gmsh::model::add(scratch_model);
        gmsh::model::setCurrent(scratch_model);
        active = true;
    }

    void restore_options() const
    {
        gmsh::option::setNumber("General.AbortOnError", abort_err);
        gmsh::option::setNumber("Mesh.Algorithm3D", alg3d);
        gmsh::option::setNumber("Mesh.Optimize", opt);
        gmsh::option::setNumber("Mesh.OptimizeNetgen", optnet);
        gmsh::option::setNumber("Mesh.MeshSizeMax", maxsz);
        gmsh::option::setNumber("Mesh.MeshSizeMin", minsz);
        gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", extbnd);
        gmsh::option::setNumber("Mesh.Renumber", renum);
    }

    ~ScopedGmshPlcState()
    {
        if (active)
        {
            try
            {
                restore_options();
                gmsh::model::setCurrent(scratch_model);
                gmsh::model::remove();
                if (!prev_model.empty() && prev_model != scratch_model)
                {
                    gmsh::model::setCurrent(prev_model);
                }
            }
            catch (...)
            {
            }
        }
    }
};

struct SortedTriHash
{
    std::size_t operator()(const std::array<Index, 3>& k) const noexcept
    {
        std::size_t h = static_cast<std::size_t>(k[0]);
        h ^= static_cast<std::size_t>(k[1]) + 0x9e3779b9u + (h << 6) + (h >> 2);
        h ^= static_cast<std::size_t>(k[2]) + 0x9e3779b9u + (h << 6) + (h >> 2);
        return h;
    }
};

inline double point_triangle_dist_sq(const Point3& p,
                                     const Point3& a,
                                     const Point3& b,
                                     const Point3& c) noexcept
{
    const double ab[3] = { b[0] - a[0], b[1] - a[1], b[2] - a[2] };
    const double ac[3] = { c[0] - a[0], c[1] - a[1], c[2] - a[2] };
    const double ap[3] = { p[0] - a[0], p[1] - a[1], p[2] - a[2] };
    const double d1    = ab[0] * ap[0] + ab[1] * ap[1] + ab[2] * ap[2];
    const double d2    = ac[0] * ap[0] + ac[1] * ap[1] + ac[2] * ap[2];
    if (d1 <= 0.0 && d2 <= 0.0)
    {
        return ap[0] * ap[0] + ap[1] * ap[1] + ap[2] * ap[2];
    }

    const double bp[3] = { p[0] - b[0], p[1] - b[1], p[2] - b[2] };
    const double d3    = ab[0] * bp[0] + ab[1] * bp[1] + ab[2] * bp[2];
    const double d4    = ac[0] * bp[0] + ac[1] * bp[1] + ac[2] * bp[2];
    if (d3 >= 0.0 && d4 <= d3)
    {
        return bp[0] * bp[0] + bp[1] * bp[1] + bp[2] * bp[2];
    }

    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
    {
        const double v  = d1 / (d1 - d3);
        const double dx = ap[0] - v * ab[0];
        const double dy = ap[1] - v * ab[1];
        const double dz = ap[2] - v * ab[2];
        return dx * dx + dy * dy + dz * dz;
    }

    const double cp[3] = { p[0] - c[0], p[1] - c[1], p[2] - c[2] };
    const double d5    = ab[0] * cp[0] + ab[1] * cp[1] + ab[2] * cp[2];
    const double d6    = ac[0] * cp[0] + ac[1] * cp[1] + ac[2] * cp[2];
    if (d6 >= 0.0 && d5 <= d6)
    {
        return cp[0] * cp[0] + cp[1] * cp[1] + cp[2] * cp[2];
    }

    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
    {
        const double w  = d2 / (d2 - d6);
        const double dx = ap[0] - w * ac[0];
        const double dy = ap[1] - w * ac[1];
        const double dz = ap[2] - w * ac[2];
        return dx * dx + dy * dy + dz * dz;
    }

    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
    {
        const double w  = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        const double qx = b[0] + w * (c[0] - b[0]) - p[0];
        const double qy = b[1] + w * (c[1] - b[1]) - p[1];
        const double qz = b[2] + w * (c[2] - b[2]) - p[2];
        return qx * qx + qy * qy + qz * qz;
    }

    const double denom = 1.0 / (va + vb + vc);
    const double v     = vb * denom;
    const double w     = vc * denom;
    const double dx    = ap[0] - v * ab[0] - w * ac[0];
    const double dy    = ap[1] - v * ab[1] - w * ac[1];
    const double dz    = ap[2] - v * ab[2] - w * ac[2];
    return dx * dx + dy * dy + dz * dz;
}

inline int next_entity_tag(int dim)
{
    std::vector<std::pair<int, int>> ents;
    gmsh::model::getEntities(ents, dim);
    int max_tag = 0;
    for (std::size_t i = 0; i < ents.size(); ++i)
    {
        if (ents[i].second > max_tag)
        {
            max_tag = ents[i].second;
        }
    }
    return max_tag + 1;
}

inline std::size_t next_node_tag()
{
    std::size_t max_tag = 0;
    gmsh::model::mesh::getMaxNodeTag(max_tag);
    return max_tag + 1;
}

inline void orient_closed_shell_outward(InputSurface& s)
{
    const std::size_t n_tris = s.triangles.size();
    if (n_tris == 0)
    {
        return;
    }

    struct EdgeAdj
    {
        std::size_t tri;
        int dir; // +1 if lo->hi, -1 if hi->lo
    };
    std::unordered_map<std::uint64_t, std::vector<EdgeAdj>> edge_map;
    edge_map.reserve(n_tris * 2);
    for (std::size_t t = 0; t < n_tris; ++t)
    {
        const auto& tr = s.triangles[t];
        for (int e = 0; e < 3; ++e)
        {
            const Index u          = tr[e];
            const Index v          = tr[(e + 1) % 3];
            const Index lo         = std::min(u, v);
            const Index hi         = std::max(u, v);
            const std::uint64_t key =
              (static_cast<std::uint64_t>(lo) << 32) | static_cast<std::uint32_t>(hi);
            edge_map[key].push_back({ t, (u == lo) ? +1 : -1 });
        }
    }

    std::vector<int> flip_state(n_tris, 0); // 0 = unvisited, +1 = keep, -1 = flip
    std::vector<std::size_t> queue;
    std::vector<std::size_t> comp;
    for (std::size_t seed = 0; seed < n_tris; ++seed)
    {
        if (flip_state[seed] != 0)
        {
            continue;
        }
        flip_state[seed] = +1;
        queue.clear();
        comp.clear();
        queue.push_back(seed);
        std::size_t head = 0;
        while (head < queue.size())
        {
            const std::size_t cur = queue[head++];
            comp.push_back(cur);
            const auto& tr = s.triangles[cur];
            for (int e = 0; e < 3; ++e)
            {
                const Index u          = tr[e];
                const Index v          = tr[(e + 1) % 3];
                const Index lo         = std::min(u, v);
                const Index hi         = std::max(u, v);
                const int cur_dir      = (u == lo) ? +1 : -1;
                const std::uint64_t key =
                  (static_cast<std::uint64_t>(lo) << 32) | static_cast<std::uint32_t>(hi);
                const auto& adj = edge_map[key];
                if (adj.size() != 2)
                {
                    continue;
                }
                for (std::size_t k = 0; k < 2; ++k)
                {
                    const std::size_t nb = adj[k].tri;
                    if (nb == cur || flip_state[nb] != 0)
                    {
                        continue;
                    }
                    // Consistent orientation across shared edge requires opposite traversal:
                    // (cur_dir * flip_state[cur]) == -(adj[k].dir * flip_state[nb])
                    flip_state[nb] = -(cur_dir * flip_state[cur] * adj[k].dir);
                    queue.push_back(nb);
                }
            }
        }

        double vol6 = 0.0;
        for (std::size_t i = 0; i < comp.size(); ++i)
        {
            const std::size_t ti = comp[i];
            auto tr              = s.triangles[ti];
            if (flip_state[ti] < 0)
            {
                std::swap(tr[1], tr[2]);
            }
            const Point3& a = s.vertices[tr[0]];
            const Point3& b = s.vertices[tr[1]];
            const Point3& c = s.vertices[tr[2]];
            vol6 += a[0] * (b[1] * c[2] - b[2] * c[1]) +
                    a[1] * (b[2] * c[0] - b[0] * c[2]) +
                    a[2] * (b[0] * c[1] - b[1] * c[0]);
        }
        const int global_sign = (vol6 < 0.0) ? -1 : +1;
        for (std::size_t i = 0; i < comp.size(); ++i)
        {
            const std::size_t ti = comp[i];
            if (flip_state[ti] * global_sign < 0)
            {
                std::swap(s.triangles[ti][1], s.triangles[ti][2]);
            }
        }
    }
}

} // namespace detail

/// @brief Native Gmsh `PlcTetrahedralizerFn` replacing TetGen (`tetgen::tetrahedralize("pYQ")`).
///        Encloses the non-manifold `plc` inside a padded outer bounding-box `Volume` with
///        `gmsh::model::mesh::embed(2, {plc_surf}, 3, box_vol)` and recovers all constrained
///        facets (including any Steiner-subdivided sub-triangles) via Gmsh's 3D Delaunay
///        boundary recovery (`Mesh.Algorithm3D = 1`, `Mesh.Renumber = 0`).
inline PlcTetOutput plc_tetrahedralizer(const Triangulation& plc)
{
    if (plc.vertices.empty() || plc.triangles.empty())
    {
        throw std::invalid_argument("gmsh_sdk::plc_tetrahedralizer received an empty PLC.");
    }

    detail::ScopedGmshPlcState guard("__mesh_arrange_plc_scratch__");

    Point3 bmin = plc.vertices.front();
    Point3 bmax = plc.vertices.front();
    for (std::size_t i = 0; i < plc.vertices.size(); ++i)
    {
        for (int k = 0; k < 3; ++k)
        {
            bmin[k] = std::min(bmin[k], plc.vertices[i][k]);
            bmax[k] = std::max(bmax[k], plc.vertices[i][k]);
        }
    }
    double diag = 0.0;
    for (int k = 0; k < 3; ++k)
    {
        const double d = bmax[k] - bmin[k];
        diag += d * d;
    }
    diag             = std::max(std::sqrt(diag), 1.0);
    const double pad = 0.35 * diag;
    for (int k = 0; k < 3; ++k)
    {
        bmin[k] -= pad;
        bmax[k] += pad;
    }

    // Collapse sub-tolerance ULP edges (< 1e-9 * diag, e.g. floating-point intersections at
    // -0.7 + 0.2 vs -0.5) so Gmsh's initial 3D Delaunay (ToleranceInitialDelaunay = 1e-12)
    // never drops coincident nodes during boundary recovery.
    const std::size_t n_raw   = plc.vertices.size();
    const double weld_tol2    = (1e-9 * diag) * (1e-9 * diag);
    std::vector<Index> uf_par(n_raw);
    for (std::size_t i = 0; i < n_raw; ++i)
    {
        uf_par[i] = static_cast<Index>(i);
    }
    auto uf_find = [&](Index x) -> Index
    {
        while (uf_par[x] != x)
        {
            uf_par[x] = uf_par[uf_par[x]];
            x         = uf_par[x];
        }
        return x;
    };

    bool any_collapsed = false;
    for (std::size_t i = 0; i < plc.triangles.size(); ++i)
    {
        const auto& tr = plc.triangles[i];
        for (int e = 0; e < 3; ++e)
        {
            const Index ra = uf_find(tr[e]);
            const Index rb = uf_find(tr[(e + 1) % 3]);
            if (ra != rb)
            {
                const double dx = plc.vertices[ra][0] - plc.vertices[rb][0];
                const double dy = plc.vertices[ra][1] - plc.vertices[rb][1];
                const double dz = plc.vertices[ra][2] - plc.vertices[rb][2];
                if (dx * dx + dy * dy + dz * dz <= weld_tol2)
                {
                    const Index lo  = std::min(ra, rb);
                    const Index hi  = std::max(ra, rb);
                    uf_par[hi]      = lo;
                    any_collapsed   = true;
                }
            }
        }
    }

    std::vector<Point3> clean_vertices;
    std::vector<std::array<Index, 3>> clean_triangles;
    std::vector<Index> clean_parent_ids;
    if (!any_collapsed)
    {
        clean_vertices  = plc.vertices;
        clean_triangles = plc.triangles;
        clean_parent_ids.resize(plc.triangles.size());
        for (std::size_t i = 0; i < plc.triangles.size(); ++i)
        {
            clean_parent_ids[i] = static_cast<Index>(i);
        }
    }
    else
    {
        std::vector<Index> raw_to_clean(n_raw, NOT_SET_Index);
        clean_vertices.reserve(n_raw);
        for (std::size_t i = 0; i < n_raw; ++i)
        {
            const Index r = uf_find(static_cast<Index>(i));
            if (raw_to_clean[r] == NOT_SET_Index)
            {
                raw_to_clean[r] = static_cast<Index>(clean_vertices.size());
                clean_vertices.push_back(plc.vertices[r]);
            }
            raw_to_clean[i] = raw_to_clean[r];
        }

        std::unordered_map<std::array<Index, 3>, bool, detail::SortedTriHash> seen_tris;
        seen_tris.reserve(plc.triangles.size());
        const double min_area2 = (1e-13 * diag * diag) * (1e-13 * diag * diag);
        for (std::size_t p = 0; p < plc.triangles.size(); ++p)
        {
            const Index c0 = raw_to_clean[plc.triangles[p][0]];
            const Index c1 = raw_to_clean[plc.triangles[p][1]];
            const Index c2 = raw_to_clean[plc.triangles[p][2]];
            if (c0 == c1 || c1 == c2 || c2 == c0)
            {
                continue;
            }
            const Point3& A  = clean_vertices[c0];
            const Point3& B  = clean_vertices[c1];
            const Point3& C  = clean_vertices[c2];
            const double nx  = (B[1] - A[1]) * (C[2] - A[2]) - (B[2] - A[2]) * (C[1] - A[1]);
            const double ny  = (B[2] - A[2]) * (C[0] - A[0]) - (B[0] - A[0]) * (C[2] - A[2]);
            const double nz  = (B[0] - A[0]) * (C[1] - A[1]) - (B[1] - A[1]) * (C[0] - A[0]);
            if (nx * nx + ny * ny + nz * nz <= min_area2)
            {
                continue;
            }
            std::array<Index, 3> ks = { c0, c1, c2 };
            std::sort(ks.begin(), ks.end());
            if (seen_tris.emplace(ks, true).second)
            {
                clean_triangles.push_back({ c0, c1, c2 });
                clean_parent_ids.push_back(static_cast<Index>(p));
            }
        }
    }

    const int box_surf = gmsh::model::addDiscreteEntity(2, 1);
    const int plc_surf = gmsh::model::addDiscreteEntity(2, 2);
    const int box_sl   = gmsh::model::geo::addSurfaceLoop({ box_surf });
    const int box_vol  = gmsh::model::geo::addVolume({ box_sl }, 1);
    gmsh::model::geo::synchronize();

    const std::size_t n_plc = clean_vertices.size();
    std::vector<std::size_t> plc_ntags(n_plc);
    std::vector<double> plc_coords;
    plc_coords.reserve(n_plc * 3);
    for (std::size_t i = 0; i < n_plc; ++i)
    {
        plc_ntags[i] = i + 1;
        plc_coords.push_back(clean_vertices[i][0]);
        plc_coords.push_back(clean_vertices[i][1]);
        plc_coords.push_back(clean_vertices[i][2]);
    }
    gmsh::model::mesh::addNodes(2, plc_surf, plc_ntags, plc_coords);

    std::vector<std::size_t> plc_tri_nodes;
    plc_tri_nodes.reserve(clean_triangles.size() * 3);
    for (std::size_t i = 0; i < clean_triangles.size(); ++i)
    {
        plc_tri_nodes.push_back(1 + static_cast<std::size_t>(clean_triangles[i][0]));
        plc_tri_nodes.push_back(1 + static_cast<std::size_t>(clean_triangles[i][1]));
        plc_tri_nodes.push_back(1 + static_cast<std::size_t>(clean_triangles[i][2]));
    }
    gmsh::model::mesh::addElementsByType(plc_surf, 2, {}, plc_tri_nodes);

    const std::array<Point3, 8> box_pts = {{
        { bmin[0], bmin[1], bmin[2] }, { bmax[0], bmin[1], bmin[2] },
        { bmax[0], bmax[1], bmin[2] }, { bmin[0], bmax[1], bmin[2] },
        { bmin[0], bmin[1], bmax[2] }, { bmax[0], bmin[1], bmax[2] },
        { bmax[0], bmax[1], bmax[2] }, { bmin[0], bmax[1], bmax[2] }
    }};
    const std::size_t b0 = n_plc + 1;
    std::vector<std::size_t> box_ntags(8);
    std::vector<double> box_coords;
    box_coords.reserve(24);
    for (std::size_t i = 0; i < 8; ++i)
    {
        box_ntags[i] = b0 + i;
        box_coords.push_back(box_pts[i][0]);
        box_coords.push_back(box_pts[i][1]);
        box_coords.push_back(box_pts[i][2]);
    }
    gmsh::model::mesh::addNodes(2, box_surf, box_ntags, box_coords);

    const std::array<std::array<std::size_t, 3>, 12> box_tris = {{
        { b0 + 0, b0 + 2, b0 + 1 }, { b0 + 0, b0 + 3, b0 + 2 },
        { b0 + 4, b0 + 5, b0 + 6 }, { b0 + 4, b0 + 6, b0 + 7 },
        { b0 + 0, b0 + 1, b0 + 5 }, { b0 + 0, b0 + 5, b0 + 4 },
        { b0 + 2, b0 + 3, b0 + 7 }, { b0 + 2, b0 + 7, b0 + 6 },
        { b0 + 1, b0 + 2, b0 + 6 }, { b0 + 1, b0 + 6, b0 + 5 },
        { b0 + 3, b0 + 0, b0 + 4 }, { b0 + 3, b0 + 4, b0 + 7 }
    }};
    std::vector<std::size_t> box_tri_nodes;
    box_tri_nodes.reserve(36);
    for (std::size_t i = 0; i < box_tris.size(); ++i)
    {
        box_tri_nodes.insert(box_tri_nodes.end(), box_tris[i].begin(), box_tris[i].end());
    }
    gmsh::model::mesh::addElementsByType(box_surf, 2, {}, box_tri_nodes);
    gmsh::model::mesh::embed(2, { plc_surf }, 3, box_vol);

    gmsh::option::setNumber("General.AbortOnError", 0);
    gmsh::option::setNumber("Mesh.Algorithm3D", 1);
    gmsh::option::setNumber("Mesh.Optimize", 0);
    gmsh::option::setNumber("Mesh.OptimizeNetgen", 0);
    gmsh::option::setNumber("Mesh.MeshSizeMax", 10.0 * diag);
    gmsh::option::setNumber("Mesh.MeshSizeMin", 10.0 * diag);
    gmsh::option::setNumber("Mesh.MeshSizeExtendFromBoundary", 0);
    gmsh::option::setNumber("Mesh.Renumber", 0);

    gmsh::model::mesh::generate(3);
    guard.restore_options();

    PlcTetOutput out;
    out.vertices = clean_vertices;
    std::unordered_map<std::size_t, Index> gmsh_to_local;
    gmsh_to_local.reserve(n_plc + 64);
    for (std::size_t i = 0; i < n_plc; ++i)
    {
        gmsh_to_local[i + 1] = static_cast<Index>(i);
    }

    std::vector<std::size_t> all_ntags;
    std::vector<double> all_coords;
    std::vector<double> all_params;
    gmsh::model::mesh::getNodes(all_ntags, all_coords, all_params);
    for (std::size_t i = 0; i < all_ntags.size(); ++i)
    {
        const std::size_t nt = all_ntags[i];
        if (gmsh_to_local.find(nt) == gmsh_to_local.end())
        {
            const Index idx   = static_cast<Index>(out.vertices.size());
            gmsh_to_local[nt] = idx;
            out.vertices.push_back(
              { all_coords[3 * i + 0], all_coords[3 * i + 1], all_coords[3 * i + 2] });
        }
    }

    std::vector<std::size_t> tet_etags;
    std::vector<std::size_t> tet_enodes;
    gmsh::model::mesh::getElementsByType(4, tet_etags, tet_enodes, box_vol);
    if (tet_etags.empty())
    {
        throw std::runtime_error("gmsh_sdk::plc_tetrahedralizer: 3D boundary recovery failed.");
    }
    out.tetrahedra.reserve(tet_etags.size());
    for (std::size_t i = 0; i < tet_etags.size(); ++i)
    {
        out.tetrahedra.push_back({ gmsh_to_local.at(tet_enodes[4 * i + 0]),
                                   gmsh_to_local.at(tet_enodes[4 * i + 1]),
                                   gmsh_to_local.at(tet_enodes[4 * i + 2]),
                                   gmsh_to_local.at(tet_enodes[4 * i + 3]) });
    }

    std::unordered_map<std::array<Index, 3>, bool, detail::SortedTriHash> tet_faces;
    std::vector<std::array<Index, 3>> unique_tet_faces;
    unique_tet_faces.reserve(out.tetrahedra.size() * 2);
    for (std::size_t i = 0; i < out.tetrahedra.size(); ++i)
    {
        const auto& tet             = out.tetrahedra[i];
        static const int F[4][3]    = { { 0, 2, 1 }, { 0, 3, 2 }, { 1, 2, 3 }, { 0, 1, 3 } };
        for (int f = 0; f < 4; ++f)
        {
            std::array<Index, 3> k  = { tet[F[f][0]], tet[F[f][1]], tet[F[f][2]] };
            std::array<Index, 3> ks = k;
            std::sort(ks.begin(), ks.end());
            if (tet_faces.emplace(ks, true).second)
            {
                unique_tet_faces.push_back(k);
            }
        }
    }

    bool any_split    = false;
    const double tol2 = (1e-6 * diag) * (1e-6 * diag);
    for (std::size_t p = 0; p < clean_triangles.size(); ++p)
    {
        const auto& pt          = clean_triangles[p];
        const Index orig_pid    = clean_parent_ids[p];
        std::array<Index, 3> ks = pt;
        std::sort(ks.begin(), ks.end());
        if (tet_faces.find(ks) != tet_faces.end())
        {
            out.surface_triangles.push_back(pt);
            out.surface_tri_parent_id.push_back(orig_pid);
        }
        else
        {
            any_split       = true;
            const Point3& A = out.vertices[pt[0]];
            const Point3& B = out.vertices[pt[1]];
            const Point3& C = out.vertices[pt[2]];
            const double np[3] = {
                (B[1] - A[1]) * (C[2] - A[2]) - (B[2] - A[2]) * (C[1] - A[1]),
                (B[2] - A[2]) * (C[0] - A[0]) - (B[0] - A[0]) * (C[2] - A[2]),
                (B[0] - A[0]) * (C[1] - A[1]) - (B[1] - A[1]) * (C[0] - A[0])
            };
            std::size_t sub_cnt = 0;
            for (std::size_t j = 0; j < unique_tet_faces.size(); ++j)
            {
                const auto& tf   = unique_tet_faces[j];
                const Point3& P0 = out.vertices[tf[0]];
                const Point3& P1 = out.vertices[tf[1]];
                const Point3& P2 = out.vertices[tf[2]];
                if (detail::point_triangle_dist_sq(P0, A, B, C) <= tol2 &&
                    detail::point_triangle_dist_sq(P1, A, B, C) <= tol2 &&
                    detail::point_triangle_dist_sq(P2, A, B, C) <= tol2)
                {
                    const double nf[3] = {
                        (P1[1] - P0[1]) * (P2[2] - P0[2]) - (P1[2] - P0[2]) * (P2[1] - P0[1]),
                        (P1[2] - P0[2]) * (P2[0] - P0[0]) - (P1[0] - P0[0]) * (P2[2] - P0[2]),
                        (P1[0] - P0[0]) * (P2[1] - P0[1]) - (P1[1] - P0[1]) * (P2[0] - P0[0])
                    };
                    std::array<Index, 3> oriented = tf;
                    if (np[0] * nf[0] + np[1] * nf[1] + np[2] * nf[2] < 0.0)
                    {
                        std::swap(oriented[1], oriented[2]);
                    }
                    out.surface_triangles.push_back(oriented);
                    out.surface_tri_parent_id.push_back(orig_pid);
                    ++sub_cnt;
                }
            }
            if (sub_cnt == 0)
            {
                throw std::runtime_error(
                  "gmsh_sdk::plc_tetrahedralizer failed to recover sub-triangles for PLC facet.");
            }
        }
    }

    if (!any_collapsed && !any_split)
    {
        out.surface_triangles.clear();
        out.surface_tri_parent_id.clear();
    }

    return out;
}

/// @brief Extract a closed `InputSurface` from a 3D volume `(3, tag)` or 2D closed shell `(2, tag)`
///        in the active Gmsh model, preserving CAD surface tags in `triangle_facet_ids`.
inline InputSurface extract_input_surface(int dim,
                                          int tag,
                                          ObjectType type = ObjectType::SOLID,
                                          int priority    = 0)
{
    std::vector<std::pair<int, int>> signed_surfaces;
    if (dim == 3)
    {
        gmsh::model::getBoundary({ { 3, tag } }, signed_surfaces, false, true, false);
    }
    else if (dim == 2)
    {
        signed_surfaces.push_back({ 2, tag });
    }
    else
    {
        throw std::invalid_argument(
          "gmsh_sdk::extract_input_surface expects dim == 2 (closed shell) or dim == 3 (volume).");
    }

    if (signed_surfaces.empty())
    {
        throw std::runtime_error("gmsh_sdk::extract_input_surface found no bounding surfaces.");
    }

    InputSurface surf;
    surf.type     = type;
    surf.priority = priority;

    std::unordered_map<std::size_t, Index> node_tag_to_vid;
    for (std::size_t s = 0; s < signed_surfaces.size(); ++s)
    {
        const int ftag = std::abs(signed_surfaces[s].second);

        std::vector<std::size_t> ntags;
        std::vector<double> coords;
        std::vector<double> params;
        gmsh::model::mesh::getNodes(ntags, coords, params, 2, ftag, true, false);
        for (std::size_t i = 0; i < ntags.size(); ++i)
        {
            if (node_tag_to_vid.find(ntags[i]) == node_tag_to_vid.end())
            {
                const Index vid           = static_cast<Index>(surf.vertices.size());
                node_tag_to_vid[ntags[i]] = vid;
                surf.vertices.push_back(
                  { coords[3 * i + 0], coords[3 * i + 1], coords[3 * i + 2] });
            }
        }

        std::vector<std::size_t> etags;
        std::vector<std::size_t> enodes;
        gmsh::model::mesh::getElementsByType(2, etags, enodes, ftag);
        for (std::size_t i = 0; i < etags.size(); ++i)
        {
            const Index v0 = node_tag_to_vid.at(enodes[3 * i + 0]);
            const Index v1 = node_tag_to_vid.at(enodes[3 * i + 1]);
            const Index v2 = node_tag_to_vid.at(enodes[3 * i + 2]);
            surf.triangles.push_back({ v0, v1, v2 });
            surf.triangle_facet_ids.push_back(static_cast<Marker>(ftag));
        }
    }

    if (surf.triangles.empty())
    {
        throw std::runtime_error(
          "gmsh_sdk::extract_input_surface found 0 2D triangle elements on input entity (run "
          "gmsh::model::mesh::generate(2) first).");
    }

    detail::orient_closed_shell_outward(surf);

    return surf;
}

/// @brief Populate the active Gmsh model with the arranged 0D/1D/2D discrete topology and
///        3D `GEO` `Volume` entities so `gmsh::model::mesh::generate(3)` meshes all arranged
///        regions directly without requiring TetGen.
/// @return Vector of `(3, volumeTag)` pairs corresponding 1-to-1 with `result.brep_volumes`.
inline std::vector<std::pair<int, int>> populate_arranged_model(const ArrangeResult& result,
                                                                bool createGeometry = false)
{
    std::vector<std::pair<int, int>> volume_dim_tags;
    if (result.brep_volumes.empty())
    {
        return volume_dim_tags;
    }

    const int base_vtag             = detail::next_entity_tag(0);
    const int base_ctag             = detail::next_entity_tag(1);
    const int base_ftag             = detail::next_entity_tag(2);
    const int base_rtag             = detail::next_entity_tag(3);
    const std::size_t base_node_tag = detail::next_node_tag();

    const std::size_t n_verts = result.surface_mesh.vertices.size();
    std::vector<int> vertex_entity_dim(n_verts, -1);
    std::vector<int> vertex_entity_tag(n_verts, -1);

    std::vector<int> point_tags(result.brep_vertices.size(), -1);
    for (std::size_t i = 0; i < result.brep_vertices.size(); ++i)
    {
        const int vtag          = gmsh::model::addDiscreteEntity(0, base_vtag + static_cast<int>(i));
        point_tags[i]           = vtag;
        const Index vid         = result.brep_vertices[i].global_vertex_id;
        vertex_entity_dim[vid]  = 0;
        vertex_entity_tag[vid]  = vtag;
    }

    std::vector<int> curve_tags(result.brep_curves.size(), -1);
    for (std::size_t i = 0; i < result.brep_curves.size(); ++i)
    {
        const auto& c = result.brep_curves[i];
        std::vector<int> bnd;
        if (c.start_brep_vertex != NOT_SET_Index && c.end_brep_vertex != NOT_SET_Index)
        {
            bnd.push_back(point_tags.at(c.start_brep_vertex));
            bnd.push_back(point_tags.at(c.end_brep_vertex));
        }
        const int ctag = gmsh::model::addDiscreteEntity(1, base_ctag + static_cast<int>(i), bnd);
        curve_tags[i]  = ctag;
        for (std::size_t k = 0; k < c.polyline_vertices.size(); ++k)
        {
            const Index vid = c.polyline_vertices[k];
            if (vertex_entity_dim[vid] < 0)
            {
                vertex_entity_dim[vid] = 1;
                vertex_entity_tag[vid] = ctag;
            }
        }
    }

    std::vector<int> face_tags(result.brep_faces.size(), -1);
    for (std::size_t i = 0; i < result.brep_faces.size(); ++i)
    {
        const auto& f = result.brep_faces[i];
        std::vector<int> bnd;
        bnd.reserve(f.signed_bounding_curves.size());
        for (std::size_t k = 0; k < f.signed_bounding_curves.size(); ++k)
        {
            const int ctag = curve_tags.at(f.signed_bounding_curves[k].first);
            bnd.push_back(f.signed_bounding_curves[k].second >= 0 ? ctag : -ctag);
        }
        const int ftag = gmsh::model::addDiscreteEntity(2, base_ftag + static_cast<int>(i), bnd);
        face_tags[i]   = ftag;
        for (std::size_t t = 0; t < f.triangles.size(); ++t)
        {
            for (int c = 0; c < 3; ++c)
            {
                const Index vid = f.triangles[t][c];
                if (vertex_entity_dim[vid] < 0)
                {
                    vertex_entity_dim[vid] = 2;
                    vertex_entity_tag[vid] = ftag;
                }
            }
        }
    }

    std::unordered_map<long long, std::pair<std::vector<std::size_t>, std::vector<double>>>
      entity_nodes;
    for (std::size_t vid = 0; vid < n_verts; ++vid)
    {
        const int dim = vertex_entity_dim[vid];
        const int tag = vertex_entity_tag[vid];
        if (dim < 0 || tag <= 0)
        {
            continue;
        }
        const long long key =
          (static_cast<long long>(dim) << 32) | static_cast<unsigned int>(tag);
        auto& bucket = entity_nodes[key];
        bucket.first.push_back(base_node_tag + vid);
        bucket.second.push_back(result.surface_mesh.vertices[vid][0]);
        bucket.second.push_back(result.surface_mesh.vertices[vid][1]);
        bucket.second.push_back(result.surface_mesh.vertices[vid][2]);
    }
    for (auto it = entity_nodes.begin(); it != entity_nodes.end(); ++it)
    {
        const int dim = static_cast<int>(it->first >> 32);
        const int tag = static_cast<int>(it->first & 0xffffffffLL);
        gmsh::model::mesh::addNodes(dim, tag, it->second.first, it->second.second);
    }

    for (std::size_t i = 0; i < result.brep_vertices.size(); ++i)
    {
        const std::size_t ntag =
          base_node_tag + static_cast<std::size_t>(result.brep_vertices[i].global_vertex_id);
        gmsh::model::mesh::addElementsByType(point_tags[i], 15, {}, { ntag });
    }
    for (std::size_t i = 0; i < result.brep_curves.size(); ++i)
    {
        const auto& poly = result.brep_curves[i].polyline_vertices;
        std::vector<std::size_t> enodes;
        if (poly.size() >= 2)
        {
            enodes.reserve((poly.size() - 1) * 2);
            for (std::size_t k = 0; k + 1 < poly.size(); ++k)
            {
                enodes.push_back(base_node_tag + static_cast<std::size_t>(poly[k]));
                enodes.push_back(base_node_tag + static_cast<std::size_t>(poly[k + 1]));
            }
        }
        gmsh::model::mesh::addElementsByType(curve_tags[i], 1, {}, enodes);
    }
    for (std::size_t i = 0; i < result.brep_faces.size(); ++i)
    {
        const auto& f = result.brep_faces[i];
        std::vector<std::size_t> enodes;
        enodes.reserve(f.triangles.size() * 3);
        for (std::size_t t = 0; t < f.triangles.size(); ++t)
        {
            enodes.push_back(base_node_tag + static_cast<std::size_t>(f.triangles[t][0]));
            enodes.push_back(base_node_tag + static_cast<std::size_t>(f.triangles[t][1]));
            enodes.push_back(base_node_tag + static_cast<std::size_t>(f.triangles[t][2]));
        }
        gmsh::model::mesh::addElementsByType(face_tags[i], 2, {}, enodes);

        if (!f.embedded_curves.empty())
        {
            std::vector<int> emb_c;
            emb_c.reserve(f.embedded_curves.size());
            for (std::size_t k = 0; k < f.embedded_curves.size(); ++k)
            {
                emb_c.push_back(curve_tags.at(f.embedded_curves[k]));
            }
            gmsh::model::mesh::embed(1, emb_c, 2, face_tags[i]);
        }
        if (!f.embedded_vertices.empty())
        {
            std::vector<int> emb_v;
            emb_v.reserve(f.embedded_vertices.size());
            for (std::size_t k = 0; k < f.embedded_vertices.size(); ++k)
            {
                emb_v.push_back(point_tags.at(f.embedded_vertices[k]));
            }
            gmsh::model::mesh::embed(0, emb_v, 2, face_tags[i]);
        }
    }

    if (createGeometry)
    {
        std::vector<std::pair<int, int>> disc_surfs;
        disc_surfs.reserve(face_tags.size());
        for (std::size_t i = 0; i < face_tags.size(); ++i)
        {
            disc_surfs.push_back({ 2, face_tags[i] });
        }
        gmsh::model::mesh::createGeometry(disc_surfs);
    }

    volume_dim_tags.reserve(result.brep_volumes.size());
    for (std::size_t i = 0; i < result.brep_volumes.size(); ++i)
    {
        const auto& bvol = result.brep_volumes[i];
        std::vector<int> ftags;
        ftags.reserve(bvol.signed_bounding_faces.size());
        for (std::size_t k = 0; k < bvol.signed_bounding_faces.size(); ++k)
        {
            ftags.push_back(face_tags.at(bvol.signed_bounding_faces[k].first));
        }
        const int sl   = gmsh::model::geo::addSurfaceLoop(ftags);
        const int vtag = gmsh::model::geo::addVolume({ sl }, base_rtag + static_cast<int>(i));
        volume_dim_tags.push_back({ 3, vtag });
    }
    gmsh::model::geo::synchronize();

    return volume_dim_tags;
}

/// @brief Public Gmsh SDK counterpart of `gmsh::model::mesh::arrangeSurfaces(...)`.
///        Extracts discrete closed surfaces/volumes from `objectDimTags` and `toolDimTags`,
///        executes `mesh_arrange::arrange_surfaces` using Gmsh's native `plc_tetrahedralizer`
///        (zero TetGen dependency), populates the active Gmsh model with meshable 3D `Volume`
///        entities, and returns `outDimTags` + `outDimTagsMap` provenance matching
///        `gmsh::model::occ::fragment`.
inline ArrangeResult arrange_surfaces(
  const std::vector<std::pair<int, int>>& objectDimTags,
  const std::vector<std::pair<int, int>>& toolDimTags,
  std::vector<std::pair<int, int>>& outDimTags,
  std::vector<std::vector<std::pair<int, int>>>& outDimTagsMap,
  const std::vector<int>& objectPriorities = std::vector<int>(),
  const std::vector<int>& toolTypes        = std::vector<int>(),
  const std::vector<int>& toolPriorities   = std::vector<int>(),
  bool removeOverlaps                      = true,
  bool optimizeSurfaces                    = false,
  double maxEdgeLength                     = 0.0,
  bool createGeometry                      = false,
  bool removeObject                        = true,
  bool removeTool                          = true)
{
    if (!objectPriorities.empty() && objectPriorities.size() != objectDimTags.size())
    {
        throw std::invalid_argument(
          "gmsh_sdk::arrange_surfaces: objectPriorities must be empty or match objectDimTags size.");
    }
    if (!toolTypes.empty() && toolTypes.size() != toolDimTags.size())
    {
        throw std::invalid_argument(
          "gmsh_sdk::arrange_surfaces: toolTypes must be empty or match toolDimTags size.");
    }
    if (!toolPriorities.empty() && toolPriorities.size() != toolDimTags.size())
    {
        throw std::invalid_argument(
          "gmsh_sdk::arrange_surfaces: toolPriorities must be empty or match toolDimTags size.");
    }

    const std::size_t n_obj  = objectDimTags.size();
    const std::size_t n_tool = toolDimTags.size();
    std::vector<InputSurface> inputs;
    inputs.reserve(n_obj + n_tool);

    for (std::size_t i = 0; i < n_obj; ++i)
    {
        const int prio = objectPriorities.empty() ? 0 : objectPriorities[i];
        inputs.push_back(extract_input_surface(
          objectDimTags[i].first, objectDimTags[i].second, ObjectType::SOLID, prio));
    }
    for (std::size_t i = 0; i < n_tool; ++i)
    {
        const int tcode       = toolTypes.empty() ? 1 : toolTypes[i];
        const ObjectType type = (tcode >= 0) ? ObjectType::CLIP_INSIDE : ObjectType::CLIP_OUTSIDE;
        const int prio        = toolPriorities.empty() ? 0 : toolPriorities[i];
        inputs.push_back(
          extract_input_surface(toolDimTags[i].first, toolDimTags[i].second, type, prio));
    }

    ArrangeOptions opts;
    opts.remove_overlaps     = removeOverlaps;
    opts.optimize_surfaces   = optimizeSurfaces;
    opts.max_edge_length     = maxEdgeLength;
    opts.plc_tetrahedralizer = plc_tetrahedralizer;
    opts.log_callback        = [](LogLevel level, const std::string& msg)
    {
        const char* lvl = "info";
        if (level == LogLevel::WARNING)
        {
            lvl = "warning";
        }
        else if (level == LogLevel::ERROR)
        {
            lvl = "error";
        }
        else if (level == LogLevel::DEBUG)
        {
            lvl = "debug";
        }
        gmsh::logger::write(msg, lvl);
    };

    ArrangeResult res = ::mesh_arrange::arrange_surfaces(inputs, opts);

    std::vector<std::pair<int, int>> to_remove;
    if (removeObject)
    {
        to_remove.insert(to_remove.end(), objectDimTags.begin(), objectDimTags.end());
    }
    if (removeTool)
    {
        to_remove.insert(to_remove.end(), toolDimTags.begin(), toolDimTags.end());
    }
    if (!to_remove.empty())
    {
        try
        {
            gmsh::model::geo::remove(to_remove, true);
            gmsh::model::geo::synchronize();
        }
        catch (...)
        {
        }
        try
        {
            gmsh::model::occ::remove(to_remove, true);
            gmsh::model::occ::synchronize();
        }
        catch (...)
        {
        }
        try
        {
            gmsh::model::removeEntities(to_remove, true);
        }
        catch (...)
        {
        }
    }

    outDimTags = populate_arranged_model(res, createGeometry);
    outDimTagsMap.assign(n_obj + n_tool, std::vector<std::pair<int, int>>());
    for (std::size_t i = 0; i < res.brep_volumes.size(); ++i)
    {
        const auto& owners = res.brep_volumes[i].owner_objects;
        for (std::size_t k = 0; k < owners.size(); ++k)
        {
            const std::size_t oid = static_cast<std::size_t>(owners[k]);
            if (oid < outDimTagsMap.size())
            {
                outDimTagsMap[oid].push_back(outDimTags[i]);
            }
        }
    }

    return res;
}

} // namespace gmsh_sdk
} // namespace mesh_arrange
