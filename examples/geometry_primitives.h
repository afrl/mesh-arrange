#pragma once

#include "mesh_arrange/MeshArrange.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mesh_arrange::examples
{

constexpr double kPi = 3.14159265358979323846264338327950288;

inline double signed_volume(const InputSurface& s)
{
    double vol = 0.0;
    for (const auto& t : s.triangles)
    {
        const auto& a = s.vertices[t[0]];
        const auto& b = s.vertices[t[1]];
        const auto& c = s.vertices[t[2]];
        vol += (a[0] * (b[1] * c[2] - b[2] * c[1]) +
                a[1] * (b[2] * c[0] - b[0] * c[2]) +
                a[2] * (b[0] * c[1] - b[1] * c[0])) /
               6.0;
    }
    return vol;
}

inline void ensure_outward_orientation(InputSurface& s)
{
    if (signed_volume(s) < 0.0)
    {
        for (auto& t : s.triangles)
        {
            std::swap(t[1], t[2]);
        }
    }
}

/// @brief Watertight subdivided axis-aligned box with 6 planar CAD facets (FacetID 0..5).
inline InputSurface make_subdivided_box(const std::array<double, 3>& bmin,
                                        const std::array<double, 3>& bmax,
                                        int nx                = 6,
                                        int ny                = 6,
                                        int nz                = 6,
                                        ObjectType type       = ObjectType::SOLID,
                                        PriorityType priority = 0)
{
    InputSurface s;
    s.type     = type;
    s.priority = priority;

    std::map<std::array<int, 3>, Index> grid_vid;
    auto get_vid = [&](int ix, int iy, int iz) -> Index
    {
        const std::array<int, 3> key = { ix, iy, iz };
        auto it                      = grid_vid.find(key);
        if (it != grid_vid.end())
        {
            return it->second;
        }
        const double fx = static_cast<double>(ix) / static_cast<double>(nx);
        const double fy = static_cast<double>(iy) / static_cast<double>(ny);
        const double fz = static_cast<double>(iz) / static_cast<double>(nz);
        const Index id  = static_cast<Index>(s.vertices.size());
        s.vertices.push_back({ bmin[0] + fx * (bmax[0] - bmin[0]),
                               bmin[1] + fy * (bmax[1] - bmin[1]),
                               bmin[2] + fz * (bmax[2] - bmin[2]) });
        grid_vid.emplace(key, id);
        return id;
    };

    auto add_quad = [&](Index v00, Index v10, Index v11, Index v01, FacetID fid, bool flip)
    {
        if (!flip)
        {
            s.triangles.push_back({ v00, v10, v11 });
            s.triangles.push_back({ v00, v11, v01 });
        }
        else
        {
            s.triangles.push_back({ v00, v11, v10 });
            s.triangles.push_back({ v00, v01, v11 });
        }
        s.triangle_facet_ids.push_back(fid);
        s.triangle_facet_ids.push_back(fid);
    };

    // -Z (fid 0) and +Z (fid 1)
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            add_quad(get_vid(i, j, 0), get_vid(i + 1, j, 0), get_vid(i + 1, j + 1, 0), get_vid(i, j + 1, 0), 0, true);
            add_quad(get_vid(i, j, nz), get_vid(i + 1, j, nz), get_vid(i + 1, j + 1, nz), get_vid(i, j + 1, nz), 1, false);
        }
    }
    // -Y (fid 2) and +Y (fid 3)
    for (int i = 0; i < nx; ++i)
    {
        for (int k = 0; k < nz; ++k)
        {
            add_quad(get_vid(i, 0, k), get_vid(i + 1, 0, k), get_vid(i + 1, 0, k + 1), get_vid(i, 0, k + 1), 2, false);
            add_quad(get_vid(i, ny, k), get_vid(i + 1, ny, k), get_vid(i + 1, ny, k + 1), get_vid(i, ny, k + 1), 3, true);
        }
    }
    // -X (fid 4) and +X (fid 5)
    for (int j = 0; j < ny; ++j)
    {
        for (int k = 0; k < nz; ++k)
        {
            add_quad(get_vid(0, j, k), get_vid(0, j + 1, k), get_vid(0, j + 1, k + 1), get_vid(0, j, k + 1), 4, true);
            add_quad(get_vid(nx, j, k), get_vid(nx, j + 1, k), get_vid(nx, j + 1, k + 1), get_vid(nx, j, k + 1), 5, false);
        }
    }

    ensure_outward_orientation(s);
    return s;
}

/// @brief Smooth watertight icosphere (subdivision level 2 -> 320 triangles, level 3 -> 1280 triangles).
inline InputSurface make_icosphere(const std::array<double, 3>& center,
                                   double radius,
                                   int subdivisions      = 2,
                                   ObjectType type       = ObjectType::SOLID,
                                   PriorityType priority = 0)
{
    const double t = (1.0 + std::sqrt(5.0)) * 0.5;
    std::vector<std::array<double, 3>> verts = {
        { -1,  t,  0 }, {  1,  t,  0 }, { -1, -t,  0 }, {  1, -t,  0 },
        {  0, -1,  t }, {  0,  1,  t }, {  0, -1, -t }, {  0,  1, -t },
        {  t,  0, -1 }, {  t,  0,  1 }, { -t,  0, -1 }, { -t,  0,  1 }
    };
    auto normalize = [](std::array<double, 3>& v)
    {
        const double len = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        v[0] /= len;
        v[1] /= len;
        v[2] /= len;
    };
    for (auto& v : verts)
    {
        normalize(v);
    }

    std::vector<std::array<Index, 3>> tris = {
        { 0, 11, 5 }, { 0, 5, 1 },  { 0, 1, 7 },   { 0, 7, 10 }, { 0, 10, 11 },
        { 1, 5, 9 },  { 5, 11, 4 }, { 11, 10, 2 }, { 10, 7, 6 }, { 7, 1, 8 },
        { 3, 9, 4 },  { 3, 4, 2 },  { 3, 2, 6 },   { 3, 6, 8 },  { 3, 8, 9 },
        { 4, 9, 5 },  { 2, 4, 11 }, { 6, 2, 10 },  { 8, 6, 7 },  { 9, 8, 1 }
    };

    for (int sub = 0; sub < subdivisions; ++sub)
    {
        std::map<std::pair<Index, Index>, Index> midpoint_cache;
        auto get_mid = [&](Index a, Index b) -> Index
        {
            const auto key = std::minmax(a, b);
            auto it        = midpoint_cache.find(key);
            if (it != midpoint_cache.end())
            {
                return it->second;
            }
            std::array<double, 3> m = { 0.5 * (verts[a][0] + verts[b][0]),
                                        0.5 * (verts[a][1] + verts[b][1]),
                                        0.5 * (verts[a][2] + verts[b][2]) };
            normalize(m);
            const Index idx = static_cast<Index>(verts.size());
            verts.push_back(m);
            midpoint_cache.emplace(key, idx);
            return idx;
        };

        std::vector<std::array<Index, 3>> next_tris;
        next_tris.reserve(tris.size() * 4);
        for (const auto& f : tris)
        {
            const Index a  = f[0];
            const Index b  = f[1];
            const Index c  = f[2];
            const Index ab = get_mid(a, b);
            const Index bc = get_mid(b, c);
            const Index ca = get_mid(c, a);
            next_tris.push_back({ a, ab, ca });
            next_tris.push_back({ b, bc, ab });
            next_tris.push_back({ c, ca, bc });
            next_tris.push_back({ ab, bc, ca });
        }
        tris.swap(next_tris);
    }

    InputSurface s;
    s.type     = type;
    s.priority = priority;
    s.vertices.reserve(verts.size());
    for (const auto& v : verts)
    {
        s.vertices.push_back({ center[0] + radius * v[0],
                               center[1] + radius * v[1],
                               center[2] + radius * v[2] });
    }
    s.triangles = std::move(tris);
    ensure_outward_orientation(s);
    return s;
}

/// @brief Smooth closed torus rotated around X and Y axes.
inline InputSurface make_tilted_torus(const std::array<double, 3>& center,
                                      double major_r,
                                      double minor_r,
                                      double rx_rad         = 0.0,
                                      double ry_rad         = 0.0,
                                      int n_major           = 40,
                                      int n_minor           = 20,
                                      ObjectType type       = ObjectType::SOLID,
                                      PriorityType priority = 0)
{
    InputSurface s;
    s.type     = type;
    s.priority = priority;

    const double cx = std::cos(rx_rad), sx = std::sin(rx_rad);
    const double cy = std::cos(ry_rad), sy = std::sin(ry_rad);

    auto rotate_and_shift = [&](double x, double y, double z) -> std::array<double, 3>
    {
        // Rotate around X
        const double y1 = cx * y - sx * z;
        const double z1 = sx * y + cx * z;
        // Rotate around Y
        const double x2 = cy * x + sy * z1;
        const double z2 = -sy * x + cy * z1;
        return { center[0] + x2, center[1] + y1, center[2] + z2 };
    };

    for (int i = 0; i < n_major; ++i)
    {
        const double u  = (2.0 * kPi * i) / static_cast<double>(n_major);
        const double cu = std::cos(u), su = std::sin(u);
        for (int j = 0; j < n_minor; ++j)
        {
            const double v  = (2.0 * kPi * j) / static_cast<double>(n_minor);
            const double cv = std::cos(v), sv = std::sin(v);
            const double r  = major_r + minor_r * cv;
            s.vertices.push_back(rotate_and_shift(r * cu, r * su, minor_r * sv));
        }
    }

    for (int i = 0; i < n_major; ++i)
    {
        const int ni = (i + 1) % n_major;
        for (int j = 0; j < n_minor; ++j)
        {
            const int nj     = (j + 1) % n_minor;
            const Index v00  = static_cast<Index>(i * n_minor + j);
            const Index v10  = static_cast<Index>(ni * n_minor + j);
            const Index v11  = static_cast<Index>(ni * n_minor + nj);
            const Index v01  = static_cast<Index>(i * n_minor + nj);
            s.triangles.push_back({ v00, v10, v11 });
            s.triangles.push_back({ v00, v11, v01 });
        }
    }

    ensure_outward_orientation(s);
    return s;
}

/// @brief Watertight cylinder between p0 and p1 with concentric cap rings (FacetID 0=barrel, 1=p0 cap, 2=p1 cap).
inline InputSurface make_cylinder(const std::array<double, 3>& p0,
                                  const std::array<double, 3>& p1,
                                  double radius,
                                  int n_radial          = 28,
                                  int n_axial           = 12,
                                  int n_cap_rings       = 2,
                                  ObjectType type       = ObjectType::SOLID,
                                  PriorityType priority = 0)
{
    std::array<double, 3> axis = { p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2] };
    const double len           = std::sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
    axis[0] /= len;
    axis[1] /= len;
    axis[2] /= len;

    // Orthonormal basis (u, v, axis)
    std::array<double, 3> ref = (std::abs(axis[2]) < 0.85) ? std::array<double, 3>{ 0, 0, 1 }
                                                           : std::array<double, 3>{ 1, 0, 0 };
    std::array<double, 3> u   = { axis[1] * ref[2] - axis[2] * ref[1],
                                  axis[2] * ref[0] - axis[0] * ref[2],
                                  axis[0] * ref[1] - axis[1] * ref[0] };
    const double ulen         = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
    u[0] /= ulen;
    u[1] /= ulen;
    u[2] /= ulen;
    std::array<double, 3> v = { axis[1] * u[2] - axis[2] * u[1],
                                axis[2] * u[0] - axis[0] * u[2],
                                axis[0] * u[1] - axis[1] * u[0] };

    InputSurface s;
    s.type     = type;
    s.priority = priority;

    auto eval_pt = [&](double s_frac, double r_frac, double theta) -> std::array<double, 3>
    {
        const double ca = std::cos(theta), sa = std::sin(theta);
        const double rr = radius * r_frac;
        return { p0[0] + s_frac * len * axis[0] + rr * (ca * u[0] + sa * v[0]),
                 p0[1] + s_frac * len * axis[1] + rr * (ca * u[1] + sa * v[1]),
                 p0[2] + s_frac * len * axis[2] + rr * (ca * u[2] + sa * v[2]) };
    };

    // Lateral rings: s = 0 .. n_axial
    for (int ia = 0; ia <= n_axial; ++ia)
    {
        const double sf = static_cast<double>(ia) / static_cast<double>(n_axial);
        for (int ir = 0; ir < n_radial; ++ir)
        {
            const double th = (2.0 * kPi * ir) / static_cast<double>(n_radial);
            s.vertices.push_back(eval_pt(sf, 1.0, th));
        }
    }

    for (int ia = 0; ia < n_axial; ++ia)
    {
        const int b0 = ia * n_radial;
        const int b1 = (ia + 1) * n_radial;
        for (int ir = 0; ir < n_radial; ++ir)
        {
            const int nr = (ir + 1) % n_radial;
            s.triangles.push_back({ static_cast<Index>(b0 + ir),
                                    static_cast<Index>(b0 + nr),
                                    static_cast<Index>(b1 + nr) });
            s.triangle_facet_ids.push_back(0);
            s.triangles.push_back({ static_cast<Index>(b0 + ir),
                                    static_cast<Index>(b1 + nr),
                                    static_cast<Index>(b1 + ir) });
            s.triangle_facet_ids.push_back(0);
        }
    }

    // Helper to build concentric cap
    auto build_cap = [&](double sf, int outer_ring_base, FacetID fid, bool flip)
    {
        int prev_ring_base = outer_ring_base;
        for (int ring = n_cap_rings - 1; ring >= 1; --ring)
        {
            const double rf    = static_cast<double>(ring) / static_cast<double>(n_cap_rings);
            const int cur_base = static_cast<int>(s.vertices.size());
            for (int ir = 0; ir < n_radial; ++ir)
            {
                const double th = (2.0 * kPi * ir) / static_cast<double>(n_radial);
                s.vertices.push_back(eval_pt(sf, rf, th));
            }
            for (int ir = 0; ir < n_radial; ++ir)
            {
                const int nr    = (ir + 1) % n_radial;
                const Index o0  = static_cast<Index>(prev_ring_base + ir);
                const Index o1  = static_cast<Index>(prev_ring_base + nr);
                const Index i0  = static_cast<Index>(cur_base + ir);
                const Index i1  = static_cast<Index>(cur_base + nr);
                if (!flip)
                {
                    s.triangles.push_back({ o0, o1, i1 });
                    s.triangles.push_back({ o0, i1, i0 });
                }
                else
                {
                    s.triangles.push_back({ o0, i1, o1 });
                    s.triangles.push_back({ o0, i0, i1 });
                }
                s.triangle_facet_ids.push_back(fid);
                s.triangle_facet_ids.push_back(fid);
            }
            prev_ring_base = cur_base;
        }
        const Index center_vid = static_cast<Index>(s.vertices.size());
        s.vertices.push_back(eval_pt(sf, 0.0, 0.0));
        for (int ir = 0; ir < n_radial; ++ir)
        {
            const int nr   = (ir + 1) % n_radial;
            const Index i0 = static_cast<Index>(prev_ring_base + ir);
            const Index i1 = static_cast<Index>(prev_ring_base + nr);
            if (!flip)
            {
                s.triangles.push_back({ i0, i1, center_vid });
            }
            else
            {
                s.triangles.push_back({ i0, center_vid, i1 });
            }
            s.triangle_facet_ids.push_back(fid);
        }
    };

    build_cap(0.0, 0, 1, true);
    build_cap(1.0, n_axial * n_radial, 2, false);

    ensure_outward_orientation(s);
    return s;
}

/// @brief Smooth sinusoidal elliptical composite tow swept along X (`along_x=true`) or Y (`along_x=false`).
inline InputSurface make_sinusoidal_tow(bool along_x,
                                        double half_length,
                                        double half_width,
                                        double half_thickness,
                                        double z_center,
                                        double wave_amp,
                                        int n_axial           = 26,
                                        int n_profile         = 24,
                                        ObjectType type       = ObjectType::SOLID,
                                        PriorityType priority = 1)
{
    InputSurface s;
    s.type     = type;
    s.priority = priority;

    auto eval_pt = [&](double t_frac, double r_frac, double phi) -> std::array<double, 3>
    {
        const double s_coord = -half_length + 2.0 * half_length * t_frac;
        const double zc      = z_center + wave_amp * std::cos(kPi * s_coord / (2.0 * half_length));
        const double w       = r_frac * half_width * std::cos(phi);
        const double h       = r_frac * half_thickness * std::sin(phi);
        if (along_x)
        {
            return { s_coord, w, zc + h };
        }
        return { w, s_coord, zc + h };
    };

    for (int ia = 0; ia <= n_axial; ++ia)
    {
        const double tf = static_cast<double>(ia) / static_cast<double>(n_axial);
        for (int ip = 0; ip < n_profile; ++ip)
        {
            const double phi = (2.0 * kPi * ip) / static_cast<double>(n_profile);
            s.vertices.push_back(eval_pt(tf, 1.0, phi));
        }
    }

    for (int ia = 0; ia < n_axial; ++ia)
    {
        const int b0 = ia * n_profile;
        const int b1 = (ia + 1) * n_profile;
        for (int ip = 0; ip < n_profile; ++ip)
        {
            const int np = (ip + 1) % n_profile;
            s.triangles.push_back({ static_cast<Index>(b0 + ip),
                                    static_cast<Index>(b1 + np),
                                    static_cast<Index>(b0 + np) });
            s.triangle_facet_ids.push_back(0);
            s.triangles.push_back({ static_cast<Index>(b0 + ip),
                                    static_cast<Index>(b1 + ip),
                                    static_cast<Index>(b1 + np) });
            s.triangle_facet_ids.push_back(0);
        }
    }

    auto add_elliptical_cap = [&](double tf, int outer_base, FacetID fid, bool flip)
    {
        const int mid_base = static_cast<int>(s.vertices.size());
        for (int ip = 0; ip < n_profile; ++ip)
        {
            const double phi = (2.0 * kPi * ip) / static_cast<double>(n_profile);
            s.vertices.push_back(eval_pt(tf, 0.5, phi));
        }
        for (int ip = 0; ip < n_profile; ++ip)
        {
            const int np   = (ip + 1) % n_profile;
            const Index o0 = static_cast<Index>(outer_base + ip);
            const Index o1 = static_cast<Index>(outer_base + np);
            const Index m0 = static_cast<Index>(mid_base + ip);
            const Index m1 = static_cast<Index>(mid_base + np);
            if (!flip)
            {
                s.triangles.push_back({ o0, o1, m1 });
                s.triangles.push_back({ o0, m1, m0 });
            }
            else
            {
                s.triangles.push_back({ o0, m1, o1 });
                s.triangles.push_back({ o0, m0, m1 });
            }
            s.triangle_facet_ids.push_back(fid);
            s.triangle_facet_ids.push_back(fid);
        }
        const Index c_vid = static_cast<Index>(s.vertices.size());
        s.vertices.push_back(eval_pt(tf, 0.0, 0.0));
        for (int ip = 0; ip < n_profile; ++ip)
        {
            const int np   = (ip + 1) % n_profile;
            const Index m0 = static_cast<Index>(mid_base + ip);
            const Index m1 = static_cast<Index>(mid_base + np);
            if (!flip)
            {
                s.triangles.push_back({ m0, m1, c_vid });
            }
            else
            {
                s.triangles.push_back({ m0, c_vid, m1 });
            }
            s.triangle_facet_ids.push_back(fid);
        }
    };

    add_elliptical_cap(0.0, 0, 1, false);
    add_elliptical_cap(1.0, n_axial * n_profile, 2, true);

    ensure_outward_orientation(s);
    return s;
}

/// @brief Feature-preserving tangential smoothing of `ArrangeResult::surface_mesh` without global subdivision.
inline void smooth_arrangement_surfaces(ArrangeResult& res, int n_iterations = 12)
{
    const auto seams = segment_seams(res.surface_mesh);
    optimize_tangential_smoothing(res.surface_mesh, seams, n_iterations, true);
    for (auto& bv : res.brep_vertices)
    {
        if (bv.global_vertex_id < res.surface_mesh.vertices.size())
        {
            bv.point = res.surface_mesh.vertices[bv.global_vertex_id];
        }
    }
}

inline void write_triangulation_json(std::ostream& os, const Triangulation& tri)
{
    os << "{\"vertices\":[";
    for (std::size_t i = 0; i < tri.vertices.size(); ++i)
    {
        if (i > 0) os << ",";
        os << "[" << tri.vertices[i][0] << "," << tri.vertices[i][1] << "," << tri.vertices[i][2] << "]";
    }
    os << "],\"triangles\":[";
    for (std::size_t i = 0; i < tri.triangles.size(); ++i)
    {
        if (i > 0) os << ",";
        os << "[" << tri.triangles[i][0] << "," << tri.triangles[i][1] << "," << tri.triangles[i][2] << "]";
    }
    os << "],\"markers\":[";
    for (std::size_t i = 0; i < tri.tri_to_marker_map.size(); ++i)
    {
        if (i > 0) os << ",";
        os << tri.tri_to_marker_map[i];
    }
    os << "]}";
}

inline void write_arrange_result_json(std::ostream& os, const ArrangeResult& res)
{
    os << "{\n  \"surface_mesh\": ";
    write_triangulation_json(os, res.surface_mesh);

    os << ",\n  \"brep_vertices\": [";
    for (std::size_t i = 0; i < res.brep_vertices.size(); ++i)
    {
        if (i > 0) os << ",";
        const auto& bv = res.brep_vertices[i];
        os << "{\"global_vid\":" << bv.global_vertex_id
           << ",\"point\":[" << bv.point[0] << "," << bv.point[1] << "," << bv.point[2] << "]}";
    }
    os << "]";

    os << ",\n  \"brep_curves\": [";
    for (std::size_t i = 0; i < res.brep_curves.size(); ++i)
    {
        if (i > 0) os << ",";
        const auto& bc = res.brep_curves[i];
        os << "{\"start_bv\":" << static_cast<int64_t>(bc.start_brep_vertex == NOT_SET_Index ? -1 : static_cast<int64_t>(bc.start_brep_vertex))
           << ",\"end_bv\":" << static_cast<int64_t>(bc.end_brep_vertex == NOT_SET_Index ? -1 : static_cast<int64_t>(bc.end_brep_vertex))
           << ",\"is_closed_loop\":" << (bc.is_closed_loop ? "true" : "false")
           << ",\"polyline\":[";
        for (std::size_t k = 0; k < bc.polyline_vertices.size(); ++k)
        {
            if (k > 0) os << ",";
            os << bc.polyline_vertices[k];
        }
        os << "],\"incident_faces\":[";
        for (std::size_t k = 0; k < bc.incident_faces.size(); ++k)
        {
            if (k > 0) os << ",";
            os << bc.incident_faces[k];
        }
        os << "]}";
    }
    os << "]";

    os << ",\n  \"brep_faces\": [";
    for (std::size_t i = 0; i < res.brep_faces.size(); ++i)
    {
        if (i > 0) os << ",";
        const auto& bf = res.brep_faces[i];
        os << "{\"boundary_id\":" << bf.boundary_id
           << ",\"bounded_regions\":["
           << (bf.bounded_regions.first == NOT_SET_RegionID ? -1 : static_cast<int>(bf.bounded_regions.first)) << ","
           << (bf.bounded_regions.second == NOT_SET_RegionID ? -1 : static_cast<int>(bf.bounded_regions.second)) << "]"
           << ",\"origin_objects\":[";
        for (std::size_t k = 0; k < bf.origin_objects.size(); ++k)
        {
            if (k > 0) os << ",";
            os << bf.origin_objects[k];
        }
        os << "],\"signed_bounding_curves\":[";
        for (std::size_t k = 0; k < bf.signed_bounding_curves.size(); ++k)
        {
            if (k > 0) os << ",";
            os << "[" << bf.signed_bounding_curves[k].first << "," << bf.signed_bounding_curves[k].second << "]";
        }
        os << "],\"embedded_curves\":[";
        for (std::size_t k = 0; k < bf.embedded_curves.size(); ++k)
        {
            if (k > 0) os << ",";
            os << bf.embedded_curves[k];
        }
        os << "],\"embedded_vertices\":[";
        for (std::size_t k = 0; k < bf.embedded_vertices.size(); ++k)
        {
            if (k > 0) os << ",";
            os << bf.embedded_vertices[k];
        }
        os << "],\"triangles\":[";
        for (std::size_t k = 0; k < bf.triangles.size(); ++k)
        {
            if (k > 0) os << ",";
            os << "[" << bf.triangles[k][0] << "," << bf.triangles[k][1] << "," << bf.triangles[k][2] << "]";
        }
        os << "]}";
    }
    os << "]";

    os << ",\n  \"brep_volumes\": [";
    for (std::size_t i = 0; i < res.brep_volumes.size(); ++i)
    {
        if (i > 0) os << ",";
        const auto& bv = res.brep_volumes[i];
        os << "{\"region_id\":" << bv.region_id
           << ",\"is_overlap\":" << (bv.is_overlap ? "true" : "false")
           << ",\"owner_objects\":[";
        for (std::size_t k = 0; k < bv.owner_objects.size(); ++k)
        {
            if (k > 0) os << ",";
            os << bv.owner_objects[k];
        }
        os << "],\"signed_bounding_faces\":[";
        for (std::size_t k = 0; k < bv.signed_bounding_faces.size(); ++k)
        {
            if (k > 0) os << ",";
            os << "[" << bv.signed_bounding_faces[k].first << "," << bv.signed_bounding_faces[k].second << "]";
        }
        os << "]}";
    }
    os << "]";

    os << ",\n  \"constituent_object_ids\": [";
    for (std::size_t i = 0; i < res.constituent_object_ids.size(); ++i)
    {
        if (i > 0) os << ",";
        os << res.constituent_object_ids[i];
    }
    os << "],\n  \"constituent_surfaces\": [";
    for (std::size_t i = 0; i < res.constituent_surfaces.size(); ++i)
    {
        if (i > 0) os << ",";
        write_triangulation_json(os, res.constituent_surfaces[i]);
    }
    os << "]\n}";
}

inline void export_demo_bundle(const std::string& json_path,
                               const std::vector<InputSurface>& inputs,
                               const std::vector<std::string>& input_names,
                               const ArrangeResult& primary_result,
                               const ArrangeResult* before_result = nullptr)
{
    std::ofstream ofs(json_path);
    if (!ofs)
    {
        throw std::runtime_error("Failed to open output JSON path: " + json_path);
    }
    ofs << std::setprecision(12);
    ofs << "{\n  \"inputs\": [";
    for (std::size_t i = 0; i < inputs.size(); ++i)
    {
        if (i > 0) ofs << ",";
        const auto& s          = inputs[i];
        const std::string name = (i < input_names.size()) ? input_names[i] : ("object_" + std::to_string(i));
        const char* tstr       = (s.type == ObjectType::SOLID)
                                   ? "SOLID"
                                   : (s.type == ObjectType::CLIP_INSIDE ? "CLIP_INSIDE" : "CLIP_OUTSIDE");
        ofs << "\n    {\"name\":\"" << name << "\",\"type\":\"" << tstr << "\",\"priority\":" << s.priority
            << ",\"vertices\":[";
        for (std::size_t k = 0; k < s.vertices.size(); ++k)
        {
            if (k > 0) ofs << ",";
            ofs << "[" << s.vertices[k][0] << "," << s.vertices[k][1] << "," << s.vertices[k][2] << "]";
        }
        ofs << "],\"triangles\":[";
        for (std::size_t k = 0; k < s.triangles.size(); ++k)
        {
            if (k > 0) ofs << ",";
            ofs << "[" << s.triangles[k][0] << "," << s.triangles[k][1] << "," << s.triangles[k][2] << "]";
        }
        ofs << "]}";
    }
    ofs << "\n  ],\n  \"result\": ";
    write_arrange_result_json(ofs, primary_result);
    if (before_result != nullptr)
    {
        ofs << ",\n  \"before_result\": ";
        write_arrange_result_json(ofs, *before_result);
    }
    ofs << "\n}\n";
}

namespace detail
{

inline std::string base64_encode_vtk_block(const void* payload, std::size_t byte_len)
{
    static constexpr char kTable[] =
      "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

    const uint64_t hdr = static_cast<uint64_t>(byte_len);
    std::vector<uint8_t> raw(sizeof(uint64_t) + byte_len);
    std::memcpy(raw.data(), &hdr, sizeof(uint64_t));
    if (byte_len > 0 && payload != nullptr)
    {
        std::memcpy(raw.data() + sizeof(uint64_t), payload, byte_len);
    }

    const std::size_t n = raw.size();
    std::string out;
    out.reserve(((n + 2) / 3) * 4);

    std::size_t i = 0;
    for (; i + 2 < n; i += 3)
    {
        const uint32_t triple = (static_cast<uint32_t>(raw[i]) << 16) |
                                (static_cast<uint32_t>(raw[i + 1]) << 8) |
                                static_cast<uint32_t>(raw[i + 2]);
        out.push_back(kTable[(triple >> 18) & 0x3F]);
        out.push_back(kTable[(triple >> 12) & 0x3F]);
        out.push_back(kTable[(triple >> 6) & 0x3F]);
        out.push_back(kTable[triple & 0x3F]);
    }
    if (i < n)
    {
        uint32_t triple = static_cast<uint32_t>(raw[i]) << 16;
        if (i + 1 < n)
        {
            triple |= static_cast<uint32_t>(raw[i + 1]) << 8;
        }
        out.push_back(kTable[(triple >> 18) & 0x3F]);
        out.push_back(kTable[(triple >> 12) & 0x3F]);
        out.push_back((i + 1 < n) ? kTable[(triple >> 6) & 0x3F] : '=');
        out.push_back('=');
    }
    return out;
}

struct VtuInt32CellArray
{
    std::string name;
    std::vector<int32_t> values;
};

inline void write_vtu_unstructured_grid(const std::string& path,
                                        const std::vector<std::array<double, 3>>& points,
                                        const std::vector<int64_t>& connectivity,
                                        const std::vector<int64_t>& offsets,
                                        const std::vector<uint8_t>& cell_types,
                                        const std::vector<VtuInt32CellArray>& cell_arrays)
{
    const std::filesystem::path fpath(path);
    if (fpath.has_parent_path())
    {
        std::filesystem::create_directories(fpath.parent_path());
    }

    std::ofstream ofs(path, std::ios::binary);
    if (!ofs)
    {
        throw std::runtime_error("Failed to open output VTU path: " + path);
    }

    std::vector<double> flat_pts(points.size() * 3);
    for (std::size_t i = 0; i < points.size(); ++i)
    {
        flat_pts[3 * i + 0] = points[i][0];
        flat_pts[3 * i + 1] = points[i][1];
        flat_pts[3 * i + 2] = points[i][2];
    }

    ofs << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n"
        << "  <UnstructuredGrid>\n"
        << "    <Piece NumberOfPoints=\"" << points.size() << "\" NumberOfCells=\"" << cell_types.size() << "\">\n"
        << "      <Points>\n"
        << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"binary\">\n"
        << "          " << base64_encode_vtk_block(flat_pts.data(), flat_pts.size() * sizeof(double)) << "\n"
        << "        </DataArray>\n"
        << "      </Points>\n"
        << "      <Cells>\n"
        << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"binary\">\n"
        << "          " << base64_encode_vtk_block(connectivity.data(), connectivity.size() * sizeof(int64_t)) << "\n"
        << "        </DataArray>\n"
        << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"binary\">\n"
        << "          " << base64_encode_vtk_block(offsets.data(), offsets.size() * sizeof(int64_t)) << "\n"
        << "        </DataArray>\n"
        << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"binary\">\n"
        << "          " << base64_encode_vtk_block(cell_types.data(), cell_types.size() * sizeof(uint8_t)) << "\n"
        << "        </DataArray>\n"
        << "      </Cells>\n"
        << "      <CellData>\n";

    for (const auto& arr : cell_arrays)
    {
        ofs << "        <DataArray type=\"Int32\" Name=\"" << arr.name << "\" format=\"binary\">\n"
            << "          " << base64_encode_vtk_block(arr.values.data(), arr.values.size() * sizeof(int32_t)) << "\n"
            << "        </DataArray>\n";
    }

    ofs << "      </CellData>\n"
        << "    </Piece>\n"
        << "  </UnstructuredGrid>\n"
        << "</VTKFile>\n";
}

} // namespace detail

/// @brief Writes all input `InputSurface` triangle meshes into a single `.vtu` file.
inline void write_inputs_vtu(const std::string& path, const std::vector<InputSurface>& inputs)
{
    std::vector<std::array<double, 3>> points;
    std::vector<int64_t> connectivity;
    std::vector<int64_t> offsets;
    std::vector<uint8_t> types;

    std::vector<int32_t> obj_id;
    std::vector<int32_t> obj_type;
    std::vector<int32_t> priority;
    std::vector<int32_t> facet_id;

    for (std::size_t i = 0; i < inputs.size(); ++i)
    {
        const auto& s          = inputs[i];
        const int64_t v_offset = static_cast<int64_t>(points.size());
        points.insert(points.end(), s.vertices.begin(), s.vertices.end());

        const int32_t tcode = (s.type == ObjectType::SOLID)
                                ? 0
                                : (s.type == ObjectType::CLIP_INSIDE ? 1 : 2);

        for (std::size_t t = 0; t < s.triangles.size(); ++t)
        {
            connectivity.push_back(v_offset + static_cast<int64_t>(s.triangles[t][0]));
            connectivity.push_back(v_offset + static_cast<int64_t>(s.triangles[t][1]));
            connectivity.push_back(v_offset + static_cast<int64_t>(s.triangles[t][2]));
            offsets.push_back(static_cast<int64_t>(connectivity.size()));
            types.push_back(5); // VTK_TRIANGLE

            obj_id.push_back(static_cast<int32_t>(i));
            obj_type.push_back(tcode);
            priority.push_back(static_cast<int32_t>(s.priority));
            facet_id.push_back(t < s.triangle_facet_ids.size()
                                 ? static_cast<int32_t>(s.triangle_facet_ids[t])
                                 : 0);
        }
    }

    detail::write_vtu_unstructured_grid(
      path,
      points,
      connectivity,
      offsets,
      types,
      {
        { "InputObjectID", std::move(obj_id) },
        { "ObjectType", std::move(obj_type) },
        { "Priority", std::move(priority) },
        { "InputFacetID", std::move(facet_id) },
      });
}

/// @brief Writes the complete non-manifold B-Rep (`BRepFace` triangles, `BRepCurve` polylines,
///        and `BRepVertex` corners) into a single `.vtu` file sharing `surface_mesh.vertices`.
inline void write_brep_vtu(const std::string& path, const ArrangeResult& res)
{
    const auto& points = res.surface_mesh.vertices;
    std::vector<int64_t> connectivity;
    std::vector<int64_t> offsets;
    std::vector<uint8_t> types;

    std::vector<int32_t> entity_dim;
    std::vector<int32_t> entity_tag;
    std::vector<int32_t> region_a;
    std::vector<int32_t> region_b;
    std::vector<int32_t> owner_obj;
    std::vector<int32_t> is_medial;
    std::vector<int32_t> valence;

    // 1. 2D BRepFace triangles (VTK_TRIANGLE = 5)
    for (std::size_t fi = 0; fi < res.brep_faces.size(); ++fi)
    {
        const auto& bf     = res.brep_faces[fi];
        const int32_t ra   = (bf.bounded_regions.first == NOT_SET_RegionID)
                               ? -1
                               : static_cast<int32_t>(bf.bounded_regions.first);
        const int32_t rb   = (bf.bounded_regions.second == NOT_SET_RegionID)
                               ? -1
                               : static_cast<int32_t>(bf.bounded_regions.second);
        const bool flip_ex = (ra < 0 && rb >= 0);
        const int32_t med  = (bf.origin_objects.empty() && ra >= 0 && rb >= 0) ? 1 : 0;
        const int32_t val  = (ra >= 0 ? 1 : 0) + (rb >= 0 ? 1 : 0);

        int32_t own = -1;
        if (bf.origin_objects.size() == 1)
        {
            own = static_cast<int32_t>(bf.origin_objects[0]);
        }
        else if ((ra >= 0) != (rb >= 0))
        {
            const int32_t r_only = (ra >= 0) ? ra : rb;
            if (static_cast<std::size_t>(r_only) < res.brep_volumes.size() &&
                res.brep_volumes[r_only].owner_objects.size() == 1)
            {
                own = static_cast<int32_t>(res.brep_volumes[r_only].owner_objects[0]);
            }
        }

        for (const auto& tri : bf.triangles)
        {
            connectivity.push_back(static_cast<int64_t>(tri[0]));
            connectivity.push_back(static_cast<int64_t>(flip_ex ? tri[2] : tri[1]));
            connectivity.push_back(static_cast<int64_t>(flip_ex ? tri[1] : tri[2]));
            offsets.push_back(static_cast<int64_t>(connectivity.size()));
            types.push_back(5); // VTK_TRIANGLE

            entity_dim.push_back(2);
            entity_tag.push_back(static_cast<int32_t>(bf.boundary_id));
            region_a.push_back(ra);
            region_b.push_back(rb);
            owner_obj.push_back(own);
            is_medial.push_back(med);
            valence.push_back(val);
        }
    }

    // 2. 1D BRepCurve polyline segments (VTK_LINE = 3)
    for (std::size_t ci = 0; ci < res.brep_curves.size(); ++ci)
    {
        const auto& bc    = res.brep_curves[ci];
        const int32_t val = static_cast<int32_t>(bc.incident_faces.size());
        for (std::size_t k = 0; k + 1 < bc.polyline_vertices.size(); ++k)
        {
            connectivity.push_back(static_cast<int64_t>(bc.polyline_vertices[k]));
            connectivity.push_back(static_cast<int64_t>(bc.polyline_vertices[k + 1]));
            offsets.push_back(static_cast<int64_t>(connectivity.size()));
            types.push_back(3); // VTK_LINE

            entity_dim.push_back(1);
            entity_tag.push_back(static_cast<int32_t>(ci));
            region_a.push_back(-1);
            region_b.push_back(-1);
            owner_obj.push_back(-1);
            is_medial.push_back(0);
            valence.push_back(val);
        }
    }

    // 3. 0D BRepVertex corners / loop anchors (VTK_VERTEX = 1)
    for (std::size_t vi = 0; vi < res.brep_vertices.size(); ++vi)
    {
        const auto& bv = res.brep_vertices[vi];
        if (bv.global_vertex_id >= points.size())
        {
            continue;
        }
        connectivity.push_back(static_cast<int64_t>(bv.global_vertex_id));
        offsets.push_back(static_cast<int64_t>(connectivity.size()));
        types.push_back(1); // VTK_VERTEX

        entity_dim.push_back(0);
        entity_tag.push_back(static_cast<int32_t>(vi));
        region_a.push_back(-1);
        region_b.push_back(-1);
        owner_obj.push_back(-1);
        is_medial.push_back(0);
        valence.push_back(static_cast<int32_t>(bv.incident_curves.size()));
    }

    detail::write_vtu_unstructured_grid(
      path,
      points,
      connectivity,
      offsets,
      types,
      {
        { "EntityDim", std::move(entity_dim) },
        { "EntityTag", std::move(entity_tag) },
        { "BoundedRegionA", std::move(region_a) },
        { "BoundedRegionB", std::move(region_b) },
        { "OwnerObjectID", std::move(owner_obj) },
        { "IsMedialSheet", std::move(is_medial) },
        { "Valence", std::move(valence) },
      });
}

/// @brief Writes the segmented 3D tetrahedral volumes (`VTK_TETRA = 10`) into a `.vtu` file.
inline void write_volumes_vtu(const std::string& path, const ArrangeResult& res)
{
    const auto& tet_mesh = res.segmented_tet.tetrahedralization;
    const auto& t2r      = res.segmented_tet.region_info.tet_to_region_map;

    std::vector<int64_t> old_to_new(tet_mesh.vertices.size(), -1);
    std::vector<std::array<double, 3>> points;
    std::vector<int64_t> connectivity;
    std::vector<int64_t> offsets;
    std::vector<uint8_t> types;

    std::vector<int32_t> vol_id;
    std::vector<int32_t> owner_obj;
    std::vector<int32_t> num_owners;
    std::vector<int32_t> is_overlap;

    auto map_vid = [&](Index v) -> int64_t
    {
        if (old_to_new[v] < 0)
        {
            old_to_new[v] = static_cast<int64_t>(points.size());
            points.push_back(tet_mesh.vertices[v]);
        }
        return old_to_new[v];
    };

    for (std::size_t tid = 0; tid < tet_mesh.tetrahedra.size(); ++tid)
    {
        if (tid >= t2r.size())
        {
            break;
        }
        const RegionID r = t2r[tid];
        if (r == NOT_SET_RegionID || static_cast<std::size_t>(r) >= res.brep_volumes.size())
        {
            continue;
        }

        auto tet       = tet_mesh.tetrahedra[tid];
        const auto& p0 = tet_mesh.vertices[tet[0]];
        const auto& p1 = tet_mesh.vertices[tet[1]];
        const auto& p2 = tet_mesh.vertices[tet[2]];
        const auto& p3 = tet_mesh.vertices[tet[3]];

        // Ensure positive VTK_TETRA Jacobian: ((p1 - p0) x (p2 - p0)) . (p3 - p0) > 0
        const double ux = p1[0] - p0[0], uy = p1[1] - p0[1], uz = p1[2] - p0[2];
        const double vx = p2[0] - p0[0], vy = p2[1] - p0[1], vz = p2[2] - p0[2];
        const double wx = p3[0] - p0[0], wy = p3[1] - p0[1], wz = p3[2] - p0[2];
        const double det =
          (uy * vz - uz * vy) * wx + (uz * vx - ux * vz) * wy + (ux * vy - uy * vx) * wz;
        if (det < 0.0)
        {
            std::swap(tet[1], tet[2]);
        }

        connectivity.push_back(map_vid(tet[0]));
        connectivity.push_back(map_vid(tet[1]));
        connectivity.push_back(map_vid(tet[2]));
        connectivity.push_back(map_vid(tet[3]));
        offsets.push_back(static_cast<int64_t>(connectivity.size()));
        types.push_back(10); // VTK_TETRA

        const auto& bvol = res.brep_volumes[r];
        vol_id.push_back(static_cast<int32_t>(r));
        owner_obj.push_back(bvol.owner_objects.size() == 1
                              ? static_cast<int32_t>(bvol.owner_objects[0])
                              : -1);
        num_owners.push_back(static_cast<int32_t>(bvol.owner_objects.size()));
        is_overlap.push_back(bvol.is_overlap ? 1 : 0);
    }

    detail::write_vtu_unstructured_grid(
      path,
      points,
      connectivity,
      offsets,
      types,
      {
        { "BRepVolumeID", std::move(vol_id) },
        { "OwnerObjectID", std::move(owner_obj) },
        { "NumOwners", std::move(num_owners) },
        { "IsOverlap", std::move(is_overlap) },
      });
}

/// @brief Exports standard `.vtu` meshes (`_inputs.vtu`, `_brep.vtu`, `_volumes.vtu`, and optional
///        `_before_volumes.vtu`) for a demonstration into `out_dir`.
inline void export_demo_meshes_vtu(const std::string& out_dir,
                                   const std::string& stem,
                                   const std::vector<InputSurface>& inputs,
                                   const ArrangeResult& primary_result,
                                   const ArrangeResult* before_result = nullptr)
{
    if (out_dir.empty())
    {
        return;
    }
    const std::filesystem::path dir(out_dir);
    std::filesystem::create_directories(dir);

    write_inputs_vtu((dir / (stem + "_inputs.vtu")).string(), inputs);
    write_brep_vtu((dir / (stem + "_brep.vtu")).string(), primary_result);
    write_volumes_vtu((dir / (stem + "_volumes.vtu")).string(), primary_result);
    if (before_result != nullptr)
    {
        write_volumes_vtu((dir / (stem + "_before_volumes.vtu")).string(), *before_result);
    }
}

struct DemoCliArgs
{
    std::string vtu_dir = ".";
    std::string json_path;
};

/// @brief Parses `[vtu_output_dir] [json_bundle_path]` (or legacy `[json_bundle_path] [vtu_output_dir]`).
inline DemoCliArgs parse_demo_cli(int argc, char** argv)
{
    DemoCliArgs args;
    auto ends_with_json = [](const std::string& s)
    {
        return s.size() >= 5 && s.compare(s.size() - 5, 5, ".json") == 0;
    };

    if (argc >= 2)
    {
        const std::string a1 = argv[1];
        if (ends_with_json(a1))
        {
            args.json_path = a1;
            args.vtu_dir   = (argc >= 3) ? argv[2] : "";
        }
        else
        {
            args.vtu_dir   = a1;
            args.json_path = (argc >= 3) ? argv[2] : "";
        }
    }
    return args;
}

} // namespace mesh_arrange::examples
