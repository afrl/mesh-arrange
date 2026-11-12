#include "mesh_arrange/MeshArrange.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <utility>
#include <vector>

#define TEST_CHECK(cond)                                                                           \
    do                                                                                             \
    {                                                                                              \
        if (!(cond))                                                                               \
        {                                                                                          \
            std::cerr << "TEST_CHECK failed at line " << __LINE__ << ": " << #cond << std::endl;   \
            return 1;                                                                              \
        }                                                                                          \
    } while (0)

using namespace mesh_arrange;

namespace
{

constexpr double kPi = 3.14159265358979323846264338327950288;

InputSurface make_box(const std::array<double, 3>& neg_bound,
                      const std::array<double, 3>& pos_bound,
                      PriorityType priority = 1)
{
    InputSurface box;
    box.type     = ObjectType::SOLID;
    box.priority = priority;
    box.vertices = {
        neg_bound,
        { pos_bound[0], neg_bound[1], neg_bound[2] },
        { pos_bound[0], pos_bound[1], neg_bound[2] },
        { neg_bound[0], pos_bound[1], neg_bound[2] },
        { neg_bound[0], neg_bound[1], pos_bound[2] },
        { pos_bound[0], neg_bound[1], pos_bound[2] },
        pos_bound,
        { neg_bound[0], pos_bound[1], pos_bound[2] }
    };
    box.triangles = {
        { 2, 1, 0 }, { 2, 0, 3 },
        { 6, 5, 1 }, { 6, 1, 2 },
        { 0, 4, 7 }, { 0, 7, 3 },
        { 1, 5, 0 }, { 0, 5, 4 },
        { 3, 7, 6 }, { 3, 6, 2 },
        { 6, 7, 4 }, { 6, 4, 5 }
    };
    box.triangle_facet_ids = { 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5 };
    return box;
}

InputSurface create_u_shape_mesh(PriorityType priority = 1)
{
    const std::array<std::array<double, 2>, 8> poly = { { { -1.0, -1.0 },
                                                          { 1.0, -1.0 },
                                                          { 1.0, 0.8 },
                                                          { 0.4, 0.8 },
                                                          { 0.4, -0.4 },
                                                          { -0.4, -0.4 },
                                                          { -0.4, 0.8 },
                                                          { -1.0, 0.8 } } };
    const double z0 = -0.35;
    const double z1 = +0.35;

    InputSurface m;
    m.type     = ObjectType::SOLID;
    m.priority = priority;
    for (int i = 0; i < 8; ++i)
    {
        m.vertices.push_back({ poly[i][0], poly[i][1], z0 });
    }
    for (int i = 0; i < 8; ++i)
    {
        m.vertices.push_back({ poly[i][0], poly[i][1], z1 });
    }

    const std::array<std::array<Index, 3>, 6> cap_tris = {
        { { 0, 5, 7 }, { 5, 6, 7 }, { 0, 1, 4 }, { 0, 4, 5 }, { 4, 1, 2 }, { 4, 2, 3 } }
    };

    for (const auto& t : cap_tris)
    {
        m.triangles.push_back({ t[0], t[2], t[1] });
        m.triangle_facet_ids.push_back(0);
    }
    for (const auto& t : cap_tris)
    {
        m.triangles.push_back({ t[0] + 8, t[1] + 8, t[2] + 8 });
        m.triangle_facet_ids.push_back(1);
    }
    for (Index i = 0; i < 8; ++i)
    {
        const Index j     = (i + 1) % 8;
        const FacetID cid = static_cast<FacetID>(2 + i);
        m.triangles.push_back({ i, j, j + 8 });
        m.triangle_facet_ids.push_back(cid);
        m.triangles.push_back({ i, j + 8, i + 8 });
        m.triangle_facet_ids.push_back(cid);
    }
    return m;
}

InputSurface create_oriented_cylinder(const double theta,
                                      const double radius,
                                      const double half_len,
                                      const int n_radial     = 12,
                                      const int n_axial      = 5,
                                      const double phase_rad = 0.0)
{
    const std::array<double, 3> axis = { std::cos(theta), std::sin(theta), 0.0 };
    const std::array<double, 3> u    = { -std::sin(theta), std::cos(theta), 0.0 };
    const std::array<double, 3> v    = { 0.0, 0.0, 1.0 };

    InputSurface m;
    m.type     = ObjectType::SOLID;
    m.priority = 1;

    for (int s = 0; s <= n_axial; ++s)
    {
        const double z_pos =
          -half_len + (2.0 * half_len * static_cast<double>(s)) / static_cast<double>(n_axial);
        for (int k = 0; k < n_radial; ++k)
        {
            const double ang = phase_rad + (2.0 * kPi * k) / static_cast<double>(n_radial);
            const double ca  = std::cos(ang);
            const double sa  = std::sin(ang);
            m.vertices.push_back({ z_pos * axis[0] + radius * (ca * u[0] + sa * v[0]),
                                   z_pos * axis[1] + radius * (ca * u[1] + sa * v[1]),
                                   z_pos * axis[2] + radius * (ca * u[2] + sa * v[2]) });
        }
    }
    const Index c0 = static_cast<Index>(m.vertices.size());
    m.vertices.push_back({ -half_len * axis[0], -half_len * axis[1], -half_len * axis[2] });
    const Index c1 = static_cast<Index>(m.vertices.size());
    m.vertices.push_back({ +half_len * axis[0], +half_len * axis[1], +half_len * axis[2] });

    for (int s = 0; s < n_axial; ++s)
    {
        const int base0 = s * n_radial;
        const int base1 = (s + 1) * n_radial;
        for (int k = 0; k < n_radial; ++k)
        {
            const Index a0 = static_cast<Index>(base0 + k);
            const Index a1 = static_cast<Index>(base0 + ((k + 1) % n_radial));
            const Index b0 = static_cast<Index>(base1 + k);
            const Index b1 = static_cast<Index>(base1 + ((k + 1) % n_radial));

            m.triangles.push_back({ a0, a1, b1 });
            m.triangle_facet_ids.push_back(0);
            m.triangles.push_back({ a0, b1, b0 });
            m.triangle_facet_ids.push_back(0);
        }
    }

    const int last_base = n_axial * n_radial;
    for (int k = 0; k < n_radial; ++k)
    {
        const Index a0 = static_cast<Index>(k);
        const Index a1 = static_cast<Index>((k + 1) % n_radial);
        const Index b0 = static_cast<Index>(last_base + k);
        const Index b1 = static_cast<Index>(last_base + ((k + 1) % n_radial));

        m.triangles.push_back({ c0, a1, a0 });
        m.triangle_facet_ids.push_back(1);
        m.triangles.push_back({ c1, b0, b1 });
        m.triangle_facet_ids.push_back(2);
    }
    return m;
}

double tet_vol(const Point3& p0, const Point3& p1, const Point3& p2, const Point3& p3)
{
    const double ax = p1[0] - p0[0], ay = p1[1] - p0[1], az = p1[2] - p0[2];
    const double bx = p2[0] - p0[0], by = p2[1] - p0[1], bz = p2[2] - p0[2];
    const double cx = p3[0] - p0[0], cy = p3[1] - p0[1], cz = p3[2] - p0[2];
    return (ax * (by * cz - bz * cy) + ay * (bz * cx - bx * cz) + az * (bx * cy - by * cx)) / 6.0;
}

double compute_tet_volume(const Tetrahedralization& t, double& min_tet_vol_out)
{
    double total    = 0.0;
    min_tet_vol_out = std::numeric_limits<double>::max();
    for (const auto& tet : t.tetrahedra)
    {
        const double v = tet_vol(t.vertices[tet[0]],
                                 t.vertices[tet[1]],
                                 t.vertices[tet[2]],
                                 t.vertices[tet[3]]);
        min_tet_vol_out = std::min(min_tet_vol_out, v);
        total += v;
    }
    return total;
}

bool is_watertight_manifold(const Triangulation& tri, double& enclosed_vol_out)
{
    if (tri.triangles.empty())
    {
        enclosed_vol_out = 0.0;
        return false;
    }
    std::map<std::pair<Index, Index>, int> undirected;
    std::map<std::pair<Index, Index>, int> directed;
    enclosed_vol_out = 0.0;

    for (const auto& f : tri.triangles)
    {
        const auto& a = tri.vertices[f[0]];
        const auto& b = tri.vertices[f[1]];
        const auto& c = tri.vertices[f[2]];
        enclosed_vol_out += (a[0] * (b[1] * c[2] - b[2] * c[1]) +
                             a[1] * (b[2] * c[0] - b[0] * c[2]) +
                             a[2] * (b[0] * c[1] - b[1] * c[0])) /
                            6.0;

        for (int e = 0; e < 3; ++e)
        {
            const Index u = f[e];
            const Index v = f[(e + 1) % 3];
            undirected[{ std::min(u, v), std::max(u, v) }]++;
            directed[{ u, v }]++;
        }
    }
    enclosed_vol_out = std::abs(enclosed_vol_out);

    for (const auto& kv : undirected)
    {
        if (kv.second != 2)
        {
            return false;
        }
    }
    for (const auto& kv : directed)
    {
        if (kv.second != 1)
        {
            return false;
        }
    }
    return true;
}

} // namespace

int main()
{
    // Case A: Convex Lens Overlap (K = 2)
    {
        std::cout << "=== Case A: Convex Lens Overlap (K = 2) ===" << std::endl;
        const auto box_a = make_box({ -1.0, -0.6, -0.6 }, { 0.3, 0.6, 0.6 });
        const auto box_b = make_box({ -0.3, -0.5, -0.5 }, { 1.0, 0.5, 0.5 });

        const auto res = arrange_surfaces({ box_a, box_b });
        TEST_CHECK(res.segmented_tet.object_region_info.regions_of_overlap.empty());
        TEST_CHECK(check_tri_edge_relations(res.segmented_tet.tetrahedralization));
        TEST_CHECK(check_tet_edge_relations(res.segmented_tet.tetrahedralization));
        TEST_CHECK(check_tet_tri_relations(res.segmented_tet.tetrahedralization));

        double min_v       = 0.0;
        const double vol   = compute_tet_volume(res.segmented_tet.tetrahedralization, min_v);
        const double exp_v = (1.3 * 1.2 * 1.2) + (1.3 * 1.0 * 1.0) - (0.6 * 1.0 * 1.0);
        TEST_CHECK(min_v > 0.0);
        TEST_CHECK(std::abs(vol - exp_v) < 1e-6);
        TEST_CHECK(res.constituent_surfaces.size() == 2);

        double v0 = 0.0, v1 = 0.0;
        TEST_CHECK(is_watertight_manifold(res.constituent_surfaces[0], v0));
        TEST_CHECK(is_watertight_manifold(res.constituent_surfaces[1], v1));
        TEST_CHECK(std::abs((v0 + v1) - vol) < 1e-6);
        TEST_CHECK(v0 > 1.27 && v1 > 0.70);
    }

    // Case B: Non-Convex U-Shape Overlap Bridging Central Matrix Gap (K = 2 + background matrix)
    {
        std::cout << "=== Case B: Non-Convex U-Shape Overlap with Central Gap ===" << std::endl;
        const auto u_solid = create_u_shape_mesh(1);
        const auto bar     = make_box({ -1.1, 0.2, -0.25 }, { 1.1, 0.75, 0.25 }, 1);
        const auto matrix  = make_box({ -1.3, -1.2, -0.55 }, { 1.3, 1.05, 0.55 }, 0);

        const auto combined = combine_surfaces({ u_solid, bar, matrix });
        auto consistent     = create_consistent_surface(combined);
        auto st             = identify_features(combined, consistent);
        TEST_CHECK(!st.object_region_info.regions_of_overlap.empty());

        double matrix_vol_before = 0.0;
        for (const RegionID r : st.object_region_info.object_to_regions_map[2])
        {
            for (const Index tid : st.region_info.region_to_tets_map[r])
            {
                const auto& tc = st.tetrahedralization.tetrahedra[tid];
                matrix_vol_before += tet_vol(st.tetrahedralization.vertices[tc[0]],
                                             st.tetrahedralization.vertices[tc[1]],
                                             st.tetrahedralization.vertices[tc[2]],
                                             st.tetrahedralization.vertices[tc[3]]);
            }
        }

        remove_interpenetrations(combined, consistent, st);

        double min_v         = 0.0;
        const double tot_vol = compute_tet_volume(st.tetrahedralization, min_v);
        const double exp_box = 2.6 * 2.25 * 1.10;
        TEST_CHECK(min_v > 0.0);
        TEST_CHECK(std::abs(tot_vol - exp_box) < 1e-6);

        double matrix_vol_after = 0.0;
        for (const RegionID r : st.object_region_info.object_to_regions_map[2])
        {
            for (const Index tid : st.region_info.region_to_tets_map[r])
            {
                const auto& tc = st.tetrahedralization.tetrahedra[tid];
                matrix_vol_after += tet_vol(st.tetrahedralization.vertices[tc[0]],
                                            st.tetrahedralization.vertices[tc[1]],
                                            st.tetrahedralization.vertices[tc[2]],
                                            st.tetrahedralization.vertices[tc[3]]);
            }
        }
        TEST_CHECK(std::abs(matrix_vol_after - matrix_vol_before) < 1e-8);

        auto rt           = recover_triangulation(st);
        auto constituents = separate_objects(rt);
        TEST_CHECK(constituents.size() == 3);
        double sum_c_vol = 0.0;
        for (const auto& c : constituents)
        {
            double cv = 0.0;
            TEST_CHECK(is_watertight_manifold(c, cv));
            TEST_CHECK(cv > 0.0);
            sum_c_vol += cv;
        }
        TEST_CHECK(std::abs(sum_c_vol - exp_box) < 1e-5);
    }

    // Case C: Complete Engulfment (Omega_A subset Omega_B, equal priority)
    {
        std::cout << "=== Case C: Complete Engulfment ===" << std::endl;
        const auto inner_a = make_box({ -0.25, -0.25, -0.25 }, { 0.25, 0.25, 0.25 }, 1);
        const auto outer_b = make_box({ -0.80, -0.80, -0.80 }, { 0.80, 0.80, 0.80 }, 1);

        const auto res     = arrange_surfaces({ inner_a, outer_b });
        const double exp_v = 1.6 * 1.6 * 1.6;
        TEST_CHECK(res.constituent_surfaces.size() == 1);
        TEST_CHECK(res.constituent_object_ids.size() == 1);
        TEST_CHECK(res.constituent_object_ids[0] == 1);

        double cv = 0.0;
        TEST_CHECK(is_watertight_manifold(res.constituent_surfaces[0], cv));
        TEST_CHECK(std::abs(cv - exp_v) < 1e-5);
    }

    // Case D1: Triple Tow Junction (K = 3 cylinders at 0, 120, 240 degrees)
    {
        std::cout << "=== Case D1: Triple Tow Junction (K = 3) ===" << std::endl;
        const auto cyl0 = create_oriented_cylinder(0.0, 0.36, 1.0, 12, 5, 0.03);
        const auto cyl1 = create_oriented_cylinder(2.0 * kPi / 3.0, 0.35, 1.0, 12, 5, 0.17);
        const auto cyl2 = create_oriented_cylinder(4.0 * kPi / 3.0, 0.34, 1.0, 12, 5, 0.31);

        ArrangeOptions opts_no_rem;
        opts_no_rem.remove_overlaps = false;
        const auto res_before       = arrange_surfaces({ cyl0, cyl1, cyl2 }, opts_no_rem);
        double min_b                = 0.0;
        const double union_vol      = compute_tet_volume(res_before.segmented_tet.tetrahedralization, min_b);

        const auto res = arrange_surfaces({ cyl0, cyl1, cyl2 });
        TEST_CHECK(res.constituent_surfaces.size() == 3);
        double sum_c = 0.0;
        for (int k = 0; k < 3; ++k)
        {
            double cv = 0.0;
            TEST_CHECK(is_watertight_manifold(res.constituent_surfaces[k], cv));
            TEST_CHECK(std::abs(cv - union_vol / 3.0) < 0.12 * union_vol);
            sum_c += cv;
        }
        TEST_CHECK(std::abs(sum_c - union_vol) < 1e-5);
    }

    // Case D2: Quadruple Overlap Junction (K = 4 mutually overlapping boxes)
    {
        std::cout << "=== Case D2: Quadruple Overlap Junction (K = 4) ===" << std::endl;
        const auto b0 = make_box({ -0.80, -0.76, -0.72 }, { 0.25, 0.22, 0.20 }, 1);
        const auto b1 = make_box({ -0.24, -0.79, -0.75 }, { 0.80, 0.24, 0.23 }, 1);
        const auto b2 = make_box({ -0.48, -0.23, -0.78 }, { 0.48, 0.80, 0.26 }, 1);
        const auto b3 = make_box({ -0.44, -0.46, -0.21 }, { 0.46, 0.44, 0.80 }, 1);

        ArrangeOptions opts_no_rem;
        opts_no_rem.remove_overlaps = false;
        const auto res_before       = arrange_surfaces({ b0, b1, b2, b3 }, opts_no_rem);
        double min_b                = 0.0;
        const double union_vol      = compute_tet_volume(res_before.segmented_tet.tetrahedralization, min_b);

        const auto res = arrange_surfaces({ b0, b1, b2, b3 });
        TEST_CHECK(res.constituent_surfaces.size() == 4);
        double sum_c = 0.0;
        for (int k = 0; k < 4; ++k)
        {
            double cv = 0.0;
            TEST_CHECK(is_watertight_manifold(res.constituent_surfaces[k], cv));
            TEST_CHECK(cv > 0.0);
            sum_c += cv;
        }
        TEST_CHECK(std::abs(sum_c - union_vol) < 1e-5);
    }

    std::cout << "test_interpenetration_cases passed." << std::endl;
    return 0;
}
