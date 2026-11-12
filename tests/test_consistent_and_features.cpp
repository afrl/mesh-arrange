#include "mesh_arrange/MeshArrange.h"

#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
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

InputSurface make_box_surface(const std::array<double, 3>& neg,
                              const std::array<double, 3>& pos,
                              ObjectType type       = ObjectType::SOLID,
                              PriorityType priority = 0,
                              bool six_facets       = true)
{
    InputSurface s;
    s.type     = type;
    s.priority = priority;
    s.vertices = {
        neg,
        { pos[0], neg[1], neg[2] },
        { pos[0], pos[1], neg[2] },
        { neg[0], pos[1], neg[2] },
        { neg[0], neg[1], pos[2] },
        { pos[0], neg[1], pos[2] },
        pos,
        { neg[0], pos[1], pos[2] }
    };
    s.triangles = {
        { 2, 1, 0 }, { 2, 0, 3 },
        { 6, 5, 1 }, { 6, 1, 2 },
        { 0, 4, 7 }, { 0, 7, 3 },
        { 1, 5, 0 }, { 0, 5, 4 },
        { 3, 7, 6 }, { 3, 6, 2 },
        { 6, 7, 4 }, { 6, 4, 5 }
    };
    if (six_facets)
    {
        s.triangle_facet_ids = { 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5 };
    }
    return s;
}

/// Verify that a `BRepVolume` assembled from signed `BRepFace` shells is a closed, outward-oriented
/// watertight 2-manifold and return its signed enclosed volume.
bool check_brep_volume_shell(const ArrangeResult& res,
                             const BRepVolume& bvol,
                             double& signed_vol_out)
{
    signed_vol_out = 0.0;
    std::map<std::pair<Index, Index>, int> undirected;
    std::map<std::pair<Index, Index>, int> directed;

    const auto& verts = res.surface_mesh.vertices;
    for (const auto& [fid, dir] : bvol.signed_bounding_faces)
    {
        if (fid >= res.brep_faces.size())
        {
            return false;
        }
        for (auto tri : res.brep_faces[fid].triangles)
        {
            if (dir < 0)
            {
                std::swap(tri[1], tri[2]);
            }
            const auto& a = verts[tri[0]];
            const auto& b = verts[tri[1]];
            const auto& c = verts[tri[2]];
            const double cx = b[1] * c[2] - b[2] * c[1];
            const double cy = b[2] * c[0] - b[0] * c[2];
            const double cz = b[0] * c[1] - b[1] * c[0];
            signed_vol_out += (a[0] * cx + a[1] * cy + a[2] * cz) / 6.0;

            for (int e = 0; e < 3; ++e)
            {
                const Index u = tri[e];
                const Index v = tri[(e + 1) % 3];
                undirected[{ std::min(u, v), std::max(u, v) }]++;
                directed[{ u, v }]++;
            }
        }
    }

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
    return signed_vol_out > 0.0;
}

} // namespace

int main()
{
    // 1. Pure arrangement WITHOUT overlap removal (remove_overlaps = false):
    //    Two equal-priority boxes decompose into 3 disjoint B-Rep volumes:
    //      - A \ B (vol = 1.3*1.2*1.2 - 0.6*1.0*1.0 = 1.272)
    //      - A n B (vol = 0.6*1.0*1.0 = 0.600, is_overlap == true, owner_objects == {0, 1})
    //      - B \ A (vol = 1.3*1.0*1.0 - 0.6*1.0*1.0 = 0.700)
    {
        const auto box_a = make_box_surface({ -1.0, -0.6, -0.6 }, { 0.3, 0.6, 0.6 }, ObjectType::SOLID, 0);
        const auto box_b = make_box_surface({ -0.3, -0.5, -0.5 }, { 1.0, 0.5, 0.5 }, ObjectType::SOLID, 0);

        ArrangeOptions opts;
        opts.remove_overlaps = false;
        opts.check_topology  = true;

        const auto res = arrange_surfaces({ box_a, box_b }, opts);
        TEST_CHECK(res.brep_volumes.size() == 3);
        TEST_CHECK(!res.brep_faces.empty());
        TEST_CHECK(!res.brep_curves.empty());
        TEST_CHECK(!res.brep_vertices.empty());

        int overlap_count = 0;
        double total_vol  = 0.0;
        for (const auto& bvol : res.brep_volumes)
        {
            double sv = 0.0;
            TEST_CHECK(check_brep_volume_shell(res, bvol, sv));
            total_vol += sv;
            if (bvol.is_overlap)
            {
                ++overlap_count;
                TEST_CHECK(bvol.owner_objects.size() == 2);
                TEST_CHECK(std::abs(sv - 0.6) < 1e-6);
            }
            else
            {
                TEST_CHECK(bvol.owner_objects.size() == 1);
            }
        }
        TEST_CHECK(overlap_count == 1);
        TEST_CHECK(std::abs(total_vol - 2.572) < 1e-6);
    }

    // 2. Priority cutting + CLIP_INSIDE hole + CLIP_OUTSIDE bounding box (including negative priority!)
    {
        // Inclusion (priority +5): [-0.8, 0.8] x [-0.3, 0.3] x [-0.3, 0.3]
        // Matrix (priority -2):    [-1.0, 1.0] x [-1.0, 1.0] x [-1.0, 1.0]
        // Hole (CLIP_INSIDE, +5):  [-0.2, 0.2] x [-0.2, 0.2] x [-1.2, 1.2] (cuts through both!)
        // Domain (CLIP_OUTSIDE,+5):[-0.5, 0.5] x [-0.5, 0.5] x [-0.5, 0.5] (clips domain to 1x1x1 cube!)
        // Net domain volume = 1.0^3 - (0.4 * 0.4 * 1.0) = 1.0 - 0.16 = 0.84!
        // Inclusion inside domain minus hole = (1.0 * 0.6 * 0.6) - (0.4 * 0.4 * 0.6) = 0.36 - 0.096 = 0.264!
        // Matrix inside domain minus inclusion & hole = 0.84 - 0.264 = 0.576!
        const auto inclusion = make_box_surface({ -0.8, -0.3, -0.3 }, { 0.8, 0.3, 0.3 }, ObjectType::SOLID, 5);
        const auto matrix    = make_box_surface({ -1.0, -1.0, -1.0 }, { 1.0, 1.0, 1.0 }, ObjectType::SOLID, -2);
        const auto hole      = make_box_surface({ -0.2, -0.2, -1.2 }, { 0.2, 0.2, 1.2 }, ObjectType::CLIP_INSIDE, 5);
        const auto domain    = make_box_surface({ -0.5, -0.5, -0.5 }, { 0.5, 0.5, 0.5 }, ObjectType::CLIP_OUTSIDE, 5);

        ArrangeOptions opts;
        opts.remove_overlaps = true;
        opts.check_topology  = true;

        const auto res = arrange_surfaces({ inclusion, matrix, hole, domain }, opts);
        TEST_CHECK(res.constituent_surfaces.size() == 2);
        TEST_CHECK(res.constituent_object_ids.size() == 2);
        TEST_CHECK(res.constituent_object_ids[0] == 0);
        TEST_CHECK(res.constituent_object_ids[1] == 1);

        double inc_vol = 0.0;
        double mat_vol = 0.0;
        for (const auto& bvol : res.brep_volumes)
        {
            TEST_CHECK(!bvol.is_overlap);
            TEST_CHECK(bvol.owner_objects.size() == 1);
            double sv = 0.0;
            TEST_CHECK(check_brep_volume_shell(res, bvol, sv));
            if (bvol.owner_objects[0] == 0)
            {
                inc_vol += sv;
            }
            else if (bvol.owner_objects[0] == 1)
            {
                mat_vol += sv;
            }
        }
        TEST_CHECK(std::abs(inc_vol - 0.264) < 1e-6);
        TEST_CHECK(std::abs(mat_vol - 0.576) < 1e-6);
    }

    std::cout << "test_consistent_and_features passed." << std::endl;
    return 0;
}
