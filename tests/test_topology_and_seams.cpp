#include "mesh_arrange/MeshArrange.h"
#include "mesh_arrange/VertexMap.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <string>
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

double tet_signed_volume_stl(const Tetrahedralization& t, const std::array<Index, 4>& tet)
{
    const auto& p0 = t.vertices[tet[0]];
    const auto& p1 = t.vertices[tet[1]];
    const auto& p2 = t.vertices[tet[2]];
    const auto& p3 = t.vertices[tet[3]];

    const double ax = p1[0] - p0[0], ay = p1[1] - p0[1], az = p1[2] - p0[2];
    const double bx = p2[0] - p0[0], by = p2[1] - p0[1], bz = p2[2] - p0[2];
    const double cx = p3[0] - p0[0], cy = p3[1] - p0[1], cz = p3[2] - p0[2];

    const double crx = by * cz - bz * cy;
    const double cry = bz * cx - bx * cz;
    const double crz = bx * cy - by * cx;
    return (ax * crx + ay * cry + az * crz) / 6.0;
}

Triangulation make_unit_cube_triangulation()
{
    Triangulation t;
    t.vertices = {
        { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 1.0, 1.0, 0.0 }, { 0.0, 1.0, 0.0 },
        { 0.0, 0.0, 1.0 }, { 1.0, 0.0, 1.0 }, { 1.0, 1.0, 1.0 }, { 0.0, 1.0, 1.0 }
    };
    t.triangles = {
        { 0, 2, 1 }, { 0, 3, 2 },
        { 4, 5, 6 }, { 4, 6, 7 },
        { 0, 1, 5 }, { 0, 5, 4 },
        { 2, 3, 7 }, { 2, 7, 6 },
        { 1, 2, 6 }, { 1, 6, 5 },
        { 3, 0, 4 }, { 3, 4, 7 }
    };
    t.tri_to_marker_map = { 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5 };
    update_edges_and_maps(t);
    return t;
}

Tetrahedralization make_kuhn_cube_tet()
{
    Tetrahedralization t;
    t.vertices = {
        { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0 }, { 1.0, 1.0, 0.0 },
        { 0.0, 0.0, 1.0 }, { 1.0, 0.0, 1.0 }, { 0.0, 1.0, 1.0 }, { 1.0, 1.0, 1.0 }
    };
    t.tetrahedra = {
        { 0, 1, 3, 7 },
        { 0, 3, 2, 7 },
        { 0, 2, 6, 7 },
        { 0, 6, 4, 7 },
        { 0, 4, 5, 7 },
        { 0, 5, 1, 7 }
    };
    t.triangles = {
        { 0, 3, 1 }, { 0, 2, 3 }, { 0, 7, 3 }
    };
    update_edges_and_maps(t);
    return t;
}

} // namespace

int main()
{
    // 1. Symmetric hash_linear_edge & Welford edge length stats
    {
        TEST_CHECK(hash_linear_edge(3, 17) == hash_linear_edge(17, 3));
        TEST_CHECK(hash_linear_edge(3, 17) != hash_linear_edge(4, 16));
        const Index big_a = static_cast<Index>(0x7fffffffU);
        const Index big_b = static_cast<Index>(0x80000001U);
        TEST_CHECK(hash_linear_edge(big_a, big_b) == hash_linear_edge(big_b, big_a));

        Triangulation empty_t;
        const auto empty_stats = query_edge_length_stats(empty_t);
        TEST_CHECK(empty_stats.min == 0.0 && empty_stats.max == 0.0 && empty_stats.mean == 0.0 &&
                   empty_stats.variance == 0.0 && empty_stats.std_dev == 0.0);

        Triangulation right_tri;
        right_tri.vertices          = { { 0.0, 0.0, 0.0 }, { 3.0, 0.0, 0.0 }, { 0.0, 4.0, 0.0 } };
        right_tri.triangles         = { { 0, 1, 2 } };
        right_tri.tri_to_marker_map = { 0 };
        update_edges_and_maps(right_tri);
        TEST_CHECK(check_tri_edge_relations(right_tri));

        const auto stats = query_edge_length_stats(right_tri);
        TEST_CHECK(std::abs(stats.min - 3.0) < 1e-12);
        TEST_CHECK(std::abs(stats.max - 5.0) < 1e-12);
        TEST_CHECK(std::abs(stats.mean - 4.0) < 1e-12);
        TEST_CHECK(std::abs(stats.variance - 1.0) < 1e-12);
        TEST_CHECK(std::abs(stats.std_dev - 1.0) < 1e-12);
    }

    // 2. Cube topology, StaticVertexMap, and seam segmentation
    {
        auto cube = make_unit_cube_triangulation();
        TEST_CHECK(check_tri_edge_relations(cube));
        TEST_CHECK(cube.edges.size() == 18);
        TEST_CHECK(cube.n_markers == 6);

        const auto v2t               = build_vert_to_tri_map(cube);
        std::size_t total_incidences = 0;
        for (Index v = 0; v < cube.vertices.size(); ++v)
        {
            for (const Index tri_id : v2t.get_indices(v))
            {
                TEST_CHECK(tri_id < cube.triangles.size());
                ++total_incidences;
            }
        }
        TEST_CHECK(total_incidences == 3 * cube.triangles.size());

        const auto seams = segment_seams(cube);
        TEST_CHECK(seams.edges_on_seams.size() == 12);
        TEST_CHECK(seams.n_seams == 12);

        const auto seams_inferred =
          segment_seams(cube.edges, cube.tri_to_edges_map, 0, cube.tri_to_marker_map);
        TEST_CHECK(seams_inferred.edges_on_seams.size() == 12);
        TEST_CHECK(seams_inferred.n_seams == 12);

        const auto seams_empty = segment_seams({}, {}, 0, {});
        TEST_CHECK(seams_empty.edges_on_seams.empty() && seams_empty.n_seams == 0);
    }

    // 3. Tetrahedralization relations & edge splits (2D & 3D)
    {
        auto tet = make_kuhn_cube_tet();
        TEST_CHECK(check_tri_edge_relations(tet));
        TEST_CHECK(check_tet_edge_relations(tet));
        TEST_CHECK(check_tet_tri_relations(tet));

        double initial_vol = 0.0;
        for (const auto& c : tet.tetrahedra)
        {
            const double v = tet_signed_volume_stl(tet, c);
            TEST_CHECK(v > 1e-12);
            initial_vol += v;
        }
        TEST_CHECK(std::abs(initial_vol - 1.0) < 1e-12);

        const std::vector<Index> edges_to_split = { 0, 1, 2, 3, 4 };
        const auto res                          = split_edges(edges_to_split, tet);
        TEST_CHECK(check_tri_edge_relations(tet));
        TEST_CHECK(check_tet_edge_relations(tet));
        TEST_CHECK(check_tet_tri_relations(tet));

        double post_vol = 0.0;
        for (const auto& c : tet.tetrahedra)
        {
            const double v = tet_signed_volume_stl(tet, c);
            TEST_CHECK(v > 1e-12);
            post_vol += v;
        }
        TEST_CHECK(std::abs(post_vol - 1.0) < 1e-12);

        std::size_t n_subdivided = 0;
        for (const auto& ue : res.updated_edges)
        {
            if (ue.update_type == UpdatedEdgeInfo::SUBDIVIDED_EDGE)
            {
                ++n_subdivided;
            }
        }
        TEST_CHECK(n_subdivided == edges_to_split.size());
    }

    // 4. Seam-preserving surface split_large_edges, optimize_edge_flips, and tangential smoothing
    {
        Triangulation oct;
        oct.vertices = {
            { 1.0, 0.0, 0.0 },
            { 0.0, 1.0, 0.0 },
            { -1.0, 0.0, 0.0 },
            { 0.0, -1.0, 0.0 },
            { 0.0, 0.0, 1.5 },
            { 0.0, 0.0, -1.5 }
        };
        oct.triangles = {
            { 0, 1, 4 }, { 1, 2, 4 }, { 2, 3, 4 }, { 3, 0, 4 },
            { 1, 0, 5 }, { 2, 1, 5 }, { 3, 2, 5 }, { 0, 3, 5 }
        };
        oct.tri_to_marker_map = { 0, 0, 0, 0, 1, 1, 1, 1 };
        update_edges_and_maps(oct);

        auto seams = segment_seams(oct);
        TEST_CHECK(seams.edges_on_seams.size() == 4);
        TEST_CHECK(seams.n_seams == 1);

        split_large_edges(0.6, oct, seams);
        TEST_CHECK(check_tri_edge_relations(oct));
        TEST_CHECK(query_edge_length_stats(oct).max <= 0.6 + 1e-12);

        double seam_perim = 0.0;
        for (const Index eid : seams.edges_on_seams)
        {
            TEST_CHECK(seams.edge_id_to_seam_id_map[eid] == 0);
            const auto& e   = oct.edges[eid];
            const double dx = oct.vertices[e[0]][0] - oct.vertices[e[1]][0];
            const double dy = oct.vertices[e[0]][1] - oct.vertices[e[1]][1];
            const double dz = oct.vertices[e[0]][2] - oct.vertices[e[1]][2];
            seam_perim += std::sqrt(dx * dx + dy * dy + dz * dz);
        }
        TEST_CHECK(std::abs(seam_perim - 4.0 * std::sqrt(2.0)) < 1e-10);

        optimize_edge_flips(oct, 3);
        TEST_CHECK(check_tri_edge_relations(oct));

        optimize_tangential_smoothing(oct, seams, 3, true);
        TEST_CHECK(check_tri_edge_relations(oct));
    }

    // 5. Equidistant closed loop tangential smoothing termination check
    {
        Triangulation hex_drum;
        const int N = 6;
        for (int k = 0; k < N; ++k)
        {
            const double theta = 2.0 * kPi * static_cast<double>(k) / static_cast<double>(N);
            hex_drum.vertices.push_back({ std::cos(theta), std::sin(theta), 0.0 });
        }
        hex_drum.vertices.push_back({ 0.15, -0.1, 0.5 });
        hex_drum.vertices.push_back({ -0.1, 0.15, -0.5 });
        for (int k = 0; k < N; ++k)
        {
            const Index next = static_cast<Index>((k + 1) % N);
            hex_drum.triangles.push_back({ static_cast<Index>(k), next, 6 });
            hex_drum.tri_to_marker_map.push_back(0);
            hex_drum.triangles.push_back({ next, static_cast<Index>(k), 7 });
            hex_drum.tri_to_marker_map.push_back(1);
        }
        update_edges_and_maps(hex_drum);
        auto hex_seams = segment_seams(hex_drum);
        split_large_edges(0.6, hex_drum, hex_seams);
        optimize_tangential_smoothing(hex_drum, hex_seams, 3, true);
        TEST_CHECK(check_tri_edge_relations(hex_drum));
    }

    std::cout << "test_topology_and_seams passed." << std::endl;
    return 0;
}
