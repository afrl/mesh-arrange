#include "mesh_arrange/GmshModelBuilder.h"
#include "mesh_arrange/MeshArrange.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
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

// =============================================================================
// Lightweight API-compatible mock of Gmsh's GModel & discrete B-Rep classes
// =============================================================================
struct MockGEntity
{
    int tag = 0;
    virtual ~MockGEntity() = default;
};

struct MockMVertex
{
    double x = 0.0, y = 0.0, z = 0.0;
    MockGEntity* on_what = nullptr;

    MockMVertex(double x_, double y_, double z_, MockGEntity* g)
      : x(x_), y(y_), z(z_), on_what(g)
    {
    }
};

struct MockMLine
{
    MockMVertex* v0 = nullptr;
    MockMVertex* v1 = nullptr;
    MockMLine(MockMVertex* a, MockMVertex* b) : v0(a), v1(b) {}
};

struct MockMTriangle
{
    MockMVertex* v[3] = { nullptr, nullptr, nullptr };
    MockMTriangle(MockMVertex* a, MockMVertex* b, MockMVertex* c) : v{ a, b, c } {}
};

struct MockGModel;

struct MockDiscreteVertex : MockGEntity
{
    std::vector<MockMVertex*> mesh_vertices;
    MockDiscreteVertex(MockGModel*, int t) { tag = t; }
    ~MockDiscreteVertex() override
    {
        for (auto* p : mesh_vertices)
        {
            delete p;
        }
    }
};

struct MockDiscreteEdge : MockGEntity
{
    MockDiscreteVertex* v_beg = nullptr;
    MockDiscreteVertex* v_end = nullptr;
    std::vector<MockMVertex*> mesh_vertices;
    std::vector<MockMLine*> lines;

    MockDiscreteEdge(MockGModel*, int t, MockDiscreteVertex* b, MockDiscreteVertex* e)
      : v_beg(b), v_end(e)
    {
        tag = t;
    }
    ~MockDiscreteEdge() override
    {
        for (auto* p : lines)
        {
            delete p;
        }
        for (auto* p : mesh_vertices)
        {
            delete p;
        }
    }
};

struct MockDiscreteFace : MockGEntity
{
    std::vector<int> bound_edge_tags;
    std::vector<int> bound_edge_dirs;
    std::vector<MockDiscreteEdge*> embedded_edges;
    std::vector<MockDiscreteVertex*> embedded_vertices;
    std::vector<MockMVertex*> mesh_vertices;
    std::vector<MockMTriangle*> triangles;

    MockDiscreteFace(MockGModel*, int t) { tag = t; }
    void setBoundEdges(const std::vector<int>& tags, const std::vector<int>& dirs)
    {
        bound_edge_tags = tags;
        bound_edge_dirs = dirs;
    }
    void addEmbeddedEdge(MockDiscreteEdge* ge) { embedded_edges.push_back(ge); }
    void addEmbeddedVertex(MockDiscreteVertex* gv) { embedded_vertices.push_back(gv); }

    ~MockDiscreteFace() override
    {
        for (auto* p : triangles)
        {
            delete p;
        }
        for (auto* p : mesh_vertices)
        {
            delete p;
        }
    }
};

struct MockDiscreteRegion : MockGEntity
{
    std::vector<int> bound_face_tags;
    std::vector<int> bound_face_dirs;

    MockDiscreteRegion(MockGModel*, int t) { tag = t; }
    void setBoundFaces(const std::vector<int>& tags, const std::vector<int>& dirs)
    {
        bound_face_tags = tags;
        bound_face_dirs = dirs;
    }
};

struct MockGModel
{
    std::vector<std::unique_ptr<MockDiscreteVertex>> vertices;
    std::vector<std::unique_ptr<MockDiscreteEdge>> edges;
    std::vector<std::unique_ptr<MockDiscreteFace>> faces;
    std::vector<std::unique_ptr<MockDiscreteRegion>> regions;

    void add(MockDiscreteVertex* v) { vertices.emplace_back(v); }
    void add(MockDiscreteEdge* e) { edges.emplace_back(e); }
    void add(MockDiscreteFace* f) { faces.emplace_back(f); }
    void add(MockDiscreteRegion* r) { regions.emplace_back(r); }
};

InputSurface make_box_surface(const std::array<double, 3>& neg, const std::array<double, 3>& pos)
{
    InputSurface s;
    s.type     = ObjectType::SOLID;
    s.priority = 0;

    // 8 corners + 6 face centers = 14 vertices, 24 outward-oriented triangles
    const double mx = 0.5 * (neg[0] + pos[0]);
    const double my = 0.5 * (neg[1] + pos[1]);
    const double mz = 0.5 * (neg[2] + pos[2]);

    s.vertices = {
        neg,                          // 0
        { pos[0], neg[1], neg[2] },   // 1
        { pos[0], pos[1], neg[2] },   // 2
        { neg[0], pos[1], neg[2] },   // 3
        { neg[0], neg[1], pos[2] },   // 4
        { pos[0], neg[1], pos[2] },   // 5
        pos,                          // 6
        { neg[0], pos[1], pos[2] },   // 7
        { mx, my, neg[2] },           // 8:  -Z center
        { pos[0], my, mz },           // 9:  +X center
        { mx, my, pos[2] },           // 10: +Z center
        { neg[0], my, mz },           // 11: -X center
        { mx, neg[1], mz },           // 12: -Y center
        { mx, pos[1], mz }            // 13: +Y center
    };

    auto add_quad_fan = [&](Index a, Index b, Index c, Index d, Index center, Marker facet)
    {
        s.triangles.push_back({ a, b, center });
        s.triangles.push_back({ b, c, center });
        s.triangles.push_back({ c, d, center });
        s.triangles.push_back({ d, a, center });
        for (int i = 0; i < 4; ++i)
        {
            s.triangle_facet_ids.push_back(facet);
        }
    };

    // Outward-oriented quads (a, b, c, d) seen from outside:
    add_quad_fan(0, 3, 2, 1, 8, 0);  // -Z
    add_quad_fan(1, 2, 6, 5, 9, 1);  // +X
    add_quad_fan(4, 5, 6, 7, 10, 2); // +Z
    add_quad_fan(0, 4, 7, 3, 11, 3); // -X
    add_quad_fan(0, 1, 5, 4, 12, 4); // -Y
    add_quad_fan(3, 7, 6, 2, 13, 5); // +Y
    return s;
}

/// Custom TetGen-free PLC tetrahedralizer for the 3-slab arrangement of
/// box_a = [-1.0, 0.5] x [-0.5, 0.5]^2 and box_b = [-0.5, 1.0] x [-0.5, 0.5]^2.
/// Places one interior Steiner vertex at the center of each convex slab
/// [-1.0, -0.5], [-0.5, 0.5], [0.5, 1.0] and cones every boundary triangle of each slab!
PlcTetOutput slab_star_plc_tetrahedralizer(const Triangulation& plc_surface)
{
    PlcTetOutput out;
    out.vertices = plc_surface.vertices;

    static constexpr double kSlabBounds[3][2] = {
        { -1.0, -0.5 },
        { -0.5, +0.5 },
        { +0.5, +1.0 }
    };

    Index slab_centers[3];
    for (int s = 0; s < 3; ++s)
    {
        slab_centers[s] = static_cast<Index>(out.vertices.size());
        out.vertices.push_back({ 0.5 * (kSlabBounds[s][0] + kSlabBounds[s][1]), 0.0, 0.0 });
    }

    for (const auto& tri : plc_surface.triangles)
    {
        const double mx = (plc_surface.vertices[tri[0]][0] + plc_surface.vertices[tri[1]][0] +
                           plc_surface.vertices[tri[2]][0]) /
                          3.0;
        for (int s = 0; s < 3; ++s)
        {
            if (mx >= kSlabBounds[s][0] - 1e-9 && mx <= kSlabBounds[s][1] + 1e-9)
            {
                out.tetrahedra.push_back({ tri[0], tri[1], tri[2], slab_centers[s] });
            }
        }
    }

    return out;
}

} // namespace

int main()
{
    const auto box_a = make_box_surface({ -1.0, -0.5, -0.5 }, { +0.5, +0.5, +0.5 });
    const auto box_b = make_box_surface({ -0.5, -0.5, -0.5 }, { +1.0, +0.5, +0.5 });

    // 1. If built without TetGen, verify calling arrange_surfaces without a custom callback
    //    raises a clear diagnostic error.
    if (!has_builtin_tetgen_backend())
    {
        bool caught = false;
        try
        {
            (void)arrange_surfaces({ box_a, box_b });
        }
        catch (const std::runtime_error&)
        {
            caught = true;
        }
        TEST_CHECK(caught);
    }

    // 2. Run full overlap-removing arrangement using ONLY the custom PlcTetrahedralizerFn
    //    and verify LogCallback captures messages cleanly without touching stdout.
    std::vector<std::pair<LogLevel, std::string>> captured_logs;
    ArrangeOptions opts;
    opts.remove_overlaps     = true;
    opts.verbose             = true;
    opts.check_topology      = true;
    opts.plc_tetrahedralizer = slab_star_plc_tetrahedralizer;
    opts.log_callback        = [&](LogLevel lvl, const std::string& msg)
    {
        captured_logs.emplace_back(lvl, msg);
    };

    const auto res = arrange_surfaces({ box_a, box_b }, opts);
    TEST_CHECK(!captured_logs.empty());
    TEST_CHECK(res.constituent_surfaces.size() == 2);
    TEST_CHECK(res.brep_volumes.size() == 2);
    TEST_CHECK(!res.brep_faces.empty());
    TEST_CHECK(!res.brep_curves.empty());
    TEST_CHECK(!res.brep_vertices.empty());

    // Verify thread_local GlobalOptions isolation
    GlobalOptions::get().output_dir = "/tmp/main_thread_dir";
    std::string worker_default_dir;
    std::thread worker([&]() { worker_default_dir = GlobalOptions::get().output_dir; });
    worker.join();
    TEST_CHECK(worker_default_dir == ".");
    GlobalOptions::get().output_dir = ".";

    // 3. Populate mock Gmsh discrete GModel via populate_gmsh_discrete_model
    MockGModel model;
    const auto mapping = populate_gmsh_discrete_model<MockGModel,
                                                      MockDiscreteVertex,
                                                      MockDiscreteEdge,
                                                      MockDiscreteFace,
                                                      MockDiscreteRegion,
                                                      MockMVertex,
                                                      MockMLine,
                                                      MockMTriangle>(&model, res);

    TEST_CHECK(mapping.vertices.size() == res.brep_vertices.size());
    TEST_CHECK(mapping.edges.size() == res.brep_curves.size());
    TEST_CHECK(mapping.faces.size() == res.brep_faces.size());
    TEST_CHECK(mapping.regions.size() == res.brep_volumes.size());
    TEST_CHECK(model.regions.size() == 2);

    // Verify 2-argument setBoundEdges(tags, dirs): positive tags, +-1 signs, and closed 1-cycles (d(dF) = 0)
    for (const auto& gf : model.faces)
    {
        TEST_CHECK(gf->bound_edge_tags.size() == gf->bound_edge_dirs.size());
        std::map<const MockDiscreteVertex*, int> vertex_boundary_flux;
        for (std::size_t i = 0; i < gf->bound_edge_tags.size(); ++i)
        {
            const int etag = gf->bound_edge_tags[i];
            const int edir = gf->bound_edge_dirs[i];
            TEST_CHECK(etag >= 1 && static_cast<std::size_t>(etag) <= mapping.edges.size());
            TEST_CHECK(edir == 1 || edir == -1);
            const auto* ge = mapping.edges.at(static_cast<std::size_t>(etag - 1));
            if (ge->v_beg && ge->v_end)
            {
                vertex_boundary_flux[ge->v_beg] -= edir;
                vertex_boundary_flux[ge->v_end] += edir;
            }
        }
        for (const auto& vkv : vertex_boundary_flux)
        {
            TEST_CHECK(vkv.second == 0);
        }
    }

    // Verify each MockDiscreteRegion's oriented MockDiscreteFace shells enclose exact volume 1.0 each!
    double total_gmsh_vol = 0.0;
    for (const auto& gr : model.regions)
    {
        TEST_CHECK(gr->bound_face_tags.size() == gr->bound_face_dirs.size());
        TEST_CHECK(!gr->bound_face_tags.empty());

        std::map<std::pair<const MockMVertex*, const MockMVertex*>, int> directed_he;
        double region_vol = 0.0;

        for (std::size_t i = 0; i < gr->bound_face_tags.size(); ++i)
        {
            const int ftag = gr->bound_face_tags[i];
            const int fdir = gr->bound_face_dirs[i];
            const auto* gf = mapping.faces.at(static_cast<std::size_t>(ftag - 1));

            for (const auto* tri : gf->triangles)
            {
                const MockMVertex* a = tri->v[0];
                const MockMVertex* b = (fdir >= 0) ? tri->v[1] : tri->v[2];
                const MockMVertex* c = (fdir >= 0) ? tri->v[2] : tri->v[1];

                region_vol += (a->x * (b->y * c->z - b->z * c->y) +
                               a->y * (b->z * c->x - b->x * c->z) +
                               a->z * (b->x * c->y - b->y * c->x)) /
                              6.0;

                directed_he[{ a, b }]++;
                directed_he[{ b, c }]++;
                directed_he[{ c, a }]++;
            }
        }

        // Every directed half-edge occurs once and its opposite half-edge occurs once!
        for (const auto& kv : directed_he)
        {
            TEST_CHECK(kv.second == 1);
            const auto opp = directed_he.find({ kv.first.second, kv.first.first });
            TEST_CHECK(opp != directed_he.end() && opp->second == 1);
        }

        TEST_CHECK(std::abs(region_vol - 1.0) < 0.05);
        total_gmsh_vol += region_vol;
    }
    TEST_CHECK(std::abs(total_gmsh_vol - 2.0) < 1e-9);

    std::cout << "test_custom_plc_and_gmsh_builder passed." << std::endl;
    return 0;
}
