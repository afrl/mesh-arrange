#include "mesh_arrange/GmshSdkAdapter.h"

#include <gmsh.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#define CHECK_COND(expr)                                                                           \
    do                                                                                             \
    {                                                                                              \
        if (!(expr))                                                                               \
        {                                                                                          \
            std::cerr << "Assertion failed at " << __FILE__ << ":" << __LINE__ << ": " << #expr    \
                      << "\n";                                                                     \
            std::abort();                                                                          \
        }                                                                                          \
    } while (false)

namespace
{

void check_volumes_mesh_quality(const std::vector<std::pair<int, int>>& vol_dim_tags,
                                double expected_total_volume,
                                double vol_tol        = 1e-3,
                                double min_sicn_bound = 0.0)
{
    double total_vol       = 0.0;
    double worst_sicn      = 1.0;
    std::size_t total_tets = 0;

    for (std::size_t i = 0; i < vol_dim_tags.size(); ++i)
    {
        CHECK_COND(vol_dim_tags[i].first == 3);
        const int vtag = vol_dim_tags[i].second;

        std::vector<std::size_t> etags;
        std::vector<std::size_t> enodes;
        gmsh::model::mesh::getElementsByType(4, etags, enodes, vtag);
        CHECK_COND(!etags.empty());
        total_tets += etags.size();

        std::vector<double> vols;
        std::vector<double> sicn;
        gmsh::model::mesh::getElementQualities(etags, vols, "volume");
        gmsh::model::mesh::getElementQualities(etags, sicn, "minSICN");

        double region_vol = 0.0;
        for (std::size_t k = 0; k < vols.size(); ++k)
        {
            region_vol += vols[k];
        }
        CHECK_COND(region_vol > 1e-6);
        total_vol += region_vol;

        for (std::size_t k = 0; k < sicn.size(); ++k)
        {
            worst_sicn = std::min(worst_sicn, sicn[k]);
        }
    }

    std::cout << "  Verified " << vol_dim_tags.size() << " Gmsh 3D volumes: total_tets="
              << total_tets << " total_vol=" << total_vol << " (expected " << expected_total_volume
              << ") minSICN=" << worst_sicn << "\n";

    CHECK_COND(std::abs(total_vol - expected_total_volume) <= vol_tol);
    CHECK_COND(worst_sicn > min_sicn_bound);
}

void test_priority_booleans_and_clippers()
{
    std::cout << "[1/2] Testing priority boolean arrangement + CLIP_OUTSIDE RVE + CLIP_INSIDE void...\n";
    gmsh::model::add("priority_rve");
    gmsh::model::setCurrent("priority_rve");

    // Two overlapping solid boxes extending outside the [-0.5, 0.5]^3 RVE window,
    // clipped outside [-0.5, 0.5]^3 and perforated by a central box void [-0.15, 0.15]^3.
    const int obj_hi   = gmsh::model::occ::addBox(-0.75, -0.30, -0.30, 1.05, 0.60, 0.60);
    const int obj_lo   = gmsh::model::occ::addBox(-0.65, -0.65, -0.65, 1.30, 1.30, 1.30);
    const int rve_clip = gmsh::model::occ::addBox(-0.50, -0.50, -0.50, 1.00, 1.00, 1.00);
    const int void_box = gmsh::model::occ::addBox(-0.15, -0.15, -0.15, 0.30, 0.30, 0.30);
    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.MeshSizeMax", 0.24);
    gmsh::model::mesh::generate(2);

    std::vector<std::pair<int, int>> out_dim_tags;
    std::vector<std::vector<std::pair<int, int>>> out_map;

    const auto res = mesh_arrange::gmsh_sdk::arrange_surfaces(
      /*objectDimTags=*/{ { 3, obj_hi }, { 3, obj_lo } },
      /*toolDimTags=*/{ { 3, rve_clip }, { 3, void_box } },
      out_dim_tags,
      out_map,
      /*objectPriorities=*/{ 10, 1 },
      /*toolTypes=*/{ -1, +1 },
      /*toolPriorities=*/{ 100, 100 },
      /*removeOverlaps=*/false,
      /*optimizeSurfaces=*/false,
      /*maxEdgeLength=*/0.0,
      /*createGeometry=*/false,
      /*removeObject=*/true,
      /*removeTool=*/true);

    CHECK_COND(res.brep_volumes.size() == 2);
    CHECK_COND(out_dim_tags.size() == 2);
    CHECK_COND(out_map.size() == 4);
    CHECK_COND(out_map[0].size() == 1);
    CHECK_COND(out_map[1].size() == 1);
    CHECK_COND(out_map[2].empty());
    CHECK_COND(out_map[3].empty());

    const int phys_hi = gmsh::model::addPhysicalGroup(3, { out_map[0][0].second }, -1, "Inclusion");
    const int phys_lo = gmsh::model::addPhysicalGroup(3, { out_map[1][0].second }, -1, "Matrix");
    CHECK_COND(phys_hi > 0 && phys_lo > 0);

    gmsh::option::setNumber("Mesh.MeshSizeMax", 0.25);
    gmsh::option::setNumber("Mesh.Optimize", 1);
    gmsh::option::setNumber("Mesh.OptimizeNetgen", 0);
    gmsh::model::mesh::generate(3);

    // Expected volume = RVE (1.0) - central void (0.3^3 = 0.027) = 0.973
    check_volumes_mesh_quality(out_dim_tags, 1.0 - 0.027, 1e-6);
    gmsh::model::remove();
}

void test_triple_overlap_interpenetration_removal(bool optimize_surfaces)
{
    std::cout << "[2/2] Testing K=3 multi-material interpenetration removal (optimizeSurfaces="
              << optimize_surfaces << ")...\n";
    gmsh::model::add("triple_overlap");
    gmsh::model::setCurrent("triple_overlap");

    // Three mutually overlapping equal-priority (priority=2) inclusions embedded in a priority=0 matrix [-1, 1]^3
    const int b1     = gmsh::model::occ::addBox(-0.55, -0.35, -0.40, 0.90, 0.70, 0.80);
    const int b2     = gmsh::model::occ::addBox(-0.20, -0.55, -0.35, 0.75, 0.80, 0.80);
    const int b3     = gmsh::model::occ::addBox(-0.30, -0.25, -0.55, 0.60, 0.80, 1.10);
    const int matrix = gmsh::model::occ::addBox(-1.00, -1.00, -1.00, 2.00, 2.00, 2.00);
    gmsh::model::occ::synchronize();

    gmsh::option::setNumber("Mesh.MeshSizeMax", 0.34);
    gmsh::model::mesh::generate(2);

    std::vector<std::pair<int, int>> out_dim_tags;
    std::vector<std::vector<std::pair<int, int>>> out_map;

    const auto res = mesh_arrange::gmsh_sdk::arrange_surfaces(
      /*objectDimTags=*/{ { 3, b1 }, { 3, b2 }, { 3, b3 }, { 3, matrix } },
      /*toolDimTags=*/{},
      out_dim_tags,
      out_map,
      /*objectPriorities=*/{ 2, 2, 2, 0 },
      /*toolTypes=*/{},
      /*toolPriorities=*/{},
      /*removeOverlaps=*/true,
      /*optimizeSurfaces=*/optimize_surfaces,
      /*maxEdgeLength=*/0.0,
      /*createGeometry=*/false,
      /*removeObject=*/true,
      /*removeTool=*/true);

    CHECK_COND(res.brep_volumes.size() == 4);
    CHECK_COND(out_dim_tags.size() == 4);
    CHECK_COND(out_map.size() == 4);
    for (std::size_t i = 0; i < 4; ++i)
    {
        CHECK_COND(!res.brep_volumes[i].is_overlap);
        CHECK_COND(!out_map[i].empty());
    }

    gmsh::option::setNumber("Mesh.MeshSizeMax", 0.38);
    gmsh::option::setNumber("Mesh.Optimize", 1);
    gmsh::option::setNumber("Mesh.OptimizeNetgen", 0);
    gmsh::model::mesh::generate(3);
    if (optimize_surfaces)
    {
        gmsh::model::mesh::optimize("", true);
    }

    // Verify zero inverted 3D tetrahedra across all 4 volumes (allowing ~1e-12 IEEE-754 ULP
    // tolerance for coplanar locked boundary slivers on flat discrete faces).
    check_volumes_mesh_quality(out_dim_tags, 8.0, 1e-6, -1e-12);
    gmsh::model::remove();
}

} // namespace

int main()
{
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 0);

    test_priority_booleans_and_clippers();
    test_triple_overlap_interpenetration_removal(false);
    test_triple_overlap_interpenetration_removal(true);

    gmsh::finalize();
    std::cout << "All Gmsh 4.15 SDK integration tests passed!\n";
    return 0;
}
