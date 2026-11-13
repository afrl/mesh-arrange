// Demo 4: Simultaneous Multi-Material Overlap Partitioning (K = 4 Tetrahedral Cluster)
// Resolves 4 mutually overlapping equal-priority spheres (15 raw arrangement subregions
// including a K=4 core overlap) into 4 watertight Voronoi-like cells meeting at a quadruple vertex.

#include "geometry_primitives.h"
#include "mesh_arrange/MeshArrange.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv)
{
    using namespace mesh_arrange;
    using namespace mesh_arrange::examples;

    const DemoCliArgs cli = parse_demo_cli(argc, argv);

    const double s          = 0.42 / std::sqrt(3.0);
    const double radius     = 0.64;
    const InputSurface sph0 = make_icosphere({ +s, +s, +s }, radius, 2, ObjectType::SOLID, 1);
    const InputSurface sph1 = make_icosphere({ +s, -s, -s }, radius, 2, ObjectType::SOLID, 1);
    const InputSurface sph2 = make_icosphere({ -s, +s, -s }, radius, 2, ObjectType::SOLID, 1);
    const InputSurface sph3 = make_icosphere({ -s, -s, +s }, radius, 2, ObjectType::SOLID, 1);

    const std::vector<InputSurface> inputs = { sph0, sph1, sph2, sph3 };

    ArrangeOptions opts_before;
    opts_before.remove_overlaps = false;
    const ArrangeResult before  = arrange_surfaces(inputs, opts_before);

    ArrangeOptions opts_after;
    opts_after.remove_overlaps = true;
    opts_after.check_topology  = true;
    ArrangeResult after        = arrange_surfaces(inputs, opts_after);
    smooth_arrangement_surfaces(after, 14);

    std::cout << "[Demo 4: K=4 Tetrahedral Sphere Cluster]\n"
              << "  Raw subregions before partitioning : " << before.brep_volumes.size() << "\n"
              << "  Disjoint volumes after partitioning: " << after.brep_volumes.size() << "\n"
              << "  BRepFaces / BRepCurves / BRepVerts : " << after.brep_faces.size() << " / "
              << after.brep_curves.size() << " / " << after.brep_vertices.size() << "\n";

    if (!cli.vtu_dir.empty())
    {
        export_demo_meshes_vtu(cli.vtu_dir, "demo4_multiway_junction", inputs, after, &before);
        std::cout << "  Wrote VTU                          : " << cli.vtu_dir
                  << "/demo4_multiway_junction_{inputs,before_volumes,brep,volumes}.vtu\n";
    }
    if (!cli.json_path.empty())
    {
        export_demo_bundle(cli.json_path, inputs, { "Sphere 0", "Sphere 1", "Sphere 2", "Sphere 3" }, after, &before);
    }
    return 0;
}
