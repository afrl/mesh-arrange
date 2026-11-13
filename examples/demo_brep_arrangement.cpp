// Demo 1: Exact Surface Arrangement & Non-Manifold B-Rep Extraction
// Intersects a smooth icosphere with a tilted torus (remove_overlaps = false)
// to extract 3D BRepVolumes (A\B, A∩B, B\A), 2D BRepFaces, 1D BRepCurves, and 0D BRepVertices.

#include "geometry_primitives.h"
#include "mesh_arrange/MeshArrange.h"

#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv)
{
    using namespace mesh_arrange;
    using namespace mesh_arrange::examples;

    const DemoCliArgs cli = parse_demo_cli(argc, argv);

    const InputSurface torus  = make_tilted_torus({ -0.14, 0.0, 0.0 }, 0.62, 0.24, 0.46, -0.22, 48, 24);
    const InputSurface sphere = make_icosphere({ +0.46, 0.0, 0.0 }, 0.52, 3);

    ArrangeOptions opts;
    opts.remove_overlaps = false; // Preserve overlap lens A ∩ B as its own BRepVolume
    opts.check_topology  = true;

    const std::vector<InputSurface> inputs = { torus, sphere };
    const ArrangeResult result             = arrange_surfaces(inputs, opts);

    std::cout << "[Demo 1: B-Rep Arrangement]\n"
              << "  BRepVolumes : " << result.brep_volumes.size() << "\n"
              << "  BRepFaces   : " << result.brep_faces.size() << "\n"
              << "  BRepCurves  : " << result.brep_curves.size() << "\n"
              << "  BRepVertices: " << result.brep_vertices.size() << "\n";

    if (!cli.vtu_dir.empty())
    {
        export_demo_meshes_vtu(cli.vtu_dir, "demo1_brep_arrangement", inputs, result);
        std::cout << "  Wrote VTU   : " << cli.vtu_dir << "/demo1_brep_arrangement_{inputs,brep,volumes}.vtu\n";
    }
    if (!cli.json_path.empty())
    {
        export_demo_bundle(cli.json_path, inputs, { "Tilted Torus", "Sphere" }, result);
    }
    return 0;
}
