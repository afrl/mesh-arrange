// Demo 5: End-to-End Woven Composite RVE -> Gmsh Discrete GModel -> 3D Tetrahedral Mesh
// Combines equal-priority interpenetrating tows (p=+2), a background resin matrix (p=0),
// and an exterior RVE unit-cell clipper (CLIP_OUTSIDE, p=+5), exporting the complete B-Rep
// for direct discrete GModel population and 3D conforming tetrahedral meshing in Gmsh.

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

    const InputSurface warp_tow = make_sinusoidal_tow(true,  0.92, 0.40, 0.20, -0.08, +0.18, 26, 24, ObjectType::SOLID, +2);
    const InputSurface weft_tow = make_sinusoidal_tow(false, 0.92, 0.40, 0.20, +0.08, -0.18, 26, 24, ObjectType::SOLID, +2);
    const InputSurface resin    = make_subdivided_box({ -0.88, -0.88, -0.52 }, { +0.88, +0.88, +0.52 }, 8, 8, 6, ObjectType::SOLID, 0);
    const InputSurface rve_box  = make_subdivided_box({ -0.68, -0.68, -0.38 }, { +0.68, +0.68, +0.38 }, 8, 8, 5, ObjectType::CLIP_OUTSIDE, +5);

    ArrangeOptions opts;
    opts.remove_overlaps              = true;
    opts.split_constrained_seam_edges = true;
    opts.check_topology               = true;

    const std::vector<InputSurface> inputs = { warp_tow, weft_tow, resin, rve_box };
    const ArrangeResult result             = arrange_surfaces(inputs, opts);

    std::cout << "[Demo 5: Woven RVE for Gmsh 3D Meshing]\n"
              << "  BRepVolumes (GRegion): " << result.brep_volumes.size() << "\n"
              << "  BRepFaces   (GFace)  : " << result.brep_faces.size() << "\n"
              << "  BRepCurves  (GEdge)  : " << result.brep_curves.size() << "\n"
              << "  BRepVertices(GVertex): " << result.brep_vertices.size() << "\n";

    if (!cli.vtu_dir.empty())
    {
        export_demo_meshes_vtu(cli.vtu_dir, "demo5_gmsh_rve", inputs, result);
        std::cout << "  Wrote VTU            : " << cli.vtu_dir << "/demo5_gmsh_rve_{inputs,brep,volumes}.vtu\n";
    }
    if (!cli.json_path.empty())
    {
        export_demo_bundle(cli.json_path, inputs, { "Warp Tow", "Weft Tow", "Resin Matrix", "RVE Unit Cell" }, result);
    }
    return 0;
}
