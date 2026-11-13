// Demo 2: Priority Hierarchy & Boolean Domain Clipping
// Combines high-priority inclusions (SOLID, p=+2), a lower-priority matrix (SOLID, p=0),
// a central cylindrical bore (CLIP_INSIDE, p=+5), and a unit-cell cutter (CLIP_OUTSIDE, p=+5).

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

    const InputSurface rod_x  = make_cylinder({ -0.90, 0.0, +0.20 }, { +0.90, 0.0, +0.20 }, 0.27, 32, 14, 2, ObjectType::SOLID, +2);
    const InputSurface rod_y  = make_cylinder({ 0.0, -0.90, -0.20 }, { 0.0, +0.90, -0.20 }, 0.27, 32, 14, 2, ObjectType::SOLID, +2);
    const InputSurface matrix = make_subdivided_box({ -0.82, -0.82, -0.68 }, { +0.82, +0.82, +0.68 }, 8, 8, 7, ObjectType::SOLID, 0);
    const InputSurface bore   = make_cylinder({ 0.0, 0.0, -0.85 }, { 0.0, 0.0, +0.85 }, 0.14, 28, 12, 2, ObjectType::CLIP_INSIDE, +5);
    const InputSurface domain = make_subdivided_box({ -0.60, -0.60, -0.52 }, { +0.60, +0.60, +0.52 }, 7, 7, 6, ObjectType::CLIP_OUTSIDE, +5);

    ArrangeOptions opts;
    opts.remove_overlaps = true;
    opts.check_topology  = true;

    const std::vector<InputSurface> inputs = { rod_x, rod_y, matrix, bore, domain };
    const ArrangeResult result             = arrange_surfaces(inputs, opts);

    std::cout << "[Demo 2: Priority & Boolean Clipping]\n"
              << "  Surviving Solids : " << result.constituent_surfaces.size() << "\n"
              << "  BRepVolumes      : " << result.brep_volumes.size() << "\n"
              << "  BRepFaces        : " << result.brep_faces.size() << "\n"
              << "  BRepCurves       : " << result.brep_curves.size() << "\n";

    if (!cli.vtu_dir.empty())
    {
        export_demo_meshes_vtu(cli.vtu_dir, "demo2_priority_clipping", inputs, result);
        std::cout << "  Wrote VTU        : " << cli.vtu_dir << "/demo2_priority_clipping_{inputs,brep,volumes}.vtu\n";
    }
    if (!cli.json_path.empty())
    {
        export_demo_bundle(
          cli.json_path,
          inputs,
          { "Inclusion X (p=+2)", "Inclusion Y (p=+2)", "Matrix (p=0)", "Bore (CLIP_INSIDE)", "Unit Cell (CLIP_OUTSIDE)" },
          result);
    }
    return 0;
}
