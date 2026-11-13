// Demo 3: Pairwise Harmonic Medial Interpenetration Removal (K = 2)
// Resolves interpenetrating curved composite tows along a smooth harmonic zero-isosurface
// saddle contact sheet while conserving total volume and preserving 2-manifold topology.

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

    const InputSurface warp_tow = make_sinusoidal_tow(true,  0.95, 0.44, 0.22, -0.09, +0.20, 28, 26, ObjectType::SOLID, 1);
    const InputSurface weft_tow = make_sinusoidal_tow(false, 0.95, 0.44, 0.22, +0.09, -0.20, 28, 26, ObjectType::SOLID, 1);

    const std::vector<InputSurface> inputs = { warp_tow, weft_tow };

    ArrangeOptions opts_before;
    opts_before.remove_overlaps = false;
    const ArrangeResult before  = arrange_surfaces(inputs, opts_before);

    ArrangeOptions opts_after;
    opts_after.remove_overlaps = true;
    opts_after.check_topology  = true;
    ArrangeResult after        = arrange_surfaces(inputs, opts_after);
    smooth_arrangement_surfaces(after, 14);

    std::cout << "[Demo 3: Harmonic Medial Interpenetration Removal]\n"
              << "  Before (remove_overlaps=false): " << before.brep_volumes.size() << " BRepVolumes (including overlap lens)\n"
              << "  After  (remove_overlaps=true) : " << after.brep_volumes.size() << " disjoint BRepVolumes, "
              << after.brep_faces.size() << " BRepFaces, " << after.brep_curves.size() << " BRepCurves\n";

    if (!cli.vtu_dir.empty())
    {
        export_demo_meshes_vtu(cli.vtu_dir, "demo3_harmonic_medial", inputs, after, &before);
        std::cout << "  Wrote VTU                     : " << cli.vtu_dir
                  << "/demo3_harmonic_medial_{inputs,before_volumes,brep,volumes}.vtu\n";
    }
    if (!cli.json_path.empty())
    {
        export_demo_bundle(cli.json_path, inputs, { "Warp Tow", "Weft Tow" }, after, &before);
    }
    return 0;
}
