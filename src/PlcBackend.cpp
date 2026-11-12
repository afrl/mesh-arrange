#include "mesh_arrange/PlcBackend.h"

#include "mesh_arrange/GlobalOptions.h"

#include <stdexcept>

#ifdef MESHARRANGE_ENABLE_TETGEN
#include "tetgen.h"
#endif

namespace mesh_arrange
{

bool has_builtin_tetgen_backend() noexcept
{
#ifdef MESHARRANGE_ENABLE_TETGEN
    return true;
#else
    return false;
#endif
}

PlcTetOutput default_plc_tetrahedralizer(const Triangulation& plc_surface)
{
#ifdef MESHARRANGE_ENABLE_TETGEN
    using namespace mesh_arrange::tetgen_internal;

    PlcTetOutput result;
    if (plc_surface.vertices.empty() || plc_surface.triangles.empty())
    {
        return result;
    }

    tetgenio in;
    in.firstnumber    = 0;
    in.numberofpoints = static_cast<int>(plc_surface.vertices.size());
    in.pointlist      = new REAL[in.numberofpoints * 3];
    for (int i = 0; i < in.numberofpoints; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            in.pointlist[3 * i + j] = plc_surface.vertices[static_cast<std::size_t>(i)][j];
        }
    }

    in.numberoffacets  = static_cast<int>(plc_surface.triangles.size());
    in.facetlist       = new tetgenio::facet[in.numberoffacets];
    in.facetmarkerlist = new int[in.numberoffacets];

    for (int i = 0; i < in.numberoffacets; ++i)
    {
        in.facetmarkerlist[i] = i + 1;
        tetgenio::facet* f    = &in.facetlist[i];
        f->numberofpolygons   = 1;
        f->polygonlist        = new tetgenio::polygon[1];
        f->numberofholes      = 0;
        f->holelist           = nullptr;
        tetgenio::polygon* p  = &f->polygonlist[0];
        p->numberofvertices   = 3;
        p->vertexlist         = new int[3];
        for (int j = 0; j < 3; ++j)
        {
            p->vertexlist[j] =
              static_cast<int>(plc_surface.triangles[static_cast<std::size_t>(i)][j]);
        }
    }

    tetgenbehavior flags;
    flags.plc            = 1;
    flags.object         = tetgenbehavior::POLY;
    flags.nobisect       = 1;
    flags.nomergefacet   = 1;
    flags.nomergevertex  = 1;
    const bool direct_stdout = GlobalOptions::get().verbose && !GlobalOptions::get().log_callback;
    flags.quiet              = direct_stdout ? 0 : 1;
    flags.verbose            = direct_stdout ? 1 : 0;
    flags.quality        = 0;
    flags.smooth_maxiter = 0;
    flags.noimprove      = 1;
    flags.neighout       = 2;
    flags.facesout       = 1;
    flags.edgesout       = 1;
    flags.steinerleft    = 1000;
    flags.epsilon        = 1e-14;

    tetgenio out;
    tetrahedralize(&flags, &in, &out);

    result.vertices.reserve(static_cast<std::size_t>(out.numberofpoints));
    for (int i = 0; i < out.numberofpoints; ++i)
    {
        result.vertices.push_back(
          { out.pointlist[3 * i], out.pointlist[3 * i + 1], out.pointlist[3 * i + 2] });
    }

    result.tetrahedra.reserve(static_cast<std::size_t>(out.numberoftetrahedra));
    for (int i = 0; i < out.numberoftetrahedra; ++i)
    {
        result.tetrahedra.push_back({ static_cast<Index>(out.tetrahedronlist[4 * i]),
                                      static_cast<Index>(out.tetrahedronlist[4 * i + 1]),
                                      static_cast<Index>(out.tetrahedronlist[4 * i + 2]),
                                      static_cast<Index>(out.tetrahedronlist[4 * i + 3]) });
    }

    for (int i = 0; i < out.numberoftrifaces; ++i)
    {
        if (out.trifacemarkerlist[i] > 0)
        {
            result.surface_triangles.push_back({ static_cast<Index>(out.trifacelist[3 * i]),
                                                 static_cast<Index>(out.trifacelist[3 * i + 1]),
                                                 static_cast<Index>(out.trifacelist[3 * i + 2]) });
            result.surface_tri_parent_id.push_back(
              static_cast<Index>(out.trifacemarkerlist[i] - 1));
        }
    }

    return result;
#else
    (void)plc_surface;
    throw std::runtime_error(
      "MeshArrange was built with MESHARRANGE_WITH_TETGEN=OFF; please provide a custom "
      "PlcTetrahedralizerFn in ArrangeOptions::plc_tetrahedralizer (or segment_regions/identify_features).");
#endif
}

} // namespace mesh_arrange
