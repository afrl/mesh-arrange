#pragma once

#include "mesh_arrange/MeshArrange.h"

#include <list>
#include <unordered_map>
#include <vector>

namespace mesh_arrange
{

/// @brief Header-only adapter that populates a Gmsh `GModel` (or API-compatible discrete B-Rep model)
///        from an `ArrangeResult` using Gmsh's native discrete entity classes:
///        `discreteVertex`, `discreteEdge`, `discreteFace`, `discreteRegion`,
///        `MVertex`, `MLine`, and `MTriangle`.
///
/// Example usage inside Gmsh (`#include "GModel.h"`, `"discreteVertex.h"`, `"discreteEdge.h"`,
/// `"discreteFace.h"`, `"discreteRegion.h"`, `"MVertex.h"`, `"MLine.h"`, `"MTriangle.h"`):
/// @code
///   auto res = mesh_arrange::arrange_surfaces(inputs, opts);
///   mesh_arrange::populate_gmsh_discrete_model<
///       GModel, discreteVertex, discreteEdge, discreteFace, discreteRegion,
///       MVertex, MLine, MTriangle>(model, res);
/// @endcode
template<typename GModelT,
         typename DiscreteVertexT,
         typename DiscreteEdgeT,
         typename DiscreteFaceT,
         typename DiscreteRegionT,
         typename MVertexT,
         typename MLineT,
         typename MTriangleT>
struct GmshDiscreteModelMapping
{
    std::vector<DiscreteVertexT*> vertices;
    std::vector<DiscreteEdgeT*> edges;
    std::vector<DiscreteFaceT*> faces;
    std::vector<DiscreteRegionT*> regions;
};

template<typename GModelT,
         typename DiscreteVertexT,
         typename DiscreteEdgeT,
         typename DiscreteFaceT,
         typename DiscreteRegionT,
         typename MVertexT,
         typename MLineT,
         typename MTriangleT>
GmshDiscreteModelMapping<GModelT,
                         DiscreteVertexT,
                         DiscreteEdgeT,
                         DiscreteFaceT,
                         DiscreteRegionT,
                         MVertexT,
                         MLineT,
                         MTriangleT>
populate_gmsh_discrete_model(GModelT* model,
                             const ArrangeResult& result,
                             int base_vertex_tag = 1,
                             int base_edge_tag   = 1,
                             int base_face_tag   = 1,
                             int base_region_tag = 1)
{
    GmshDiscreteModelMapping<GModelT,
                             DiscreteVertexT,
                             DiscreteEdgeT,
                             DiscreteFaceT,
                             DiscreteRegionT,
                             MVertexT,
                             MLineT,
                             MTriangleT>
      mapping;

    if (!model)
    {
        return mapping;
    }

    const auto& coords = result.surface_mesh.vertices;
    std::unordered_map<Index, MVertexT*> global_vid_to_mvertex;
    global_vid_to_mvertex.reserve(coords.size());

    // 1. Create 0D discreteVertex entities for all B-Rep corner / junction vertices
    mapping.vertices.resize(result.brep_vertices.size(), nullptr);
    for (std::size_t i = 0; i < result.brep_vertices.size(); ++i)
    {
        const auto& bv = result.brep_vertices[i];
        auto* gv       = new DiscreteVertexT(model, base_vertex_tag + static_cast<int>(i));
        auto* mv       = new MVertexT(bv.point[0], bv.point[1], bv.point[2], gv);
        gv->mesh_vertices.push_back(mv);
        model->add(gv);
        mapping.vertices[i]                      = gv;
        global_vid_to_mvertex[bv.global_vertex_id] = mv;
    }

    // 2. Create 1D discreteEdge entities for all B-Rep seam curves
    mapping.edges.resize(result.brep_curves.size(), nullptr);
    for (std::size_t i = 0; i < result.brep_curves.size(); ++i)
    {
        const auto& bc = result.brep_curves[i];
        DiscreteVertexT* v_beg =
          (bc.start_brep_vertex != NOT_SET_Index) ? mapping.vertices[bc.start_brep_vertex] : nullptr;
        DiscreteVertexT* v_end =
          (bc.end_brep_vertex != NOT_SET_Index) ? mapping.vertices[bc.end_brep_vertex] : nullptr;

        auto* ge = new DiscreteEdgeT(model, base_edge_tag + static_cast<int>(i), v_beg, v_end);
        if (bc.polyline_vertices.size() >= 2)
        {
            const std::size_t last_interior = bc.polyline_vertices.size() - 1;
            for (std::size_t k = 1; k < last_interior; ++k)
            {
                const Index vid = bc.polyline_vertices[k];
                auto it         = global_vid_to_mvertex.find(vid);
                if (it == global_vid_to_mvertex.end())
                {
                    const auto& p  = coords[vid];
                    auto* mv       = new MVertexT(p[0], p[1], p[2], ge);
                    ge->mesh_vertices.push_back(mv);
                    global_vid_to_mvertex.emplace(vid, mv);
                }
            }
            for (std::size_t k = 0; k + 1 < bc.polyline_vertices.size(); ++k)
            {
                MVertexT* mv0 = global_vid_to_mvertex.at(bc.polyline_vertices[k]);
                MVertexT* mv1 = global_vid_to_mvertex.at(bc.polyline_vertices[k + 1]);
                ge->lines.push_back(new MLineT(mv0, mv1));
            }
        }
        model->add(ge);
        mapping.edges[i] = ge;
    }

    // 3. Create 2D discreteFace entities for all B-Rep region-boundary patches
    mapping.faces.resize(result.brep_faces.size(), nullptr);
    for (std::size_t i = 0; i < result.brep_faces.size(); ++i)
    {
        const auto& bf = result.brep_faces[i];
        auto* gf       = new DiscreteFaceT(model, base_face_tag + static_cast<int>(i));

        std::vector<int> edge_tags;
        edge_tags.reserve(bf.signed_bounding_curves.size());
        for (const auto& [cid, dir] : bf.signed_bounding_curves)
        {
            const int tag = base_edge_tag + static_cast<int>(cid);
            edge_tags.push_back(dir >= 0 ? tag : -tag);
        }
        if (!edge_tags.empty())
        {
            gf->setBoundEdges(edge_tags);
        }

        for (const auto& tri : bf.triangles)
        {
            MVertexT* mvs[3] = { nullptr, nullptr, nullptr };
            for (int c = 0; c < 3; ++c)
            {
                const Index vid = tri[c];
                auto it         = global_vid_to_mvertex.find(vid);
                if (it == global_vid_to_mvertex.end())
                {
                    const auto& p = coords[vid];
                    auto* mv      = new MVertexT(p[0], p[1], p[2], gf);
                    gf->mesh_vertices.push_back(mv);
                    it = global_vid_to_mvertex.emplace(vid, mv).first;
                }
                mvs[c] = it->second;
            }
            gf->triangles.push_back(new MTriangleT(mvs[0], mvs[1], mvs[2]));
        }

        model->add(gf);
        mapping.faces[i] = gf;
    }

    // 4. Create 3D discreteRegion entities for all B-Rep volumes
    mapping.regions.resize(result.brep_volumes.size(), nullptr);
    for (std::size_t i = 0; i < result.brep_volumes.size(); ++i)
    {
        const auto& bvol = result.brep_volumes[i];
        auto* gr         = new DiscreteRegionT(model, base_region_tag + static_cast<int>(i));

        std::vector<int> face_tags;
        std::vector<int> face_dirs;
        face_tags.reserve(bvol.signed_bounding_faces.size());
        face_dirs.reserve(bvol.signed_bounding_faces.size());
        for (const auto& [fid, dir] : bvol.signed_bounding_faces)
        {
            face_tags.push_back(base_face_tag + static_cast<int>(fid));
            face_dirs.push_back(dir);
        }
        if (!face_tags.empty())
        {
            gr->setBoundFaces(face_tags, face_dirs);
        }

        model->add(gr);
        mapping.regions[i] = gr;
    }

    return mapping;
}

} // namespace mesh_arrange
