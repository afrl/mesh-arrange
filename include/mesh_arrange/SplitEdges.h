#pragma once

#include "mesh_arrange/Seams.h"
#include "mesh_arrange/Tetrahedralization.h"
#include "mesh_arrange/Triangulation.h"

#include <vector>

namespace mesh_arrange
{

struct UpdatedEdgeInfo
{
    enum Type
    {
        NEW_EDGE,
        SUBDIVIDED_EDGE
    };

    Type update_type   = NEW_EDGE;
    Index orig_edge_id = NOT_SET_Index;
    Index new_edge_id  = NOT_SET_Index;
    Index new_node_id  = NOT_SET_Index;
};

struct UpdatedTriangleInfo
{
    Index orig_tri_id = NOT_SET_Index;
    std::vector<Index> new_tri_ids;
};

struct UpdatedTetrahedronInfo
{
    Index orig_tet_id = NOT_SET_Index;
    std::vector<Index> new_tet_ids;
};

struct UpdateResult
{
    std::vector<UpdatedEdgeInfo> updated_edges;
    std::vector<UpdatedTriangleInfo> updated_triangles;
    std::vector<UpdatedTetrahedronInfo> updated_tetrahedra;
};

UpdateResult split_edges(const std::vector<Index>& edges_to_split, Tetrahedralization& t);

UpdateResult split_edges(const std::vector<Index>& edges_to_split,
                         const std::vector<Index>& existing_vertices_for_split_edges,
                         Tetrahedralization& t);

UpdateResult split_edges(const std::vector<Index>& edges_to_split, Triangulation& t);

UpdateResult split_edges(const std::vector<Index>& edges_to_split,
                         const std::vector<Index>& existing_vertices_for_split_edges,
                         Triangulation& t);

/// @brief Splits interior chord edges whose endpoints both lie on seams.
UpdateResult split_constrained_edges(Triangulation& t, SeamInfo& seams);

/// @brief Splits edges exceeding `max_edge_length` while preserving `SeamInfo`.
void split_large_edges(const double max_edge_length, Triangulation& t, SeamInfo& seams);

} // namespace mesh_arrange
