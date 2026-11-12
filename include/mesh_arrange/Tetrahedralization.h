#pragma once

#include "mesh_arrange/Triangulation.h"

#include <array>
#include <string>
#include <utility>
#include <vector>

namespace mesh_arrange
{

/// @brief 3D tetrahedral mesh with tracked boundary/interface triangles and full incidence maps.
struct Tetrahedralization
{
    using LocalizedMap = std::vector<std::vector<std::pair<Index, ConnIndex>>>;

    /// Coordinates of all vertices.
    std::vector<std::array<double, 3>> vertices;
    /// Tetrahedron connectivities (4 vertex indices per tet).
    std::vector<std::array<Index, 4>> tetrahedra;
    /// Tracked surface/interface triangles (3 vertex indices per triangle).
    std::vector<std::array<Index, 3>> triangles;
    /// Unique undirected edges across all tetrahedra and tracked triangles.
    std::vector<std::array<Index, 2>> edges;

    /// Edge-triangle incidence: { triangle_id, local_edge_pos (0..2) }.
    LocalizedMap edge_to_tris_map;
    std::vector<std::array<Index, 3>> tri_to_edges_map;

    /// Edge-tetrahedron incidence: { tet_id, local_edge_pos (0..5) }.
    LocalizedMap edge_to_tets_map;
    std::vector<std::array<Index, 6>> tet_to_edges_map;

    /// Face-tetrahedron incidence for tracked triangles (`NOT_SET_Index` when untracked / exterior).
    std::vector<std::array<Index, 4>> tet_to_tri_map;
    std::vector<std::array<Index, 2>> tri_to_tet_map;
};

/// @brief Rebuilds `edges`, `edge_to_*`, `*_to_edges_map`, `tet_to_tri_map`, and `tri_to_tet_map`.
void update_edges_and_maps(Tetrahedralization& t);

/// @brief Writes `<filestem>_triangles.vtk` and `<filestem>_tetrahedra.vtk` for debugging.
void write_debug_mesh_state_to_vtk(const Tetrahedralization& t, const std::string& filestem);

bool check_tri_edge_relations(const Tetrahedralization& t);
bool check_tet_edge_relations(const Tetrahedralization& t);
bool check_tet_tri_relations(const Tetrahedralization& t);

} // namespace mesh_arrange
