#pragma once

#include "mesh_arrange/Common.h"

#include <array>
#include <utility>
#include <vector>

namespace mesh_arrange
{

/// @brief Surface triangulation with patch markers and edge-triangle incidence maps.
struct Triangulation
{
    /// Each entry stores pairs of { triangle_index, local_edge_index_in_triangle (0..2) }.
    using LocalizedMap = std::vector<std::vector<std::pair<Index, ConnIndex>>>;

    /// Coordinates of all vertices.
    std::vector<std::array<double, 3>> vertices;
    /// Triangle vertex connectivities (0-based indices into `vertices`).
    std::vector<std::array<Index, 3>> triangles;
    /// Total number of markers, equal to max(tri_to_marker_map) + 1 when non-empty.
    Marker n_markers = 0;
    /// Patch / facet / boundary marker per triangle.
    std::vector<Marker> tri_to_marker_map;

    /// Unique undirected edges (0-based indices into `vertices`).
    std::vector<std::array<Index, 2>> edges;
    /// For each edge: list of { triangle_id, local_edge_position (0..2) }.
    LocalizedMap edge_to_tris_map;
    /// For each triangle: indices of its 3 edges ({v0,v1}, {v1,v2}, {v2,v0}).
    std::vector<std::array<Index, 3>> tri_to_edges_map;
};

/// @brief Rebuilds `n_markers`, `edges`, `edge_to_tris_map`, and `tri_to_edges_map`.
void update_edges_and_maps(Triangulation& t);

struct EdgeLengthStats
{
    double min      = 0.0;
    double max      = 0.0;
    double mean     = 0.0;
    double variance = 0.0;
    double std_dev  = 0.0;
};

/// @brief Computes Welford edge-length statistics over `t.edges`.
EdgeLengthStats query_edge_length_stats(const Triangulation& t);

/// @brief Verifies bidirectional consistency between triangles and edges.
bool check_tri_edge_relations(const Triangulation& t);

} // namespace mesh_arrange
