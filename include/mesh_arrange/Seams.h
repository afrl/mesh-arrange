#pragma once

#include "mesh_arrange/Common.h"
#include "mesh_arrange/Triangulation.h"

#include <array>
#include <vector>

namespace mesh_arrange
{

/// @brief Topological seam curves where two or more distinct surface markers meet.
struct SeamInfo
{
    Marker n_seams = 0;
    /// Edge indices in `edges` that lie on a seam.
    std::vector<Index> edges_on_seams;
    /// Map from `edge_id` to `SeamID` (`NOT_SET_SeamID` for non-seam edges).
    std::vector<Marker> edge_id_to_seam_id_map;
};

/// @brief Segments edges shared by 2 or more distinct triangle markers into unique `SeamID`s.
SeamInfo segment_seams(const std::vector<std::array<Index, 2>>& edges,
                       const std::vector<std::array<Index, 3>>& tri_to_edges_map,
                       const Marker& n_markers,
                       const std::vector<Marker>& triangle_to_marker_map);

SeamInfo segment_seams(const Triangulation& t);

} // namespace mesh_arrange
