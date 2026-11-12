#include "mesh_arrange/Seams.h"

#include "detail/Containers.h"

#include <algorithm>
#include <map>

namespace mesh_arrange
{

SeamInfo segment_seams(const std::vector<std::array<Index, 2>>& edges,
                       const std::vector<std::array<Index, 3>>& tri_to_edges_map,
                       const Marker& n_markers,
                       const std::vector<Marker>& triangle_to_marker_map)
{
    Marker effective_n_markers = n_markers;
    if (!triangle_to_marker_map.empty())
    {
        const Marker max_marker =
          *std::max_element(triangle_to_marker_map.begin(), triangle_to_marker_map.end());
        effective_n_markers = std::max<Marker>(effective_n_markers, max_marker + 1);
    }

    std::vector<detail::DynamicBitset> edge_id_to_interface_bitset_map(
      edges.size(),
      detail::DynamicBitset(effective_n_markers));

    for (std::size_t tri_id = 0; tri_id < tri_to_edges_map.size(); ++tri_id)
    {
        const Marker marker = triangle_to_marker_map[tri_id];
        for (const Index edge_id : tri_to_edges_map[tri_id])
        {
            edge_id_to_interface_bitset_map[edge_id].set(marker, true);
        }
    }

    std::vector<Index> edges_on_seams;
    std::vector<Marker> edge_id_to_seam_id_map(edges.size(), NOT_SET_SeamID);
    std::map<detail::DynamicBitset, SeamID> bitset_to_seam_id_map;

    for (std::size_t edge_id = 0; edge_id < edge_id_to_interface_bitset_map.size(); ++edge_id)
    {
        const auto& bitset = edge_id_to_interface_bitset_map[edge_id];
        if (bitset.count() > 1)
        {
            edges_on_seams.push_back(static_cast<Index>(edge_id));
            auto it = bitset_to_seam_id_map.find(bitset);
            if (it == bitset_to_seam_id_map.end())
            {
                const auto next_seam_id         = static_cast<SeamID>(bitset_to_seam_id_map.size());
                bitset_to_seam_id_map[bitset]   = next_seam_id;
                edge_id_to_seam_id_map[edge_id] = next_seam_id;
            }
            else
            {
                edge_id_to_seam_id_map[edge_id] = it->second;
            }
        }
    }

    return { static_cast<SeamID>(bitset_to_seam_id_map.size()),
             edges_on_seams,
             edge_id_to_seam_id_map };
}

SeamInfo segment_seams(const Triangulation& t)
{
    return segment_seams(t.edges, t.tri_to_edges_map, t.n_markers, t.tri_to_marker_map);
}

} // namespace mesh_arrange
