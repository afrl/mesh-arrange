#include "mesh_arrange/Triangulation.h"

#include "detail/ElementTables.h"
#include "detail/MathHelpers.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <unordered_map>

namespace mesh_arrange
{

void update_edges_and_maps(Triangulation& t)
{
    if (!t.tri_to_marker_map.empty())
    {
        t.n_markers = *std::max_element(t.tri_to_marker_map.begin(), t.tri_to_marker_map.end()) + 1;
    }
    else
    {
        t.n_markers = 0;
    }

    t.edges.clear();
    t.edge_to_tris_map.clear();
    t.tri_to_edges_map.resize(t.triangles.size());

    std::unordered_map<uint64_t, Index> edge_hash_map;
    edge_hash_map.reserve(t.triangles.size() * 2);

    for (std::size_t tri_id = 0; tri_id < t.triangles.size(); ++tri_id)
    {
        const auto& tri_conn = t.triangles[tri_id];
        for (ConnIndex i = 0; i < 3; ++i)
        {
            const Index v1    = tri_conn[detail::kTriEdges[i][0]];
            const Index v2    = tri_conn[detail::kTriEdges[i][1]];
            const uint64_t h  = hash_linear_edge(v1, v2);
            const auto it     = edge_hash_map.find(h);
            if (it == edge_hash_map.end())
            {
                const Index edge_id = static_cast<Index>(t.edges.size());
                t.edges.push_back({ v1, v2 });
                edge_hash_map.emplace(h, edge_id);
                t.edge_to_tris_map.push_back({ { static_cast<Index>(tri_id), i } });
                t.tri_to_edges_map[tri_id][i] = edge_id;
            }
            else
            {
                const Index edge_id = it->second;
                t.edge_to_tris_map[edge_id].push_back({ static_cast<Index>(tri_id), i });
                t.tri_to_edges_map[tri_id][i] = edge_id;
            }
        }
    }
}

EdgeLengthStats query_edge_length_stats(const Triangulation& t)
{
    if (t.edges.empty())
    {
        return { 0.0, 0.0, 0.0, 0.0, 0.0 };
    }

    double mean = 0.0;
    double m2   = 0.0;
    double min  = std::numeric_limits<double>::max();
    double max  = 0.0;
    for (std::size_t k = 0; k < t.edges.size(); ++k)
    {
        const auto& edge       = t.edges[k];
        const double length    = detail::distance(t.vertices[edge[0]], t.vertices[edge[1]]);
        const double prev_mean = mean;
        mean += (length - mean) / static_cast<double>(k + 1);
        m2 += (length - prev_mean) * (length - mean);
        min = std::min(min, length);
        max = std::max(max, length);
    }
    const double variance =
      (t.edges.size() > 1) ? (m2 / static_cast<double>(t.edges.size() - 1)) : 0.0;
    return { min, max, mean, variance, std::sqrt(variance) };
}

} // namespace mesh_arrange
