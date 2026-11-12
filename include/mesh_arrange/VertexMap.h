#pragma once

#include "mesh_arrange/Common.h"
#include "mesh_arrange/Triangulation.h"

#include <cstddef>
#include <vector>

namespace mesh_arrange
{

struct StaticVertexMapAccessor;

/// @brief Compressed-row static vertex-to-triangle incidence table.
struct StaticVertexMap
{
    std::vector<Index> id_list;
    std::vector<std::size_t> offsets;

    StaticVertexMapAccessor get_indices(const Index& vertex_id) const;
};

struct StaticVertexMapAccessor
{
    auto begin() const
    {
        return map.id_list.begin() + static_cast<std::ptrdiff_t>(map.offsets[vertex_id]);
    }
    auto end() const
    {
        return map.id_list.begin() + static_cast<std::ptrdiff_t>(map.offsets[vertex_id + 1]);
    }

    const StaticVertexMap& map;
    const Index vertex_id;
};

StaticVertexMap build_vert_to_tri_map(const Triangulation& t);

} // namespace mesh_arrange
