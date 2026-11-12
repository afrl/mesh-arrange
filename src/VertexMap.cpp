#include "mesh_arrange/VertexMap.h"

namespace mesh_arrange
{

StaticVertexMapAccessor StaticVertexMap::get_indices(const Index& vertex_id) const
{
    return { *this, vertex_id };
}

StaticVertexMap build_vert_to_tri_map(const Triangulation& t)
{
    std::vector<std::size_t> n_tris_per_vertex(t.vertices.size(), 0);
    for (const auto& tri : t.triangles)
    {
        for (Index v_id : tri)
        {
            n_tris_per_vertex[v_id]++;
        }
    }

    StaticVertexMap m;
    m.offsets.resize(t.vertices.size() + 1, 0);
    for (std::size_t v_id = 0; v_id < t.vertices.size(); ++v_id)
    {
        m.offsets[v_id + 1] = m.offsets[v_id] + n_tris_per_vertex[v_id];
    }

    m.id_list.resize(m.offsets.back());
    for (std::size_t tri_id = 0; tri_id < t.triangles.size(); ++tri_id)
    {
        for (Index v_id : t.triangles[tri_id])
        {
            const std::size_t position = m.offsets[v_id + 1] - n_tris_per_vertex[v_id];
            m.id_list[position]        = static_cast<Index>(tri_id);
            n_tris_per_vertex[v_id]--;
        }
    }

    return m;
}

} // namespace mesh_arrange
