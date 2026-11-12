#include "mesh_arrange/Tetrahedralization.h"

#include "detail/DebugVTK.h"
#include "detail/ElementTables.h"

#include <algorithm>
#include <stdexcept>
#include <unordered_map>

namespace mesh_arrange
{

void update_edges_and_maps(Tetrahedralization& t)
{
    t.edges.clear();
    t.edge_to_tets_map.clear();
    t.edge_to_tris_map.clear();
    t.tet_to_edges_map.assign(
      t.tetrahedra.size(),
      { NOT_SET_Index, NOT_SET_Index, NOT_SET_Index, NOT_SET_Index, NOT_SET_Index, NOT_SET_Index });
    t.tri_to_edges_map.assign(t.triangles.size(), { NOT_SET_Index, NOT_SET_Index, NOT_SET_Index });

    std::unordered_map<uint64_t, Index> edge_hash_map;
    edge_hash_map.reserve(t.tetrahedra.size() * 2 + t.triangles.size() * 2);

    auto get_or_add_edge = [&](Index v1, Index v2) -> Index
    {
        const uint64_t h = hash_linear_edge(v1, v2);
        const auto it    = edge_hash_map.find(h);
        if (it != edge_hash_map.end())
        {
            return it->second;
        }
        const Index edge_id = static_cast<Index>(t.edges.size());
        t.edges.push_back({ v1, v2 });
        t.edge_to_tets_map.emplace_back();
        t.edge_to_tris_map.emplace_back();
        edge_hash_map.emplace(h, edge_id);
        return edge_id;
    };

    for (std::size_t tet_id = 0; tet_id < t.tetrahedra.size(); ++tet_id)
    {
        const auto& tet_conn = t.tetrahedra[tet_id];
        for (ConnIndex i = 0; i < 6; ++i)
        {
            const Index v1      = tet_conn[detail::kTetEdges[i][0]];
            const Index v2      = tet_conn[detail::kTetEdges[i][1]];
            const Index edge_id = get_or_add_edge(v1, v2);
            t.tet_to_edges_map[tet_id][i] = edge_id;
            t.edge_to_tets_map[edge_id].push_back({ static_cast<Index>(tet_id), i });
        }
    }

    for (std::size_t tri_id = 0; tri_id < t.triangles.size(); ++tri_id)
    {
        const auto& tri_conn = t.triangles[tri_id];
        for (ConnIndex i = 0; i < 3; ++i)
        {
            const Index v1      = tri_conn[detail::kTriEdges[i][0]];
            const Index v2      = tri_conn[detail::kTriEdges[i][1]];
            const Index edge_id = get_or_add_edge(v1, v2);
            t.tri_to_edges_map[tri_id][i] = edge_id;
            t.edge_to_tris_map[edge_id].push_back({ static_cast<Index>(tri_id), i });
        }
    }

    struct SortedTriHash
    {
        std::size_t operator()(const std::array<Index, 3>& k) const noexcept
        {
            std::size_t seed = 0;
            for (int m = 0; m < 3; ++m)
            {
                seed ^= std::hash<Index>{}(k[m]) + 0x9e3779b9ULL + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };

    std::unordered_map<std::array<Index, 3>, Index, SortedTriHash> tri_lookup;
    tri_lookup.reserve(t.triangles.size());
    for (std::size_t tri_id = 0; tri_id < t.triangles.size(); ++tri_id)
    {
        auto key = t.triangles[tri_id];
        std::sort(key.begin(), key.end());
        tri_lookup.emplace(key, static_cast<Index>(tri_id));
    }

    t.tet_to_tri_map.assign(
      t.tetrahedra.size(),
      { NOT_SET_Index, NOT_SET_Index, NOT_SET_Index, NOT_SET_Index });
    t.tri_to_tet_map.assign(t.triangles.size(), { NOT_SET_Index, NOT_SET_Index });

    for (std::size_t tet_id = 0; tet_id < t.tetrahedra.size(); ++tet_id)
    {
        const auto& tet_conn = t.tetrahedra[tet_id];
        for (ConnIndex f = 0; f < 4; ++f)
        {
            std::array<Index, 3> key = { tet_conn[detail::kTetFaces[f][0]],
                                         tet_conn[detail::kTetFaces[f][1]],
                                         tet_conn[detail::kTetFaces[f][2]] };
            std::sort(key.begin(), key.end());
            const auto it = tri_lookup.find(key);
            if (it != tri_lookup.end())
            {
                const Index tri_id          = it->second;
                t.tet_to_tri_map[tet_id][f] = tri_id;
                if (t.tri_to_tet_map[tri_id][0] == NOT_SET_Index)
                {
                    t.tri_to_tet_map[tri_id][0] = static_cast<Index>(tet_id);
                }
                else if (t.tri_to_tet_map[tri_id][1] == NOT_SET_Index)
                {
                    t.tri_to_tet_map[tri_id][1] = static_cast<Index>(tet_id);
                }
                else
                {
                    throw std::runtime_error("Triangle seems to exist in more than 2 tetrahedra.");
                }
            }
        }
    }
}

void write_debug_mesh_state_to_vtk(const Tetrahedralization& t, const std::string& filestem)
{
    detail::write_vtk_triangles(filestem + "_triangles.vtk", t.vertices, t.triangles);
    detail::write_vtk_tetrahedra(filestem + "_tetrahedra.vtk", t.vertices, t.tetrahedra);
}

namespace
{
template<typename T>
bool check_tri_edge_relations_kernel(const T& t)
{
    auto require_check = [](bool cond, const std::string& msg)
    {
        if (!cond)
        {
            throw std::runtime_error(msg);
        }
    };

    require_check(t.tri_to_edges_map.size() == t.triangles.size(),
                  "tri_to_edges_map size (" + std::to_string(t.tri_to_edges_map.size()) +
                    ") does not match triangles size (" + std::to_string(t.triangles.size()) + ")");
    require_check(t.edge_to_tris_map.size() == t.edges.size(),
                  "edge_to_tris_map size (" + std::to_string(t.edge_to_tris_map.size()) +
                    ") does not match edges size (" + std::to_string(t.edges.size()) + ")");

    for (std::size_t tri_id = 0; tri_id < t.triangles.size(); ++tri_id)
    {
        const auto& tri_conn = t.triangles[tri_id];
        for (int i = 0; i < 3; ++i)
        {
            require_check(tri_conn[i] < t.vertices.size(),
                          "Triangle " + std::to_string(tri_id) + " references invalid vertex " +
                            std::to_string(tri_conn[i]));
        }

        for (ConnIndex i = 0; i < 3; ++i)
        {
            const Index edge_id = t.tri_to_edges_map[tri_id][i];
            require_check(edge_id < t.edges.size(),
                          "Triangle " + std::to_string(tri_id) + " references invalid edge " +
                            std::to_string(edge_id));

            for (Index v : t.edges[edge_id])
            {
                require_check(std::find(tri_conn.begin(), tri_conn.end(), v) != tri_conn.end(),
                              "The vertices of the edge do not belong to the triangle.");
            }

            require_check(hash_linear_edge(tri_conn[detail::kTriEdges[i][0]],
                                           tri_conn[detail::kTriEdges[i][1]]) ==
                            hash_linear_edge(t.edges[edge_id][0], t.edges[edge_id][1]),
                          "The edge is in the wrong position of tri_to_edges_map.");

            bool tri_found = false;
            for (const auto& entry : t.edge_to_tris_map[edge_id])
            {
                if (entry.first == tri_id)
                {
                    tri_found = true;
                }
            }
            require_check(tri_found, "The triangle does not appear in edge_to_tris_map.");
        }
    }

    for (std::size_t edge_id = 0; edge_id < t.edges.size(); ++edge_id)
    {
        for (const auto& entry : t.edge_to_tris_map[edge_id])
        {
            require_check(entry.first < t.triangles.size() && entry.second < 3,
                          "edge_to_tris_map contains out-of-range triangle entry.");
            require_check(t.tri_to_edges_map[entry.first][entry.second] == edge_id,
                          "edge_to_tris_map disagrees with tri_to_edges_map.");
        }
    }

    return true;
}
} // namespace

bool check_tri_edge_relations(const Triangulation& t)
{
    return check_tri_edge_relations_kernel(t);
}

bool check_tri_edge_relations(const Tetrahedralization& t)
{
    return check_tri_edge_relations_kernel(t);
}

bool check_tet_edge_relations(const Tetrahedralization& t)
{
    auto require_check = [](bool cond, const std::string& msg)
    {
        if (!cond)
        {
            throw std::runtime_error(msg);
        }
    };

    require_check(t.tet_to_edges_map.size() == t.tetrahedra.size(),
                  "tet_to_edges_map size does not match tetrahedra size.");
    require_check(t.edge_to_tets_map.size() == t.edges.size(),
                  "edge_to_tets_map size does not match edges size.");

    for (std::size_t tet_id = 0; tet_id < t.tetrahedra.size(); ++tet_id)
    {
        const auto& tet_conn = t.tetrahedra[tet_id];
        for (ConnIndex i = 0; i < 6; ++i)
        {
            const Index edge_id = t.tet_to_edges_map[tet_id][i];
            require_check(edge_id < t.edges.size(),
                          "tet_to_edges_map references invalid edge " + std::to_string(edge_id));

            for (Index v : t.edges[edge_id])
            {
                require_check(std::find(tet_conn.begin(), tet_conn.end(), v) != tet_conn.end(),
                              "The vertices of the edge do not belong to the tetrahedron.");
            }

            require_check(hash_linear_edge(tet_conn[detail::kTetEdges[i][0]],
                                           tet_conn[detail::kTetEdges[i][1]]) ==
                            hash_linear_edge(t.edges[edge_id][0], t.edges[edge_id][1]),
                          "The edge is in the wrong position of tet_to_edges_map.");

            bool tet_found = false;
            for (const auto& entry : t.edge_to_tets_map[edge_id])
            {
                if (entry.first == tet_id)
                {
                    tet_found = true;
                }
            }
            require_check(tet_found, "The tetrahedron does not appear in edge_to_tets_map.");
        }
    }

    for (std::size_t edge_id = 0; edge_id < t.edges.size(); ++edge_id)
    {
        for (const auto& entry : t.edge_to_tets_map[edge_id])
        {
            require_check(entry.first < t.tetrahedra.size() && entry.second < 6,
                          "edge_to_tets_map contains out-of-range tetrahedron entry.");
            require_check(t.tet_to_edges_map[entry.first][entry.second] == edge_id,
                          "edge_to_tets_map disagrees with tet_to_edges_map.");
        }
    }

    return true;
}

bool check_tet_tri_relations(const Tetrahedralization& t)
{
    auto require_check = [](bool cond, const std::string& msg)
    {
        if (!cond)
        {
            throw std::runtime_error(msg);
        }
    };

    require_check(t.tet_to_tri_map.size() == t.tetrahedra.size(),
                  "tet_to_tri_map size does not match tetrahedra size.");
    require_check(t.tri_to_tet_map.size() == t.triangles.size(),
                  "tri_to_tet_map size does not match triangles size.");

    for (std::size_t tet_id = 0; tet_id < t.tetrahedra.size(); ++tet_id)
    {
        const auto& tet_conn = t.tetrahedra[tet_id];
        for (ConnIndex i = 0; i < 4; ++i)
        {
            const Index tri_id = t.tet_to_tri_map[tet_id][i];
            if (tri_id == NOT_SET_Index)
            {
                continue;
            }
            require_check(tri_id < t.triangles.size(),
                          "tet_to_tri_map references out-of-range triangle.");

            std::array<Index, 3> sorted_1 = { tet_conn[detail::kTetFaces[i][0]],
                                              tet_conn[detail::kTetFaces[i][1]],
                                              tet_conn[detail::kTetFaces[i][2]] };
            std::sort(sorted_1.begin(), sorted_1.end());
            std::array<Index, 3> sorted_2 = t.triangles[tri_id];
            std::sort(sorted_2.begin(), sorted_2.end());

            require_check(sorted_1 == sorted_2,
                          "The triangle is in the wrong position of tet_to_tri_map.");
            require_check(t.tri_to_tet_map[tri_id][0] == tet_id ||
                            t.tri_to_tet_map[tri_id][1] == tet_id,
                          "The triangle does not appear in tri_to_tet_map.");
        }
    }

    for (std::size_t tri_id = 0; tri_id < t.triangles.size(); ++tri_id)
    {
        for (const Index tet_id : t.tri_to_tet_map[tri_id])
        {
            if (tet_id == NOT_SET_Index)
            {
                continue;
            }
            require_check(tet_id < t.tetrahedra.size(),
                          "tri_to_tet_map references out-of-range tetrahedron.");
            require_check(std::find(t.tet_to_tri_map[tet_id].begin(),
                                    t.tet_to_tri_map[tet_id].end(),
                                    tri_id) != t.tet_to_tri_map[tet_id].end(),
                          "The triangle does not appear in tet_to_tri_map.");
        }
    }

    return true;
}

} // namespace mesh_arrange
