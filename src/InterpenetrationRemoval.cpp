#include "mesh_arrange/InterpenetrationRemoval.h"

#include "detail/DebugVTK.h"
#include "detail/ElementTables.h"
#include "detail/MathHelpers.h"
#include "mesh_arrange/GlobalOptions.h"
#include "mesh_arrange/Refinement.h"

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>
#include <deque>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace mesh_arrange
{

namespace
{

struct OverlapRegionGroup
{
    RegionID group_id = 0;
    std::vector<RegionID> region_ids;
    std::set<RegionBoundaryID> internal_boundary_ids;
    std::set<RegionBoundaryID> external_boundary_ids;
};

std::vector<OverlapRegionGroup> find_overlap_region_groups(
  const RegionBoundaryInfo& region_boundary_info,
  const std::vector<RegionID>& all_overlap_regions)
{
    const auto& region_id_to_boundary_ids_map =
      region_boundary_info.region_id_to_region_boundary_ids_map;
    const auto& boundary_id_to_region_ids_map =
      region_boundary_info.region_boundary_id_to_region_ids_map;

    const std::set<RegionID> all_overlap_set(all_overlap_regions.begin(),
                                             all_overlap_regions.end());
    std::set<RegionID> remaining_overlap_regions = all_overlap_set;

    RegionID curr_group_id = 0;
    std::vector<OverlapRegionGroup> overlap_region_groups;
    while (!remaining_overlap_regions.empty())
    {
        std::vector<RegionID> overlap_regions_group;
        std::set<RegionBoundaryID> internal_boundary_ids;
        std::set<RegionBoundaryID> external_boundary_ids;

        std::set<RegionID> adjacency_to_search_queue = { *remaining_overlap_regions.begin() };
        while (!adjacency_to_search_queue.empty())
        {
            const RegionID curr_region_id = *adjacency_to_search_queue.begin();
            overlap_regions_group.push_back(curr_region_id);
            remaining_overlap_regions.erase(curr_region_id);
            adjacency_to_search_queue.erase(adjacency_to_search_queue.begin());

            for (const RegionBoundaryID boundary_id : region_id_to_boundary_ids_map[curr_region_id])
            {
                bool boundary_connects_to_overlap = false;
                for (const RegionID adjacent_region_id : boundary_id_to_region_ids_map[boundary_id])
                {
                    if (adjacent_region_id != curr_region_id &&
                        all_overlap_set.count(adjacent_region_id) > 0)
                    {
                        boundary_connects_to_overlap = true;
                        if (remaining_overlap_regions.count(adjacent_region_id) > 0)
                        {
                            adjacency_to_search_queue.insert(adjacent_region_id);
                        }
                    }
                }

                if (boundary_connects_to_overlap)
                {
                    internal_boundary_ids.insert(boundary_id);
                }
                else
                {
                    external_boundary_ids.insert(boundary_id);
                }
            }
        }

        overlap_region_groups.push_back({ curr_group_id++,
                                          std::move(overlap_regions_group),
                                          std::move(internal_boundary_ids),
                                          std::move(external_boundary_ids) });
    }
    return overlap_region_groups;
}

struct SortedTriKeyHash
{
    std::size_t operator()(const std::array<Index, 3>& k) const noexcept
    {
        std::size_t seed = 0;
        for (int i = 0; i < 3; ++i)
        {
            seed ^= std::hash<Index>{}(k[i]) + 0x9e3779b9ULL + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

struct GroupHarmonicField
{
    RegionID group_id = 0;
    std::vector<Index> tets_in_overlap;
    std::vector<Index> overlap_to_global_vertex_map;
    std::unordered_map<Index, Index> global_to_overlap_vertex_map;
    std::vector<RegionID> local_id_to_object_id_map;
    std::unordered_map<RegionID, RegionID> object_id_to_local_id_map;
    std::vector<Eigen::VectorXd> weights;
    std::vector<std::vector<RegionID>> vertex_label_sets;
    std::set<RegionBoundaryID> removable_boundary_ids;
};

GroupHarmonicField compute_group_harmonic_field(const SegmentedTetrahedralization& sf,
                                                const OverlapRegionGroup& overlap_group)
{
    static std::atomic<int> overlap_count{ 0 };
    const int overlap_index = overlap_count.fetch_add(1, std::memory_order_relaxed);

    GroupHarmonicField out;
    out.group_id = overlap_group.group_id;

    const auto& overlap_regions_group = overlap_group.region_ids;
    const auto& internal_boundary_ids = overlap_group.internal_boundary_ids;
    const auto& external_boundary_ids = overlap_group.external_boundary_ids;

    for (const RegionID region_id : overlap_regions_group)
    {
        for (const Index tet_id : sf.region_info.region_to_tets_map[region_id])
        {
            out.tets_in_overlap.push_back(tet_id);
        }
    }

    std::set<Index> vertices_in_overlap_set;
    for (const Index tet_index : out.tets_in_overlap)
    {
        for (const Index n_id : sf.tetrahedralization.tetrahedra[tet_index])
        {
            vertices_in_overlap_set.insert(n_id);
        }
    }

    out.overlap_to_global_vertex_map.reserve(vertices_in_overlap_set.size());
    out.global_to_overlap_vertex_map.reserve(vertices_in_overlap_set.size());
    Index local_v_index = 0;
    for (const Index v_id : vertices_in_overlap_set)
    {
        out.overlap_to_global_vertex_map.push_back(v_id);
        out.global_to_overlap_vertex_map.emplace(v_id, local_v_index++);
    }

    const int n_verts = static_cast<int>(out.overlap_to_global_vertex_map.size());
    const int n_tets  = static_cast<int>(out.tets_in_overlap.size());

    Eigen::MatrixXd V_ov(n_verts, 3);
    for (int i = 0; i < n_verts; ++i)
    {
        const auto& p = sf.tetrahedralization.vertices[out.overlap_to_global_vertex_map[i]];
        V_ov(i, 0)    = p[0];
        V_ov(i, 1)    = p[1];
        V_ov(i, 2)    = p[2];
    }

    Eigen::MatrixXi T_ov(n_tets, 4);
    for (int t = 0; t < n_tets; ++t)
    {
        const auto& conn = sf.tetrahedralization.tetrahedra[out.tets_in_overlap[t]];
        for (int c = 0; c < 4; ++c)
        {
            T_ov(t, c) = static_cast<int>(out.global_to_overlap_vertex_map.at(conn[c]));
        }
    }

    std::set<RegionID> overlapping_objects;
    for (const RegionID region_id : overlap_regions_group)
    {
        for (const RegionID object_id : sf.object_region_info.region_to_objects_map[region_id])
        {
            overlapping_objects.insert(object_id);
        }
    }
    out.local_id_to_object_id_map.assign(overlapping_objects.begin(), overlapping_objects.end());
    for (RegionID k = 0; k < static_cast<RegionID>(out.local_id_to_object_id_map.size()); ++k)
    {
        out.object_id_to_local_id_map.emplace(out.local_id_to_object_id_map[k], k);
    }
    const int K = static_cast<int>(out.local_id_to_object_id_map.size());

    const std::set<RegionID> group_region_set(overlap_regions_group.begin(),
                                              overlap_regions_group.end());
    std::unordered_map<RegionID, std::vector<bool>> region_has_local_obj;
    for (const RegionID region_id : overlap_regions_group)
    {
        std::vector<bool> mask(K, false);
        for (const RegionID obj_id : sf.object_region_info.region_to_objects_map[region_id])
        {
            mask[out.object_id_to_local_id_map.at(obj_id)] = true;
        }
        region_has_local_obj.emplace(region_id, std::move(mask));
    }

    std::vector<std::vector<bool>> vertex_present_objects(n_verts, std::vector<bool>(K, false));
    for (const RegionID region_id : overlap_regions_group)
    {
        const auto& r_mask = region_has_local_obj.at(region_id);
        for (const Index tet_id : sf.region_info.region_to_tets_map[region_id])
        {
            for (const Index g_vid : sf.tetrahedralization.tetrahedra[tet_id])
            {
                const int lv = static_cast<int>(out.global_to_overlap_vertex_map.at(g_vid));
                for (int k = 0; k < K; ++k)
                {
                    if (r_mask[k])
                    {
                        vertex_present_objects[lv][k] = true;
                    }
                }
            }
        }
    }

    std::map<RegionBoundaryID, std::vector<int>> boundary_to_local_owners;
    for (const RegionBoundaryID b_id : internal_boundary_ids)
    {
        boundary_to_local_owners[b_id] = {};
    }
    for (const RegionBoundaryID b_id : external_boundary_ids)
    {
        boundary_to_local_owners[b_id] = {};
    }
    for (int k = 0; k < K; ++k)
    {
        const RegionID obj_id = out.local_id_to_object_id_map[k];
        for (const RegionBoundaryID b_id :
             sf.object_region_boundary_info.object_to_external_boundary_ids_map[obj_id])
        {
            auto it = boundary_to_local_owners.find(b_id);
            if (it != boundary_to_local_owners.end())
            {
                it->second.push_back(k);
            }
        }
    }
    for (const auto& kv : boundary_to_local_owners)
    {
        if (kv.second.size() == 1)
        {
            out.removable_boundary_ids.insert(kv.first);
        }
    }

    std::vector<bool> has_active_cap(n_verts, false);
    std::vector<std::vector<bool>> terminates_at_v(n_verts, std::vector<bool>(K, false));
    std::vector<std::vector<bool>> continues_at_v(n_verts, std::vector<bool>(K, false));

    for (const auto& kv : boundary_to_local_owners)
    {
        const RegionBoundaryID b_id = kv.first;
        const auto& owners_on_b     = kv.second;
        if (owners_on_b.empty())
        {
            continue;
        }

        std::vector<bool> is_owner_on_b(K, false);
        for (int k : owners_on_b)
        {
            is_owner_on_b[k] = true;
        }

        for (const Index tri_id : sf.region_boundaries_info.region_boundary_to_triangles_map[b_id])
        {
            const auto& br = sf.region_info.tri_to_bounded_regions_map[tri_id];
            std::vector<bool> c_in(K, false);
            if (br.first != NOT_SET_RegionID && group_region_set.count(br.first) > 0)
            {
                const auto& m = region_has_local_obj.at(br.first);
                for (int k = 0; k < K; ++k)
                {
                    c_in[k] = c_in[k] || m[k];
                }
            }
            if (br.second != NOT_SET_RegionID && group_region_set.count(br.second) > 0)
            {
                const auto& m = region_has_local_obj.at(br.second);
                for (int k = 0; k < K; ++k)
                {
                    c_in[k] = c_in[k] || m[k];
                }
            }

            std::vector<int> term_list;
            std::vector<int> surv_list;
            for (int k = 0; k < K; ++k)
            {
                if (c_in[k])
                {
                    if (is_owner_on_b[k])
                    {
                        term_list.push_back(k);
                    }
                    else
                    {
                        surv_list.push_back(k);
                    }
                }
            }

            if (!term_list.empty() && !surv_list.empty())
            {
                for (const Index g_vid : sf.tetrahedralization.triangles[tri_id])
                {
                    const auto it_lv = out.global_to_overlap_vertex_map.find(g_vid);
                    if (it_lv == out.global_to_overlap_vertex_map.end())
                    {
                        continue;
                    }
                    const int lv       = static_cast<int>(it_lv->second);
                    has_active_cap[lv] = true;
                    for (int k : term_list)
                    {
                        terminates_at_v[lv][k] = true;
                    }
                    for (int k : surv_list)
                    {
                        continues_at_v[lv][k] = true;
                    }
                }
            }
        }
    }

    struct PairHash
    {
        std::size_t operator()(const std::pair<int, int>& p) const noexcept
        {
            return (static_cast<std::size_t>(static_cast<uint32_t>(p.first)) << 32) |
                   static_cast<std::size_t>(static_cast<uint32_t>(p.second));
        }
    };
    std::unordered_map<std::pair<int, int>, double, PairHash> edge_cot_weight;
    edge_cot_weight.reserve(static_cast<std::size_t>(n_tets) * 4);

    for (int t = 0; t < n_tets; ++t)
    {
        int vid[4] = { T_ov(t, 0), T_ov(t, 1), T_ov(t, 2), T_ov(t, 3) };
        for (int e = 0; e < 6; ++e)
        {
            int u = vid[detail::kTetEdges[e][0]];
            int v = vid[detail::kTetEdges[e][1]];
            if (u > v)
            {
                std::swap(u, v);
            }
            edge_cot_weight.emplace(std::make_pair(u, v), 0.0);
        }

        Eigen::Vector3d q[4] = { V_ov.row(vid[0]).transpose(),
                                 V_ov.row(vid[1]).transpose(),
                                 V_ov.row(vid[2]).transpose(),
                                 V_ov.row(vid[3]).transpose() };
        double six_v         = (q[1] - q[0]).dot((q[2] - q[0]).cross(q[3] - q[0]));
        if (six_v < 0.0)
        {
            std::swap(vid[1], vid[2]);
            std::swap(q[1], q[2]);
            six_v = -six_v;
        }
        if (six_v <= 0.0)
        {
            continue;
        }

        // Outward face double-area vectors M_a opposite vertex q[a] (||M_a|| = 2 * Area_a):
        // -int_{tet} (grad phi_a . grad phi_b) dV = -(M_a . M_b) / (36 * V_tet) = (1/6) l_opp cot(theta_opp)
        const Eigen::Vector3d M[4] = { (q[2] - q[1]).cross(q[3] - q[1]),
                                       (q[3] - q[0]).cross(q[2] - q[0]),
                                       (q[1] - q[0]).cross(q[3] - q[0]),
                                       (q[2] - q[0]).cross(q[1] - q[0]) };
        const double inv_36_vol    = 1.0 / (6.0 * six_v);

        for (int a = 0; a < 4; ++a)
        {
            for (int b = a + 1; b < 4; ++b)
            {
                int u = vid[a];
                int v = vid[b];
                if (u > v)
                {
                    std::swap(u, v);
                }
                edge_cot_weight[std::make_pair(u, v)] += -M[a].dot(M[b]) * inv_36_vol;
            }
        }
    }

    std::vector<double> finite_pos_cots;
    finite_pos_cots.reserve(edge_cot_weight.size());
    for (auto& kv : edge_cot_weight)
    {
        const double val   = kv.second;
        const double w_pos = (std::isfinite(val) && val > 0.0) ? val : 0.0;
        kv.second          = w_pos;
        if (w_pos > 0.0)
        {
            finite_pos_cots.push_back(w_pos);
        }
    }

    double ref_cot_scale = 1.0;
    if (!finite_pos_cots.empty())
    {
        const std::size_t mid = finite_pos_cots.size() / 2;
        std::nth_element(finite_pos_cots.begin(),
                         finite_pos_cots.begin() + static_cast<std::ptrdiff_t>(mid),
                         finite_pos_cots.end());
        ref_cot_scale = std::max(1e-12, finite_pos_cots[mid]);
    }
    const double max_cot_cap = 100.0 * ref_cot_scale;

    std::vector<Eigen::Triplet<double>> a_triplets;
    a_triplets.reserve(edge_cot_weight.size() * 4 + n_verts);
    std::vector<std::vector<int>> adj(n_verts);
    for (const auto& kv : edge_cot_weight)
    {
        const int u          = kv.first.first;
        const int v          = kv.first.second;
        const double w_clamp = std::min(kv.second, max_cot_cap);
        const double w       = 0.95 * w_clamp + 0.05 * ref_cot_scale;
        adj[u].push_back(v);
        adj[v].push_back(u);
        a_triplets.emplace_back(u, v, -w);
        a_triplets.emplace_back(v, u, -w);
        a_triplets.emplace_back(u, u, +w);
        a_triplets.emplace_back(v, v, +w);
    }

    Eigen::SparseMatrix<double> A(n_verts, n_verts);
    A.setFromTriplets(a_triplets.begin(), a_triplets.end());

    std::vector<std::vector<bool>> is_dirichlet(K, std::vector<bool>(n_verts, false));
    std::vector<Eigen::VectorXd> bc_values(K, Eigen::VectorXd::Zero(n_verts));
    std::vector<std::vector<int>> seam_tied_local_objs(n_verts);

    for (int v = 0; v < n_verts; ++v)
    {
        for (int k = 0; k < K; ++k)
        {
            if (!vertex_present_objects[v][k])
            {
                is_dirichlet[k][v] = true;
                bc_values[k](v)    = 0.0;
            }
        }

        if (!has_active_cap[v])
        {
            continue;
        }

        std::vector<int> strict_survivors;
        for (int k = 0; k < K; ++k)
        {
            if (vertex_present_objects[v][k] && continues_at_v[v][k] && !terminates_at_v[v][k])
            {
                strict_survivors.push_back(k);
            }
        }

        if (strict_survivors.size() == 1)
        {
            const int winner = strict_survivors[0];
            for (int k = 0; k < K; ++k)
            {
                is_dirichlet[k][v] = true;
                bc_values[k](v)    = (k == winner) ? 1.0 : 0.0;
            }
        }
        else if (strict_survivors.size() >= 2)
        {
            std::vector<bool> is_surv(K, false);
            for (int k : strict_survivors)
            {
                is_surv[k] = true;
            }
            for (int k = 0; k < K; ++k)
            {
                if (!is_surv[k])
                {
                    is_dirichlet[k][v] = true;
                    bc_values[k](v)    = 0.0;
                }
            }
        }
        else
        {
            std::vector<int> tied_seam;
            for (int k = 0; k < K; ++k)
            {
                if (vertex_present_objects[v][k] && (continues_at_v[v][k] || terminates_at_v[v][k]))
                {
                    tied_seam.push_back(k);
                }
            }
            if (!tied_seam.empty())
            {
                const double eq_val = 1.0 / static_cast<double>(tied_seam.size());
                std::vector<bool> in_seam(K, false);
                for (int k : tied_seam)
                {
                    in_seam[k]         = true;
                    is_dirichlet[k][v] = true;
                    bc_values[k](v)    = eq_val;
                }
                for (int k = 0; k < K; ++k)
                {
                    if (!in_seam[k])
                    {
                        is_dirichlet[k][v] = true;
                        bc_values[k](v)    = 0.0;
                    }
                }
                seam_tied_local_objs[v] = std::move(tied_seam);
            }
        }
    }

    std::vector<Eigen::VectorXd> raw_u(K, Eigen::VectorXd::Zero(n_verts));
    for (int k = 0; k < K; ++k)
    {
        std::vector<bool> visited(n_verts, false);
        for (int start = 0; start < n_verts; ++start)
        {
            if (visited[start] || is_dirichlet[k][start])
            {
                continue;
            }
            std::vector<int> comp;
            bool touches_dirichlet = false;
            std::deque<int> q      = { start };
            visited[start]         = true;
            while (!q.empty())
            {
                const int cur = q.front();
                q.pop_front();
                comp.push_back(cur);
                for (int nb : adj[cur])
                {
                    if (is_dirichlet[k][nb])
                    {
                        touches_dirichlet = true;
                    }
                    else if (!visited[nb])
                    {
                        visited[nb] = true;
                        q.push_back(nb);
                    }
                }
            }
            if (!touches_dirichlet)
            {
                for (int v : comp)
                {
                    int first_present = 0;
                    while (first_present < K && !vertex_present_objects[v][first_present])
                    {
                        ++first_present;
                    }
                    is_dirichlet[k][v] = true;
                    bc_values[k](v)    = (k == first_present) ? 1.0 : 0.0;
                }
            }
        }

        std::vector<int> free_to_vert;
        std::vector<int> vert_to_free(n_verts, -1);
        for (int v = 0; v < n_verts; ++v)
        {
            if (!is_dirichlet[k][v])
            {
                vert_to_free[v] = static_cast<int>(free_to_vert.size());
                free_to_vert.push_back(v);
            }
            else
            {
                raw_u[k](v) = bc_values[k](v);
            }
        }

        const int n_free = static_cast<int>(free_to_vert.size());
        if (n_free == 0)
        {
            continue;
        }

        std::vector<Eigen::Triplet<double>> sub_trips;
        sub_trips.reserve(static_cast<std::size_t>(A.nonZeros()));
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n_free);

        for (int col = 0; col < A.outerSize(); ++col)
        {
            const int fc = vert_to_free[col];
            for (Eigen::SparseMatrix<double>::InnerIterator it(A, col); it; ++it)
            {
                const int row = static_cast<int>(it.row());
                const int fr  = vert_to_free[row];
                if (fr >= 0)
                {
                    if (fc >= 0)
                    {
                        sub_trips.emplace_back(fr, fc, it.value());
                    }
                    else
                    {
                        rhs(fr) -= it.value() * bc_values[k](col);
                    }
                }
            }
        }

        if (rhs.squaredNorm() <= 1e-30)
        {
            for (int f = 0; f < n_free; ++f)
            {
                raw_u[k](free_to_vert[f]) = 0.0;
            }
            continue;
        }

        Eigen::SparseMatrix<double> A_in(n_free, n_free);
        A_in.setFromTriplets(sub_trips.begin(), sub_trips.end());

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
        ldlt.compute(A_in);
        Eigen::VectorXd u_in = Eigen::VectorXd::Zero(n_free);
        bool solved          = false;
        if (ldlt.info() == Eigen::Success)
        {
            u_in   = ldlt.solve(rhs);
            solved = (ldlt.info() == Eigen::Success && u_in.allFinite());
        }
        if (!solved)
        {
            Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
            lu.compute(A_in);
            if (lu.info() == Eigen::Success)
            {
                u_in   = lu.solve(rhs);
                solved = (lu.info() == Eigen::Success && u_in.allFinite());
            }
        }
        if (!solved)
        {
            Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> cg;
            cg.compute(A_in);
            if (cg.info() == Eigen::Success)
            {
                u_in = cg.solve(rhs);
            }
        }

        for (int f = 0; f < n_free; ++f)
        {
            raw_u[k](free_to_vert[f]) = std::clamp(u_in(f), 0.0, 1.0);
        }
    }

    out.weights.assign(K, Eigen::VectorXd::Zero(n_verts));
    for (int v = 0; v < n_verts; ++v)
    {
        double sum_u = 0.0;
        for (int k = 0; k < K; ++k)
        {
            if (vertex_present_objects[v][k])
            {
                sum_u += std::max(0.0, raw_u[k](v));
            }
        }
        if (sum_u > 1e-14)
        {
            for (int k = 0; k < K; ++k)
            {
                if (vertex_present_objects[v][k])
                {
                    out.weights[k](v) = std::max(0.0, raw_u[k](v)) / sum_u;
                }
            }
        }
        else
        {
            for (int k = 0; k < K; ++k)
            {
                if (vertex_present_objects[v][k])
                {
                    out.weights[k](v) = 1.0;
                    break;
                }
            }
        }
    }

    out.vertex_label_sets.resize(n_verts);
    auto pick_single_winner = [&](int v) -> RegionID
    {
        int best_k      = -1;
        double best_val = -1.0;
        for (int k = 0; k < K; ++k)
        {
            if (!vertex_present_objects[v][k])
            {
                continue;
            }
            const double val = out.weights[k](v);
            if (best_k < 0 || val > best_val + 1e-14)
            {
                best_val = val;
                best_k   = k;
            }
        }
        return out.local_id_to_object_id_map[std::max(0, best_k)];
    };

    for (int v = 0; v < n_verts; ++v)
    {
        if (!seam_tied_local_objs[v].empty())
        {
            for (int k : seam_tied_local_objs[v])
            {
                out.vertex_label_sets[v].push_back(out.local_id_to_object_id_map[k]);
            }
            std::sort(out.vertex_label_sets[v].begin(), out.vertex_label_sets[v].end());
        }
        else
        {
            out.vertex_label_sets[v] = { pick_single_winner(v) };
        }
    }

    static constexpr int kLocalTetFaces[4][3] = {
        { 0, 1, 2 }, { 0, 1, 3 }, { 0, 2, 3 }, { 1, 2, 3 }
    };
    auto sets_intersect = [](const std::vector<RegionID>& a, const std::vector<RegionID>& b)
    {
        std::size_t i = 0, j = 0;
        while (i < a.size() && j < b.size())
        {
            if (a[i] == b[j])
            {
                return true;
            }
            if (a[i] < b[j])
            {
                ++i;
            }
            else
            {
                ++j;
            }
        }
        return false;
    };
    auto sets_intersect3 = [](const std::vector<RegionID>& a,
                              const std::vector<RegionID>& b,
                              const std::vector<RegionID>& c)
    {
        for (const RegionID x : a)
        {
            if (std::binary_search(b.begin(), b.end(), x) &&
                std::binary_search(c.begin(), c.end(), x))
            {
                return true;
            }
        }
        return false;
    };

    for (int t = 0; t < n_tets; ++t)
    {
        for (int f = 0; f < 4; ++f)
        {
            const int u    = T_ov(t, kLocalTetFaces[f][0]);
            const int v    = T_ov(t, kLocalTetFaces[f][1]);
            const int w    = T_ov(t, kLocalTetFaces[f][2]);
            const auto& lu = out.vertex_label_sets[u];
            const auto& lv = out.vertex_label_sets[v];
            const auto& lw = out.vertex_label_sets[w];
            if (sets_intersect(lu, lv) && sets_intersect(lv, lw) && sets_intersect(lu, lw) &&
                !sets_intersect3(lu, lv, lw))
            {
                out.vertex_label_sets[u] = { pick_single_winner(u) };
                out.vertex_label_sets[v] = { pick_single_winner(v) };
                out.vertex_label_sets[w] = { pick_single_winner(w) };
            }
        }
    }

    if (GlobalOptions::get().write_debug_meshes)
    {
        std::vector<Point3> v_pts(n_verts);
        for (int i = 0; i < n_verts; ++i)
        {
            v_pts[i] = { V_ov(i, 0), V_ov(i, 1), V_ov(i, 2) };
        }
        std::vector<std::array<Index, 4>> t_els(n_tets);
        for (int t = 0; t < n_tets; ++t)
        {
            t_els[t] = { static_cast<Index>(T_ov(t, 0)),
                         static_cast<Index>(T_ov(t, 1)),
                         static_cast<Index>(T_ov(t, 2)),
                         static_cast<Index>(T_ov(t, 3)) };
        }
        std::vector<detail::VtkScalarFieldF64> pt_fields;
        pt_fields.reserve(K);
        for (int k = 0; k < K; ++k)
        {
            std::vector<double> vals(out.weights[k].data(),
                                     out.weights[k].data() + out.weights[k].size());
            pt_fields.push_back({ "harmonic_weight_" + std::to_string(k), std::move(vals) });
        }
        detail::write_vtk_tetrahedra(GlobalOptions::get().output_dir + "/harmonic_weights_" +
                                       std::to_string(overlap_index) + ".vtk",
                                     v_pts,
                                     t_els,
                                     {},
                                     pt_fields);
    }

    return out;
}

std::array<Index, 4> make_positive_tet(const std::vector<Point3>& vertices,
                                       Index v0,
                                       Index v1,
                                       Index v2,
                                       Index v3)
{
    if (detail::tetrahedron_signed_volume(vertices[v0], vertices[v1], vertices[v2], vertices[v3]) <
        0.0)
    {
        return { v0, v2, v1, v3 };
    }
    return { v0, v1, v2, v3 };
}

void subdivide_face_canonically(
  const std::array<Index, 3>& orig_face,
  const std::unordered_map<std::size_t, Index>& edge_cut_map,
  std::unordered_map<std::array<Index, 3>, Index, SortedTriKeyHash>& face_steiner_map,
  std::vector<Point3>& vertices,
  std::vector<std::array<Index, 3>>& out_sub_tris)
{
    std::array<Index, 3> s = orig_face;
    bool flipped           = false;
    if (s[0] > s[1])
    {
        std::swap(s[0], s[1]);
        flipped = !flipped;
    }
    if (s[1] > s[2])
    {
        std::swap(s[1], s[2]);
        flipped = !flipped;
    }
    if (s[0] > s[1])
    {
        std::swap(s[0], s[1]);
        flipped = !flipped;
    }

    auto find_cut = [&](Index a, Index b) -> Index
    {
        const auto it = edge_cut_map.find(hash_linear_edge(a, b));
        return (it != edge_cut_map.end()) ? it->second : NOT_SET_Index;
    };

    const Index p01 = find_cut(s[0], s[1]);
    const Index p12 = find_cut(s[1], s[2]);
    const Index p02 = find_cut(s[0], s[2]);

    const int n_cuts = (p01 != NOT_SET_Index ? 1 : 0) + (p12 != NOT_SET_Index ? 1 : 0) +
                       (p02 != NOT_SET_Index ? 1 : 0);

    auto emit = [&](Index a, Index b, Index c)
    {
        if (flipped)
        {
            out_sub_tris.push_back({ a, c, b });
        }
        else
        {
            out_sub_tris.push_back({ a, b, c });
        }
    };

    if (n_cuts == 0)
    {
        out_sub_tris.push_back(orig_face);
        return;
    }

    if (n_cuts == 1)
    {
        if (p01 != NOT_SET_Index)
        {
            emit(s[0], p01, s[2]);
            emit(p01, s[1], s[2]);
        }
        else if (p12 != NOT_SET_Index)
        {
            emit(s[1], p12, s[0]);
            emit(p12, s[2], s[0]);
        }
        else
        {
            emit(s[0], s[1], p02);
            emit(p02, s[1], s[2]);
        }
        return;
    }

    if (n_cuts == 2)
    {
        Index c         = NOT_SET_Index;
        Index a         = NOT_SET_Index;
        Index b         = NOT_SET_Index;
        Index pca       = NOT_SET_Index;
        Index pcb       = NOT_SET_Index;
        bool local_flip = false;

        if (p12 == NOT_SET_Index)
        {
            c          = s[0];
            a          = s[1];
            b          = s[2];
            pca        = p01;
            pcb        = p02;
            local_flip = false;
        }
        else if (p02 == NOT_SET_Index)
        {
            c          = s[1];
            a          = s[0];
            b          = s[2];
            pca        = p01;
            pcb        = p12;
            local_flip = true;
        }
        else
        {
            c          = s[2];
            a          = s[0];
            b          = s[1];
            pca        = p02;
            pcb        = p12;
            local_flip = false;
        }

        auto emit_oriented = [&](Index x, Index y, Index z)
        {
            if (local_flip)
            {
                emit(x, z, y);
            }
            else
            {
                emit(x, y, z);
            }
        };

        emit_oriented(c, pca, pcb);
        emit_oriented(a, b, pcb);
        emit_oriented(a, pcb, pca);
        return;
    }

    Index pf        = NOT_SET_Index;
    const auto it_f = face_steiner_map.find(s);
    if (it_f != face_steiner_map.end())
    {
        pf = it_f->second;
    }
    else
    {
        pf             = static_cast<Index>(vertices.size());
        const auto& a0 = vertices[p01];
        const auto& a1 = vertices[p12];
        const auto& a2 = vertices[p02];
        vertices.push_back(
          { (a0[0] + a1[0] + a2[0]) / 3.0,
            (a0[1] + a1[1] + a2[1]) / 3.0,
            (a0[2] + a1[2] + a2[2]) / 3.0 });
        face_steiner_map.emplace(s, pf);
    }

    emit(s[0], p01, pf);
    emit(s[0], pf, p02);
    emit(s[1], p12, pf);
    emit(s[1], pf, p01);
    emit(s[2], p02, pf);
    emit(s[2], pf, p12);
}

RegionID find_max_priority_object_for_region(const std::vector<RegionID>& region_objects,
                                             const std::vector<PriorityType>& object_priorities)
{
    PriorityType best_p = std::numeric_limits<PriorityType>::min();
    RegionID best_obj   = NOT_SET_RegionID;
    for (const RegionID obj_id : region_objects)
    {
        const PriorityType p = object_priorities[obj_id];
        if (best_obj == NOT_SET_RegionID || p > best_p || (p == best_p && obj_id < best_obj))
        {
            best_p   = p;
            best_obj = obj_id;
        }
    }
    return best_obj;
}

void subdivide_and_partition_overlaps(const CombinedSurfaces& inconsistent_combined_surfaces,
                                      const ConsistentSurfaces& consistent_surfaces,
                                      const std::vector<OverlapRegionGroup>& /*overlap_groups*/,
                                      const std::vector<GroupHarmonicField>& group_fields,
                                      SegmentedTetrahedralization& st)
{
    const std::size_t n_orig_verts = st.tetrahedralization.vertices.size();
    const std::size_t n_orig_tets  = st.tetrahedralization.tetrahedra.size();
    const std::size_t n_orig_tris  = st.tetrahedralization.triangles.size();

    std::vector<int> tet_to_group_idx(n_orig_tets, -1);
    std::set<RegionBoundaryID> all_removable_boundaries;
    for (std::size_t g = 0; g < group_fields.size(); ++g)
    {
        for (const Index tet_id : group_fields[g].tets_in_overlap)
        {
            tet_to_group_idx[tet_id] = static_cast<int>(g);
        }
        all_removable_boundaries.insert(group_fields[g].removable_boundary_ids.begin(),
                                        group_fields[g].removable_boundary_ids.end());
    }

    struct EdgeCutAccum
    {
        Index v_min  = 0;
        Index v_max  = 0;
        double t_sum = 0.0;
        int count    = 0;
    };
    std::unordered_map<std::size_t, EdgeCutAccum> cut_accum_map;

    auto label_sets_disjoint = [](const std::vector<RegionID>& a, const std::vector<RegionID>& b)
    {
        std::size_t i = 0, j = 0;
        while (i < a.size() && j < b.size())
        {
            if (a[i] == b[j])
            {
                return false;
            }
            if (a[i] < b[j])
            {
                ++i;
            }
            else
            {
                ++j;
            }
        }
        return true;
    };

    for (const auto& gf : group_fields)
    {
        std::unordered_set<std::size_t> visited_group_edges;
        for (const Index tet_id : gf.tets_in_overlap)
        {
            const auto& conn = st.tetrahedralization.tetrahedra[tet_id];
            for (int e = 0; e < 6; ++e)
            {
                Index ga = conn[detail::kTetEdges[e][0]];
                Index gb = conn[detail::kTetEdges[e][1]];
                if (ga > gb)
                {
                    std::swap(ga, gb);
                }
                const std::size_t h = hash_linear_edge(ga, gb);
                if (!visited_group_edges.insert(h).second)
                {
                    continue;
                }

                const int la      = static_cast<int>(gf.global_to_overlap_vertex_map.at(ga));
                const int lb      = static_cast<int>(gf.global_to_overlap_vertex_map.at(gb));
                const auto& set_a = gf.vertex_label_sets[la];
                const auto& set_b = gf.vertex_label_sets[lb];
                if (!label_sets_disjoint(set_a, set_b))
                {
                    continue;
                }

                auto eval_max_group_w = [&](const std::vector<RegionID>& s, int lv) -> double
                {
                    double m = 0.0;
                    for (const RegionID obj : s)
                    {
                        const RegionID k = gf.object_id_to_local_id_map.at(obj);
                        m                = std::max(m, gf.weights[k](lv));
                    }
                    return m;
                };

                const double da    = eval_max_group_w(set_a, la) - eval_max_group_w(set_b, la);
                const double db    = eval_max_group_w(set_a, lb) - eval_max_group_w(set_b, lb);
                const double denom = da - db;
                const double t =
                  (std::abs(denom) > 1e-14) ? std::clamp(da / denom, 0.05, 0.95) : 0.5;

                auto& acc = cut_accum_map[h];
                acc.v_min = ga;
                acc.v_max = gb;
                acc.t_sum += t;
                acc.count += 1;
            }
        }
    }

    std::vector<EdgeCutAccum> sorted_cuts;
    sorted_cuts.reserve(cut_accum_map.size());
    for (const auto& kv : cut_accum_map)
    {
        sorted_cuts.push_back(kv.second);
    }
    std::sort(sorted_cuts.begin(),
              sorted_cuts.end(),
              [](const EdgeCutAccum& a, const EdgeCutAccum& b)
              { return (a.v_min != b.v_min) ? (a.v_min < b.v_min) : (a.v_max < b.v_max); });

    std::unordered_map<std::size_t, Index> edge_cut_map;
    edge_cut_map.reserve(sorted_cuts.size());
    struct EdgeSteinerOrigin
    {
        Index v0 = 0;
        Index v1 = 0;
        double t = 0.5;
    };
    std::unordered_map<Index, EdgeSteinerOrigin> edge_steiner_origin;
    edge_steiner_origin.reserve(sorted_cuts.size());

    for (const auto& c : sorted_cuts)
    {
        const double t   = c.t_sum / static_cast<double>(c.count);
        const auto& p0   = st.tetrahedralization.vertices[c.v_min];
        const auto& p1   = st.tetrahedralization.vertices[c.v_max];
        const Index vid  = static_cast<Index>(st.tetrahedralization.vertices.size());
        st.tetrahedralization.vertices.push_back({ (1.0 - t) * p0[0] + t * p1[0],
                                                   (1.0 - t) * p0[1] + t * p1[1],
                                                   (1.0 - t) * p0[2] + t * p1[2] });
        edge_cut_map.emplace(hash_linear_edge(c.v_min, c.v_max), vid);
        edge_steiner_origin.emplace(vid, EdgeSteinerOrigin{ c.v_min, c.v_max, t });
    }

    std::unordered_map<std::array<Index, 3>, Index, SortedTriKeyHash> face_steiner_map;

    std::vector<RegionID> orig_region_owner(st.region_info.n_regions, NOT_SET_RegionID);
    for (RegionID r = 0; r < st.region_info.n_regions; ++r)
    {
        orig_region_owner[r] = find_max_priority_object_for_region(
          st.object_region_info.region_to_objects_map[r],
          inconsistent_combined_surfaces.object_priorities);
    }

    auto eval_weight_at_vert = [&](const GroupHarmonicField& gf,
                                   int k,
                                   Index vid,
                                   const std::array<double, 4>& tet_corner_w,
                                   Index tet_steiner_vid) -> double
    {
        if (vid < n_orig_verts)
        {
            const auto it = gf.global_to_overlap_vertex_map.find(vid);
            return (it != gf.global_to_overlap_vertex_map.end()) ? gf.weights[k](it->second) : 0.0;
        }
        if (vid == tet_steiner_vid)
        {
            return 0.25 * (tet_corner_w[0] + tet_corner_w[1] + tet_corner_w[2] + tet_corner_w[3]);
        }
        const auto it_e = edge_steiner_origin.find(vid);
        if (it_e != edge_steiner_origin.end())
        {
            const auto it0  = gf.global_to_overlap_vertex_map.find(it_e->second.v0);
            const auto it1  = gf.global_to_overlap_vertex_map.find(it_e->second.v1);
            const double w0 =
              (it0 != gf.global_to_overlap_vertex_map.end()) ? gf.weights[k](it0->second) : 0.0;
            const double w1 =
              (it1 != gf.global_to_overlap_vertex_map.end()) ? gf.weights[k](it1->second) : 0.0;
            return (1.0 - it_e->second.t) * w0 + it_e->second.t * w1;
        }
        return 0.25 * (tet_corner_w[0] + tet_corner_w[1] + tet_corner_w[2] + tet_corner_w[3]);
    };

    struct SubTetMeta
    {
        bool is_overlap      = false;
        RegionID group_id    = 0;
        RegionID owner       = NOT_SET_RegionID;
        RegionID orig_region = 0;
    };

    std::vector<std::array<Index, 4>> new_tetrahedra;
    std::vector<SubTetMeta> new_tet_meta;
    new_tetrahedra.reserve(n_orig_tets * 2);
    new_tet_meta.reserve(n_orig_tets * 2);

    for (std::size_t tid = 0; tid < n_orig_tets; ++tid)
    {
        const auto& conn      = st.tetrahedralization.tetrahedra[tid];
        const RegionID orig_r = st.region_info.tet_to_region_map[tid];
        const int g_idx       = tet_to_group_idx[tid];
        const bool is_ov      = (g_idx >= 0);

        std::vector<std::pair<int, Index>> tet_cuts;
        for (int e = 0; e < 6; ++e)
        {
            const Index va = conn[detail::kTetEdges[e][0]];
            const Index vb = conn[detail::kTetEdges[e][1]];
            const auto it  = edge_cut_map.find(hash_linear_edge(va, vb));
            if (it != edge_cut_map.end())
            {
                tet_cuts.emplace_back(e, it->second);
            }
        }

        const std::size_t sub_begin = new_tetrahedra.size();
        Index tet_steiner_vid       = NOT_SET_Index;

        if (tet_cuts.empty())
        {
            new_tetrahedra.push_back(conn);
        }
        else if (tet_cuts.size() == 1)
        {
            const int e    = tet_cuts[0].first;
            const Index pe = tet_cuts[0].second;
            const int ia   = detail::kTetEdges[e][0];
            const int ib   = detail::kTetEdges[e][1];
            int ic = -1, id = -1;
            for (int c = 0; c < 4; ++c)
            {
                if (c != ia && c != ib)
                {
                    if (ic < 0)
                    {
                        ic = c;
                    }
                    else
                    {
                        id = c;
                    }
                }
            }
            new_tetrahedra.push_back(
              make_positive_tet(st.tetrahedralization.vertices, conn[ia], pe, conn[ic], conn[id]));
            new_tetrahedra.push_back(
              make_positive_tet(st.tetrahedralization.vertices, pe, conn[ib], conn[ic], conn[id]));
        }
        else
        {
            int single_face_apex = -1;
            for (int apex_c = 0; apex_c < 4; ++apex_c)
            {
                bool any_incident_cut = false;
                for (const auto& ec : tet_cuts)
                {
                    if (detail::kTetEdges[ec.first][0] == apex_c ||
                        detail::kTetEdges[ec.first][1] == apex_c)
                    {
                        any_incident_cut = true;
                        break;
                    }
                }
                if (!any_incident_cut)
                {
                    single_face_apex = apex_c;
                    break;
                }
            }

            if (single_face_apex >= 0)
            {
                int i0 = -1, i1 = -1, i2 = -1;
                for (int c = 0; c < 4; ++c)
                {
                    if (c != single_face_apex)
                    {
                        if (i0 < 0)
                        {
                            i0 = c;
                        }
                        else if (i1 < 0)
                        {
                            i1 = c;
                        }
                        else
                        {
                            i2 = c;
                        }
                    }
                }
                std::vector<std::array<Index, 3>> base_sub_tris;
                subdivide_face_canonically({ conn[i0], conn[i1], conn[i2] },
                                           edge_cut_map,
                                           face_steiner_map,
                                           st.tetrahedralization.vertices,
                                           base_sub_tris);
                const Index v_apex = conn[single_face_apex];
                for (const auto& st_tri : base_sub_tris)
                {
                    new_tetrahedra.push_back(make_positive_tet(st.tetrahedralization.vertices,
                                                               st_tri[0],
                                                               st_tri[1],
                                                               st_tri[2],
                                                               v_apex));
                }
            }
            else
            {
                Point3 p_center = { 0.0, 0.0, 0.0 };
                for (const auto& ec : tet_cuts)
                {
                    const auto& pt = st.tetrahedralization.vertices[ec.second];
                    p_center[0] += pt[0];
                    p_center[1] += pt[1];
                    p_center[2] += pt[2];
                }
                const double inv_n = 1.0 / static_cast<double>(tet_cuts.size());
                p_center[0] *= inv_n;
                p_center[1] *= inv_n;
                p_center[2] *= inv_n;

                tet_steiner_vid = static_cast<Index>(st.tetrahedralization.vertices.size());
                st.tetrahedralization.vertices.push_back(p_center);

                std::vector<std::array<Index, 3>> face_sub_tris;
                for (int f = 0; f < 4; ++f)
                {
                    const auto& fc = detail::kTetFaces[f];
                    face_sub_tris.clear();
                    subdivide_face_canonically({ conn[fc[0]], conn[fc[1]], conn[fc[2]] },
                                               edge_cut_map,
                                               face_steiner_map,
                                               st.tetrahedralization.vertices,
                                               face_sub_tris);
                    for (const auto& st_tri : face_sub_tris)
                    {
                        new_tetrahedra.push_back(make_positive_tet(st.tetrahedralization.vertices,
                                                                   st_tri[0],
                                                                   st_tri[1],
                                                                   st_tri[2],
                                                                   tet_steiner_vid));
                    }
                }
            }
        }

        const std::size_t sub_end = new_tetrahedra.size();
        if (!is_ov)
        {
            const RegionID owner = orig_region_owner[orig_r];
            for (std::size_t s = sub_begin; s < sub_end; ++s)
            {
                new_tet_meta.push_back(SubTetMeta{ false, 0, owner, orig_r });
            }
        }
        else
        {
            const auto& gf             = group_fields[g_idx];
            const auto& allowed_global = st.object_region_info.region_to_objects_map[orig_r];

            for (std::size_t s = sub_begin; s < sub_end; ++s)
            {
                const auto& sc = new_tetrahedra[s];

                std::vector<RegionID> cand = allowed_global;
                std::sort(cand.begin(), cand.end());
                bool has_orig_corner = false;
                for (int c = 0; c < 4; ++c)
                {
                    if (sc[c] < n_orig_verts)
                    {
                        has_orig_corner  = true;
                        const int lv     = static_cast<int>(gf.global_to_overlap_vertex_map.at(sc[c]));
                        const auto& lset = gf.vertex_label_sets[lv];
                        std::vector<RegionID> next_cand;
                        std::set_intersection(cand.begin(),
                                              cand.end(),
                                              lset.begin(),
                                              lset.end(),
                                              std::back_inserter(next_cand));
                        cand = std::move(next_cand);
                    }
                }

                RegionID chosen_owner = NOT_SET_RegionID;
                if (has_orig_corner && cand.size() == 1)
                {
                    chosen_owner = cand[0];
                }
                else
                {
                    const auto& search_objs = (!cand.empty()) ? cand : allowed_global;
                    double best_score       = -1.0;
                    for (const RegionID obj_id : search_objs)
                    {
                        const int k = static_cast<int>(gf.object_id_to_local_id_map.at(obj_id));
                        std::array<double, 4> corner_w;
                        for (int c = 0; c < 4; ++c)
                        {
                            const int lv =
                              static_cast<int>(gf.global_to_overlap_vertex_map.at(conn[c]));
                            corner_w[c] = gf.weights[k](lv);
                        }
                        double score = 0.0;
                        for (int c = 0; c < 4; ++c)
                        {
                            score += eval_weight_at_vert(gf, k, sc[c], corner_w, tet_steiner_vid);
                        }
                        if (chosen_owner == NOT_SET_RegionID || score > best_score + 1e-14 ||
                            (std::abs(score - best_score) <= 1e-14 && obj_id < chosen_owner))
                        {
                            best_score   = score;
                            chosen_owner = obj_id;
                        }
                    }
                }

                new_tet_meta.push_back(SubTetMeta{ true, gf.group_id, chosen_owner, orig_r });
            }
        }
    }

    // 2-manifold edge-star and vertex-star repair
    {
        const std::size_t n_sub_tets = new_tetrahedra.size();
        std::vector<std::array<Index, 4>> sub_tet_adj(
          n_sub_tets,
          { NOT_SET_Index, NOT_SET_Index, NOT_SET_Index, NOT_SET_Index });
        std::vector<double> sub_tet_vol(n_sub_tets, 0.0);

        {
            std::unordered_map<std::array<Index, 3>, std::pair<Index, ConnIndex>, SortedTriKeyHash>
              first_on_face;
            first_on_face.reserve(n_sub_tets * 3);
            for (std::size_t tid = 0; tid < n_sub_tets; ++tid)
            {
                const auto& tc   = new_tetrahedra[tid];
                sub_tet_vol[tid] = std::max(
                  0.0,
                  detail::tetrahedron_signed_volume(st.tetrahedralization.vertices[tc[0]],
                                                    st.tetrahedralization.vertices[tc[1]],
                                                    st.tetrahedralization.vertices[tc[2]],
                                                    st.tetrahedralization.vertices[tc[3]]));

                for (ConnIndex f = 0; f < 4; ++f)
                {
                    const auto& fc           = detail::kTetFaces[f];
                    std::array<Index, 3> key = { tc[fc[0]], tc[fc[1]], tc[fc[2]] };
                    std::sort(key.begin(), key.end());
                    auto [it, ins] =
                      first_on_face.emplace(key, std::make_pair(static_cast<Index>(tid), f));
                    if (!ins)
                    {
                        sub_tet_adj[tid][f] = it->second.first;
                        sub_tet_adj[it->second.first][it->second.second] =
                          static_cast<Index>(tid);
                    }
                }
            }
        }

        std::unordered_map<std::size_t, std::vector<Index>> ov_edge_stars;
        std::unordered_map<Index, std::vector<Index>> ov_vert_stars;
        for (std::size_t tid = 0; tid < n_sub_tets; ++tid)
        {
            if (!new_tet_meta[tid].is_overlap)
            {
                continue;
            }
            const auto& tc = new_tetrahedra[tid];
            for (int e = 0; e < 6; ++e)
            {
                ov_edge_stars[hash_linear_edge(tc[detail::kTetEdges[e][0]],
                                               tc[detail::kTetEdges[e][1]])];
            }
            for (int c = 0; c < 4; ++c)
            {
                ov_vert_stars[tc[c]];
            }
        }
        for (std::size_t tid = 0; tid < n_sub_tets; ++tid)
        {
            const auto& tc = new_tetrahedra[tid];
            for (int e = 0; e < 6; ++e)
            {
                auto it = ov_edge_stars.find(
                  hash_linear_edge(tc[detail::kTetEdges[e][0]], tc[detail::kTetEdges[e][1]]));
                if (it != ov_edge_stars.end())
                {
                    it->second.push_back(static_cast<Index>(tid));
                }
            }
            for (int c = 0; c < 4; ++c)
            {
                auto it = ov_vert_stars.find(tc[c]);
                if (it != ov_vert_stars.end())
                {
                    it->second.push_back(static_cast<Index>(tid));
                }
            }
        }

        auto repair_star = [&](const std::vector<Index>& star_tets) -> bool
        {
            if (star_tets.size() < 2)
            {
                return false;
            }
            std::unordered_set<Index> star_set(star_tets.begin(), star_tets.end());
            std::set<RegionID> owners_present;
            for (const Index t : star_tets)
            {
                if (new_tet_meta[t].owner != NOT_SET_RegionID)
                {
                    owners_present.insert(new_tet_meta[t].owner);
                }
            }
            if (owners_present.size() <= 1)
            {
                return false;
            }

            bool modified = false;
            for (const RegionID owner_o : owners_present)
            {
                std::vector<std::vector<Index>> comps;
                std::unordered_set<Index> visited;
                for (const Index seed : star_tets)
                {
                    if (new_tet_meta[seed].owner != owner_o || visited.count(seed) > 0)
                    {
                        continue;
                    }
                    std::vector<Index> comp;
                    std::deque<Index> q = { seed };
                    visited.insert(seed);
                    while (!q.empty())
                    {
                        const Index cur = q.front();
                        q.pop_front();
                        comp.push_back(cur);
                        for (int f = 0; f < 4; ++f)
                        {
                            const Index nb = sub_tet_adj[cur][f];
                            if (nb != NOT_SET_Index && star_set.count(nb) > 0 &&
                                new_tet_meta[nb].owner == owner_o && visited.count(nb) == 0)
                            {
                                visited.insert(nb);
                                q.push_back(nb);
                            }
                        }
                    }
                    comps.push_back(std::move(comp));
                }

                if (comps.size() <= 1)
                {
                    continue;
                }

                int keep_idx    = -1;
                bool keep_fixed = false;
                double keep_vol = -1.0;
                for (std::size_t ci = 0; ci < comps.size(); ++ci)
                {
                    bool has_fixed = false;
                    double c_vol   = 0.0;
                    for (const Index t : comps[ci])
                    {
                        if (!new_tet_meta[t].is_overlap)
                        {
                            has_fixed = true;
                        }
                        c_vol += sub_tet_vol[t];
                    }
                    if (keep_idx < 0 || (has_fixed && !keep_fixed) ||
                        (has_fixed == keep_fixed && c_vol > keep_vol))
                    {
                        keep_idx   = static_cast<int>(ci);
                        keep_fixed = has_fixed;
                        keep_vol   = c_vol;
                    }
                }

                for (std::size_t ci = 0; ci < comps.size(); ++ci)
                {
                    if (static_cast<int>(ci) == keep_idx)
                    {
                        continue;
                    }
                    bool all_ov = true;
                    std::vector<RegionID> common_allowed;
                    for (std::size_t idx = 0; idx < comps[ci].size(); ++idx)
                    {
                        const Index t = comps[ci][idx];
                        if (!new_tet_meta[t].is_overlap)
                        {
                            all_ov = false;
                            break;
                        }
                        std::vector<RegionID> r_objs =
                          st.object_region_info.region_to_objects_map[new_tet_meta[t].orig_region];
                        std::sort(r_objs.begin(), r_objs.end());
                        if (idx == 0)
                        {
                            common_allowed = std::move(r_objs);
                        }
                        else
                        {
                            std::vector<RegionID> inter;
                            std::set_intersection(common_allowed.begin(),
                                                  common_allowed.end(),
                                                  r_objs.begin(),
                                                  r_objs.end(),
                                                  std::back_inserter(inter));
                            common_allowed = std::move(inter);
                        }
                    }
                    if (!all_ov || common_allowed.empty())
                    {
                        continue;
                    }

                    std::map<RegionID, int> adj_counts;
                    for (const Index t : comps[ci])
                    {
                        for (int f = 0; f < 4; ++f)
                        {
                            const Index nb = sub_tet_adj[t][f];
                            if (nb != NOT_SET_Index && star_set.count(nb) > 0)
                            {
                                const RegionID nb_o = new_tet_meta[nb].owner;
                                if (nb_o != owner_o && nb_o != NOT_SET_RegionID &&
                                    std::binary_search(common_allowed.begin(),
                                                       common_allowed.end(),
                                                       nb_o))
                                {
                                    adj_counts[nb_o]++;
                                }
                            }
                        }
                    }

                    RegionID best_alt = NOT_SET_RegionID;
                    int best_cnt      = -1;
                    for (const auto& kv : adj_counts)
                    {
                        if (kv.second > best_cnt)
                        {
                            best_cnt = kv.second;
                            best_alt = kv.first;
                        }
                    }
                    if (best_alt != NOT_SET_RegionID)
                    {
                        for (const Index t : comps[ci])
                        {
                            new_tet_meta[t].owner = best_alt;
                        }
                        modified = true;
                    }
                }
            }
            return modified;
        };

        for (int pass = 0; pass < 24; ++pass)
        {
            bool any_change = false;
            for (const auto& kv : ov_edge_stars)
            {
                if (repair_star(kv.second))
                {
                    any_change = true;
                }
            }
            for (const auto& kv : ov_vert_stars)
            {
                if (repair_star(kv.second))
                {
                    any_change = true;
                }
            }
            if (!any_change)
            {
                break;
            }
        }
    }

    struct OrigSubTriInfo
    {
        std::array<Index, 3> oriented_conn;
        FacetID facet_id             = 0;
        RegionBoundaryID boundary_id = 0;
    };
    std::unordered_map<std::array<Index, 3>, OrigSubTriInfo, SortedTriKeyHash> orig_sub_tri_lookup;
    orig_sub_tri_lookup.reserve(n_orig_tris * 2);

    std::vector<std::array<Index, 3>> sub_tris_buf;
    for (std::size_t tr = 0; tr < n_orig_tris; ++tr)
    {
        const auto& tri_conn       = st.tetrahedralization.triangles[tr];
        const FacetID fid          = st.region_info.tri_to_facet_map[tr];
        const RegionBoundaryID bid = st.region_boundaries_info.triangle_to_boundary_id_map[tr];

        sub_tris_buf.clear();
        subdivide_face_canonically(tri_conn,
                                   edge_cut_map,
                                   face_steiner_map,
                                   st.tetrahedralization.vertices,
                                   sub_tris_buf);
        for (const auto& st_conn : sub_tris_buf)
        {
            auto key = st_conn;
            std::sort(key.begin(), key.end());
            orig_sub_tri_lookup.emplace(key, OrigSubTriInfo{ st_conn, fid, bid });
        }
    }

    struct TetFaceIncidence
    {
        Index tet0                = NOT_SET_Index;
        Index tet1                = NOT_SET_Index;
        std::array<Index, 3> conn = { 0, 0, 0 };
    };
    std::unordered_map<std::array<Index, 3>, TetFaceIncidence, SortedTriKeyHash> all_tet_faces;
    all_tet_faces.reserve(new_tetrahedra.size() * 3);
    std::vector<std::array<Index, 3>> face_keys_in_order;
    face_keys_in_order.reserve(new_tetrahedra.size() * 3);

    for (std::size_t tid = 0; tid < new_tetrahedra.size(); ++tid)
    {
        const auto& tc = new_tetrahedra[tid];
        for (int f = 0; f < 4; ++f)
        {
            const auto& fc             = detail::kTetFaces[f];
            std::array<Index, 3> f_conn = { tc[fc[0]], tc[fc[1]], tc[fc[2]] };
            auto key                   = f_conn;
            std::sort(key.begin(), key.end());
            auto [it, inserted] = all_tet_faces.emplace(
              key,
              TetFaceIncidence{ static_cast<Index>(tid), NOT_SET_Index, f_conn });
            if (inserted)
            {
                face_keys_in_order.push_back(key);
            }
            else
            {
                it->second.tet1 = static_cast<Index>(tid);
            }
        }
    }

    const FacetID n_consistent_facets = static_cast<FacetID>(
      consistent_surfaces.consistent_facet_to_inconsistent_facets_map.size());
    FacetID max_existing_facet = 0;
    if (!st.region_info.tri_to_facet_map.empty())
    {
        max_existing_facet =
          *std::max_element(st.region_info.tri_to_facet_map.begin(),
                            st.region_info.tri_to_facet_map.end()) +
          1;
    }
    const FacetID medial_facet_base = std::max(n_consistent_facets, max_existing_facet);
    FacetID next_medial_facet       = medial_facet_base;
    std::map<std::tuple<RegionID, RegionID, RegionID>, FacetID> medial_pair_to_facet;

    std::vector<std::array<Index, 3>> final_triangles;
    std::vector<FacetID> final_tri_to_facet_map;
    final_triangles.reserve(orig_sub_tri_lookup.size() + 256);
    final_tri_to_facet_map.reserve(orig_sub_tri_lookup.size() + 256);

    for (const auto& key : face_keys_in_order)
    {
        const auto& inc     = all_tet_faces.at(key);
        const auto it_orig  = orig_sub_tri_lookup.find(key);
        const bool has_orig = (it_orig != orig_sub_tri_lookup.end());

        if (inc.tet1 == NOT_SET_Index)
        {
            if (has_orig)
            {
                final_triangles.push_back(it_orig->second.oriented_conn);
                final_tri_to_facet_map.push_back(it_orig->second.facet_id);
            }
            else
            {
                final_triangles.push_back(inc.conn);
                final_tri_to_facet_map.push_back(0);
            }
            continue;
        }

        const auto& m0 = new_tet_meta[inc.tet0];
        const auto& m1 = new_tet_meta[inc.tet1];

        if (!m0.is_overlap && !m1.is_overlap)
        {
            if (has_orig)
            {
                final_triangles.push_back(it_orig->second.oriented_conn);
                final_tri_to_facet_map.push_back(it_orig->second.facet_id);
            }
            continue;
        }

        if (m0.owner == m1.owner)
        {
            continue;
        }

        if (has_orig && all_removable_boundaries.count(it_orig->second.boundary_id) == 0 &&
            (!m0.is_overlap || !m1.is_overlap))
        {
            final_triangles.push_back(it_orig->second.oriented_conn);
            final_tri_to_facet_map.push_back(it_orig->second.facet_id);
        }
        else
        {
            const RegionID grp   = m0.is_overlap ? m0.group_id : m1.group_id;
            const RegionID obj_a = std::min(m0.owner, m1.owner);
            const RegionID obj_b = std::max(m0.owner, m1.owner);
            const auto mkey      = std::make_tuple(grp, obj_a, obj_b);
            auto [mit, mins]     = medial_pair_to_facet.emplace(mkey, next_medial_facet);
            if (mins)
            {
                ++next_medial_facet;
            }
            final_triangles.push_back(has_orig ? it_orig->second.oriented_conn : inc.conn);
            final_tri_to_facet_map.push_back(mit->second);
        }
    }

    st.tetrahedralization.tetrahedra = std::move(new_tetrahedra);
    st.tetrahedralization.triangles  = std::move(final_triangles);
    update_edges_and_maps(st.tetrahedralization);

    const std::size_t n_final_tets = st.tetrahedralization.tetrahedra.size();
    std::vector<std::array<Index, 4>> tet_adj(
      n_final_tets,
      { NOT_SET_Index, NOT_SET_Index, NOT_SET_Index, NOT_SET_Index });

    {
        std::unordered_map<std::array<Index, 3>, std::pair<Index, ConnIndex>, SortedTriKeyHash>
          first_tet_on_face;
        first_tet_on_face.reserve(n_final_tets * 3);
        for (std::size_t tid = 0; tid < n_final_tets; ++tid)
        {
            const auto& tc = st.tetrahedralization.tetrahedra[tid];
            for (ConnIndex f = 0; f < 4; ++f)
            {
                const auto& fc           = detail::kTetFaces[f];
                std::array<Index, 3> key = { tc[fc[0]], tc[fc[1]], tc[fc[2]] };
                std::sort(key.begin(), key.end());
                auto [it, inserted] =
                  first_tet_on_face.emplace(key, std::make_pair(static_cast<Index>(tid), f));
                if (!inserted)
                {
                    const Index other_t       = it->second.first;
                    const ConnIndex other_f   = it->second.second;
                    tet_adj[tid][f]           = other_t;
                    tet_adj[other_t][other_f] = static_cast<Index>(tid);
                }
            }
        }
    }

    std::vector<RegionID> new_tet_to_region(n_final_tets, NOT_SET_RegionID);
    std::vector<std::vector<Index>> new_region_to_tets;
    std::vector<RegionID> new_region_owner;

    for (std::size_t seed = 0; seed < n_final_tets; ++seed)
    {
        if (new_tet_to_region[seed] != NOT_SET_RegionID)
        {
            continue;
        }
        const RegionID r_id       = static_cast<RegionID>(new_region_to_tets.size());
        const RegionID seed_owner = new_tet_meta[seed].owner;
        new_region_to_tets.emplace_back();
        new_region_owner.push_back(seed_owner);

        std::deque<Index> bfs   = { static_cast<Index>(seed) };
        new_tet_to_region[seed] = r_id;
        while (!bfs.empty())
        {
            const Index cur = bfs.front();
            bfs.pop_front();
            new_region_to_tets[r_id].push_back(cur);

            for (ConnIndex f = 0; f < 4; ++f)
            {
                if (st.tetrahedralization.tet_to_tri_map[cur][f] == NOT_SET_Index)
                {
                    const Index nb = tet_adj[cur][f];
                    if (nb != NOT_SET_Index && new_tet_to_region[nb] == NOT_SET_RegionID)
                    {
                        assert(new_tet_meta[nb].owner == seed_owner);
                        new_tet_to_region[nb] = r_id;
                        bfs.push_back(nb);
                    }
                }
            }
        }
        std::sort(new_region_to_tets[r_id].begin(), new_region_to_tets[r_id].end());
    }

    const RegionID n_new_regions = static_cast<RegionID>(new_region_to_tets.size());

    st.region_info.n_regions          = n_new_regions;
    st.region_info.tet_to_region_map  = std::move(new_tet_to_region);
    st.region_info.region_to_tets_map = std::move(new_region_to_tets);
    st.region_info.tri_to_facet_map   = std::move(final_tri_to_facet_map);

    const std::size_t n_final_tris = st.tetrahedralization.triangles.size();
    st.region_info.tri_to_bounded_regions_map.assign(
      n_final_tris,
      { NOT_SET_RegionID, NOT_SET_RegionID });
    st.region_info.region_to_bounding_triangles_map.assign(n_new_regions, {});

    for (std::size_t tr = 0; tr < n_final_tris; ++tr)
    {
        const auto& tt = st.tetrahedralization.tri_to_tet_map[tr];
        auto& br       = st.region_info.tri_to_bounded_regions_map[tr];
        if (tt[0] != NOT_SET_Index)
        {
            br.first = st.region_info.tet_to_region_map[tt[0]];
            st.region_info.region_to_bounding_triangles_map[br.first].push_back(
              static_cast<Index>(tr));
        }
        if (tt[1] != NOT_SET_Index)
        {
            br.second = st.region_info.tet_to_region_map[tt[1]];
            if (br.second != br.first)
            {
                st.region_info.region_to_bounding_triangles_map[br.second].push_back(
                  static_cast<Index>(tr));
            }
        }
        if (br.first > br.second)
        {
            std::swap(br.first, br.second);
        }
    }

    const std::size_t n_objects = st.object_region_info.object_to_regions_map.size();
    st.object_region_info.regions_of_overlap.clear();
    st.object_region_info.region_types.assign(n_new_regions, RegionType::KEPT);
    st.object_region_info.region_to_objects_map.assign(n_new_regions, {});
    st.object_region_info.object_to_regions_map.assign(n_objects, {});

    for (RegionID r = 0; r < n_new_regions; ++r)
    {
        const RegionID owner                           = new_region_owner[r];
        st.object_region_info.region_to_objects_map[r] = { owner };
        st.object_region_info.object_to_regions_map[owner].push_back(r);
    }

    std::vector<std::set<FacetID>> object_pred_facet_sets(n_objects);
    for (std::size_t obj = 0; obj < n_objects; ++obj)
    {
        object_pred_facet_sets[obj].insert(
          inconsistent_combined_surfaces.object_to_facets_map[obj].begin(),
          inconsistent_combined_surfaces.object_to_facets_map[obj].end());
    }

    for (std::size_t tr = 0; tr < n_final_tris; ++tr)
    {
        FacetID& fid = st.region_info.tri_to_facet_map[tr];
        if (fid >= n_consistent_facets)
        {
            continue;
        }
        const auto& pred_facets =
          consistent_surfaces.consistent_facet_to_inconsistent_facets_map[fid];
        if (pred_facets.size() <= 1)
        {
            continue;
        }

        const auto& br = st.region_info.tri_to_bounded_regions_map[tr];
        std::set<FacetID> matched_preds;
        auto collect_for_region = [&](RegionID r)
        {
            if (r == NOT_SET_RegionID)
            {
                return;
            }
            const RegionID owner = new_region_owner[r];
            for (const FacetID pf : pred_facets)
            {
                if (object_pred_facet_sets[owner].count(pf) > 0)
                {
                    matched_preds.insert(pf);
                }
            }
        };
        collect_for_region(br.first);
        collect_for_region(br.second);

        if (matched_preds.size() == 1)
        {
            fid = *matched_preds.begin();
        }
    }

    if (!check_tri_edge_relations(st.tetrahedralization) ||
        !check_tet_edge_relations(st.tetrahedralization) ||
        !check_tet_tri_relations(st.tetrahedralization))
    {
        throw std::runtime_error(
          "Topological relation check failed after multi-material interpenetration removal.");
    }
}

} // namespace

void remove_interpenetrations(const CombinedSurfaces& inconsistent_combined_surfaces,
                              const ConsistentSurfaces& consistent_surfaces,
                              SegmentedTetrahedralization& st)
{
    if (st.object_region_info.regions_of_overlap.empty())
    {
        GlobalOptions::get().log_info("* No interpenetrations");
        return;
    }

    GlobalOptions::get().log_info("* remove_interpenetrations");
    GlobalOptions::get().log_info(" - Initial # of overlap regions: " +
                                  std::to_string(st.object_region_info.regions_of_overlap.size()));

    split_overlap_edges(st);

    const auto overlap_groups = find_overlap_region_groups(
      st.region_boundaries_info,
      st.object_region_info.regions_of_overlap);

    GlobalOptions::get().log_info(" - Discovered " + std::to_string(overlap_groups.size()) +
                                  " groups of contiguous overlap regions");

    std::vector<GroupHarmonicField> group_fields;
    group_fields.reserve(overlap_groups.size());
    for (const auto& overlap_group : overlap_groups)
    {
        group_fields.push_back(compute_group_harmonic_field(st, overlap_group));
    }

    subdivide_and_partition_overlaps(inconsistent_combined_surfaces,
                                     consistent_surfaces,
                                     overlap_groups,
                                     group_fields,
                                     st);

    auto region_boundary_result    = segment_region_boundaries(inconsistent_combined_surfaces,
                                                            consistent_surfaces,
                                                            st.tetrahedralization,
                                                            st.region_info,
                                                            st.object_region_info);
    st.region_boundaries_info      = std::move(region_boundary_result.region_boundaries_info);
    st.object_region_boundary_info = std::move(region_boundary_result.object_region_boundary_info);

    st.seams_info = segment_seams(st.tetrahedralization, st.region_boundaries_info);

    if (GlobalOptions::get().write_intermediate_meshes)
    {
        write_features_to_vtk(st, "remove_interpenetrations_");
    }
}

} // namespace mesh_arrange
