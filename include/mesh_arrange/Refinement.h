#pragma once

#include "mesh_arrange/FeatureIdentification.h"

#include <vector>

namespace mesh_arrange
{

void split_edges(const std::vector<Index>& edges_to_split, SegmentedTetrahedralization& st);

void split_edges(const std::vector<Index>& edges_to_split,
                 const std::vector<Index>& existing_vertices_for_split_edges,
                 SegmentedTetrahedralization& st);

void split_constrained_edges(SegmentedTetrahedralization& st);

void split_overlap_edges(SegmentedTetrahedralization& st);

void split_large_edges(const double max_edge_length, SegmentedTetrahedralization& st);

} // namespace mesh_arrange
