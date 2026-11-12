#pragma once

#include "mesh_arrange/Common.h"

#include <array>

namespace mesh_arrange::detail
{

/// Exact T2D3N edge connectivity table: e0={0,1}, e1={1,2}, e2={2,0}
inline constexpr std::array<std::array<ConnIndex, 2>, 3> kTriEdges = { {
    { 0, 1 },
    { 1, 2 },
    { 2, 0 },
} };

/// Exact T3D4N edge connectivity table: e0={0,1}, e1={0,2}, e2={0,3}, e3={1,2}, e4={1,3}, e5={2,3}
inline constexpr std::array<std::array<ConnIndex, 2>, 6> kTetEdges = { {
    { 0, 1 },
    { 0, 2 },
    { 0, 3 },
    { 1, 2 },
    { 1, 3 },
    { 2, 3 },
} };

/// Exact T3D4N face vertex connectivity table (outward when signed volume > 0):
/// f0={0,2,1}, f1={0,3,2}, f2={1,2,3}, f3={0,1,3}
inline constexpr std::array<std::array<ConnIndex, 3>, 4> kTetFaces = { {
    { 0, 2, 1 },
    { 0, 3, 2 },
    { 1, 2, 3 },
    { 0, 1, 3 },
} };

} // namespace mesh_arrange::detail
