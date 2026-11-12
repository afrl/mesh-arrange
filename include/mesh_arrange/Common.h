#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>

namespace mesh_arrange
{

using ConnIndex = uint16_t;
using Index     = uint32_t;
using Marker    = uint32_t;
using FacetID   = Marker;
using SeamID    = Marker;
using Point3    = std::array<double, 3>;

constexpr Index NOT_SET_Index     = std::numeric_limits<Index>::max();
constexpr Marker NOT_SET_Marker   = std::numeric_limits<Marker>::max();
constexpr FacetID NOT_SET_FacetID = std::numeric_limits<FacetID>::max();
constexpr SeamID NOT_SET_SeamID   = std::numeric_limits<SeamID>::max();

/// @brief Fast, collision-free 64-bit undirected hash for 2-node linear edges.
inline uint64_t hash_linear_edge(const Index& v1, const Index& v2) noexcept
{
    const uint64_t lo = static_cast<uint32_t>(std::min(v1, v2));
    const uint64_t hi = static_cast<uint32_t>(std::max(v1, v2));
    return lo | (hi << 32);
}

} // namespace mesh_arrange
