#pragma once

#include "mesh_arrange/Common.h"

#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <utility>
#include <vector>

namespace mesh_arrange::detail
{

inline Eigen::Vector3d point_as_vector(const Point3& p) noexcept
{
    return Eigen::Vector3d(p[0], p[1], p[2]);
}

inline Point3 to_point(const Eigen::Vector3d& v) noexcept
{
    return { v.x(), v.y(), v.z() };
}

inline Point3 add_points(const Point3& a, const Point3& b) noexcept
{
    return { a[0] + b[0], a[1] + b[1], a[2] + b[2] };
}

inline Point3 midpoint(const Point3& a, const Point3& b) noexcept
{
    return { 0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2]) };
}

inline double distance(const Point3& a, const Point3& b) noexcept
{
    const double dx = a[0] - b[0];
    const double dy = a[1] - b[1];
    const double dz = a[2] - b[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

/// Exact signed volume matching `(p2 - p1).dot((p3 - p1).cross(p4 - p1)) / 6.0`.
inline double tetrahedron_signed_volume(const Eigen::Vector3d& p1,
                                        const Eigen::Vector3d& p2,
                                        const Eigen::Vector3d& p3,
                                        const Eigen::Vector3d& p4) noexcept
{
    return (p2 - p1).dot((p3 - p1).cross(p4 - p1)) / 6.0;
}

inline double tetrahedron_signed_volume(const Point3& p1,
                                        const Point3& p2,
                                        const Point3& p3,
                                        const Point3& p4) noexcept
{
    return tetrahedron_signed_volume(point_as_vector(p1),
                                     point_as_vector(p2),
                                     point_as_vector(p3),
                                     point_as_vector(p4));
}

struct SegmentClosestPoint
{
    double distance      = 0.0;
    Point3 closest_point = { 0.0, 0.0, 0.0 };
};

inline SegmentClosestPoint closest_point_on_line_segment(const Point3& p,
                                                         const Point3& seg_p1,
                                                         const Point3& seg_p2) noexcept
{
    const Eigen::Vector3d pv   = point_as_vector(p);
    const Eigen::Vector3d p1   = point_as_vector(seg_p1);
    const Eigen::Vector3d p2   = point_as_vector(seg_p2);
    const Eigen::Vector3d diff = p2 - p1;
    const double len2          = diff.squaredNorm();
    if (len2 <= 1e-30)
    {
        return { (pv - p1).norm(), seg_p1 };
    }
    const double t = std::clamp((pv - p1).dot(diff) / len2, 0.0, 1.0);
    if (t <= 0.0)
    {
        return { (pv - p1).norm(), seg_p1 };
    }
    if (t >= 1.0)
    {
        return { (pv - p2).norm(), seg_p2 };
    }
    const Eigen::Vector3d cp = p1 + t * diff;
    return { (pv - cp).norm(), to_point(cp) };
}

inline std::size_t mix_hash(std::size_t x) noexcept
{
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    x = x ^ (x >> 31);
    return x;
}

inline std::size_t hash_combine_pair(std::size_t x, std::size_t y) noexcept
{
    const auto hx = mix_hash(x);
    return mix_hash(y) + 0x9e3779b9ULL + (hx << 6) + (hx >> 2);
}

template<typename Container>
inline std::size_t hash_container(const Container& vec) noexcept
{
    std::size_t seed = vec.size();
    for (const auto& v : vec)
    {
        seed ^= mix_hash(static_cast<std::size_t>(v)) + 0x9e3779b9ULL + (seed << 6) + (seed >> 2);
    }
    return seed;
}

template<typename Container>
inline std::size_t hash_sorted_container(const Container& vec)
{
    std::vector<std::size_t> idx(vec.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(),
                     idx.end(),
                     [&vec](std::size_t a, std::size_t b) { return vec[a] < vec[b]; });
    std::size_t seed = vec.size();
    for (std::size_t i : idx)
    {
        seed ^=
          mix_hash(static_cast<std::size_t>(vec[i])) + 0x9e3779b9ULL + (seed << 6) + (seed >> 2);
    }
    return seed;
}

template<typename T>
inline std::vector<std::size_t> find_sorted_permutation(const std::vector<T>& vec)
{
    std::vector<std::size_t> p(vec.size());
    std::iota(p.begin(), p.end(), 0);
    std::stable_sort(p.begin(), p.end(), [&](std::size_t i, std::size_t j) { return vec[i] < vec[j]; });
    return p;
}

} // namespace mesh_arrange::detail
