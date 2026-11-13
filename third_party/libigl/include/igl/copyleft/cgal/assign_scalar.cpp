// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2015 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "assign_scalar.h"
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

template <>
IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const CGAL::Epeck::FT & rhs,
  const bool & slow_and_more_precise,
  double & lhs)
{
  // Always perform canonical IEEE-754 round-to-nearest-even via binary search
  // on CGAL::to_interval(rhs.exact()) (typically 0-2 comparisons with GMP) so
  // equal rationals from distinct intersection constructions always map to bit-
  // identical doubles before remesh_intersections deduplicates vertices.
  (void)slow_and_more_precise;
  return assign_scalar(rhs,lhs);
}

template <>
IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const CGAL::Epeck::FT & rhs,
  const bool & slow_and_more_precise,
  float & lhs)
{
  double d;
  igl::copyleft::cgal::assign_scalar(rhs,slow_and_more_precise,d);
  lhs = float(d);
}

// If we haven't specialized the types then `slow_and_more_precise` doesn't make sense.
template <typename RHS, typename LHS>
IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const RHS & rhs,
  const bool & slow_and_more_precise,
  LHS & lhs)
{
  return assign_scalar(rhs,lhs);
}

IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const CGAL::Epeck::FT & cgal,
  CGAL::Epeck::FT & d)
{
  d = cgal;
}

IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const CGAL::Epeck::FT & _cgal,
  double & d)
{
  const auto double_to_ord = [](double x) -> std::uint64_t {
    std::uint64_t bits = 0;
    std::memcpy(&bits, &x, sizeof(double));
    return (bits & (1ULL << 63)) ? ~bits : (bits | (1ULL << 63));
  };
  const auto ord_to_double = [](std::uint64_t ord) -> double {
    const std::uint64_t bits =
      (ord & (1ULL << 63)) ? (ord & ~(1ULL << 63)) : ~ord;
    double x = 0.0;
    std::memcpy(&x, &bits, sizeof(double));
    return x;
  };

  // FORCE evaluation of the exact type otherwise interval might be huge.
  const auto & exact = _cgal.exact();
  if (CGAL::is_zero(exact)) {
    d = 0.0;
    return;
  }

  const auto interval = CGAL::to_interval(exact);
  const double lo_bound = std::isfinite(interval.first)
                            ? interval.first
                            : -std::numeric_limits<double>::max();
  const double hi_bound = std::isfinite(interval.second)
                            ? interval.second
                            : std::numeric_limits<double>::max();

  if (lo_bound >= hi_bound) {
    d = (lo_bound == 0.0) ? 0.0 : lo_bound;
    return;
  }

  std::uint64_t lo = double_to_ord(lo_bound);
  std::uint64_t hi = double_to_ord(hi_bound);

  // Binary search across the IEEE-754 lattice guarantees <= 64 exact comparisons
  // even when Quotient<MP_Float>::to_interval overflows to [-inf, +inf].
  while (hi - lo > 1ULL) {
    const std::uint64_t mid = lo + (hi - lo) / 2ULL;
    const double mid_d = ord_to_double(mid);
    const CGAL::Comparison_result cmp = CGAL::compare(exact, mid_d);
    if (cmp == CGAL::LARGER) {
      lo = mid;
    } else if (cmp == CGAL::SMALLER) {
      hi = mid;
    } else {
      d = (mid_d == 0.0) ? 0.0 : mid_d;
      return;
    }
  }

  const double d_lo = ord_to_double(lo);
  const double d_hi = ord_to_double(hi);
  const CGAL::Comparison_result cmp =
    CGAL::compare(CGAL::abs(exact - d_lo), CGAL::abs(exact - d_hi));
  if (cmp == CGAL::SMALLER) {
    d = d_lo;
  } else if (cmp == CGAL::LARGER) {
    d = d_hi;
  } else {
    std::uint64_t bits_lo = 0;
    std::memcpy(&bits_lo, &d_lo, sizeof(double));
    d = ((bits_lo & 1ULL) == 0ULL) ? d_lo : d_hi;
  }
  if (d == 0.0) d = 0.0;
}

IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const CGAL::Epeck::FT & _cgal,
  float& d)
{
  // FORCE evaluation of the exact type otherwise interval might be huge.
  const CGAL::Epeck::FT cgal = _cgal.exact();
  const auto interval = CGAL::to_interval(cgal);
  d = interval.first;
  do {
      const float next = nextafter(d, float(interval.second));
      if (CGAL::abs(cgal-d) < CGAL::abs(cgal-next)) break;
      d = next;
  } while (d < float(interval.second));
}

IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const double & c,
  double & d)
{
  d = c;
}

IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const float& c,
  float& d)
{
  d = c;
}

IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const float& c,
  double& d)
{
  d = c;
}

#if !defined(WIN32) && !defined(CGAL_DISABLE_GMP)
IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT & cgal,
  CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT & d)
{
  d = cgal;
}

IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT & cgal,
  double & d)
{
  const auto interval = CGAL::to_interval(cgal);
  d = interval.first;
  do {
      const double next = nextafter(d, interval.second);
      if (CGAL::abs(cgal-d) < CGAL::abs(cgal-next)) break;
      d = next;
  } while (d < interval.second);
}

IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt::FT & cgal,
  float& d)
{
  const auto interval = CGAL::to_interval(cgal);
  d = interval.first;
  do {
      const float next = nextafter(d, float(interval.second));
      if (CGAL::abs(cgal-d) < CGAL::abs(cgal-next)) break;
      d = next;
  } while (d < float(interval.second));
}

IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const CGAL::Simple_cartesian<mpq_class>::FT & cgal,
  CGAL::Simple_cartesian<mpq_class>::FT & d)
{
  d = cgal;
}

IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const CGAL::Simple_cartesian<mpq_class>::FT & cgal,
  double & d)
{
  const auto interval = CGAL::to_interval(cgal);
  d = interval.first;
  do {
      const double next = nextafter(d, interval.second);
      if (CGAL::abs(cgal-d) < CGAL::abs(cgal-next)) break;
      d = next;
  } while (d < interval.second);
}

IGL_INLINE void igl::copyleft::cgal::assign_scalar(
  const CGAL::Simple_cartesian<mpq_class>::FT & cgal,
  float& d)
{
  const auto interval = CGAL::to_interval(cgal);
  d = interval.first;
  do {
      const float next = nextafter(d, float(interval.second));
      if (CGAL::abs(cgal-d) < CGAL::abs(cgal-next)) break;
      d = next;
  } while (d < float(interval.second));
}

#endif // !WIN32 && !CGAL_DISABLE_GMP

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
// generated by autoexplicit.sh
template void igl::copyleft::cgal::assign_scalar<float, double>(float const&, bool const&, double&);
// generated by autoexplicit.sh
template void igl::copyleft::cgal::assign_scalar<float, CGAL::Epeck::FT >(float const&, bool const&, CGAL::Epeck::FT&);
template void igl::copyleft::cgal::assign_scalar<double, double>(double const&, bool const&, double&);
template void igl::copyleft::cgal::assign_scalar<double, CGAL::Epeck::FT >(double const&, bool const&, CGAL::Epeck::FT&);
template void igl::copyleft::cgal::assign_scalar<CGAL::Epeck::FT, CGAL::Epeck::FT >(CGAL::Epeck::FT const&, bool const&, CGAL::Epeck::FT&);
#endif
