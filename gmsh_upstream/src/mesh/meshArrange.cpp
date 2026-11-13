// Gmsh - Copyright (C) 1997-2026 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file in the Gmsh root directory for license information.
// Please report all issues on https://gitlab.onelab.info/gmsh/gmsh/issues.

#include "GmshConfig.h"
#include "meshArrange.h"
#include "GmshMessage.h"

#if defined(HAVE_MESH_ARRANGE) && defined(HAVE_MESH)

#include "Context.h"
#include "GmshDefines.h"
#include "GModel.h"
#include "GModelIO_GEO.h"
#include "GRegion.h"
#include "discreteVertex.h"
#include "discreteEdge.h"
#include "discreteFace.h"
#include "discreteRegion.h"
#include "MVertex.h"
#include "MPoint.h"
#include "MLine.h"
#include "MTriangle.h"
#include "MTetrahedron.h"
#include "meshGRegion.h"

#include "mesh_arrange/MeshArrange.h"
#include "mesh_arrange/PlcBackend.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

using mesh_arrange::ArrangeOptions;
using mesh_arrange::ArrangeResult;
using mesh_arrange::Index;
using mesh_arrange::InputSurface;
using mesh_arrange::LogLevel;
using mesh_arrange::Marker;
using mesh_arrange::NOT_SET_Index;
using mesh_arrange::ObjectType;
using mesh_arrange::PlcTetOutput;
using mesh_arrange::Point3;
using mesh_arrange::Triangulation;

struct SortedTriHash {
  std::size_t operator()(const std::array<Index, 3> &k) const noexcept
  {
    std::size_t h = static_cast<std::size_t>(k[0]);
    h ^= static_cast<std::size_t>(k[1]) + 0x9e3779b9u + (h << 6) + (h >> 2);
    h ^= static_cast<std::size_t>(k[2]) + 0x9e3779b9u + (h << 6) + (h >> 2);
    return h;
  }
};

double pointTriangleDistSq(const Point3 &p, const Point3 &a, const Point3 &b,
                           const Point3 &c) noexcept
{
  const double ab[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
  const double ac[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
  const double ap[3] = {p[0] - a[0], p[1] - a[1], p[2] - a[2]};
  const double d1 = ab[0] * ap[0] + ab[1] * ap[1] + ab[2] * ap[2];
  const double d2 = ac[0] * ap[0] + ac[1] * ap[1] + ac[2] * ap[2];
  if(d1 <= 0.0 && d2 <= 0.0)
    return ap[0] * ap[0] + ap[1] * ap[1] + ap[2] * ap[2];

  const double bp[3] = {p[0] - b[0], p[1] - b[1], p[2] - b[2]};
  const double d3 = ab[0] * bp[0] + ab[1] * bp[1] + ab[2] * bp[2];
  const double d4 = ac[0] * bp[0] + ac[1] * bp[1] + ac[2] * bp[2];
  if(d3 >= 0.0 && d4 <= d3)
    return bp[0] * bp[0] + bp[1] * bp[1] + bp[2] * bp[2];

  const double vc = d1 * d4 - d3 * d2;
  if(vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
    const double v = d1 / (d1 - d3);
    const double dx = ap[0] - v * ab[0];
    const double dy = ap[1] - v * ab[1];
    const double dz = ap[2] - v * ab[2];
    return dx * dx + dy * dy + dz * dz;
  }

  const double cp[3] = {p[0] - c[0], p[1] - c[1], p[2] - c[2]};
  const double d5 = ab[0] * cp[0] + ab[1] * cp[1] + ab[2] * cp[2];
  const double d6 = ac[0] * cp[0] + ac[1] * cp[1] + ac[2] * cp[2];
  if(d6 >= 0.0 && d5 <= d6)
    return cp[0] * cp[0] + cp[1] * cp[1] + cp[2] * cp[2];

  const double vb = d5 * d2 - d1 * d6;
  if(vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
    const double w = d2 / (d2 - d6);
    const double dx = ap[0] - w * ac[0];
    const double dy = ap[1] - w * ac[1];
    const double dz = ap[2] - w * ac[2];
    return dx * dx + dy * dy + dz * dz;
  }

  const double va = d3 * d6 - d5 * d4;
  if(va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
    const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    const double qx = b[0] + w * (c[0] - b[0]) - p[0];
    const double qy = b[1] + w * (c[1] - b[1]) - p[1];
    const double qz = b[2] + w * (c[2] - b[2]) - p[2];
    return qx * qx + qy * qy + qz * qz;
  }

  const double denom = 1.0 / (va + vb + vc);
  const double v = vb * denom;
  const double w = vc * denom;
  const double dx = ap[0] - v * ab[0] - w * ac[0];
  const double dy = ap[1] - v * ab[1] - w * ac[1];
  const double dz = ap[2] - v * ab[2] - w * ac[2];
  return dx * dx + dy * dy + dz * dz;
}

struct ScopedMeshPlcContext {
  int abortOnError;
  int verbosity;
  int algo3d;
  int optimize;
  int optimizeNetgen;
  double lcMax;
  double lcMin;
  int lcExtendFromBoundary;
  int renumber;
  double angleToleranceFacetOverlap;
  double toleranceInitialDelaunay;

  ScopedMeshPlcContext(double diag, int algo)
  {
    abortOnError = CTX::instance()->abortOnError;
    verbosity = Msg::GetVerbosity();
    algo3d = CTX::instance()->mesh.algo3d;
    optimize = CTX::instance()->mesh.optimize;
    optimizeNetgen = CTX::instance()->mesh.optimizeNetgen;
    lcMax = CTX::instance()->mesh.lcMax;
    lcMin = CTX::instance()->mesh.lcMin;
    lcExtendFromBoundary = CTX::instance()->mesh.lcExtendFromBoundary;
    renumber = CTX::instance()->mesh.renumber;
    angleToleranceFacetOverlap = CTX::instance()->mesh.angleToleranceFacetOverlap;
    toleranceInitialDelaunay = CTX::instance()->mesh.toleranceInitialDelaunay;

    CTX::instance()->abortOnError = 0;
    if(verbosity > 1) Msg::SetVerbosity(1);
    CTX::instance()->mesh.algo3d = algo;
    CTX::instance()->mesh.optimize = 0;
    CTX::instance()->mesh.optimizeNetgen = 0;
    CTX::instance()->mesh.lcMax = 10.0 * diag;
    CTX::instance()->mesh.lcMin = 10.0 * diag;
    CTX::instance()->mesh.lcExtendFromBoundary = 0;
    CTX::instance()->mesh.renumber = 0;
    CTX::instance()->mesh.angleToleranceFacetOverlap = 0.0;
    CTX::instance()->mesh.toleranceInitialDelaunay = 1e-14;
  }

  ~ScopedMeshPlcContext()
  {
    CTX::instance()->abortOnError = abortOnError;
    Msg::SetVerbosity(verbosity);
    CTX::instance()->mesh.algo3d = algo3d;
    CTX::instance()->mesh.optimize = optimize;
    CTX::instance()->mesh.optimizeNetgen = optimizeNetgen;
    CTX::instance()->mesh.lcMax = lcMax;
    CTX::instance()->mesh.lcMin = lcMin;
    CTX::instance()->mesh.lcExtendFromBoundary = lcExtendFromBoundary;
    CTX::instance()->mesh.renumber = renumber;
    CTX::instance()->mesh.angleToleranceFacetOverlap = angleToleranceFacetOverlap;
    CTX::instance()->mesh.toleranceInitialDelaunay = toleranceInitialDelaunay;
  }
};

struct ScopedCurrentModel {
  int prevIndex;
  bool switched;

  explicit ScopedCurrentModel(GModel *m)
    : prevIndex(GModel::getCurrentIndex()), switched(false)
  {
    if(m && GModel::current() != m) {
      GModel::setCurrent(m);
      switched = true;
    }
  }

  ~ScopedCurrentModel()
  {
    if(switched) GModel::setCurrentIndex(prevIndex);
  }
};

// Native Gmsh PLC boundary-recovery callback replacing TetGen ("pYQ").
// Encloses the non-manifold PLC in a padded bounding-box discreteRegion with an
// embedded discreteFace and runs MeshDelaunayVolume (tetgenBR).
PlcTetOutput gmshInternalPlcTetrahedralizer(const Triangulation &plc)
{
  if(plc.vertices.empty() || plc.triangles.empty())
    throw std::invalid_argument("Empty PLC passed to gmshInternalPlcTetrahedralizer.");

  Point3 bmin = plc.vertices.front();
  Point3 bmax = plc.vertices.front();
  for(std::size_t i = 0; i < plc.vertices.size(); ++i) {
    for(int k = 0; k < 3; ++k) {
      bmin[k] = std::min(bmin[k], plc.vertices[i][k]);
      bmax[k] = std::max(bmax[k], plc.vertices[i][k]);
    }
  }
  double diag = 0.0;
  for(int k = 0; k < 3; ++k) {
    const double d = bmax[k] - bmin[k];
    diag += d * d;
  }
  diag = std::max(std::sqrt(diag), 1.0);
  const double pad = 0.35 * diag;
  for(int k = 0; k < 3; ++k) {
    bmin[k] -= pad;
    bmax[k] += pad;
  }

  // Only weld exact coincident coordinates (<= 1e-15 * diag); MeshArrange's
  // weld_sub_ulp_consistent_edges already performs CGAL-verified sub-ULP edge
  // welding, and larger unverified collapses can sweep coplanar edges across
  // neighboring vertices.
  const std::size_t nRaw = plc.vertices.size();
  const double weldTol2 = (1e-15 * diag) * (1e-15 * diag);
  const double minArea2 = (1e-15 * diag * diag) * (1e-15 * diag * diag);
  std::vector<Index> ufPar(nRaw);
  for(std::size_t i = 0; i < nRaw; ++i) ufPar[i] = static_cast<Index>(i);
  auto ufFind = [&](Index x) -> Index {
    while(ufPar[x] != x) {
      ufPar[x] = ufPar[ufPar[x]];
      x = ufPar[x];
    }
    return x;
  };

  std::vector<std::vector<std::size_t> > vertTris(nRaw);
  for(std::size_t i = 0; i < plc.triangles.size(); ++i) {
    for(int c = 0; c < 3; ++c) vertTris[plc.triangles[i][c]].push_back(i);
  }

  auto canCollapseWithoutFold = [&](Index fromRep, Index toRep) -> bool {
    const Point3 &pTo = plc.vertices[toRep];
    for(std::size_t k = 0; k < vertTris[fromRep].size(); ++k) {
      const auto &tr = plc.triangles[vertTris[fromRep][k]];
      Index r[3] = {ufFind(tr[0]), ufFind(tr[1]), ufFind(tr[2])};
      if(r[0] == r[1] || r[1] == r[2] || r[2] == r[0]) continue;
      if(r[0] == toRep || r[1] == toRep || r[2] == toRep) continue;
      if(r[0] != fromRep && r[1] != fromRep && r[2] != fromRep) continue;

      const Point3 &A0 = plc.vertices[r[0]];
      const Point3 &B0 = plc.vertices[r[1]];
      const Point3 &C0 = plc.vertices[r[2]];
      const double ox = (B0[1] - A0[1]) * (C0[2] - A0[2]) - (B0[2] - A0[2]) * (C0[1] - A0[1]);
      const double oy = (B0[2] - A0[2]) * (C0[0] - A0[0]) - (B0[0] - A0[0]) * (C0[2] - A0[2]);
      const double oz = (B0[0] - A0[0]) * (C0[1] - A0[1]) - (B0[1] - A0[1]) * (C0[0] - A0[0]);
      const double oldLen2 = ox * ox + oy * oy + oz * oz;
      if(oldLen2 <= minArea2) continue;

      const Point3 &A1 = (r[0] == fromRep) ? pTo : A0;
      const Point3 &B1 = (r[1] == fromRep) ? pTo : B0;
      const Point3 &C1 = (r[2] == fromRep) ? pTo : C0;
      const double nx = (B1[1] - A1[1]) * (C1[2] - A1[2]) - (B1[2] - A1[2]) * (C1[1] - A1[1]);
      const double ny = (B1[2] - A1[2]) * (C1[0] - A1[0]) - (B1[0] - A1[0]) * (C1[2] - A1[2]);
      const double nz = (B1[0] - A1[0]) * (C1[1] - A1[1]) - (B1[1] - A1[1]) * (C1[0] - A1[0]);
      const double newLen2 = nx * nx + ny * ny + nz * nz;
      if(newLen2 <= minArea2) return false;

      const double dot = ox * nx + oy * ny + oz * nz;
      if(dot <= 0.0 || dot * dot < 0.998 * oldLen2 * newLen2) return false;
    }
    return true;
  };

  for(std::size_t i = 0; i < plc.triangles.size(); ++i) {
    const auto &tr = plc.triangles[i];
    for(int e = 0; e < 3; ++e) {
      const Index ra = ufFind(tr[e]);
      const Index rb = ufFind(tr[(e + 1) % 3]);
      if(ra != rb) {
        const double dx = plc.vertices[ra][0] - plc.vertices[rb][0];
        const double dy = plc.vertices[ra][1] - plc.vertices[rb][1];
        const double dz = plc.vertices[ra][2] - plc.vertices[rb][2];
        if(dx * dx + dy * dy + dz * dz <= weldTol2) {
          const Index keep = std::min(ra, rb);
          const Index drop = std::max(ra, rb);
          if(canCollapseWithoutFold(drop, keep)) {
            ufPar[drop] = keep;
            vertTris[keep].insert(vertTris[keep].end(), vertTris[drop].begin(),
                                  vertTris[drop].end());
          }
          else if(canCollapseWithoutFold(keep, drop)) {
            ufPar[keep] = drop;
            vertTris[drop].insert(vertTris[drop].end(), vertTris[keep].begin(),
                                  vertTris[keep].end());
          }
        }
      }
    }
  }

  std::vector<Point3> cleanVertices;
  std::vector<std::array<Index, 3> > cleanTriangles;
  std::vector<Index> cleanParentIds;
  {
    std::vector<Index> rawToClean(nRaw, NOT_SET_Index);
    cleanVertices.reserve(nRaw);
    for(std::size_t i = 0; i < nRaw; ++i) {
      const Index r = ufFind(static_cast<Index>(i));
      if(rawToClean[r] == NOT_SET_Index) {
        rawToClean[r] = static_cast<Index>(cleanVertices.size());
        cleanVertices.push_back(plc.vertices[r]);
      }
      rawToClean[i] = rawToClean[r];
    }
    std::unordered_map<std::array<Index, 3>, bool, SortedTriHash> seenTris;
    seenTris.reserve(plc.triangles.size());
    cleanTriangles.reserve(plc.triangles.size());
    cleanParentIds.reserve(plc.triangles.size());
    for(std::size_t p = 0; p < plc.triangles.size(); ++p) {
      const Index c0 = rawToClean[plc.triangles[p][0]];
      const Index c1 = rawToClean[plc.triangles[p][1]];
      const Index c2 = rawToClean[plc.triangles[p][2]];
      if(c0 == c1 || c1 == c2 || c2 == c0) continue;
      const Point3 &A = cleanVertices[c0];
      const Point3 &B = cleanVertices[c1];
      const Point3 &C = cleanVertices[c2];
      const double nx = (B[1] - A[1]) * (C[2] - A[2]) - (B[2] - A[2]) * (C[1] - A[1]);
      const double ny = (B[2] - A[2]) * (C[0] - A[0]) - (B[0] - A[0]) * (C[2] - A[2]);
      const double nz = (B[0] - A[0]) * (C[1] - A[1]) - (B[1] - A[1]) * (C[0] - A[0]);
      if(nx * nx + ny * ny + nz * nz <= minArea2) continue;
      std::array<Index, 3> ks = {c0, c1, c2};
      std::sort(ks.begin(), ks.end());
      if(seenTris.emplace(ks, true).second) {
        cleanTriangles.push_back({c0, c1, c2});
        cleanParentIds.push_back(static_cast<Index>(p));
      }
    }
  }

  auto buildAndRecover = [&](GModel &scratch, std::vector<MVertex *> &plcVerts,
                             GRegion *&grOut, int algo) {
    auto *boxFace = new discreteFace(&scratch, 1, false);
    auto *plcFace = new discreteFace(&scratch, 2, false);
    scratch.add(boxFace);
    scratch.add(plcFace);

    const std::size_t nPlc = cleanVertices.size();
    plcVerts.assign(nPlc, nullptr);
    for(std::size_t i = 0; i < nPlc; ++i) {
      const auto &p = cleanVertices[i];
      auto *mv = new MVertex(p[0], p[1], p[2], plcFace, i + 1);
      plcFace->mesh_vertices.push_back(mv);
      plcVerts[i] = mv;
    }
    for(std::size_t i = 0; i < cleanTriangles.size(); ++i) {
      const auto &t = cleanTriangles[i];
      plcFace->triangles.push_back(
        new MTriangle(plcVerts[t[0]], plcVerts[t[1]], plcVerts[t[2]]));
    }

    const std::array<Point3, 8> boxPts = {{
      {bmin[0], bmin[1], bmin[2]}, {bmax[0], bmin[1], bmin[2]},
      {bmax[0], bmax[1], bmin[2]}, {bmin[0], bmax[1], bmin[2]},
      {bmin[0], bmin[1], bmax[2]}, {bmax[0], bmin[1], bmax[2]},
      {bmax[0], bmax[1], bmax[2]}, {bmin[0], bmax[1], bmax[2]}
    }};
    std::array<MVertex *, 8> boxVerts;
    for(std::size_t i = 0; i < 8; ++i) {
      auto *mv = new MVertex(boxPts[i][0], boxPts[i][1], boxPts[i][2], boxFace,
                             nPlc + 1 + i);
      boxFace->mesh_vertices.push_back(mv);
      boxVerts[i] = mv;
    }
    static const int boxTris[12][3] = {
      {0, 2, 1}, {0, 3, 2}, {4, 5, 6}, {4, 6, 7},
      {0, 1, 5}, {0, 5, 4}, {2, 3, 7}, {2, 7, 6},
      {1, 2, 6}, {1, 6, 5}, {3, 0, 4}, {3, 4, 7}
    };
    for(int i = 0; i < 12; ++i) {
      boxFace->triangles.push_back(new MTriangle(
        boxVerts[boxTris[i][0]], boxVerts[boxTris[i][1]], boxVerts[boxTris[i][2]]));
    }

    auto *gr = new discreteRegion(&scratch, 1, false);
    scratch.add(gr);
    gr->setBoundFaces(std::vector<GFace *>{boxFace}, std::vector<int>{1});
    gr->addEmbeddedFace(plcFace);

    ScopedMeshPlcContext ctxGuard(diag, algo);
    std::vector<GRegion *> regions = {gr};
    MeshDelaunayVolume(regions);
    scratch.setAllVolumesPositive();
    grOut = gr;
  };

  const int prevErrCount = Msg::GetErrorCount();
  Msg::NoThrow noThrowGuard;

  GModel scratch("__mesh_arrange_plc__");
  std::vector<MVertex *> plcVerts;
  GRegion *gr = nullptr;
  buildAndRecover(scratch, plcVerts, gr, ALGO_3D_INITIAL_ONLY);
  if(!gr || gr->tetrahedra.empty()) {
    if(prevErrCount == 0) Msg::ResetErrorCounter();
    scratch.destroy();
    gr = nullptr;
    buildAndRecover(scratch, plcVerts, gr, ALGO_3D_DELAUNAY);
  }
  if(!gr || gr->tetrahedra.empty()) {
    if(prevErrCount == 0) Msg::ResetErrorCounter();
    throw std::runtime_error("Gmsh 3D PLC boundary recovery failed.");
  }
  if(prevErrCount == 0 && Msg::GetErrorCount() > 0)
    Msg::ResetErrorCounter();

  const std::size_t nPlc = cleanVertices.size();
  PlcTetOutput out;
  out.vertices = cleanVertices;
  std::unordered_map<MVertex *, Index> mvToLocal;
  mvToLocal.reserve(nPlc + 64);
  for(std::size_t i = 0; i < nPlc; ++i) {
    mvToLocal[plcVerts[i]] = static_cast<Index>(i);
  }

  out.tetrahedra.reserve(gr->tetrahedra.size());
  for(std::size_t i = 0; i < gr->tetrahedra.size(); ++i) {
    MTetrahedron *mt = gr->tetrahedra[i];
    std::array<Index, 4> conn;
    for(int c = 0; c < 4; ++c) {
      MVertex *mv = mt->getVertex(c);
      auto it = mvToLocal.find(mv);
      if(it == mvToLocal.end()) {
        const Index idx = static_cast<Index>(out.vertices.size());
        mvToLocal[mv] = idx;
        out.vertices.push_back({mv->x(), mv->y(), mv->z()});
        conn[c] = idx;
      }
      else {
        conn[c] = it->second;
      }
    }
    out.tetrahedra.push_back(conn);
  }

  std::unordered_map<std::array<Index, 3>, bool, SortedTriHash> tetFaces;
  std::vector<std::array<Index, 3> > uniqueTetFaces;
  uniqueTetFaces.reserve(out.tetrahedra.size() * 2);
  for(std::size_t i = 0; i < out.tetrahedra.size(); ++i) {
    const auto &tet = out.tetrahedra[i];
    static const int F[4][3] = {{0, 2, 1}, {0, 3, 2}, {1, 2, 3}, {0, 1, 3}};
    for(int f = 0; f < 4; ++f) {
      std::array<Index, 3> k = {tet[F[f][0]], tet[F[f][1]], tet[F[f][2]]};
      std::array<Index, 3> ks = k;
      std::sort(ks.begin(), ks.end());
      if(tetFaces.emplace(ks, true).second)
        uniqueTetFaces.push_back(k);
    }
  }

  bool anySplit = false;
  const double tol2 = (1e-6 * diag) * (1e-6 * diag);
  for(std::size_t p = 0; p < cleanTriangles.size(); ++p) {
    const auto &pt = cleanTriangles[p];
    const Index origPid = cleanParentIds[p];
    std::array<Index, 3> ks = pt;
    std::sort(ks.begin(), ks.end());
    if(tetFaces.find(ks) != tetFaces.end()) {
      out.surface_triangles.push_back(pt);
      out.surface_tri_parent_id.push_back(origPid);
    }
    else {
      anySplit = true;
      const Point3 &A = out.vertices[pt[0]];
      const Point3 &B = out.vertices[pt[1]];
      const Point3 &C = out.vertices[pt[2]];
      const double np[3] = {
        (B[1] - A[1]) * (C[2] - A[2]) - (B[2] - A[2]) * (C[1] - A[1]),
        (B[2] - A[2]) * (C[0] - A[0]) - (B[0] - A[0]) * (C[2] - A[2]),
        (B[0] - A[0]) * (C[1] - A[1]) - (B[1] - A[1]) * (C[0] - A[0])
      };
      std::size_t subCnt = 0;
      for(std::size_t j = 0; j < uniqueTetFaces.size(); ++j) {
        const auto &tf = uniqueTetFaces[j];
        const Point3 &P0 = out.vertices[tf[0]];
        const Point3 &P1 = out.vertices[tf[1]];
        const Point3 &P2 = out.vertices[tf[2]];
        if(pointTriangleDistSq(P0, A, B, C) <= tol2 &&
           pointTriangleDistSq(P1, A, B, C) <= tol2 &&
           pointTriangleDistSq(P2, A, B, C) <= tol2) {
          const double nf[3] = {
            (P1[1] - P0[1]) * (P2[2] - P0[2]) - (P1[2] - P0[2]) * (P2[1] - P0[1]),
            (P1[2] - P0[2]) * (P2[0] - P0[0]) - (P1[0] - P0[0]) * (P2[2] - P0[2]),
            (P1[0] - P0[0]) * (P2[1] - P0[1]) - (P1[1] - P0[1]) * (P2[0] - P0[0])
          };
          std::array<Index, 3> oriented = tf;
          if(np[0] * nf[0] + np[1] * nf[1] + np[2] * nf[2] < 0.0)
            std::swap(oriented[1], oriented[2]);
          out.surface_triangles.push_back(oriented);
          out.surface_tri_parent_id.push_back(origPid);
          ++subCnt;
        }
      }
      if(subCnt == 0)
        throw std::runtime_error("Failed to recover Steiner sub-triangles for PLC facet.");
    }
  }

  const bool anyModified = (cleanVertices.size() != plc.vertices.size() ||
                            cleanTriangles.size() != plc.triangles.size());
  if(!anyModified && !anySplit) {
    out.surface_triangles.clear();
    out.surface_tri_parent_id.clear();
  }
  return out;
}

void orientClosedShellOutward(InputSurface &s)
{
  const std::size_t nTris = s.triangles.size();
  if(nTris == 0) return;

  struct EdgeAdj {
    std::size_t tri;
    int dir;
  };
  std::unordered_map<std::uint64_t, std::vector<EdgeAdj> > edgeMap;
  edgeMap.reserve(nTris * 2);
  for(std::size_t t = 0; t < nTris; ++t) {
    const auto &tr = s.triangles[t];
    for(int e = 0; e < 3; ++e) {
      const Index u = tr[e];
      const Index v = tr[(e + 1) % 3];
      const Index lo = std::min(u, v);
      const Index hi = std::max(u, v);
      const std::uint64_t key =
        (static_cast<std::uint64_t>(lo) << 32) | static_cast<std::uint32_t>(hi);
      edgeMap[key].push_back({t, (u == lo) ? 1 : -1});
    }
  }

  std::vector<int> flipState(nTris, 0);
  std::vector<std::size_t> queue;
  std::vector<std::size_t> comp;
  for(std::size_t seed = 0; seed < nTris; ++seed) {
    if(flipState[seed] != 0) continue;
    flipState[seed] = 1;
    queue.clear();
    comp.clear();
    queue.push_back(seed);
    std::size_t head = 0;
    while(head < queue.size()) {
      const std::size_t cur = queue[head++];
      comp.push_back(cur);
      const auto &tr = s.triangles[cur];
      for(int e = 0; e < 3; ++e) {
        const Index u = tr[e];
        const Index v = tr[(e + 1) % 3];
        const Index lo = std::min(u, v);
        const Index hi = std::max(u, v);
        const int curDir = (u == lo) ? 1 : -1;
        const std::uint64_t key =
          (static_cast<std::uint64_t>(lo) << 32) | static_cast<std::uint32_t>(hi);
        const auto &adj = edgeMap[key];
        if(adj.size() != 2) continue;
        for(std::size_t k = 0; k < 2; ++k) {
          const std::size_t nb = adj[k].tri;
          if(nb == cur || flipState[nb] != 0) continue;
          flipState[nb] = -(curDir * flipState[cur] * adj[k].dir);
          queue.push_back(nb);
        }
      }
    }

    double vol6 = 0.0;
    for(std::size_t i = 0; i < comp.size(); ++i) {
      auto tr = s.triangles[comp[i]];
      if(flipState[comp[i]] < 0) std::swap(tr[1], tr[2]);
      const Point3 &a = s.vertices[tr[0]];
      const Point3 &b = s.vertices[tr[1]];
      const Point3 &c = s.vertices[tr[2]];
      vol6 += a[0] * (b[1] * c[2] - b[2] * c[1]) +
              a[1] * (b[2] * c[0] - b[0] * c[2]) +
              a[2] * (b[0] * c[1] - b[1] * c[0]);
    }
    const int globalSign = (vol6 < 0.0) ? -1 : 1;
    for(std::size_t i = 0; i < comp.size(); ++i) {
      if(flipState[comp[i]] * globalSign < 0)
        std::swap(s.triangles[comp[i]][1], s.triangles[comp[i]][2]);
    }
  }
}

bool extractInputSurface(GModel *m, int dim, int tag, ObjectType type,
                         int priority, InputSurface &out)
{
  std::vector<GFace *> faces;
  if(dim == 3) {
    GRegion *gr = m->getRegionByTag(tag);
    if(!gr) {
      Msg::Error("Volume %d does not exist in model", tag);
      return false;
    }
    faces = gr->faces();
  }
  else if(dim == 2) {
    GFace *gf = m->getFaceByTag(std::abs(tag));
    if(!gf) {
      Msg::Error("Surface %d does not exist in model", tag);
      return false;
    }
    faces.push_back(gf);
  }
  else {
    Msg::Error("arrangeSurfaces expects entities of dimension 2 or 3 (got dim=%d)", dim);
    return false;
  }

  out = InputSurface{};
  out.type = type;
  out.priority = priority;

  std::unordered_map<MVertex *, Index> mvMap;
  for(std::size_t f = 0; f < faces.size(); ++f) {
    GFace *gf = faces[f];
    const Marker facetId = static_cast<Marker>(gf->tag());
    for(std::size_t t = 0; t < gf->triangles.size(); ++t) {
      MTriangle *mt = gf->triangles[t];
      Index idx[3];
      for(int c = 0; c < 3; ++c) {
        MVertex *mv = mt->getVertex(c);
        auto it = mvMap.find(mv);
        if(it == mvMap.end()) {
          const Index vid = static_cast<Index>(out.vertices.size());
          mvMap[mv] = vid;
          out.vertices.push_back({mv->x(), mv->y(), mv->z()});
          idx[c] = vid;
        }
        else {
          idx[c] = it->second;
        }
      }
      out.triangles.push_back({idx[0], idx[1], idx[2]});
      out.triangle_facet_ids.push_back(facetId);
    }
  }

  if(out.triangles.empty()) {
    Msg::Error("Entity (%d, %d) has no 2D triangle mesh (generate 2D mesh first)", dim, tag);
    return false;
  }

  orientClosedShellOutward(out);
  return true;
}

} // namespace

int arrangeDiscreteSurfaces(
  GModel *m,
  const std::vector<std::pair<int, int> > &objectDimTags,
  const std::vector<std::pair<int, int> > &toolDimTags,
  std::vector<std::pair<int, int> > &outDimTags,
  std::vector<std::vector<std::pair<int, int> > > &outDimTagsMap,
  const std::vector<int> &objectPriorities,
  const std::vector<int> &toolTypes,
  const std::vector<int> &toolPriorities,
  bool removeOverlaps,
  bool optimizeSurfaces,
  double maxEdgeLength,
  bool createGeometry,
  bool removeObject,
  bool removeTool)
{
  if(!m) return 1;
  outDimTags.clear();
  outDimTagsMap.clear();

  if(!objectPriorities.empty() && objectPriorities.size() != objectDimTags.size()) {
    Msg::Error("objectPriorities size (%zu) must match objectDimTags size (%zu)",
               objectPriorities.size(), objectDimTags.size());
    return 1;
  }
  if(!toolTypes.empty() && toolTypes.size() != toolDimTags.size()) {
    Msg::Error("toolTypes size (%zu) must match toolDimTags size (%zu)",
               toolTypes.size(), toolDimTags.size());
    return 1;
  }
  if(!toolPriorities.empty() && toolPriorities.size() != toolDimTags.size()) {
    Msg::Error("toolPriorities size (%zu) must match toolDimTags size (%zu)",
               toolPriorities.size(), toolDimTags.size());
    return 1;
  }

  const std::size_t nObj = objectDimTags.size();
  const std::size_t nTool = toolDimTags.size();
  std::vector<InputSurface> inputs;
  inputs.reserve(nObj + nTool);

  for(std::size_t i = 0; i < nObj; ++i) {
    InputSurface s;
    const int prio = objectPriorities.empty() ? 0 : objectPriorities[i];
    if(!extractInputSurface(m, objectDimTags[i].first, objectDimTags[i].second,
                            ObjectType::SOLID, prio, s))
      return 1;
    inputs.push_back(std::move(s));
  }
  for(std::size_t i = 0; i < nTool; ++i) {
    InputSurface s;
    const int tcode = toolTypes.empty() ? 1 : toolTypes[i];
    const ObjectType type =
      (tcode >= 0) ? ObjectType::CLIP_INSIDE : ObjectType::CLIP_OUTSIDE;
    const int prio = toolPriorities.empty() ? 0 : toolPriorities[i];
    if(!extractInputSurface(m, toolDimTags[i].first, toolDimTags[i].second,
                            type, prio, s))
      return 1;
    inputs.push_back(std::move(s));
  }

  ArrangeOptions opts;
  opts.remove_overlaps = removeOverlaps;
  opts.optimize_surfaces = optimizeSurfaces;
  opts.max_edge_length = maxEdgeLength;
  opts.plc_tetrahedralizer = gmshInternalPlcTetrahedralizer;
  opts.log_callback = [](LogLevel lvl, const std::string &msg) {
    switch(lvl) {
    case LogLevel::ERROR: Msg::Error("%s", msg.c_str()); break;
    case LogLevel::WARNING: Msg::Warning("%s", msg.c_str()); break;
    case LogLevel::DEBUG: Msg::Debug("%s", msg.c_str()); break;
    case LogLevel::INFO:
    default: Msg::Info("%s", msg.c_str()); break;
    }
  };

  ArrangeResult res;
  try {
    res = mesh_arrange::arrange_surfaces(inputs, opts);
  }
  catch(const std::exception &ex) {
    Msg::Error("arrangeSurfaces failed: %s", ex.what());
    return 1;
  }

  ScopedCurrentModel currentGuard(m);

  std::vector<std::pair<int, int> > toRemove;
  if(removeObject)
    toRemove.insert(toRemove.end(), objectDimTags.begin(), objectDimTags.end());
  if(removeTool)
    toRemove.insert(toRemove.end(), toolDimTags.begin(), toolDimTags.end());
  if(!toRemove.empty()) {
    if(m->getGEOInternals()) {
      m->getGEOInternals()->remove(toRemove, true);
      m->getGEOInternals()->synchronize(m);
    }
    std::vector<GEntity *> removed;
    m->remove(toRemove, removed, true);
    for(std::size_t i = 0; i < removed.size(); ++i) delete removed[i];
  }

  const int baseVTag = m->getMaxElementaryNumber(0) + 1;
  const int baseETag = m->getMaxElementaryNumber(1) + 1;
  const int baseFTag = m->getMaxElementaryNumber(2) + 1;
  const int baseRTag = m->getMaxElementaryNumber(3) + 1;

  const auto &coords = res.surface_mesh.vertices;
  std::unordered_map<Index, MVertex *> vidToMVertex;
  vidToMVertex.reserve(coords.size());

  std::vector<discreteVertex *> dVerts(res.brep_vertices.size(), nullptr);
  for(std::size_t i = 0; i < res.brep_vertices.size(); ++i) {
    const auto &bv = res.brep_vertices[i];
    auto *gv = new discreteVertex(m, baseVTag + static_cast<int>(i),
                                  bv.point[0], bv.point[1], bv.point[2]);
    auto *mv = new MVertex(bv.point[0], bv.point[1], bv.point[2], gv);
    gv->mesh_vertices.push_back(mv);
    gv->points.push_back(new MPoint(mv));
    m->add(gv);
    dVerts[i] = gv;
    vidToMVertex[bv.global_vertex_id] = mv;
  }

  std::vector<discreteEdge *> dEdges(res.brep_curves.size(), nullptr);
  for(std::size_t i = 0; i < res.brep_curves.size(); ++i) {
    const auto &bc = res.brep_curves[i];
    GVertex *vBeg =
      (bc.start_brep_vertex != NOT_SET_Index) ? dVerts[bc.start_brep_vertex] : nullptr;
    GVertex *vEnd =
      (bc.end_brep_vertex != NOT_SET_Index) ? dVerts[bc.end_brep_vertex] : nullptr;
    auto *ge = new discreteEdge(m, baseETag + static_cast<int>(i), vBeg, vEnd);
    if(bc.polyline_vertices.size() >= 2) {
      const std::size_t lastInterior = bc.polyline_vertices.size() - 1;
      for(std::size_t k = 1; k < lastInterior; ++k) {
        const Index vid = bc.polyline_vertices[k];
        if(vidToMVertex.find(vid) == vidToMVertex.end()) {
          const auto &p = coords[vid];
          auto *mv = new MVertex(p[0], p[1], p[2], ge);
          ge->mesh_vertices.push_back(mv);
          vidToMVertex.emplace(vid, mv);
        }
      }
      for(std::size_t k = 0; k + 1 < bc.polyline_vertices.size(); ++k) {
        MVertex *mv0 = vidToMVertex.at(bc.polyline_vertices[k]);
        MVertex *mv1 = vidToMVertex.at(bc.polyline_vertices[k + 1]);
        ge->lines.push_back(new MLine(mv0, mv1));
      }
    }
    m->add(ge);
    dEdges[i] = ge;
  }

  std::vector<discreteFace *> dFaces(res.brep_faces.size(), nullptr);
  std::vector<std::pair<int, int> > discreteDimTagsForParam;
  for(std::size_t i = 0; i < res.brep_faces.size(); ++i) {
    const auto &bf = res.brep_faces[i];
    const int fTag = baseFTag + static_cast<int>(i);
    auto *gf = new discreteFace(m, fTag);

    std::vector<int> eTags;
    std::vector<int> eDirs;
    eTags.reserve(bf.signed_bounding_curves.size());
    eDirs.reserve(bf.signed_bounding_curves.size());
    for(std::size_t k = 0; k < bf.signed_bounding_curves.size(); ++k) {
      eTags.push_back(dEdges.at(bf.signed_bounding_curves[k].first)->tag());
      eDirs.push_back(bf.signed_bounding_curves[k].second >= 0 ? 1 : -1);
    }
    if(!eTags.empty())
      gf->setBoundEdges(eTags, eDirs);

    for(std::size_t k = 0; k < bf.embedded_curves.size(); ++k)
      gf->addEmbeddedEdge(dEdges.at(bf.embedded_curves[k]));
    for(std::size_t k = 0; k < bf.embedded_vertices.size(); ++k)
      gf->addEmbeddedVertex(dVerts.at(bf.embedded_vertices[k]));

    for(std::size_t t = 0; t < bf.triangles.size(); ++t) {
      const auto &tri = bf.triangles[t];
      MVertex *mvs[3] = {nullptr, nullptr, nullptr};
      for(int c = 0; c < 3; ++c) {
        const Index vid = tri[c];
        auto it = vidToMVertex.find(vid);
        if(it == vidToMVertex.end()) {
          const auto &p = coords[vid];
          auto *mv = new MVertex(p[0], p[1], p[2], gf);
          gf->mesh_vertices.push_back(mv);
          it = vidToMVertex.emplace(vid, mv).first;
        }
        mvs[c] = it->second;
      }
      gf->triangles.push_back(new MTriangle(mvs[0], mvs[1], mvs[2]));
    }
    m->add(gf);
    dFaces[i] = gf;
    discreteDimTagsForParam.emplace_back(2, fTag);
  }

  if(createGeometry && !discreteDimTagsForParam.empty())
    m->createGeometryOfDiscreteEntities(discreteDimTagsForParam);

  // Create built-in GEO Volumes bounded by the arranged discreteFaces so that
  // gmsh::model::mesh::generate(3) meshes all output regions directly.
  outDimTags.reserve(res.brep_volumes.size());
  for(std::size_t i = 0; i < res.brep_volumes.size(); ++i) {
    const auto &bvol = res.brep_volumes[i];
    std::vector<int> fTags;
    fTags.reserve(bvol.signed_bounding_faces.size());
    for(std::size_t k = 0; k < bvol.signed_bounding_faces.size(); ++k) {
      fTags.push_back(dFaces.at(bvol.signed_bounding_faces[k].first)->tag());
    }
    int slTag = -1;
    int rTag = baseRTag + static_cast<int>(i);
    m->getGEOInternals()->addSurfaceLoop(slTag, fTags);
    m->getGEOInternals()->addVolume(rTag, {slTag});
    outDimTags.emplace_back(3, rTag);
  }
  m->getGEOInternals()->synchronize(m);
  m->destroyMeshCaches();

  outDimTagsMap.assign(nObj + nTool, std::vector<std::pair<int, int> >());
  for(std::size_t i = 0; i < res.brep_volumes.size(); ++i) {
    const auto &owners = res.brep_volumes[i].owner_objects;
    for(std::size_t k = 0; k < owners.size(); ++k) {
      const std::size_t oid = static_cast<std::size_t>(owners[k]);
      if(oid < outDimTagsMap.size())
        outDimTagsMap[oid].push_back(outDimTags[i]);
    }
  }

  return 0;
}

#else

int arrangeDiscreteSurfaces(
  GModel *,
  const std::vector<std::pair<int, int> > &,
  const std::vector<std::pair<int, int> > &,
  std::vector<std::pair<int, int> > &,
  std::vector<std::vector<std::pair<int, int> > > &,
  const std::vector<int> &,
  const std::vector<int> &,
  const std::vector<int> &,
  bool, bool, double, bool, bool, bool)
{
  Msg::Error("Gmsh must be compiled with ENABLE_MESH_ARRANGE and ENABLE_MESH "
             "to arrange discrete surfaces");
  return 1;
}

#endif
