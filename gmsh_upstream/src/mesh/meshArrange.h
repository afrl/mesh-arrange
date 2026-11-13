// Gmsh - Copyright (C) 1997-2026 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file in the Gmsh root directory for license information.
// Please report all issues on https://gitlab.onelab.info/gmsh/gmsh/issues.

#ifndef MESH_ARRANGE_GMSH_H
#define MESH_ARRANGE_GMSH_H

#include <utility>
#include <vector>

class GModel;

// Intersect a collection of closed discrete volumes (dim == 3) or closed
// surface shells (dim == 2) in `objectDimTags` and optional boolean clippers
// in `toolDimTags` (`toolTypes[k] >= 0` -> clipInside, `< 0` -> clipOutside),
// segment 3D regions, 2D region boundaries, and 1D/0D non-manifold seams, and
// optionally partition equal-priority multi-material overlaps along smooth
// harmonic medial surfaces.
//
// Populates `outDimTags` with all resulting (3, volumeTag) regions and
// `outDimTagsMap` with the child (3, volumeTag) regions belonging to each
// input entity in `objectDimTags` followed by `toolDimTags` (matching the
// convention of `gmsh::model::occ::fragment`).
int arrangeDiscreteSurfaces(
  GModel *m,
  const std::vector<std::pair<int, int> > &objectDimTags,
  const std::vector<std::pair<int, int> > &toolDimTags,
  std::vector<std::pair<int, int> > &outDimTags,
  std::vector<std::vector<std::pair<int, int> > > &outDimTagsMap,
  const std::vector<int> &objectPriorities = std::vector<int>(),
  const std::vector<int> &toolTypes = std::vector<int>(),
  const std::vector<int> &toolPriorities = std::vector<int>(),
  bool removeOverlaps = true,
  bool optimizeSurfaces = false,
  double maxEdgeLength = 0.,
  bool createGeometry = false,
  bool removeObject = true,
  bool removeTool = true);

#endif
