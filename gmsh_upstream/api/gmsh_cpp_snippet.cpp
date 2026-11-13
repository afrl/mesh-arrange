// Insert #include "meshArrange.h" under #if defined(HAVE_MESH) at top of
// Gmsh src/common/gmsh.cpp, and place this implementation right after
// gmsh::model::mesh::createTopology:

GMSH_API void gmsh::model::mesh::arrangeSurfaces(
  const vectorpair &objectDimTags, const vectorpair &toolDimTags,
  vectorpair &outDimTags, std::vector<vectorpair> &outDimTagsMap,
  const std::vector<int> &objectPriorities, const std::vector<int> &toolTypes,
  const std::vector<int> &toolPriorities, const bool removeOverlaps,
  const bool optimizeSurfaces, const double maxEdgeLength,
  const bool createGeometry, const bool removeObject, const bool removeTool)
{
  if(!_checkInit()) return;
#if defined(HAVE_MESH)
  arrangeDiscreteSurfaces(GModel::current(), objectDimTags, toolDimTags,
                          outDimTags, outDimTagsMap, objectPriorities,
                          toolTypes, toolPriorities, removeOverlaps,
                          optimizeSurfaces, maxEdgeLength, createGeometry,
                          removeObject, removeTool);
#else
  Msg::Error("Arranging discrete surfaces requires the mesh module");
#endif
}
