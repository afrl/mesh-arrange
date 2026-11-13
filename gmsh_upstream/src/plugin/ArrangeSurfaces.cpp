// Gmsh - Copyright (C) 1997-2026 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file in the Gmsh root directory for license information.
// Please report all issues on https://gitlab.onelab.info/gmsh/gmsh/issues.

#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include "GmshConfig.h"
#include "GModel.h"
#include "GModelIO_GEO.h"
#include "GRegion.h"
#include "GmshMessage.h"
#include "ArrangeSurfaces.h"

#if defined(HAVE_MESH) && defined(HAVE_MESH_ARRANGE)
#include "meshArrange.h"
#endif

GMSH_ArrangeSurfacesPlugin::GMSH_ArrangeSurfacesPlugin()
  : GMSH_MeshPlugin(
      {{GMSH_FULLRC, "RemoveOverlaps", nullptr, 1.,
        "Partition equal-priority multi-material overlaps along smooth "
        "harmonic medial surfaces (1: yes, 0: no)"},
       {GMSH_FULLRC, "OptimizeSurfaces", nullptr, 0.,
        "Optimize arranged surface triangulations via edge flips and "
        "feature-preserving tangential smoothing (1: yes, 0: no)"},
       {GMSH_FULLRC, "MaxEdgeLength", nullptr, 0.,
        "Optional upper bound on arranged surface edge length (0: disabled)"},
       {GMSH_FULLRC, "CreateGeometry", nullptr, 0.,
        "Parameterize arranged discrete surfaces for subsequent 2D/3D "
        "remeshing (1: yes, 0: no)"},
       {GMSH_FULLRC, "RemoveObject", nullptr, 1.,
        "Remove input solid volumes/surfaces after arrangement (1: yes, 0: no)"},
       {GMSH_FULLRC, "RemoveTool", nullptr, 1.,
        "Remove input clipping volumes/surfaces after arrangement "
        "(1: yes, 0: no)"},
       {GMSH_FULLRC, "PropagatePhysicals", nullptr, 1.,
        "Propagate 3D Physical Groups from input volumes to arranged child "
        "volumes (1: yes, 0: no)"}},
      {{GMSH_FULLRC, "VolumeTags", nullptr, "",
        "Comma-separated 3D volume tags to arrange (empty and SurfaceTags "
        "empty => all non-tool 3D volumes in current model)"},
       {GMSH_FULLRC, "SurfaceTags", nullptr, "",
        "Comma-separated closed 2D surface shell tags to arrange"},
       {GMSH_FULLRC, "Priorities", nullptr, "",
        "Optional comma-separated integer priorities matching VolumeTags "
        "followed by SurfaceTags"},
       {GMSH_FULLRC, "ClipInsideVolumeTags", nullptr, "",
        "Comma-separated 3D tool volume tags whose interior is subtracted "
        "(voids)"},
       {GMSH_FULLRC, "ClipOutsideVolumeTags", nullptr, "",
        "Comma-separated 3D tool volume tags whose exterior is clipped "
        "(RVE window)"},
       {GMSH_FULLRC, "ClipPriorities", nullptr, "",
        "Optional comma-separated integer priorities for ClipInsideVolumeTags "
        "followed by ClipOutsideVolumeTags"}})
{
}

std::string GMSH_ArrangeSurfacesPlugin::getHelp() const
{
  return "Plugin(ArrangeSurfaces) intersects closed discrete 3D volumes or 2D "
         "surface shells, resolves boolean priorities and RVE/void clippers, "
         "segments non-manifold 3D/2D/1D/0D B-Rep topology, and optionally "
         "removes multi-material interpenetrations along smooth harmonic "
         "medial surfaces.\n\n"
         "If `VolumeTags' and `SurfaceTags' are both empty, all 3D volumes in "
         "the current model (excluding any volumes listed in "
         "`ClipInsideVolumeTags' or `ClipOutsideVolumeTags') are arranged.\n\n"
         "If `PropagatePhysicals' is enabled, 3D Physical Groups defined on "
         "input volumes are automatically transferred to their resulting "
         "arranged child volumes.\n\n"
         "For programmatic access with full parent-to-child volume provenance "
         "(outDimTagsMap), see gmsh::model::mesh::arrangeSurfaces.";
}

#if defined(HAVE_MESH) && defined(HAVE_MESH_ARRANGE)

int GMSH_ArrangeSurfacesPlugin::run()
{
  const bool removeOverlaps = (option(0) != 0.);
  const bool optimizeSurfaces = (option(1) != 0.);
  const double maxEdgeLength = option(2);
  const bool createGeometry = (option(3) != 0.);
  const bool removeObject = (option(4) != 0.);
  const bool removeTool = (option(5) != 0.);
  const bool propagatePhysicals = (option(6) != 0.);

  std::vector<int> volTags, surfTags, objPrios;
  std::vector<int> clipInVols, clipOutVols, clipPrios;
  if(!optionIntList(0, volTags) || !optionIntList(1, surfTags) ||
     !optionIntList(2, objPrios) || !optionIntList(3, clipInVols) ||
     !optionIntList(4, clipOutVols) || !optionIntList(5, clipPrios))
    return 0;

  GModel *m = GModel::current();
  if(!m) {
    Msg::Error("Plugin(ArrangeSurfaces): no current model");
    return 0;
  }

  std::set<int> toolVolSet(clipInVols.begin(), clipInVols.end());
  toolVolSet.insert(clipOutVols.begin(), clipOutVols.end());

  std::vector<std::pair<int, int> > objectDimTags;
  for(std::size_t i = 0; i < volTags.size(); ++i)
    objectDimTags.emplace_back(3, volTags[i]);
  for(std::size_t i = 0; i < surfTags.size(); ++i)
    objectDimTags.emplace_back(2, surfTags[i]);

  if(objectDimTags.empty()) {
    for(auto it = m->firstRegion(); it != m->lastRegion(); ++it) {
      int rTag = (*it)->tag();
      if(toolVolSet.find(rTag) == toolVolSet.end())
        objectDimTags.emplace_back(3, rTag);
    }
  }

  if(objectDimTags.empty()) {
    Msg::Error("Plugin(ArrangeSurfaces): no input volumes or surfaces found");
    return 0;
  }

  std::vector<std::pair<int, int> > toolDimTags;
  std::vector<int> toolTypes;
  for(std::size_t i = 0; i < clipInVols.size(); ++i) {
    toolDimTags.emplace_back(3, clipInVols[i]);
    toolTypes.push_back(1);
  }
  for(std::size_t i = 0; i < clipOutVols.size(); ++i) {
    toolDimTags.emplace_back(3, clipOutVols[i]);
    toolTypes.push_back(-1);
  }

  // Capture existing 3D Physical Groups before input volumes are removed
  struct SavedPhysGroup {
    std::string name;
    std::vector<int> origRegionTags;
  };
  std::map<int, SavedPhysGroup> savedGroups3D;
  std::map<int, std::size_t> volTagToObjIdx;
  if(propagatePhysicals) {
    for(std::size_t i = 0; i < objectDimTags.size(); ++i) {
      if(objectDimTags[i].first == 3)
        volTagToObjIdx[objectDimTags[i].second] = i;
    }
    std::map<int, std::vector<GEntity *> > groups;
    m->getPhysicalGroups(3, groups);
    for(auto &kv : groups) {
      SavedPhysGroup sg;
      sg.name = m->getPhysicalName(3, kv.first);
      for(std::size_t j = 0; j < kv.second.size(); ++j)
        sg.origRegionTags.push_back(kv.second[j]->tag());
      savedGroups3D[kv.first] = std::move(sg);
    }
    std::vector<int> emptyTags;
    for(auto &kv : savedGroups3D) {
      if(m->getGEOInternals())
        m->getGEOInternals()->modifyPhysicalGroup(3, kv.first, 2, emptyTags);
      m->removePhysicalGroup(3, kv.first);
    }
  }

  std::vector<std::pair<int, int> > outDimTags;
  std::vector<std::vector<std::pair<int, int> > > outDimTagsMap;
  int rc = arrangeDiscreteSurfaces(m, objectDimTags, toolDimTags, outDimTags,
                                   outDimTagsMap, objPrios, toolTypes,
                                   clipPrios, removeOverlaps, optimizeSurfaces,
                                   maxEdgeLength, createGeometry, removeObject,
                                   removeTool);
  if(rc != 0) {
    for(auto &kv : savedGroups3D) {
      if(m->getGEOInternals())
        m->getGEOInternals()->modifyPhysicalGroup(3, kv.first, 0,
                                                  kv.second.origRegionTags);
      m->addPhysicalGroup(3, kv.first, kv.second.origRegionTags);
      if(!kv.second.name.empty())
        m->setPhysicalName(kv.second.name, 3, kv.first);
    }
    return 0;
  }

  if(propagatePhysicals && !savedGroups3D.empty()) {
    for(auto &kv : savedGroups3D) {
      const int physTag = kv.first;
      const SavedPhysGroup &sg = kv.second;
      std::vector<int> newTags;
      std::set<int> seen;
      for(std::size_t j = 0; j < sg.origRegionTags.size(); ++j) {
        int origTag = sg.origRegionTags[j];
        auto itObj = volTagToObjIdx.find(origTag);
        if(itObj != volTagToObjIdx.end() &&
           itObj->second < outDimTagsMap.size()) {
          if(!removeObject && m->getRegionByTag(origTag) &&
             seen.insert(origTag).second) {
            newTags.push_back(origTag);
          }
          const auto &children = outDimTagsMap[itObj->second];
          for(std::size_t k = 0; k < children.size(); ++k) {
            if(children[k].first == 3 &&
               seen.insert(children[k].second).second) {
              newTags.push_back(children[k].second);
            }
          }
        }
        else if(m->getRegionByTag(origTag) && seen.insert(origTag).second) {
          newTags.push_back(origTag);
        }
      }
      std::vector<int> emptyTags;
      if(m->getGEOInternals())
        m->getGEOInternals()->modifyPhysicalGroup(3, physTag, 2, emptyTags);
      m->removePhysicalGroup(3, physTag);
      if(!newTags.empty()) {
        if(m->getGEOInternals())
          m->getGEOInternals()->modifyPhysicalGroup(3, physTag, 0, newTags);
        m->addPhysicalGroup(3, physTag, newTags);
        if(!sg.name.empty()) m->setPhysicalName(sg.name, 3, physTag);
      }
    }
  }

  return (int)outDimTags.size();
}

#else

int GMSH_ArrangeSurfacesPlugin::run()
{
  Msg::Error("Plugin(ArrangeSurfaces) requires MeshArrange and the mesh "
             "module");
  return 0;
}

#endif
