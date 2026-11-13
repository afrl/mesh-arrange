// Gmsh - Copyright (C) 1997-2026 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file in the Gmsh root directory for license information.
// Please report all issues on https://gitlab.onelab.info/gmsh/gmsh/issues.

#ifndef ARRANGE_SURFACES_PLUGIN_H
#define ARRANGE_SURFACES_PLUGIN_H

#include <string>
#include "Plugin.h"

class GMSH_ArrangeSurfacesPlugin : public GMSH_MeshPlugin {
public:
  GMSH_ArrangeSurfacesPlugin();
  std::string getName() const { return "ArrangeSurfaces"; }
  std::string getShortHelp() const
  {
    return "Non-manifold discrete surface arrangement and overlap removal";
  }
  std::string getHelp() const;
  std::string getAuthor() const
  {
    return "AFRL / Mesh-Arrange Contributors";
  }
  int run();
};

#endif
