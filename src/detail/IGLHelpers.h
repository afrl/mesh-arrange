#pragma once

#include "mesh_arrange/Common.h"
#include "mesh_arrange/Triangulation.h"

#include <Eigen/Dense>
#include "igl/AABB.h"

#include <vector>

namespace mesh_arrange::detail
{

struct IGLTriSurf
{
    Eigen::MatrixX3d vertices;
    Eigen::MatrixX3i triangles;
    igl::AABB<Eigen::MatrixX3d, 3> tree;
};

struct IGLToTriangulationResult
{
    Triangulation t;
    std::vector<std::vector<Index>> element_to_origin_tris_map;
};

IGLTriSurf triangulation_to_igl(const Triangulation& t, const Marker& marker_to_extract);

std::vector<IGLTriSurf> decompose_into_igl_facets(const Triangulation& combined_surface_mesh);

IGLToTriangulationResult igl_to_triangulation(const Eigen::MatrixXd& vertices,
                                              const Eigen::MatrixXi& triangles);

Eigen::Vector3d igl_closest_point(const IGLTriSurf& igl_surf, const Eigen::Vector3d& p);

double signed_volume(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

void orient_closed_surface(const Eigen::MatrixXd& V, Eigen::MatrixXi& F, bool positive);

void orient_volume(const Eigen::MatrixXd& V, Eigen::MatrixXi& F);

} // namespace mesh_arrange::detail
