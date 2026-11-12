#pragma once

#include "mesh_arrange/Common.h"
#include "mesh_arrange/Triangulation.h"

#include <Eigen/Dense>

#include <vector>

namespace mesh_arrange::detail
{

/// @brief Compact static 3D triangle AABB bounding-volume hierarchy for closest-point queries.
class TriangleAABB
{
  public:
    struct Node
    {
        Eigen::Vector3d bmin = Eigen::Vector3d::Zero();
        Eigen::Vector3d bmax = Eigen::Vector3d::Zero();
        int left             = -1;
        int right            = -1;
        int start            = 0;
        int count            = 0;
    };

    void init(const Eigen::MatrixX3d& V, const Eigen::MatrixX3i& F);

    Eigen::Vector3d closest_point(const Eigen::MatrixX3d& V,
                                  const Eigen::MatrixX3i& F,
                                  const Eigen::Vector3d& p) const;

    bool empty() const noexcept { return nodes_.empty(); }

  private:
    std::vector<int> tri_order_;
    std::vector<Node> nodes_;

    int build_recursive(const Eigen::MatrixX3d& V,
                        const Eigen::MatrixX3i& F,
                        const std::vector<Eigen::Vector3d>& centroids,
                        int start,
                        int count);
};

struct IGLTriSurf
{
    Eigen::MatrixX3d vertices;
    Eigen::MatrixX3i triangles;
    TriangleAABB tree;
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
