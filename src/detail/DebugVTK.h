#pragma once

#include "mesh_arrange/Common.h"
#include "mesh_arrange/Triangulation.h"

#include <array>
#include <string>
#include <utility>
#include <vector>

namespace mesh_arrange::detail
{

struct VtkScalarFieldU32
{
    std::string name;
    std::vector<uint32_t> values;
};

struct VtkScalarFieldF64
{
    std::string name;
    std::vector<double> values;
};

void write_vtk_triangles(const std::string& filename,
                         const std::vector<Point3>& vertices,
                         const std::vector<std::array<Index, 3>>& triangles,
                         const std::vector<VtkScalarFieldU32>& cell_fields  = {},
                         const std::vector<VtkScalarFieldU32>& point_fields = {});

void write_vtk_triangulation(const std::string& filename, const Triangulation& t);

void write_vtk_lines(const std::string& filename,
                     const std::vector<Point3>& vertices,
                     const std::vector<std::array<Index, 2>>& lines,
                     const std::vector<VtkScalarFieldU32>& cell_fields = {});

void write_vtk_tetrahedra(const std::string& filename,
                          const std::vector<Point3>& vertices,
                          const std::vector<std::array<Index, 4>>& tetrahedra,
                          const std::vector<VtkScalarFieldU32>& cell_fields  = {},
                          const std::vector<VtkScalarFieldF64>& point_fields = {});

} // namespace mesh_arrange::detail
