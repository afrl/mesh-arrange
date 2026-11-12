#include "detail/DebugVTK.h"

#include <fstream>
#include <iomanip>

namespace mesh_arrange::detail
{

namespace
{
void write_vtk_points(std::ofstream& out, const std::vector<Point3>& vertices)
{
    out << "# vtk DataFile Version 3.0\nMeshArrange Debug Output\nASCII\nDATASET UNSTRUCTURED_GRID\n";
    out << "POINTS " << vertices.size() << " double\n" << std::setprecision(16);
    for (const auto& p : vertices)
    {
        out << p[0] << " " << p[1] << " " << p[2] << "\n";
    }
}
} // namespace

void write_vtk_triangles(const std::string& filename,
                         const std::vector<Point3>& vertices,
                         const std::vector<std::array<Index, 3>>& triangles,
                         const std::vector<VtkScalarFieldU32>& cell_fields,
                         const std::vector<VtkScalarFieldU32>& point_fields)
{
    std::ofstream out(filename);
    if (!out)
    {
        return;
    }
    write_vtk_points(out, vertices);
    out << "CELLS " << triangles.size() << " " << triangles.size() * 4 << "\n";
    for (const auto& tri : triangles)
    {
        out << "3 " << tri[0] << " " << tri[1] << " " << tri[2] << "\n";
    }
    out << "CELL_TYPES " << triangles.size() << "\n";
    for (std::size_t i = 0; i < triangles.size(); ++i)
    {
        out << "5\n";
    }

    if (!point_fields.empty() && !vertices.empty())
    {
        out << "POINT_DATA " << vertices.size() << "\n";
        for (const auto& f : point_fields)
        {
            if (f.values.size() != vertices.size())
            {
                continue;
            }
            out << "SCALARS " << f.name << " unsigned_int 1\nLOOKUP_TABLE default\n";
            for (uint32_t v : f.values)
            {
                out << v << "\n";
            }
        }
    }

    if (!cell_fields.empty() && !triangles.empty())
    {
        out << "CELL_DATA " << triangles.size() << "\n";
        for (const auto& f : cell_fields)
        {
            if (f.values.size() != triangles.size())
            {
                continue;
            }
            out << "SCALARS " << f.name << " unsigned_int 1\nLOOKUP_TABLE default\n";
            for (uint32_t v : f.values)
            {
                out << v << "\n";
            }
        }
    }
}

void write_vtk_triangulation(const std::string& filename, const Triangulation& t)
{
    std::vector<VtkScalarFieldU32> fields;
    if (t.tri_to_marker_map.size() == t.triangles.size())
    {
        fields.push_back({ "marker", t.tri_to_marker_map });
    }
    write_vtk_triangles(filename, t.vertices, t.triangles, fields);
}

void write_vtk_lines(const std::string& filename,
                     const std::vector<Point3>& vertices,
                     const std::vector<std::array<Index, 2>>& lines,
                     const std::vector<VtkScalarFieldU32>& cell_fields)
{
    std::ofstream out(filename);
    if (!out)
    {
        return;
    }
    write_vtk_points(out, vertices);
    out << "CELLS " << lines.size() << " " << lines.size() * 3 << "\n";
    for (const auto& l : lines)
    {
        out << "2 " << l[0] << " " << l[1] << "\n";
    }
    out << "CELL_TYPES " << lines.size() << "\n";
    for (std::size_t i = 0; i < lines.size(); ++i)
    {
        out << "3\n";
    }
    if (!cell_fields.empty() && !lines.empty())
    {
        out << "CELL_DATA " << lines.size() << "\n";
        for (const auto& f : cell_fields)
        {
            if (f.values.size() != lines.size())
            {
                continue;
            }
            out << "SCALARS " << f.name << " unsigned_int 1\nLOOKUP_TABLE default\n";
            for (uint32_t v : f.values)
            {
                out << v << "\n";
            }
        }
    }
}

void write_vtk_tetrahedra(const std::string& filename,
                          const std::vector<Point3>& vertices,
                          const std::vector<std::array<Index, 4>>& tetrahedra,
                          const std::vector<VtkScalarFieldU32>& cell_fields,
                          const std::vector<VtkScalarFieldF64>& point_fields)
{
    std::ofstream out(filename);
    if (!out)
    {
        return;
    }
    write_vtk_points(out, vertices);
    out << "CELLS " << tetrahedra.size() << " " << tetrahedra.size() * 5 << "\n";
    for (const auto& tet : tetrahedra)
    {
        out << "4 " << tet[0] << " " << tet[1] << " " << tet[2] << " " << tet[3] << "\n";
    }
    out << "CELL_TYPES " << tetrahedra.size() << "\n";
    for (std::size_t i = 0; i < tetrahedra.size(); ++i)
    {
        out << "10\n";
    }

    if (!point_fields.empty() && !vertices.empty())
    {
        out << "POINT_DATA " << vertices.size() << "\n";
        for (const auto& f : point_fields)
        {
            if (f.values.size() != vertices.size())
            {
                continue;
            }
            out << "SCALARS " << f.name << " double 1\nLOOKUP_TABLE default\n";
            for (double v : f.values)
            {
                out << v << "\n";
            }
        }
    }

    if (!cell_fields.empty() && !tetrahedra.empty())
    {
        out << "CELL_DATA " << tetrahedra.size() << "\n";
        for (const auto& f : cell_fields)
        {
            if (f.values.size() != tetrahedra.size())
            {
                continue;
            }
            out << "SCALARS " << f.name << " unsigned_int 1\nLOOKUP_TABLE default\n";
            for (uint32_t v : f.values)
            {
                out << v << "\n";
            }
        }
    }
}

} // namespace mesh_arrange::detail
