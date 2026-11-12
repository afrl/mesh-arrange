#pragma once

#include <string>

namespace mesh_arrange
{

/// @brief Runtime diagnostic and logging controls for mesh_arrange pipelines.
struct GlobalOptions
{
    static GlobalOptions& get()
    {
        static GlobalOptions opts;
        return opts;
    }

    bool write_intermediate_meshes   = false;
    bool write_debug_meshes          = false;
    bool check_common_errors         = false;
    bool verbose                     = false;
    std::string curr_status_padding  = "";
    std::string output_dir           = ".";

  private:
    GlobalOptions() = default;
};

} // namespace mesh_arrange
