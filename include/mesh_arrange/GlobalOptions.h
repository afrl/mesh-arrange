#pragma once

#include <functional>
#include <iostream>
#include <string>

namespace mesh_arrange
{

/// @brief Severity level passed to `LogCallback`.
enum class LogLevel
{
    INFO    = 0,
    DEBUG   = 1,
    WARNING = 2,
    ERROR   = 3
};

using LogCallback = std::function<void(LogLevel level, const std::string& message)>;

/// @brief Thread-local diagnostic and logging controls for mesh_arrange pipelines.
struct GlobalOptions
{
    static GlobalOptions& get()
    {
        thread_local GlobalOptions opts;
        return opts;
    }

    void log(LogLevel level, const std::string& message) const
    {
        if (log_callback)
        {
            if (verbose || level != LogLevel::DEBUG)
            {
                log_callback(level, message);
            }
        }
        else if (verbose)
        {
            std::cout << message << std::endl;
        }
    }

    void log_info(const std::string& message) const
    {
        if (log_callback)
        {
            log_callback(verbose ? LogLevel::INFO : LogLevel::DEBUG, message);
        }
        else if (verbose)
        {
            std::cout << message << std::endl;
        }
    }

    bool write_intermediate_meshes  = false;
    bool write_debug_meshes         = false;
    bool check_common_errors        = false;
    bool verbose                    = false;
    std::string curr_status_padding = "";
    std::string output_dir          = ".";
    LogCallback log_callback        = {};

  private:
    GlobalOptions() = default;
};

} // namespace mesh_arrange
