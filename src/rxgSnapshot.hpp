#pragma once

#include <string>

#include "rxgProjectArea.hpp"

namespace rxgaming {
    constexpr int PROJECTAREA_SNAPSHOT_SCHEMA_VERSION = 1;

    void save_project_area(const RxGamingProjectArea& projectArea, const std::string& path);
    RxGamingProjectArea load_project_area(const std::string& path);
}
