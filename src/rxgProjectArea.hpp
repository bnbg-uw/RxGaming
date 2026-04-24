/*
Copyright(C) 2024  University of Washington
This program is free software : you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.If not, see https ://www.gnu.org/licenses/.

Bryce Bartl - Geller
University of Washington Forest Resilience Lab
2 / 22 / 2026

rxgprojectarea.hpp
*/

#pragma once
#include <functional>
#include <string>
#include "projectSettings.hpp"
#include "rxgRxUnit.hpp"
#include "readProcessedFolder.hpp"

namespace rxgaming {

    struct ProgressEvent {
        std::string stage;
        std::string message;
        int unitIndex = -1;
        std::string unitName;
        int completed = -1;
        int total = -1;
    };

    using ProgressCallback = std::function<void(const ProgressEvent&)>;

    class RxGamingProjectArea {
    public:
        std::vector<RxGamingRxUnit> rxUnits;
        
        RxGamingProjectArea(const ProjectSettings& ps, const ProgressCallback& progressCallback = ProgressCallback());
        
    };
}
