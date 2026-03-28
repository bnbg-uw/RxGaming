/*
Copyright(C) 2024  University of Washington
This program is free software : you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.If not, see https ://www.gnu.org/licenses/.

Bryce Bartl - Geller
University of Washington Forest Resilience Lab
2 / 22 / 2026

rxgprojectarea.cpp
*/

#include "rxgProjectArea.hpp"

namespace rxgaming {

    RxGamingProjectArea::RxGamingProjectArea(const ProjectSettings& ps) {
        auto unitPolygon = lapis::MultiPolygon(ps.unitPolyPath);

        auto lidar = processedfolder::readProcessedFolder(ps.lidarPath);

        //clean and reproject unitpolygon

        //allometry on convex hull

        //do preprocessing
    }
}