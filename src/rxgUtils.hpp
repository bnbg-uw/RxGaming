/*
Copyright(C) 2024  University of Washington
This program is free software : you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.If not, see https ://www.gnu.org/licenses/.

Bryce Bartl - Geller
University of Washington Forest Resilience Lab
2 / 22 / 2026

rxgUtils.cpp
*/

#pragma once
#include <random>
#include "LapisGis.hpp"

namespace rxgaming {

    void set_proj_db_path(const std::string& path) {
        std::cout << "Setting PROJ data directory to: " << path << "\n";
        lapis::lapisGisInit(path);

        if(!lapis::isProjDirectorySet()) {
            std::cout << "Failed to set PROJ data directory. Is the path correct?\n";
        }
        if(lapis::projDbExists()) {
            std::cout << "PROJ database found at: " << lapis::getProjDirectory() << "\n";
        }
        else {
            std::cout << "Failed to find PROJ database. Is the PROJ data directory correct?\n";
        }
    }

    inline std::mt19937& globalRng() {
        static std::mt19937 rng{std::random_device{}()};
        return rng;
    }

    inline void set_seed(uint32_t seed) {
        globalRng().seed(seed);
    }
}