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
#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <initializer_list>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <utility>
#include "LapisGis.hpp"
#include "taolist.hpp"
#include <pybind11/numpy.h>

namespace rxgaming {
    namespace py = pybind11;

    inline void set_proj_db_path(const std::string& path) {
        std::cout << "Setting PROJ data directory to: " << path << "\n";
        lapis::setProjDefaultDirectory(path);
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

    inline bool timing_enabled() {
        static const bool enabled = [] {
            std::string value;
#ifdef _MSC_VER
            char* raw = nullptr;
            size_t len = 0;
            if (_dupenv_s(&raw, &len, "RXGAMING_TIMING") != 0 || raw == nullptr) {
                return false;
            }
            value.assign(raw);
            free(raw);
#else
            if (const char* raw = std::getenv("RXGAMING_TIMING")) {
                value.assign(raw);
            }
            else {
                return false;
            }
#endif
            std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
                return static_cast<char>(std::tolower(c));
            });
            return value == "1" || value == "true" || value == "yes" || value == "on" || value == "debug";
        }();
        return enabled;
    }

    class ScopedTimer {
    public:
        ScopedTimer(std::string name, std::initializer_list<std::pair<std::string, std::string>> fields = {})
            : _name(std::move(name)), _start(std::chrono::steady_clock::now()) {
            if (!timing_enabled()) {
                return;
            }
            bool first = true;
            for (const auto& [key, value] : fields) {
                if (!first) {
                    _context << ", ";
                }
                _context << key << "=" << value;
                first = false;
            }
        }

        ~ScopedTimer() {
            if (!timing_enabled()) {
                return;
            }
            const auto end = std::chrono::steady_clock::now();
            const auto elapsed = std::chrono::duration<double, std::milli>(end - _start).count();
            std::cout << "[rxgaming timing] " << _name << ": " << std::fixed << std::setprecision(2) << elapsed << " ms";
            if (!_context.str().empty()) {
                std::cout << " (" << _context.str() << ")";
            }
            std::cout << "\n" << std::flush;
        }

    private:
        std::string _name;
        std::chrono::steady_clock::time_point _start;
        std::ostringstream _context;
    };
}
