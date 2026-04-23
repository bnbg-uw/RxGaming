/*
Copyright(C) 2024  University of Washington
This program is free software : you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.If not, see https ://www.gnu.org/licenses/.

Bryce Bartl - Geller
University of Washington Forest Resilience Lab
2 / 22 / 2026

rxgRxUnit.cpp
*/

#include "rxgRxUnit.hpp"

namespace rxgaming {

    RxGamingRxUnit::RxGamingRxUnit(std::string name, lapis::Raster<lapis::cell_t> mask, lapis::Raster<double> chm, lapis::Raster<int> basinMap,
            rxtools::TaoList taos) : RxUnit(mask, taos) {
            this->name = name;
            this->chm = chm;
            this->basinMap = basinMap;
            hillshade = computeHillshade(chm);
            targetStructure = currentStructure;
            treatedStructure = currentStructure;
        };

    py::array_t<lapis::cell_t> RxGamingRxUnit::get_mask() const {
        return raster_to_numpy(unitMask);
    }

    py::array_t<double> RxGamingRxUnit::get_chm() const {
        return raster_to_numpy(chm);
    }

    py::array_t<lapis::cell_t> RxGamingRxUnit::get_basin() const {
        return raster_to_numpy(basinMap);
    }

    py::array_t<double> RxGamingRxUnit::get_hillshade() const {
        return raster_to_numpy(hillshade);
    }

    py::array_t<lapis::cell_t> RxGamingRxUnit::get_clump_map() const {
        try {    
            std::unordered_map<int, int> taoIds;
            lapis::Raster<int> clumpMap((lapis::Alignment)basinMap);
            auto rawClumps = getRawClumps(taos);
            for (size_t i = 0; i < taos.size(); ++i) {
                auto e = clumpMap.extract(taos.x(i), taos.y(i), lapis::ExtractMethod::near);
                if (e.has_value() && e.value() != 1) {
                    taoIds.emplace(std::make_pair(e.value(), (int)rawClumps[i]));
                }
            }
    
            for (lapis::cell_t j = 0; j < clumpMap.ncell(); ++j) {
                if (clumpMap[j].has_value()) {
                    auto x = taoIds.find(clumpMap[j].value());
                    if (x != taoIds.end()) {
                        clumpMap[j].value() = x->second;
                    }
                    else {
                        clumpMap[j].value() = 0;
                    }
                }
            }
            return raster_to_numpy(clumpMap);
        }
        catch (std::exception e) {
            std::cout << e.what() << "\n";
            std::abort();
        }
    }

    py::array_t<double> RxGamingRxUnit::get_taos() const {
        return taolist_to_numpy(taos);
    }

    py::array_t<double> RxGamingRxUnit::get_treat_chm() const {
        std::unordered_set<int> basinIds;
        for (size_t i = 0; i < treatedTaos.size(); ++i) {
            auto v = basinMap.extract(treatedTaos.x(i), treatedTaos.y(i), lapis::ExtractMethod::near);
            if (v.has_value()) {
                basinIds.emplace(v.value());
            }
        }

        auto thisChm = chm;
        for (lapis::cell_t i = 0; i < thisChm.ncell(); ++i) {
            if (thisChm[i].has_value()) {
                if (basinIds.find(basinMap[i].value()) == basinIds.end()) {
                    thisChm[i].value() = 0;
                }
            }
        }

        return(raster_to_numpy(thisChm));
    }

    py::array_t<lapis::cell_t> RxGamingRxUnit::get_treat_basin() const {
        return(raster_to_numpy(getTreatBasin()));
    }

    py::array_t<double> RxGamingRxUnit::get_treat_hillshade() const {
        auto b = getTreatBasin();
        auto thisHill = hillshade;
        for (lapis::cell_t i = 0; i < thisHill.ncell(); ++i) {
            if (thisHill[i].has_value()) {
                if (b[i].has_value() && b[i].value() == 1) {
                    thisHill[i].value() = 200;
                }
            }
        }
        return(raster_to_numpy(thisHill));
    }

    py::array_t<lapis::cell_t> RxGamingRxUnit::get_treat_clump_map() const {
        try {
            auto b = getTreatBasin();
    
            std::unordered_map<int, int> taoIds;
            
            auto groupsizes = getRawClumps(treatedTaos);
            for (size_t i = 0; i < treatedTaos.size(); ++i) {
                auto e = b.extract(treatedTaos.x(i), treatedTaos.y(i), lapis::ExtractMethod::near);
                if (e.has_value() && e.value() != 1) {
                    taoIds.emplace(std::make_pair(e.value(), (int)groupsizes[i]));
                }
            }
    
            for (lapis::cell_t j = 0; j < b.ncell(); ++j) {
                if (b[j].has_value()) {
                    auto x = taoIds.find(b[j].value());
                    if (x != taoIds.end()) {
                        b[j].value() = x->second;
                    }
                    else {
                        b[j].value() = 0;
                    }
                }
            }
    
            return(raster_to_numpy(b));
        }
        catch (std::exception e) {
            std::cout << e.what() << "\n";
            std::abort();
        }
    }

    py::array_t<double> RxGamingRxUnit::get_treat_taos() const {
        return taolist_to_numpy(treatedTaos);
    }

    py::array_t<double> RxGamingRxUnit::get_cut_taos() const {
        return taolist_to_numpy(cutTaos);
    }

    std::vector<rxtools::StructureSummary> RxGamingRxUnit::get_simulated_structures(double bbDbh) const {
        try {
            auto align = lapis::Alignment((lapis::Extent)unitMask, 1, 1);
    
            std::vector<rxtools::StructureSummary> structures;
    
            std::deque<size_t> notBBbase;
            rxtools::TaoList base(taos.crs());
            for (size_t i = 0; i < taos.size(); i++) {
                if (taos.dbh(i) >= bbDbh) {
                    base.addTao(taos.xy(i),
                        taos.height(i),
                        taos.radius(i),
                        taos.area(i),
                        taos.dbh(i));
                }
                else
                {
                    notBBbase.push_back(i);
                }
            }
            auto testTaos = base;
            auto notBBidx = notBBbase;
    
            structures.push_back(rxtools::StructureSummary(testTaos, align, areaHa));
            size_t step = (taos.size() - testTaos.size()) / 10;
    
            //short to tall.
            for (int i = 1; i < 11; ++i) {
                auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
                char buf[26];
                ctime_s(buf, sizeof(buf), &t);
                std::cout << i << " " << buf << "\n";
                size_t limit = std::min<size_t>(taos.size(), testTaos.size() + step);
                while (testTaos.size() < limit && notBBidx.size()) {
                    auto idx = notBBidx.back();
                    notBBidx.pop_back();
                    testTaos.addTao(taos.xy(idx),
                        taos.height(idx),
                        taos.radius(idx),
                        taos.area(idx),
                        taos.dbh(idx));
                }
                structures.push_back(rxtools::StructureSummary(testTaos, align, areaHa));
            }
    
            //tall to short.
            testTaos = base;
            notBBidx = notBBbase;
            for (int i = 1; i < 11; ++i) {
                auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
                size_t limit = std::min<size_t>(taos.size(), testTaos.size() + step);
                while (testTaos.size() < limit && notBBidx.size()) {
                    auto idx = notBBidx.front();
                    notBBidx.pop_front();
                    testTaos.addTao(taos.xy(idx),
                        taos.height(idx),
                        taos.radius(idx),
                        taos.area(idx),
                        taos.dbh(idx));
                }
                structures.push_back(rxtools::StructureSummary(testTaos, align, areaHa));
            }
    
            //random
            testTaos = base;
            notBBidx = notBBbase;
            std::shuffle(std::begin(notBBidx), std::end(notBBidx), globalRng());
            for (int i = 1; i < 11; ++i) {
                auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
                size_t limit = std::min<size_t>(taos.size(), testTaos.size() + step);
                while (testTaos.size() < limit && notBBidx.size()) {
                    auto idx = notBBidx.back();
                    notBBidx.pop_back();
                    testTaos.addTao(taos.xy(idx),
                        taos.height(idx),
                        taos.radius(idx),
                        taos.area(idx),
                        taos.dbh(idx));
                }
                structures.push_back(rxtools::StructureSummary(testTaos, align, areaHa));
            }
            return structures;
        }
        catch (std::exception e) {
            std::cout << e.what() << "\n";
            std::abort();
        }
    }
        
    py::array_t<double> RxGamingRxUnit::taolist_to_numpy(rxtools::TaoList taos) const {
        py::array_t<double> arr({(py::ssize_t)taos.size(), (py::ssize_t)5});
        auto buf = arr.mutable_unchecked<2>();
        for (size_t i = 0; i < taos.size(); ++i) {
            buf(i, 0) = taos.xy(i).x;
            buf(i, 1) = taos.xy(i).y;
            buf(i, 2) = taos.height(i);
            buf(i, 3) = taos.radius(i);
            buf(i, 4) = taos.dbh(i);
        }
        return arr;
    }

    lapis::Raster<double> RxGamingRxUnit::computeHillshade(const lapis::Raster<double>& chm, double az, double elev) {
        az = (360-az+90) * M_PI / 180.0;
        elev = (90- elev) * M_PI / 180.0;
        
        lapis::Raster<double> hillshade((lapis::Alignment)chm);
        for (lapis::rowcol_t x = 0; x < chm.ncol(); ++x) {
            auto xl = std::max(x-1, 0);
            auto xr = std::min(x+1, chm.ncol()-1);
            for (lapis::rowcol_t y = 0; y < chm.nrow(); ++y) {
                if(!chm.atRCUnsafe(y, x).has_value()) {
                    continue;
                }

                auto yl = std::max(y-1, 0);
                auto yr = std::min(y+1, chm.nrow()-1);
                
                if(!chm.atRCUnsafe(yl, x).has_value() ||
                 !chm.atRCUnsafe(yr, x).has_value() ||
                 !chm.atRCUnsafe(y, xl).has_value() ||
                 !chm.atRCUnsafe(y, xr).has_value()
                ) {
                    continue;
                }

                auto sx = (chm.atRCUnsafe(y, xr).value() - chm.atRCUnsafe(y, xl).value()) / (2.0*chm.xres());
                auto sy = (chm.atRCUnsafe(yr, x).value() - chm.atRCUnsafe(yl, x).value()) / (2.0*chm.yres());

                auto asp_rad = std::atan2(sy, sx);
                auto s_mag_rad = std::atan(std::sqrt(sx*sx + sy*sy));
                hillshade.atRCUnsafe(y, x).value() = 255.0 * (std::cos(elev) * std::cos(s_mag_rad) + std::sin(elev) * std::sin(s_mag_rad) * std::cos(az - asp_rad));
                hillshade.atRCUnsafe(y, x).has_value() = true;
            }
        }
        return hillshade;
    }

    std::vector<size_t> RxGamingRxUnit::getRawClumps(rxtools::TaoList taos) const {
        lapis::lico::GraphLico g{ lapis::Alignment((lapis::Extent)unitMask, 1, 1) };
        for (size_t i = 0; i < taos.size(); ++i) {
            g.addTAO(taos.node(i), lapis::lico::NodeStatus::on);
        }
        std::vector<size_t> clumpSizes;
        for (size_t i = 0; i < g.nodes.size(); ++i) {
            clumpSizes.push_back(g.clumpSize(i));
        }
        return clumpSizes;
    }

    lapis::Raster<int> RxGamingRxUnit::getTreatBasin() const {
        std::unordered_set<int> basinIds;
        for (size_t i = 0; i < treatedTaos.size(); ++i) {
            auto v = basinMap.extract(treatedTaos.x(i), treatedTaos.y(i), lapis::ExtractMethod::near);
            if (v.has_value()) {
                basinIds.emplace(v.value());
            }
        }

        auto thisBasin = basinMap;
        for (lapis::cell_t i = 0; i < thisBasin.ncell(); ++i) {
            if (thisBasin[i].has_value()) {
                if (basinIds.find(basinMap[i].value()) == basinIds.end()) {
                    thisBasin[i].value() = 1;
                }
            }
        }

        return thisBasin;
    }
}
