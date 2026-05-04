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
#include <atomic>
#include <iostream>
#include <sstream>
#include <stdexcept>

namespace rxgaming {
    RxGamingProjectArea::RxGamingProjectArea(std::vector<RxGamingRxUnit> units) : rxUnits(std::move(units)) {}

    namespace {
        void emit_progress(
            const ProgressCallback& progressCallback,
            const std::string& stage,
            const std::string& message,
            int unitIndex = -1,
            const std::string& unitName = std::string(),
            int completed = -1,
            int total = -1)
        {
            if (!progressCallback) {
                return;
            }
            progressCallback(ProgressEvent{stage, message, unitIndex, unitName, completed, total});
        }
    }

    RxGamingProjectArea::RxGamingProjectArea(const ProjectSettings& ps, const ProgressCallback& progressCallback) {
        emit_progress(progressCallback, "load_units", "Loading unit polygons from " + ps.unitPolyPath);
        auto unitPolygon = lapis::VectorDataset<lapis::MultiPolygon>(ps.unitPolyPath);
        emit_progress(progressCallback, "load_lidar", "Loading lidar dataset from " + ps.lidarPath);
        auto lidar = processedfolder::readProcessedFolder(ps.lidarPath);

        emit_progress(progressCallback, "project_units", "Projecting units into lidar CRS");
        unitPolygon.projectInPlace(lidar->crs());

        emit_progress(progressCallback, "read_fia", "Reading FIA data from " + ps.fiaPath);
        auto reader = rxtools::allometry::FIAReader(ps.fiaPath);
        auto dist = lidar->units().value().convertOneToThis(10000, lapis::linearUnitPresets::meter);
        lapis::Extent e(
            unitPolygon.extent().xmin() - dist,
            unitPolygon.extent().xmax() + dist, 
            unitPolygon.extent().ymin() - dist,
            unitPolygon.extent().ymax() + dist,
            unitPolygon.crs());
        auto n = reader.limitByExtent(e);
        if (!n) {
            throw std::runtime_error("No FIA plots in buffered extent of project area.");
        }

        emit_progress(progressCallback, "build_fia_map", "Building FIA plot tree map");
        reader.makePlotTreeMap(std::vector<std::string>{ "DIA" });
        auto allTrees = reader.collapsePlotTreeMap();
        emit_progress(progressCallback, "fit_dbh_model", "Fitting DBH model");
        auto dbhModel = rxtools::allometry::DbhModel(allTrees.height, allTrees.get("DIA"), lapis::linearUnitPresets::internationalFoot, lapis::linearUnitPresets::internationalInch);

        std::pair<lapis::coord_t, lapis::coord_t> expectedRes{};
        auto units = lidar->units();
        if (lidar->type() == processedfolder::RunType::fusion) {
            if(units->name() == "metre") {
                expectedRes.first = 0.75;
            }
            else {
                expectedRes.first = 2.4606;
            }
        }
        else {
            expectedRes.first = lidar->csmAlignment()->xres();
        }
        expectedRes.second = 0;

        emit_progress(progressCallback, "load_mask", "Loading shared mask raster");
        auto l = lapis::VectorDataset<lapis::MultiPolygon>(processedfolder::stringOrThrow(lidar->tileLayoutVector()));
        auto mask = lapis::Raster<lapis::cell_t>(processedfolder::stringOrThrow(lidar->maskRaster()));
        const int totalUnits = (int)unitPolygon.nFeature();
        emit_progress(
            progressCallback,
            "process_units",
            "Processing " + std::to_string(totalUnits) + " units across " + std::to_string(lidar->nTiles()) + " tiles",
            -1,
            std::string(),
            0,
            totalUnits);

        bool nameFieldExists = false;
        if(std::find(
            unitPolygon.getAllFieldNames().begin(),
            unitPolygon.getAllFieldNames().end(),
            ps.unitName) != unitPolygon.getAllFieldNames().end()) 
        {
            nameFieldExists = true;
        }

        std::atomic<int> completedUnits(0);
        std::cout << "NUMBER OF THREADS: " << ps.nThread << "\n";
        #pragma omp parallel num_threads(ps.nThread)
        {
            #pragma omp single
            {
                std::cout
                    << "omp_in_parallel=" << omp_in_parallel()
                    << ", team_size=" << omp_get_num_threads()
                    << ", requested=" << ps.nThread
                    << '\n';
            }
            std::vector<RxGamingRxUnit> localUnits;

            #pragma omp for nowait
            for(int i = 0; i < unitPolygon.nFeature(); ++i) {
                int thisThreadNum =  omp_get_thread_num();
                auto unit = unitPolygon.getFeature(i);
                std::string name;
                if(nameFieldExists) {
                    name = unit.getStringField(ps.unitName);
                } else {
                    std::ostringstream nameFallback;
                    nameFallback << "unit_" << i;
                    name = nameFallback.str();
                }
                
                emit_progress(
                    progressCallback,
                    "unit_start",
                    "[Thread " + std::to_string(thisThreadNum) + "] Starting unit " + std::to_string(i) + " (" + name + ")",
                    i,
                    name,
                    completedUnits.load(),
                    totalUnits);
                std::vector<lapis::Raster<double>> chm;
                std::vector<lapis::Raster<int>> basinMap;
                std::unordered_set<std::string> usedTiles;
                rxtools::TaoList taos(mask.crs());
                try {
                    if (mask.dataOverlapsMultiPolygon(unit.getGeometry())) {
                        std::cout << "[ Thread" + std::to_string(thisThreadNum) + "] Unit: " + std::to_string(i) + "\n";
                        for (int j = 0; j < lidar->nTiles(); j++) {
                            auto e = lidar->extentByTile(j);
                            if (e) {
                                if (unit.getGeometry().overlapsExtent(e.value())) {
                                    auto thisBasinMap = lapis::Raster<int>(processedfolder::stringOrThrow(lidar->watershedSegmentRaster(j))) + i * 10000;
                                    //thisBasinMap.repairResolution(expectedRes.first, expectedRes.first, expectedRes.second, expectedRes.second);

                                    basinMap.push_back(thisBasinMap);

                                    auto s = processedfolder::stringOrThrow(lidar->csmRaster(j));
                                    if (std::find(usedTiles.begin(), usedTiles.end(), s) == usedTiles.end()) {
                                        auto thisChm = lapis::Raster<double>(s);
                                        //thisChm.repairResolution(expectedRes.first, expectedRes.first, expectedRes.second, expectedRes.second);
                                        chm.push_back(thisChm);
                                        usedTiles.emplace(s);
                                    }

                                    std::unordered_map<int, int> basinIdToNcells;
                                    for (lapis::cell_t c = 0; c < thisBasinMap.ncell(); ++c) {
                                        if(thisBasinMap[c].has_value()) {
                                            auto val = thisBasinMap[c].value();
                                            if(basinIdToNcells.find(val) == basinIdToNcells.end()) {
                                                basinIdToNcells[val] = 1;
                                            }
                                            else {
                                                basinIdToNcells[val]++;
                                            }
                                        }
                                    }

                                    lapis::VectorDataset<lapis::Point> pts(lidar->highPoints(j)->string());
                                    auto unitE = unit.getGeometry().boundingBox();
                                    for (int k = 0; k < pts.nFeature(); k++) {
                                        auto xy = lidar->coordGetter()(pts.getFeature(k));
                                        if (unitE.contains(xy.x, xy.y)) {
                                            if (unit.getGeometry().containsPoint(lapis::Point(xy, thisBasinMap.crs()))) {
                                                auto val = thisBasinMap.atXY(xy.x, xy.y);
                                                if (!val.has_value()) {
                                                    continue;
                                                }
                                                double a = (double)basinIdToNcells[val.value()] * thisBasinMap.xres() * thisBasinMap.yres();
                                                auto height = lidar->heightGetter()(pts.getFeature(k));
                                                //at some point the vert unit for dbh should not be hard coded.
                                                taos.addTao(xy,
                                                    height,
                                                    3,
                                                    a,
                                                    dbhModel.predict(height, lapis::linearUnitPresets::meter, lapis::linearUnitPresets::centimeter)
                                                );
                                            }
                                        }
                                    }
                                }
                            }
                        } // for(int j = 0;...)

                        if (basinMap.empty() || chm.empty()) {
                            int completed = completedUnits.fetch_add(1) + 1;
                            emit_progress(
                                progressCallback,
                                "unit_skipped",
                                "[Thread " + std::to_string(thisThreadNum) + "] Skipped unit " + std::to_string(i) + " (" + name + ") because no overlapping lidar rasters were collected",
                                i,
                                name,
                                completed,
                                totalUnits);
                            continue;
                        }

                        auto mergeVectorInt = [](std::vector<lapis::Raster<int>>& v) {
                            std::vector<lapis::Raster<int>*> toMerge;
                            for (int i = 0; i < v.size(); i++) {
                                toMerge.push_back(&v[i]);
                            }
                            return lapis::mosaicInside(toMerge);
                        };
                        auto outBasin = mergeVectorInt(basinMap);

                        auto mergeVectorDouble = [](std::vector<lapis::Raster<double>>& v) {
                            std::vector<lapis::Raster<double>*> toMerge;
                            for (int i = 0; i < v.size(); i++) {
                                toMerge.push_back(&v[i]);
                            }
                            return lapis::mosaicInside(toMerge);
                        };
                        auto outChm = mergeVectorDouble(chm);

                        auto thisMask = mask;
                        thisMask = lapis::cropRaster(thisMask, outBasin, lapis::SnapType::out);
                        thisMask.maskByMultiPolygon(unit.getGeometry());
                        thisMask = lapis::trimRaster(thisMask);
                        
                        auto tmp = outBasin;
                        for (lapis::cell_t cell = 0; cell < outBasin.ncell(); ++cell) {
                            auto x = outBasin.xFromCellUnsafe(cell);
                            auto y = outBasin.yFromCellUnsafe(cell);
                            if(((lapis::Extent)thisMask).contains(x,y)) {
                                if (!thisMask.atXY(outBasin.xFromCellUnsafe(cell),outBasin.yFromCellUnsafe(cell)).has_value()) {
                                    tmp.atCell(cell).has_value() = false;
                                }
                                else {
                                    tmp.atCell(cell).has_value() = true;
                                }
                            }
                            else {
                                tmp.atCell(cell).has_value() = false;
                            }
                        }
                        outChm.mask(tmp);
                        outBasin.mask(tmp);

                        outChm = lapis::trimRaster(outChm);
                        outBasin = lapis::trimRaster(outBasin);

                        outBasin = lapis::cropRaster(outBasin, outChm, lapis::SnapType::out);
                        outBasin = lapis::extendRaster(outBasin, outChm, lapis::SnapType::in);

                        localUnits.emplace_back(name, thisMask, outChm, outBasin, taos);
                        int completed = completedUnits.fetch_add(1) + 1;
                        emit_progress(
                            progressCallback,
                            "unit_complete",
                            "[Thread " + std::to_string(thisThreadNum) + "] Completed unit " + std::to_string(i) + " (" + name + ")",
                            i,
                            name,
                            completed,
                            totalUnits);
                    }
                    else {
                        int completed = completedUnits.fetch_add(1) + 1;
                        emit_progress(
                            progressCallback,
                            "unit_skipped",
                            "[Thread " + std::to_string(thisThreadNum) + "] Skipped unit " + std::to_string(i) + " (" + name + ") because the project mask does not overlap",
                            i,
                            name,
                            completed,
                            totalUnits);
                        continue;
                    }
                } // try:
                catch (const processedfolder::FileNotFoundException& e) {
                    int completed = completedUnits.fetch_add(1) + 1;
                    emit_progress(
                        progressCallback,
                        "unit_skipped",
                        "[Thread " + std::to_string(thisThreadNum) + "] Skipped unit " + std::to_string(i) + " (" + name + ") due to missing file: " + std::string(e.what()),
                        i,
                        name,
                        completed,
                        totalUnits);
                    continue;
                }
                catch (const std::exception& e) {
                    emit_progress(
                        progressCallback,
                        "unit_failed",
                        "[Thread " + std::to_string(thisThreadNum) + "] Unit " + std::to_string(i) + " (" + name + ") failed: " + std::string(e.what()),
                        i,
                        name,
                        completedUnits.load(),
                        totalUnits);
                    throw;
                }
                catch (...) {
                    emit_progress(
                        progressCallback,
                        "unit_failed",
                        "[Thread " + std::to_string(thisThreadNum) + "] Unit " + std::to_string(i) + " (" + name + ") failed with an unknown exception",
                        i,
                        name,
                        completedUnits.load(),
                        totalUnits);
                    throw;
                }
            } // for(lapis::ConstFeature<lapis::MultiPolygon> unit : unitPolygon)

            #pragma omp critical
            {
                rxUnits.insert(rxUnits.end(), localUnits.begin(), localUnits.end());
            }
        } // #pragma omp parallel
        emit_progress(
            progressCallback,
            "finished",
            "Finished building project area with " + std::to_string(rxUnits.size()) + " units",
            -1,
            std::string(),
            totalUnits,
            totalUnits);
    }
}
