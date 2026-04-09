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
#include <iostream>
#include <sstream>

namespace rxgaming {

    RxGamingProjectArea::RxGamingProjectArea(const ProjectSettings& ps) {
        std::cout << "RxGamingProjectArea: loading unit polygons from " << ps.unitPolyPath << "\n" << std::flush;
        auto unitPolygon = lapis::VectorDataset<lapis::MultiPolygon>(ps.unitPolyPath);
        std::cout << "RxGamingProjectArea: loading lidar dataset from " << ps.lidarPath << "\n" << std::flush;
        auto lidar = processedfolder::readProcessedFolder(ps.lidarPath);

        std::cout << "RxGamingProjectArea: projecting units into lidar CRS\n" << std::flush;
        unitPolygon.projectInPlace(lidar->crs());

        std::cout << "RxGamingProjectArea: reading FIA data from " << ps.fiaPath << "\n" << std::flush;
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
            std::cout << "no fia plots in buffered extent of project area; aborting";
            std::abort();
        }

        std::cout << "RxGamingProjectArea: building FIA plot tree map\n" << std::flush;
        reader.makePlotTreeMap(std::vector<std::string>{ "DIA" });
        auto allTrees = reader.collapsePlotTreeMap();
        std::cout << "RxGamingProjectArea: fitting DBH model\n" << std::flush;
        auto dbhModel = rxtools::allometry::DbhModel(allTrees.height, allTrees.get("DIA"), lapis::linearUnitPresets::internationalFoot, lapis::linearUnitPresets::internationalInch);

        std::pair<lapis::coord_t, lapis::coord_t> expectedRes{};
        auto units = lidar->units();
        if (lidar->type() == processedfolder::RunType::fusion) {
            std::cout << "b\n";
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

        std::cout << "RxGamingProjectArea: loading shared mask raster\n" << std::flush;
        auto l = lapis::VectorDataset<lapis::MultiPolygon>(processedfolder::stringOrThrow(lidar->tileLayoutVector()));
        auto mask = lapis::Raster<lapis::cell_t>(processedfolder::stringOrThrow(lidar->maskRaster()));
        std::cout << "RxGamingProjectArea: processing " << unitPolygon.nFeature() << " units across " << lidar->nTiles() << " tiles\n" << std::flush;

        bool nameFieldExists = false;
        if(std::find(
            unitPolygon.getAllFieldNames().begin(),
            unitPolygon.getAllFieldNames().end(),
            ps.unitName) != unitPolygon.getAllFieldNames().end()) 
        {
            nameFieldExists = true;
        }

        #pragma omp parallel
        {
            std::vector<RxGamingRxUnit> localUnits;

            #pragma omp for nowait
            for(int i = 0; i < unitPolygon.nFeature(); ++i) {
                std::cout << "RxGamingProjectArea: loop " << i << std::flush;
                auto unit = unitPolygon.getFeature(i);
                std::cout << "RxGamingProjectArea: got  feature " << i << std::flush;
                std::string name;
                if(nameFieldExists) {
                    name = unit.getStringField(ps.unitName);
                } else {
                    std::ostringstream nameFallback;
                    nameFallback << "unit_" << i;
                    name = nameFallback.str();
                }
                
                std::cout << "RxGamingProjectArea: unit " << i << " (" << name << ") start\n" << std::flush;
                std::vector<lapis::Raster<double>> chm;
                std::vector<lapis::Raster<int>> basinMap;
                std::unordered_set<std::string> usedTiles;
                rxtools::TaoList taos(mask.crs());
                try {
                    if (mask.dataOverlapsMultiPolygon(unit.getGeometry())) {
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

                        
                        
                        std::cout
                            << "RxGamingProjectArea: unit " << i
                            << " collected " << basinMap.size() << " basin rasters, "
                            << chm.size() << " CHM rasters, "
                            << taos.size() << " TAOs\n"
                            << std::flush;

                        if (basinMap.empty() || chm.empty()) {
                            std::cout
                                << "RxGamingProjectArea: unit " << i
                                << " skipped because no overlapping lidar rasters were collected\n"
                                << std::flush;
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
                        std::cout << "RxGamingProjectArea: unit " << i << " complete\n" << std::flush;
                    }
                    else {
                        std::cout << "RxGamingProjectArea: unit " << i << " skipped because project mask does not overlap\n" << std::flush;
                        continue;
                    }
                } // try:
                catch (const processedfolder::FileNotFoundException& e) {
                    std::cout << "RxGamingProjectArea: unit " << i << " skipped due to missing file: " << e.what() << "\n" << std::flush;
                    continue;
                }
                catch (const std::exception& e) {
                    std::cout << "RxGamingProjectArea: unit " << i << " failed with std::exception: " << e.what() << "\n" << std::flush;
                    throw;
                }
                catch (...) {
                    std::cout << "RxGamingProjectArea: unit " << i << " failed with unknown exception\n" << std::flush;
                    continue;
                }
            } // for(lapis::ConstFeature<lapis::MultiPolygon> unit : unitPolygon)

            #pragma omp critical
            {
                rxUnits.insert(rxUnits.end(), localUnits.begin(), localUnits.end());
            }
        } // #pragma omp parallel
        std::cout << "RxGamingProjectArea: finished with " << rxUnits.size() << " units\n" << std::flush;
    }
}
