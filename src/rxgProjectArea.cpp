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
        auto unitPolygon = lapis::VectorDataset<lapis::MultiPolygon>(ps.unitPolyPath);
        auto lidar = processedfolder::readProcessedFolder(ps.lidarPath);

        unitPolygon.projectInPlace(lidar->crs());

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

        reader.makePlotTreeMap(std::vector<std::string>{ "DIA" });
        auto allTrees = reader.collapsePlotTreeMap();
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

        auto l = lapis::VectorDataset<lapis::MultiPolygon>(processedfolder::stringOrThrow(lidar->tileLayoutVector()));
        auto mask = lapis::Raster<lapis::cell_t>(processedfolder::stringOrThrow(lidar->maskRaster()));

        #pragma omp parallel
        {
            std::vector<RxGamingRxUnit> localUnits;

            #pragma omp for nowait
            for(int i = 0; i < unitPolygon.nFeature(); ++i) {
                auto unit = unitPolygon.getFeature(i);
                auto name = unit.getStringField(ps.unitName);
                std::vector<lapis::Raster<int>> mhm;
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
                                    auto thisMhm = lapis::Raster<int>(processedfolder::stringOrThrow(lidar->maxHeightRaster(j)));
                                    //thisMhm.repairResolution(expectedRes.first, expectedRes.first, expectedRes.second, expectedRes.second);
                                    auto thisBasinMap = lapis::Raster<int>(processedfolder::stringOrThrow(lidar->watershedSegmentRaster(j))) + i * 10000;
                                    //thisBasinMap.repairResolution(expectedRes.first, expectedRes.first, expectedRes.second, expectedRes.second);

                                    mhm.push_back(thisMhm);
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
                                            if (unit.getGeometry().containsPoint(lapis::Point(xy, thisMhm.crs()))) {
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

                        auto mergeVectorInt = [](std::vector<lapis::Raster<int>>& v) {
                            std::vector<lapis::Raster<int>*> toMerge;
                            for (int i = 0; i < v.size(); i++) {
                                toMerge.push_back(&v[i]);
                            }
                            return lapis::mosaicInside(toMerge);
                        };
                        auto outMhm = mergeVectorInt(mhm);
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
                        thisMask = lapis::cropRaster(thisMask, outMhm, lapis::SnapType::out);
                        thisMask.maskByMultiPolygon(unit.getGeometry());
                        thisMask = lapis::trimRaster(thisMask);
                        
                        auto tmp = outMhm;
                        for (lapis::cell_t cell = 0; cell < outMhm.ncell(); ++cell) {
                            auto x = outMhm.xFromCellUnsafe(cell);
                            auto y = outMhm.yFromCellUnsafe(cell);
                            if(((lapis::Extent)thisMask).contains(x,y)) {
                                if (!thisMask.atXY(outMhm.xFromCellUnsafe(cell),outMhm.yFromCellUnsafe(cell)).has_value()) {
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
                        outMhm.mask(tmp);
                        outChm.mask(tmp);
                        outBasin.mask(tmp);

                        outMhm = lapis::trimRaster(outMhm);
                        outChm = lapis::trimRaster(outChm);
                        outBasin = lapis::trimRaster(outBasin);

                        outMhm = lapis::cropRaster(outMhm, outChm, lapis::SnapType::out);
                        outMhm = lapis::extendRaster(outMhm, outChm, lapis::SnapType::in);
                        outBasin = lapis::cropRaster(outBasin, outChm, lapis::SnapType::out);
                        outBasin = lapis::extendRaster(outBasin, outChm, lapis::SnapType::in);

                        localUnits.emplace_back(name, thisMask, outMhm, outChm, outBasin, taos);
                    }
                    else {
                        continue;
                    }
                } // try:
                catch (processedfolder::FileNotFoundException e) {
                    continue;
                }
            } // for(lapis::ConstFeature<lapis::MultiPolygon> unit : unitPolygon)

            #pragma omp critical
            {
                rxUnits.insert(rxUnits.end(), localUnits.begin(), localUnits.end());
            }
        } // #pragma omp parallel
    }
}