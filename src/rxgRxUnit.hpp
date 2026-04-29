/*
Copyright(C) 2024  University of Washington
This program is free software : you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.If not, see https ://www.gnu.org/licenses/.

Bryce Bartl - Geller
University of Washington Forest Resilience Lab
2 / 22 / 2026

rxgRxUnit.hpp
*/

#pragma once
#include "rxgUtils.hpp"
#include "rxunit.hpp"
#include "treatment.hpp"
#include <pybind11/numpy.h>

namespace rxgaming {
    namespace py = pybind11;

    class RxGamingRxUnit : public rxtools::RxUnit {
    public:
        std::string name;
        lapis::Raster<double> chm;
        lapis::Raster<int> basinMap;
        lapis::Raster<double> hillshade;
        rxtools::TaoList cutTaos;
        rxtools::treatmentResult result = rxtools::treatmentResult::success;

        RxGamingRxUnit(std::string name, lapis::Raster<lapis::cell_t> mask, lapis::Raster<double> chm, lapis::Raster<int> basinMap, rxtools::TaoList taos);

        RxGamingRxUnit() = default;

        py::array_t<lapis::cell_t> get_mask() const;
        py::array_t<double> get_chm() const;
        py::array_t<int> get_basin() const;
        py::array_t<double> get_hillshade() const;
        py::array_t<int> get_clump_map() const;
        py::array_t<double> get_taos() const;
        py::array_t<uint64_t> get_clump_sizes() const;

        py::array_t<double> get_treat_chm() const;
        py::array_t<int> get_treat_basin() const;
        py::array_t<double> get_treat_hillshade() const;
        py::array_t<int> get_treat_clump_map() const;
        py::array_t<double> get_treat_taos() const;
        py::array_t<uint64_t> get_treat_clump_sizes() const;

        py::array_t<double> get_cut_taos() const;

        void export_rendered_geotiff(
            const std::string& outputPath,
            const py::array_t<std::uint8_t, py::array::c_style | py::array::forcecast>& image,
            int mapLeftPx,
            int mapTopPx,
            int mapWidthPx,
            int mapHeightPx
        ) const;
        void write_chm_raster(const std::string& outputPath, bool treated) const;
        void write_basin_raster(const std::string& outputPath, bool treated) const;
        void write_clumpmap_raster(const std::string& outputPath, bool treated) const;

        std::vector<rxtools::StructureSummary> get_simulated_structures(double bbDbh) const;
        void refresh_derived_state();

    private:
        template <class T>
        py::array_t<T> raster_to_numpy(const lapis::Raster<T>& r) const {
            py::array_t<T> arr({r.nrow(), r.ncol()});
            auto buf = arr.template mutable_unchecked<2>();

            for (lapis::rowcol_t i = 0; i < r.nrow(); ++i) {
                for (lapis::rowcol_t j = 0; j < r.ncol(); ++j) {
                    auto v = r.atRCUnsafe(i, j);
                    if (v.has_value()) {
                        buf(i, j) = v.value();
                    } else {
                        if constexpr (std::is_floating_point_v<T>) {
                            buf(i, j) = std::numeric_limits<T>::quiet_NaN();
                        } else {
                            buf(i, j) = std::numeric_limits<T>::min();
                        }
                    }
                }
            }
            return arr;
        }

        py::array_t<double> taolist_to_numpy(rxtools::TaoList taos) const;
        lapis::Raster<double> computeHillshade(const lapis::Raster<double>& chm, double az = 315, double elev = 45);
        std::vector<size_t> getRawClumps(rxtools::TaoList taos) const;
        lapis::Raster<double> getTreatChmRaster() const;
        lapis::Raster<int> getClumpMapRaster(rxtools::TaoList taos, const lapis::Raster<int>& sourceBasin) const;
        lapis::Raster<int> getTreatBasin() const;
    };
} // namespace rxgaming
