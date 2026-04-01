/*
Copyright(C) 2024  University of Washington
This program is free software : you can redistribute it and /or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.If not, see https ://www.gnu.org/licenses/.

Bryce Bartl - Geller
University of Washington Forest Resilience Lab
2 / 22 / 2026

rxgTreatmentEngine.hpp
*/

#pragma once
#include "treatment.hpp"
#include "rxgRxUnit.hpp"

namespace rxgaming {
    class TreatmentEngine {
    public:
        TreatmentEngine() : treater() {};

        inline void do_treatment(RxgRxUnit& unit, double dbhMin, double dbhMax) {
            try {
                auto trt =  treater.doTreatment(unit, dbhMin, dbhMax, 10, false, "E:/dropbox/rxgaming paper/treatmentvisual/data/");
                unit.treatedTaos = std::get<0>(trt);
                unit.cutTaos = std::get<1>(trt);
                unit.result = std::get<2>(trt);
                unit.treatedStructure = rxtools::StructureSummary(
                    unit.treatedTaos,
                    lapis::Alignment((lapis::Extent)unit.unitMask, 1, 1),
                    unit.areaHa);
            }
            catch (std::exception e) {
                std::cout << e.what();
                std::abort();
            }
        }
    private:
        rxtools::Treatment treater;
    };
} // namespace rxgaming