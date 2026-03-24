#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "ProjectSettings.h"
#include "RxUnit.h"

namespace py = pybind11;

PYBIND11_MODULE(rxgaming_core, m) {
    py::class_<ProjectSettings>(m, "ProjectSettings")
        .def(py::init<std::string, std::string, std::string,
                      std::string, std::string, int>())
        .def("get_name", &ProjectSettings::getName)
        .def("get_units", &ProjectSettings::getUnits)
        .def("save", &ProjectSettings::save)
        .def("is_cache_valid", &ProjectSettings::isCacheValid)
        .def("load_from_cache", &ProjectSettings::loadCache)
        .def("recompute_and_cache", &ProjectSettings::recomputeAndCache)
        .def("get_mcs_prop", &ProjectSettings::getMcsProp);

    py::class_<RxUnit>(m, "RxUnit")
        .def("get_name", &RxUnit::getName)
        .def("get_chm", &RxUnit::getChmNumpy)        // returns np.ndarray
        .def("get_tao_points", &RxUnit::getTaoNumpy) // returns np.ndarray
        .def("get_basin_map", &RxUnit::getBasinNumpy)
        .def("get_clump_map", &RxUnit::getClumpNumpy)
        .def("get_hillshade", &RxUnit::getHillshadeNumpy)
        .def("get_current_structure", &RxUnit::getCurrentStructure)
        .def("get_target_structure", &RxUnit::getTargetStructure)
        .def("set_target_structure", &RxUnit::setTargetStructure)
        .def("do_treatment", &RxUnit::doTreatment)
        .def("get_treat_chm", &RxUnit::getTreatChmNumpy)
        .def("get_treat_taos", &RxUnit::getTreatTaosNumpy)
        .def("get_simulated_structures", &RxUnit::getSimulatedStructures);

    py::class_<StructureSummary>(m, "StructureSummary")
        .def_readonly("tpa", &StructureSummary::tpa)
        .def_readonly("ba", &StructureSummary::ba)
        .def_readonly("mcs", &StructureSummary::mcs)
        .def_readonly("cc", &StructureSummary::cc);
}