#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "projectSettings.hpp"
#include "rxgRxUnit.hpp"
#include "rxgTreatmentEngine.hpp"

namespace py = pybind11;

/*
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
        */

PYBIND11_MODULE(rxgaming_core, m) {
    m.def("set_seed", &rxgaming::set_seed, "Set the random seed");
    
    py::class_<rxtools::StructureSummary>(m, "StructureSummary")
        .def_readwrite("ba", &rxtools::StructureSummary::ba)
        .def_readwrite("tph", &rxtools::StructureSummary::tph)
        .def_readwrite("mcs", &rxtools::StructureSummary::mcs)
        .def_readwrite("cc", &rxtools::StructureSummary::cc)

        .def_readwrite("csd", &rxtools::StructureSummary::csd)
        .def_readwrite("binMins", &rxtools::StructureSummary::binMins)
        .def_readwrite("binMaxs", &rxtools::StructureSummary::binMaxs)

        .def("__getitem__", [](const rxtools::StructureSummary& s, int i) {
            if (i < 0 || i > 3) throw py::index_error();
            return s[i];
        });

    py::class_<rxgaming::ProjectSettings>(m, "ProjectSettings")
        .def(py::init<std::string, std::string, std::string, std::string, std::string,
             std::string, std::string, std::string, std::string, int>(),
             py::arg("name"), py::arg("unitPolyPath"), py::arg("refDataPath"), py::arg("mcsPropPath"), py::arg("fiaPath"),
             py::arg("projDbPath"), py::arg("lidarPath"), py::arg("unitName"), py::arg("savePath"),py::arg("nThread"))
        .def_readonly("name", &rxgaming::ProjectSettings::name)
        .def_readonly("unitPolyPath", &rxgaming::ProjectSettings::unitPolyPath)
        .def_readonly("refDataPath", &rxgaming::ProjectSettings::refDataPath)
        .def_readonly("mcsPropPath", &rxgaming::ProjectSettings::mcsPropPath)
        .def_readonly("fiaPath", &rxgaming::ProjectSettings::fiaPath)
        .def_readonly("projDbPath", &rxgaming::ProjectSettings::projDbPath)
        .def_readonly("lidarPath", &rxgaming::ProjectSettings::lidarPath)
        .def_readonly("unitName", &rxgaming::ProjectSettings::unitName)
        .def_readonly("nThread", &rxgaming::ProjectSettings::nThread);
    
    py::class_<rxgaming::RxGamingRxUnit>(m, "RxUnit")
        .def_readonly("name", &rxgaming::RxGamingRxUnit::name)
        .def_readonly("areaHa", &rxgaming::RxGamingRxUnit::areaHa)
        .def_readwrite("currentStructure", &rxgaming::RxGamingRxUnit::currentStructure)
        .def_readwrite("targetStructure", &rxgaming::RxGamingRxUnit::targetStructure)
        .def_readwrite("treatedStructure", &rxgaming::RxGamingRxUnit::treatedStructure)

        .def("get_mask", &rxgaming::RxGamingRxUnit::get_mask)
        .def("get_chm", &rxgaming::RxGamingRxUnit::get_chm)
        .def("get_basin", &rxgaming::RxGamingRxUnit::get_basin)
        .def("get_hillshade", &rxgaming::RxGamingRxUnit::get_hillshade)
        .def("get_clump_map", &rxgaming::RxGamingRxUnit::get_clump_map)
        .def("get_taos", &rxgaming::RxGamingRxUnit::get_taos)

        .def("get_treat_chm", &rxgaming::RxGamingRxUnit::get_treat_chm)
        .def("get_treat_basin", &rxgaming::RxGamingRxUnit::get_treat_basin)
        .def("get_treat_hillshade", &rxgaming::RxGamingRxUnit::get_treat_hillshade)
        .def("get_treat_clump_map", &rxgaming::RxGamingRxUnit::get_treat_clump_map)
        .def("get_treat_taos", &rxgaming::RxGamingRxUnit::get_treat_taos)
        
        .def("get_cut_taos", &rxgaming::RxGamingRxUnit::get_cut_taos)

        .def("get_simulated_structures", &rxgaming::RxGamingRxUnit::get_simulated_structures);

    py::class_<rxgaming::TreatmentEngine>(m, "TreatmentEngine")
        .def(py::init<>())
        .def("do_treatment", &rxgaming::TreatmentEngine::do_treatment);

}