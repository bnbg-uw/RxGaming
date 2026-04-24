#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <mutex>
#include "projectSettings.hpp"
#include "rxgProjectArea.hpp"
#include "rxgSnapshot.hpp"
#include "rxgTreatmentEngine.hpp"

namespace py = pybind11;

namespace {
    rxgaming::RxGamingProjectArea build_project_area_with_progress(const rxgaming::ProjectSettings& ps, py::object callback) {
        if (callback.is_none()) {
            py::gil_scoped_release release;
            return rxgaming::RxGamingProjectArea(ps);
        }

        std::mutex callbackMutex;
        rxgaming::ProgressCallback progressCallback = [callback, &callbackMutex](const rxgaming::ProgressEvent& event) {
            std::lock_guard<std::mutex> lock(callbackMutex);
            py::gil_scoped_acquire gil;
            callback(event);
        };

        py::gil_scoped_release release;
        return rxgaming::RxGamingProjectArea(ps, progressCallback);
    }
}

PYBIND11_MODULE(rxgaming_core, m) {
    m.def("set_proj_db_path", &rxgaming::set_proj_db_path, "Set the PROJ data directory path");
    m.def("set_seed", &rxgaming::set_seed, "Set the random seed");
    m.def(
        "build_project_area_with_progress",
        &build_project_area_with_progress,
        py::arg("ps"),
        py::arg("callback") = py::none(),
        "Build a ProjectArea while reporting structured progress events.");
    m.def(
        "save_project_area",
        &rxgaming::save_project_area,
        py::arg("projectArea"),
        py::arg("path"),
        "Save a fully self-contained ProjectArea snapshot.");
    m.def(
        "load_project_area",
        &rxgaming::load_project_area,
        py::arg("path"),
        "Load a ProjectArea snapshot without requiring source lidar data.");
    
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

    py::enum_<rxtools::treatmentResult>(m, "TreatmentResult")
        .value("success", rxtools::treatmentResult::success)
        .value("diameterFailure", rxtools::treatmentResult::diameterFailure)
        .value("cuttingFailure", rxtools::treatmentResult::cuttingFailure);

    py::class_<rxgaming::ProgressEvent>(m, "ProgressEvent")
        .def_readonly("stage", &rxgaming::ProgressEvent::stage)
        .def_readonly("message", &rxgaming::ProgressEvent::message)
        .def_readonly("unitIndex", &rxgaming::ProgressEvent::unitIndex)
        .def_readonly("unitName", &rxgaming::ProgressEvent::unitName)
        .def_readonly("completed", &rxgaming::ProgressEvent::completed)
        .def_readonly("total", &rxgaming::ProgressEvent::total);

    py::class_<rxgaming::ProjectSettings>(m, "ProjectSettings")
        .def(py::init<std::string, std::string, std::string, std::string,
             std::string, std::string, std::string, std::string, int>(),
             py::arg("name"), py::arg("unitPolyPath"), py::arg("refDataPath"), py::arg("mcsPropPath"), py::arg("fiaPath"),
             py::arg("lidarPath"), py::arg("unitName"), py::arg("savePath"),py::arg("nThread"))
        .def_readonly("name", &rxgaming::ProjectSettings::name)
        .def_readonly("unitPolyPath", &rxgaming::ProjectSettings::unitPolyPath)
        .def_readonly("refDataPath", &rxgaming::ProjectSettings::refDataPath)
        .def_readonly("mcsPropPath", &rxgaming::ProjectSettings::mcsPropPath)
        .def_readonly("fiaPath", &rxgaming::ProjectSettings::fiaPath)
        .def_readonly("lidarPath", &rxgaming::ProjectSettings::lidarPath)
        .def_readonly("unitName", &rxgaming::ProjectSettings::unitName)
        .def_readonly("savePath", &rxgaming::ProjectSettings::savePath)
        .def_readonly("nThread", &rxgaming::ProjectSettings::nThread);
    
    py::class_<rxgaming::RxGamingRxUnit>(m, "RxUnit")
        .def_readonly("name", &rxgaming::RxGamingRxUnit::name)
        .def_readonly("areaHa", &rxgaming::RxGamingRxUnit::areaHa)
        .def_readwrite("result", &rxgaming::RxGamingRxUnit::result)
        .def_readwrite("currentStructure", &rxgaming::RxGamingRxUnit::currentStructure)
        .def_readwrite("targetStructure", &rxgaming::RxGamingRxUnit::targetStructure)
        .def_readwrite("treatedStructure", &rxgaming::RxGamingRxUnit::treatedStructure)

        .def("get_mask", &rxgaming::RxGamingRxUnit::get_mask)
        .def("get_chm", &rxgaming::RxGamingRxUnit::get_chm)
        .def("get_basin", &rxgaming::RxGamingRxUnit::get_basin)
        .def("get_hillshade", &rxgaming::RxGamingRxUnit::get_hillshade)
        .def("get_clump_map", &rxgaming::RxGamingRxUnit::get_clump_map)
        .def("get_taos", &rxgaming::RxGamingRxUnit::get_taos)
        .def("get_clump_sizes", &rxgaming::RxGamingRxUnit::get_clump_sizes)

        .def("get_treat_chm", &rxgaming::RxGamingRxUnit::get_treat_chm)
        .def("get_treat_basin", &rxgaming::RxGamingRxUnit::get_treat_basin)
        .def("get_treat_hillshade", &rxgaming::RxGamingRxUnit::get_treat_hillshade)
        .def("get_treat_clump_map", &rxgaming::RxGamingRxUnit::get_treat_clump_map)
        .def("get_treat_taos", &rxgaming::RxGamingRxUnit::get_treat_taos)
        .def("get_treat_clump_sizes", &rxgaming::RxGamingRxUnit::get_treat_clump_sizes)
        
        .def("get_cut_taos", &rxgaming::RxGamingRxUnit::get_cut_taos)

        .def("get_simulated_structures", &rxgaming::RxGamingRxUnit::get_simulated_structures);

    py::class_<rxgaming::TreatmentEngine>(m, "TreatmentEngine")
        .def(py::init<>())
        .def("do_treatment", &rxgaming::TreatmentEngine::do_treatment);
    
    py::class_<rxgaming::RxGamingProjectArea>(m, "ProjectArea")
        .def(py::init<>())
        .def(py::init([](const rxgaming::ProjectSettings& ps) {
            py::gil_scoped_release release;
            return rxgaming::RxGamingProjectArea(ps);
        }))
        .def_readwrite("rxUnits", &rxgaming::RxGamingProjectArea::rxUnits);
}
