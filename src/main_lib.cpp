#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <thread>
#include <utility>
#include "projectSettings.hpp"
#include "rxgProjectArea.hpp"
#include "rxgSnapshot.hpp"
#include "rxgTreatmentEngine.hpp"

namespace py = pybind11;

namespace {
    struct ProjectAreaBuildSnapshot {
        std::string stage;
        std::string message;
        int completed = -1;
        int total = -1;
        std::string status = "running";
        std::string error;
    };

    struct ProjectAreaBuildState {
        std::mutex mutex;
        ProjectAreaBuildSnapshot snapshot;
        std::optional<rxgaming::RxGamingProjectArea> projectArea;
        bool finalized = false;

        void applyEvent(const rxgaming::ProgressEvent& event) {
            std::lock_guard<std::mutex> lock(mutex);
            snapshot.stage = event.stage;
            snapshot.message = event.message;
            snapshot.completed = event.completed;
            snapshot.total = event.total;
            if (event.stage == "unit_failed") {
                snapshot.status = "failed";
                snapshot.error = event.message;
            }
        }

        ProjectAreaBuildSnapshot snapshotCopy() {
            std::lock_guard<std::mutex> lock(mutex);
            return snapshot;
        }

        void storeSuccess(rxgaming::RxGamingProjectArea&& area) {
            std::lock_guard<std::mutex> lock(mutex);
            projectArea = std::move(area);
            snapshot.status = "succeeded";
            snapshot.error.clear();
        }

        void storeFailure(const std::string& errorText) {
            std::lock_guard<std::mutex> lock(mutex);
            snapshot.status = "failed";
            snapshot.error = errorText;
            if (snapshot.message.empty()) {
                snapshot.message = errorText;
            }
        }
    };

    class ProjectAreaBuildHandle {
    public:
        explicit ProjectAreaBuildHandle(const rxgaming::ProjectSettings& ps)
            : state_(std::make_shared<ProjectAreaBuildState>())
        {
            worker_ = std::thread([state = state_, ps]() {
                rxgaming::ProgressCallback progressCallback = [state](const rxgaming::ProgressEvent& event) {
                    state->applyEvent(event);
                };

                try {
                    auto projectArea = rxgaming::RxGamingProjectArea(ps, progressCallback);
                    state->storeSuccess(std::move(projectArea));
                }
                catch (const std::exception& e) {
                    state->storeFailure(e.what());
                }
                catch (...) {
                    state->storeFailure("Project build failed with an unknown exception.");
                }
            });
        }

        ProjectAreaBuildHandle(const ProjectAreaBuildHandle&) = delete;
        ProjectAreaBuildHandle& operator=(const ProjectAreaBuildHandle&) = delete;

        ~ProjectAreaBuildHandle() {
            if (worker_.joinable()) {
                worker_.join();
            }
        }

        ProjectAreaBuildSnapshot poll() {
            return state_->snapshotCopy();
        }

        rxgaming::RxGamingProjectArea finish() {
            {
                std::lock_guard<std::mutex> lock(state_->mutex);
                if (state_->snapshot.status == "running") {
                    throw std::runtime_error("Project build is still running.");
                }
                if (state_->finalized) {
                    throw std::runtime_error("Project build has already been finalized.");
                }
                state_->finalized = true;
            }

            if (worker_.joinable()) {
                py::gil_scoped_release release;
                worker_.join();
            }

            std::optional<rxgaming::RxGamingProjectArea> completedArea;
            std::string errorText;
            {
                std::lock_guard<std::mutex> lock(state_->mutex);
                if (state_->snapshot.status == "failed") {
                    errorText = state_->snapshot.error.empty() ? state_->snapshot.message : state_->snapshot.error;
                }
                else if (state_->projectArea.has_value()) {
                    completedArea = std::move(state_->projectArea);
                }
            }

            if (!errorText.empty()) {
                throw std::runtime_error(errorText);
            }
            if (!completedArea.has_value()) {
                throw std::runtime_error("Project build completed without a resulting ProjectArea.");
            }
            return std::move(completedArea.value());
        }

    private:
        std::shared_ptr<ProjectAreaBuildState> state_;
        std::thread worker_;
    };

    std::shared_ptr<ProjectAreaBuildHandle> start_project_area_build(const rxgaming::ProjectSettings& ps) {
        return std::make_shared<ProjectAreaBuildHandle>(ps);
    }

    ProjectAreaBuildSnapshot poll_project_area_build(const std::shared_ptr<ProjectAreaBuildHandle>& handle) {
        if (!handle) {
            throw std::runtime_error("Project build handle is invalid.");
        }
        return handle->poll();
    }

    rxgaming::RxGamingProjectArea finish_project_area_build(const std::shared_ptr<ProjectAreaBuildHandle>& handle) {
        if (!handle) {
            throw std::runtime_error("Project build handle is invalid.");
        }
        return handle->finish();
    }
}

PYBIND11_MODULE(rxgaming_core, m) {
    m.def("set_proj_db_path", &rxgaming::set_proj_db_path, "Set the PROJ data directory path");
    m.def("set_seed", &rxgaming::set_seed, "Set the random seed");
    m.def("start_project_area_build", &start_project_area_build, py::arg("ps"), "Start a ProjectArea build in the background.");
    m.def("poll_project_area_build", &poll_project_area_build, py::arg("handle"), "Poll the latest ProjectArea build progress.");
    m.def("finish_project_area_build", &finish_project_area_build, py::arg("handle"), "Finalize a completed ProjectArea build.");
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

    py::class_<ProjectAreaBuildSnapshot>(m, "ProjectAreaBuildSnapshot")
        .def(py::init<>())
        .def_readonly("stage", &ProjectAreaBuildSnapshot::stage)
        .def_readonly("message", &ProjectAreaBuildSnapshot::message)
        .def_readonly("completed", &ProjectAreaBuildSnapshot::completed)
        .def_readonly("total", &ProjectAreaBuildSnapshot::total)
        .def_readonly("status", &ProjectAreaBuildSnapshot::status)
        .def_readonly("error", &ProjectAreaBuildSnapshot::error);

    py::class_<ProjectAreaBuildHandle, std::shared_ptr<ProjectAreaBuildHandle>>(m, "ProjectAreaBuildHandle");

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
        .def("export_rendered_geotiff", &rxgaming::RxGamingRxUnit::export_rendered_geotiff)
        .def("write_tao_shapefile", &rxgaming::RxGamingRxUnit::write_tao_shapefile)
        .def("write_chm_raster", &rxgaming::RxGamingRxUnit::write_chm_raster)
        .def("write_basin_raster", &rxgaming::RxGamingRxUnit::write_basin_raster)
        .def("write_clumpmap_raster", &rxgaming::RxGamingRxUnit::write_clumpmap_raster)

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
