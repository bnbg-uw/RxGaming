#include "projectSettings.hpp"
#include "rxgProjectArea.hpp"
#include "rxgUtils.hpp"
#include "rxgTreatmentEngine.hpp"

#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iostream>
#include <string>
#include <string_view>

namespace {
    using rxgaming::ProjectSettings;
    namespace fs = std::filesystem;

    void printUsage() {
        std::cout
            << "Usage:\n"
            << "  rxgaming_projectarea_smoke"
            << " <proj-data-dir>"
            << " <unit-shapefile>"
            << " <mcs-prop-csv>"
            << " <fia-dir>"
            << " <lidar-dir>"
            << " <unit-name-field>"
            << " <threads>"
            << " [project-name]"
            << " [ref-data-csv]"
            << " [save-path]\n\n"
            << "Example:\n"
            << "  rxgaming_projectarea_smoke resources resources/units.shp resources/mcs_prop.csv"
            << " resources/fia data/lidar UNIT_NAME 1 DemoProject\n";
    }

    bool isHelpArg(std::string_view arg) {
        return arg == "--help" || arg == "-h";
    }

    ProjectSettings parseArgs(int argc, char** argv) {
        if (argc >= 2 && isHelpArg(argv[1])) {
            printUsage();
            std::exit(0);
        }

        if (argc < 8) {
            printUsage();
            throw std::runtime_error("Not enough arguments provided.");
        }

        ProjectSettings settings{};
        settings.name = argc >= 9 ? argv[8] : "rxgaming_projectarea_smoke";
        settings.unitPolyPath = argv[2];
        settings.refDataPath = argc >= 10 ? argv[9] : "";
        settings.mcsPropPath = argv[3];
        settings.fiaPath = argv[4];
        settings.lidarPath = argv[5];
        settings.unitName = argv[6];
        settings.savePath = argc >= 11 ? argv[10] : "";
        settings.nThread = std::stoi(argv[7]);
        return settings;
    }

    void validatePath(const fs::path& path, const char* label, bool expectDirectory) {
        if (!fs::exists(path)) {
            throw std::runtime_error(std::string(label) + " does not exist: " + path.string());
        }
        if (expectDirectory && !fs::is_directory(path)) {
            throw std::runtime_error(std::string(label) + " is not a directory: " + path.string());
        }
        if (!expectDirectory && fs::is_directory(path)) {
            throw std::runtime_error(std::string(label) + " is unexpectedly a directory: " + path.string());
        }
    }

    void printSettings(const ProjectSettings& settings, const fs::path& projDataDir) {
        std::cout << "ProjectArea smoke test configuration:\n";
        std::cout << "  PROJ data dir: " << projDataDir << "\n";
        std::cout << "  Project name: " << settings.name << "\n";
        std::cout << "  Unit shapefile: " << settings.unitPolyPath << "\n";
        std::cout << "  Reference data: " << (settings.refDataPath.empty() ? "<none>" : settings.refDataPath) << "\n";
        std::cout << "  MCS prop CSV: " << settings.mcsPropPath << "\n";
        std::cout << "  FIA dir: " << settings.fiaPath << "\n";
        std::cout << "  Lidar dir: " << settings.lidarPath << "\n";
        std::cout << "  Unit field: " << settings.unitName << "\n";
        std::cout << "  Threads: " << settings.nThread << "\n";
        std::cout << std::flush;
    }
}

int main(int argc, char** argv) {
    try {
        auto settings = parseArgs(argc, argv);
        const fs::path projDataDir = fs::path(argv[1]);

        validatePath(projDataDir, "PROJ data directory", true);
        validatePath(projDataDir / "proj.db", "PROJ database", false);
        validatePath(settings.unitPolyPath, "Unit polygon shapefile", false);
        validatePath(settings.mcsPropPath, "MCS prop CSV", false);
        validatePath(settings.fiaPath, "FIA directory", true);
        validatePath(settings.lidarPath, "Lidar directory", true);

        printSettings(settings, projDataDir);
        auto projDataDirArg = projDataDir.string();
        if (!projDataDirArg.empty() && projDataDirArg.back() != fs::path::preferred_separator) {
            projDataDirArg.push_back(static_cast<char>(fs::path::preferred_separator));
        }
        rxgaming::set_proj_db_path(projDataDirArg);

        std::cout << "Constructing RxGamingProjectArea...\n" << std::flush;
        rxgaming::RxGamingProjectArea projectArea(settings);
        std::cout << "Construction complete.\n";
        std::cout << "Units loaded: " << projectArea.rxUnits.size() << "\n";

        for (size_t i = 0; i < projectArea.rxUnits.size(); ++i) {
            const auto& unit = projectArea.rxUnits[i];
            std::cout << "  [" << i << "] " << unit.name << " areaHa=" << unit.areaHa << "\n";
        }
        auto treater = rxgaming::TreatmentEngine();
        for(int i = 0; i < 20; ++i) {
            projectArea.rxUnits[0].targetStructure.ba = projectArea.rxUnits[0].currentStructure.ba - (projectArea.rxUnits[0].currentStructure.ba / 20) * i;
            auto before = std::chrono::high_resolution_clock::now();
            treater.do_treatment(projectArea.rxUnits[0], 0, 9999);
            auto after = std::chrono::high_resolution_clock::now();
            auto dur = after-before;
            std::cout << "Run " << i << " of 20:" << projectArea.rxUnits[0].treatedStructure.ba << "\n";
        }
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "rxgaming_projectarea_smoke failed: " << e.what() << "\n";
        return 1;
    }
    catch (...) {
        std::cerr << "rxgaming_projectarea_smoke failed with an unknown exception.\n";
        return 1;
    }
}
