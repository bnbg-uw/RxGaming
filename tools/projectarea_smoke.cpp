#include "projectSettings.hpp"
#include "rxgProjectArea.hpp"
#include "rxgUtils.hpp"
#include "rxgTreatmentEngine.hpp"
#include <cmath>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>

namespace {
    using rxgaming::ProjectSettings;
    using rxgaming::RxGamingRxUnit;
    using rxgaming::TreatmentEngine;
    namespace fs = std::filesystem;

    constexpr double kBaFloor = 0.0;
    constexpr int kMaxCoverIterations = 24;

    bool nearlyEqual(double lhs, double rhs, double tolerance = 1e-6) {
        return std::abs(lhs - rhs) <= tolerance;
    }

    bool sameTreatmentOutcome(const RxGamingRxUnit& lhs, const RxGamingRxUnit& rhs) {
        return lhs.treatedTaos.size() == rhs.treatedTaos.size()
            && lhs.cutTaos.size() == rhs.cutTaos.size()
            && nearlyEqual(lhs.treatedStructure.ba, rhs.treatedStructure.ba)
            && nearlyEqual(lhs.treatedStructure.tph, rhs.treatedStructure.tph)
            && nearlyEqual(lhs.treatedStructure.cc, rhs.treatedStructure.cc);
    }

    void optimizeCoverTreatment(RxGamingRxUnit& unit, double targetCc, TreatmentEngine& treater) {
        unit.targetStructure.cc = targetCc;

        RxGamingRxUnit bestUnit = unit;
        RxGamingRxUnit anchorUnit = unit;
        const bool currentIsFeasible = unit.currentStructure.cc >= targetCc;
        double bestGap = currentIsFeasible ? unit.currentStructure.cc - targetCc : std::numeric_limits<double>::infinity();

        double anchorBa = unit.currentStructure.ba;
        double stepBa = std::max(unit.currentStructure.ba * 0.25, 0.5);
        const double minStepBa = std::max(unit.currentStructure.ba * 0.001, 0.01);
        size_t lastCutCount = anchorUnit.cutTaos.size();
        int unchangedRuns = 0;

        for (int iteration = 0; iteration < kMaxCoverIterations; ++iteration) {
            if (stepBa < minStepBa) {
                break;
            }

            const double candidateBa = std::max(kBaFloor, anchorBa - stepBa);
            if (nearlyEqual(candidateBa, anchorBa, minStepBa * 0.5)) {
                stepBa *= 0.5;
                continue;
            }

            RxGamingRxUnit candidateUnit = anchorUnit;
            candidateUnit.targetStructure.ba = candidateBa;
            treater.do_treatment(candidateUnit, candidateUnit.dbhMin, candidateUnit.dbhMax);

            if (sameTreatmentOutcome(candidateUnit, anchorUnit)) {
                ++unchangedRuns;
                stepBa *= 0.5;
                if (unchangedRuns >= 2) {
                    break;
                }
                continue;
            }

            unchangedRuns = 0;
            const double candidateCc = candidateUnit.treatedStructure.cc;
            const bool candidateIsFeasible = candidateCc >= targetCc;
            const size_t candidateCutCount = candidateUnit.cutTaos.size();
            const bool oneTreeBreakpoint =
                candidateCutCount == lastCutCount + 1 || (candidateCutCount > 0 && candidateCutCount + 1 == lastCutCount);

            if (candidateIsFeasible) {
                const double candidateGap = candidateCc - targetCc;
                if (candidateGap < bestGap
                    || (nearlyEqual(candidateGap, bestGap) && candidateUnit.treatedStructure.ba < bestUnit.treatedStructure.ba)) {
                    bestGap = candidateGap;
                    bestUnit = candidateUnit;
                }

                anchorUnit = candidateUnit;
                anchorBa = candidateBa;
                lastCutCount = candidateCutCount;

                if (oneTreeBreakpoint) {
                    stepBa *= 0.5;
                }
                continue;
            }

            stepBa *= 0.5;
        }

        unit = bestUnit;
    }

    void appendStructures(
        std::ofstream& stream,
        std::string& setting,
        std::string& prescription,
        std::string& type,
        double& presc_targ,
        rxtools::StructureSummary& current, 
        rxtools::StructureSummary& target,
        rxtools::StructureSummary& treated,
        rxtools::treatmentResult result) {
            stream << setting << "," << prescription << "," << type << "," << presc_targ << ",";
            stream << current.ba  << "," << current.mcs << "," << current.tph << "," << current.cc << ",";
            stream << target.ba  << "," << target.mcs << "," << target.tph << "," << target.cc << ",";
            stream << treated.ba  << "," << treated.mcs << "," << treated.tph << "," << treated.cc << ",";
            if(result == rxtools::treatmentResult::cuttingFailure) {
                stream << "cutting_failure";
            } else if (result == rxtools::treatmentResult::diameterFailure) {
                stream << "diameter_failure";
            } else {
                stream << "success";
            }
            stream << "\n";
        }

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

        auto projectpoly = lapis::VectorDataset<lapis::MultiPolygon>(settings.unitPolyPath);
        auto treater = rxgaming::TreatmentEngine();

        fs::path base = "F:/backbone_outputs";
        if(fs::exists(base)) {
            fs::remove_all(base) ;
        }
        fs::create_directory(base);

        std::ofstream outTable(base/"table.csv");
        outTable << "DA_SETTING" << "," << "prescription" << "," << "type" << "," << "prescription_target" << ",";
        outTable << "current_ba"  << "," << "current_mcs" << "," << "current_tph" << "," << "current_cc" << ",";
        outTable << "target_ba"  << "," << "target_mcs" << "," << "target_tph" << "," << "target_cc" << ",";
        outTable << "treated_ba"  << "," << "treated_mcs" << "," << "treated_tph" << "," << "treated_cc" << ",";
        outTable << "treatment_result"<< "\n";
        
        for (size_t i = 0; i < projectArea.rxUnits.size(); ++i) {
            std::cout << i << "\n";
            auto name = projectpoly.getStringField(i, "name");
            auto type = projectpoly.getStringField(i, "type");
            auto target = projectpoly.getRealField(i, "target");
            auto setting = projectpoly.getStringField(i, "DA_SETTING");
            std::cout << "Processing unit name=" << name
                << ", type=" << type
                << ", setting=" << setting
                << ", target=" << target << "\n";

            if(!fs::exists(base / setting)) {
                fs::create_directory(base / setting);
            }

            auto& unit = projectArea.rxUnits[i];
            unit.dbhMin = 0;
            unit.dbhMax = 76.2;

            if(type == "sdi") {
                double qmd = unit.currentStructure.ba / unit.currentStructure.tph;
                qmd = qmd / 0.00007854;
                qmd = std::sqrt(qmd);

                unit.targetStructure.tph = target / std::pow(qmd / 25, 1.605);
                unit.targetStructure.ba = unit.targetStructure.tph * 0.00007854 * std::pow(qmd, 2);
                if(name == "E1" || name == "E2" || name == "A" || name == "B") {
                    unit.dbhMax = 60.706;
                }
                treater.do_treatment(unit, unit.dbhMin, unit.dbhMax);
            } else if(type == "cover") {
                optimizeCoverTreatment(unit, target, treater);
            }

            unit.chm.writeRaster((base/setting/(setting+"_chm.tif")).string());
            unit.getTreatChmRaster().writeRaster((base/setting/(setting+"_treated_chm.tif")).string());

            unit.taos.writeShapefile(base/setting/(setting+"_taos.shp"));
            unit.treatedTaos.writeShapefile(base/setting/(setting+"_treated_taos.shp"));
            appendStructures(outTable, setting, name, type, target, unit.currentStructure, unit.targetStructure, unit.treatedStructure, unit.result);
        }
        outTable.close();
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
