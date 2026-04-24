#include "rxgSnapshot.hpp"

#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include <H5Cpp.h>

namespace rxgaming {
    namespace {
        using H5::Attribute;
        using H5::DataSet;
        using H5::DataSpace;
        using H5::Group;
        using H5::H5File;
        using H5::PredType;
        using H5::StrType;

        template <typename T>
        const PredType& hdf5_type();

        template <>
        const PredType& hdf5_type<double>() {
            return PredType::NATIVE_DOUBLE;
        }

        template <>
        const PredType& hdf5_type<int>() {
            return PredType::NATIVE_INT;
        }

        template <>
        const PredType& hdf5_type<int64_t>() {
            return PredType::NATIVE_INT64;
        }

        template <>
        const PredType& hdf5_type<uint8_t>() {
            return PredType::NATIVE_UINT8;
        }

        template <>
        const PredType& hdf5_type<uint64_t>() {
            return PredType::NATIVE_UINT64;
        }

        template <typename T>
        void write_scalar_attribute(H5::H5Object& object, const std::string& name, const T& value) {
            DataSpace scalar(H5S_SCALAR);
            Attribute attribute = object.createAttribute(name, hdf5_type<T>(), scalar);
            attribute.write(hdf5_type<T>(), &value);
        }

        template <typename T>
        T read_scalar_attribute(const H5::H5Object& object, const std::string& name) {
            T value{};
            Attribute attribute = object.openAttribute(name);
            attribute.read(hdf5_type<T>(), &value);
            return value;
        }

        void write_string_attribute(H5::H5Object& object, const std::string& name, const std::string& value) {
            DataSpace scalar(H5S_SCALAR);
            StrType stringType(PredType::C_S1, H5T_VARIABLE);
            Attribute attribute = object.createAttribute(name, stringType, scalar);
            attribute.write(stringType, value);
        }

        std::string read_string_attribute(const H5::H5Object& object, const std::string& name) {
            StrType stringType(PredType::C_S1, H5T_VARIABLE);
            Attribute attribute = object.openAttribute(name);
            std::string value;
            attribute.read(stringType, value);
            return value;
        }

        template <typename T>
        void write_vector_dataset(Group& group, const std::string& name, const std::vector<T>& values) {
            const hsize_t dimensions[] = {static_cast<hsize_t>(values.size())};
            DataSpace dataSpace(1, dimensions);
            DataSet dataSet = group.createDataSet(name, hdf5_type<T>(), dataSpace);
            if (!values.empty()) {
                dataSet.write(values.data(), hdf5_type<T>());
            }
        }

        template <typename T>
        std::vector<T> read_vector_dataset(const Group& group, const std::string& name) {
            DataSet dataSet = group.openDataSet(name);
            DataSpace dataSpace = dataSet.getSpace();
            const int rank = dataSpace.getSimpleExtentNdims();
            if (rank != 1) {
                throw std::runtime_error("Expected a one-dimensional dataset for " + name);
            }

            hsize_t dimensions[1]{};
            dataSpace.getSimpleExtentDims(dimensions);
            std::vector<T> values(static_cast<size_t>(dimensions[0]));
            if (!values.empty()) {
                dataSet.read(values.data(), hdf5_type<T>());
            }
            return values;
        }

        template <typename T>
        void write_raster(Group& parent, const std::string& name, const lapis::Raster<T>& raster) {
            Group group = parent.createGroup(name);
            write_scalar_attribute<int64_t>(group, "nrow", raster.nrow());
            write_scalar_attribute<int64_t>(group, "ncol", raster.ncol());
            write_scalar_attribute<double>(group, "xmin", raster.xmin());
            write_scalar_attribute<double>(group, "ymin", raster.ymin());
            write_scalar_attribute<double>(group, "xres", raster.xres());
            write_scalar_attribute<double>(group, "yres", raster.yres());
            write_string_attribute(group, "crs_wkt", raster.crs().getCompleteWKT());

            std::vector<T> values(static_cast<size_t>(raster.ncell()));
            std::vector<uint8_t> mask(static_cast<size_t>(raster.ncell()), 0);
            for (lapis::cell_t cell = 0; cell < raster.ncell(); ++cell) {
                const auto value = raster[cell];
                if (value.has_value()) {
                    values[static_cast<size_t>(cell)] = value.value();
                    mask[static_cast<size_t>(cell)] = 1;
                }
            }

            write_vector_dataset(group, "values", values);
            write_vector_dataset(group, "mask", mask);
        }

        template <typename T>
        lapis::Raster<T> read_raster(const Group& parent, const std::string& name) {
            Group group = parent.openGroup(name);
            const auto nrow = static_cast<lapis::rowcol_t>(read_scalar_attribute<int64_t>(group, "nrow"));
            const auto ncol = static_cast<lapis::rowcol_t>(read_scalar_attribute<int64_t>(group, "ncol"));
            const auto xmin = read_scalar_attribute<double>(group, "xmin");
            const auto ymin = read_scalar_attribute<double>(group, "ymin");
            const auto xres = read_scalar_attribute<double>(group, "xres");
            const auto yres = read_scalar_attribute<double>(group, "yres");
            const auto crsWkt = read_string_attribute(group, "crs_wkt");

            lapis::Alignment alignment(xmin, ymin, nrow, ncol, xres, yres, lapis::CoordRef(crsWkt));
            lapis::Raster<T> raster(alignment);

            const auto values = read_vector_dataset<T>(group, "values");
            const auto mask = read_vector_dataset<uint8_t>(group, "mask");
            if (values.size() != mask.size() || values.size() != static_cast<size_t>(raster.ncell())) {
                throw std::runtime_error("Raster payload sizes did not match alignment for " + name);
            }

            for (lapis::cell_t cell = 0; cell < raster.ncell(); ++cell) {
                auto value = raster[cell];
                if (mask[static_cast<size_t>(cell)] != 0) {
                    value.has_value() = true;
                    value.value() = values[static_cast<size_t>(cell)];
                }
                else {
                    value.has_value() = false;
                }
            }

            return raster;
        }

        void write_taolist(Group& parent, const std::string& name, const rxtools::TaoList& taos) {
            Group group = parent.createGroup(name);
            write_string_attribute(group, "crs_wkt", taos.crs().getCompleteWKT());

            std::vector<double> xs;
            std::vector<double> ys;
            std::vector<double> heights;
            std::vector<double> radii;
            std::vector<double> areas;
            std::vector<double> dbhs;
            xs.reserve(taos.size());
            ys.reserve(taos.size());
            heights.reserve(taos.size());
            radii.reserve(taos.size());
            areas.reserve(taos.size());
            dbhs.reserve(taos.size());
            for (size_t index = 0; index < taos.size(); ++index) {
                xs.push_back(taos.x(index));
                ys.push_back(taos.y(index));
                heights.push_back(taos.height(index));
                radii.push_back(taos.radius(index));
                areas.push_back(taos.area(index));
                dbhs.push_back(taos.dbh(index));
            }

            write_vector_dataset(group, "x", xs);
            write_vector_dataset(group, "y", ys);
            write_vector_dataset(group, "height", heights);
            write_vector_dataset(group, "radius", radii);
            write_vector_dataset(group, "area", areas);
            write_vector_dataset(group, "dbh", dbhs);
        }

        rxtools::TaoList read_taolist(const Group& parent, const std::string& name) {
            Group group = parent.openGroup(name);
            const auto crsWkt = read_string_attribute(group, "crs_wkt");
            const auto xs = read_vector_dataset<double>(group, "x");
            const auto ys = read_vector_dataset<double>(group, "y");
            const auto heights = read_vector_dataset<double>(group, "height");
            const auto radii = read_vector_dataset<double>(group, "radius");
            const auto areas = read_vector_dataset<double>(group, "area");
            const auto dbhs = read_vector_dataset<double>(group, "dbh");

            if (
                xs.size() != ys.size()
                || xs.size() != heights.size()
                || xs.size() != radii.size()
                || xs.size() != areas.size()
                || xs.size() != dbhs.size()
            ) {
                throw std::runtime_error("TAO list arrays did not share a consistent length");
            }

            rxtools::TaoList taos{lapis::CoordRef(crsWkt)};
            for (size_t index = 0; index < xs.size(); ++index) {
                taos.addTao(
                    lapis::CoordXY{xs[index], ys[index]},
                    heights[index],
                    radii[index],
                    areas[index],
                    dbhs[index]
                );
            }
            return taos;
        }

        void write_structure_summary(Group& parent, const std::string& name, const rxtools::StructureSummary& summary) {
            Group group = parent.createGroup(name);
            write_scalar_attribute<double>(group, "ba", summary.ba);
            write_scalar_attribute<double>(group, "tph", summary.tph);
            write_scalar_attribute<double>(group, "mcs", summary.mcs);
            write_scalar_attribute<double>(group, "cc", summary.cc);
            write_vector_dataset(group, "csd", summary.csd);

            std::vector<int64_t> binMins(summary.binMins.begin(), summary.binMins.end());
            std::vector<int64_t> binMaxs(summary.binMaxs.begin(), summary.binMaxs.end());
            write_vector_dataset(group, "bin_mins", binMins);
            write_vector_dataset(group, "bin_maxs", binMaxs);
        }

        rxtools::StructureSummary read_structure_summary(const Group& parent, const std::string& name) {
            Group group = parent.openGroup(name);
            const auto ba = read_scalar_attribute<double>(group, "ba");
            const auto tph = read_scalar_attribute<double>(group, "tph");
            const auto mcs = read_scalar_attribute<double>(group, "mcs");
            const auto cc = read_scalar_attribute<double>(group, "cc");
            const auto csd = read_vector_dataset<double>(group, "csd");
            const auto binMinsRaw = read_vector_dataset<int>(group, "bin_mins");
            const auto binMaxsRaw = read_vector_dataset<int>(group, "bin_maxs");

            std::vector<int> binMins(binMinsRaw.begin(), binMinsRaw.end());
            std::vector<int> binMaxs(binMaxsRaw.begin(), binMaxsRaw.end());
            return rxtools::StructureSummary(ba, tph, mcs, cc, csd, binMins, binMaxs);
        }

        void write_unit(Group& unitsGroup, size_t index, const RxGamingRxUnit& unit) {
            Group group = unitsGroup.createGroup(std::to_string(index));
            write_string_attribute(group, "name", unit.name);
            write_scalar_attribute<double>(group, "area_ha", unit.areaHa);
            write_scalar_attribute<double>(group, "dbh_min", unit.dbhMin);
            write_scalar_attribute<double>(group, "dbh_max", unit.dbhMax);
            write_scalar_attribute<uint8_t>(group, "paired", unit.paired ? 1u : 0u);
            write_scalar_attribute<uint8_t>(group, "treated", unit.treated ? 1u : 0u);
            write_scalar_attribute<int64_t>(group, "result", static_cast<int64_t>(unit.result));

            write_taolist(group, "taos", unit.taos);
            write_taolist(group, "treated_taos", unit.treatedTaos);
            write_taolist(group, "cut_taos", unit.cutTaos);
            write_raster(group, "unit_mask", unit.unitMask);
            write_raster(group, "chm", unit.chm);
            write_raster(group, "basin_map", unit.basinMap);
            write_structure_summary(group, "current_structure", unit.currentStructure);
            write_structure_summary(group, "target_structure", unit.targetStructure);
            write_structure_summary(group, "treated_structure", unit.treatedStructure);
        }

        RxGamingRxUnit read_unit(const Group& unitsGroup, size_t index) {
            Group group = unitsGroup.openGroup(std::to_string(index));
            RxGamingRxUnit unit;
            unit.name = read_string_attribute(group, "name");
            unit.areaHa = read_scalar_attribute<double>(group, "area_ha");
            unit.dbhMin = read_scalar_attribute<double>(group, "dbh_min");
            unit.dbhMax = read_scalar_attribute<double>(group, "dbh_max");
            unit.paired = read_scalar_attribute<uint8_t>(group, "paired") != 0;
            unit.treated = read_scalar_attribute<uint8_t>(group, "treated") != 0;
            unit.result = static_cast<rxtools::treatmentResult>(read_scalar_attribute<int64_t>(group, "result"));

            unit.taos = read_taolist(group, "taos");
            unit.treatedTaos = read_taolist(group, "treated_taos");
            unit.cutTaos = read_taolist(group, "cut_taos");
            unit.unitMask = read_raster<lapis::cell_t>(group, "unit_mask");
            unit.chm = read_raster<double>(group, "chm");
            unit.basinMap = read_raster<int>(group, "basin_map");
            unit.currentStructure = read_structure_summary(group, "current_structure");
            unit.targetStructure = read_structure_summary(group, "target_structure");
            unit.treatedStructure = read_structure_summary(group, "treated_structure");
            unit.refresh_derived_state();
            return unit;
        }
    }

    void save_project_area(const RxGamingProjectArea& projectArea, const std::string& path) {
        try {
            const auto parentPath = std::filesystem::path(path).parent_path();
            if (!parentPath.empty()) {
                std::filesystem::create_directories(parentPath);
            }

            H5File file(path, H5F_ACC_TRUNC);
            write_scalar_attribute<int64_t>(file, "schema_version", PROJECTAREA_SNAPSHOT_SCHEMA_VERSION);
            write_scalar_attribute<int64_t>(file, "unit_count", static_cast<int64_t>(projectArea.rxUnits.size()));

            Group unitsGroup = file.createGroup("units");
            for (size_t index = 0; index < projectArea.rxUnits.size(); ++index) {
                write_unit(unitsGroup, index, projectArea.rxUnits[index]);
            }
        }
        catch (const H5::Exception& exception) {
            throw std::runtime_error("Failed to save project snapshot: " + std::string(exception.getCDetailMsg()));
        }
    }

    RxGamingProjectArea load_project_area(const std::string& path) {
        try {
            H5File file(path, H5F_ACC_RDONLY);
            const auto schemaVersion = read_scalar_attribute<int64_t>(file, "schema_version");
            if (schemaVersion != PROJECTAREA_SNAPSHOT_SCHEMA_VERSION) {
                throw std::runtime_error("Unsupported project snapshot schema version");
            }

            const auto unitCount = read_scalar_attribute<int64_t>(file, "unit_count");
            Group unitsGroup = file.openGroup("units");

            std::vector<RxGamingRxUnit> units;
            units.reserve(static_cast<size_t>(unitCount));
            for (int64_t index = 0; index < unitCount; ++index) {
                units.push_back(read_unit(unitsGroup, static_cast<size_t>(index)));
            }
            return RxGamingProjectArea(std::move(units));
        }
        catch (const H5::Exception& exception) {
            throw std::runtime_error("Failed to load project snapshot: " + std::string(exception.getCDetailMsg()));
        }
    }
}
