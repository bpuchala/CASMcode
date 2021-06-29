#include "casm/app/import/methods/ImportConfigurations.hh"

#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/QueryHandler.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/dataformatter/DataFormatter_impl.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/casm_io/json/optional.hh"
#include "casm/clex/ConfigMapping.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/MappedPropertiesTools.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/SimpleStructureTools.hh"
#include "casm/clex/io/json/ConfigMapping_json_io.hh"
#include "casm/clex/io/json/Configuration_json_io.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/io/SimpleStructureIO.hh"
#include "casm/database/ConfigDatabaseTools.hh"
#include "casm/database/Import.hh"  // TODO: move ImportSettings?
#include "casm/database/PropertiesDatabase.hh"
#include "casm/global/definitions.hh"

namespace CASM {

std::string ImportConfigurations::desc() const {
  return "Import configurations options ...";
}

std::string ImportConfigurations::name() const {
  return traits<Configuration>::short_name;
}

namespace import_configurations_impl {

/// Read structure paths from JSON input
///
/// Read:
/// - "structures": (array of string) structure file paths
/// - "batch": (string) path to batch file containing one structure file path
///   per line
/// - "selection: (string) name of Configuration selection to read
///   properties.calc.json files from using default calctype
std::vector<fs::path> parse_structure_paths(ParentInputParser &parser,
                                            PrimClex const *primclex_ptr) {
  std::vector<fs::path> structure_paths;

  // read "structures"
  if (parser.self.contains("structures")) {
    std::vector<fs::path> array_of_paths;
    parser.optional(array_of_paths, "structures");
    for (auto const &structure_path : array_of_paths) {
      structure_paths.push_back(fs::absolute(structure_path));
    }
  }

  // read "batch" file
  if (parser.self.contains("batch")) {
    fs::path batch_path;
    parser.optional(batch_path, "batch");
    if (!batch_path.empty()) {
      try {
        fs::ifstream batchfile(batch_path);
        fs::path tpath;
        while (batchfile >> tpath) {
          structure_paths.push_back(fs::absolute(tpath));
          batchfile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
      } catch (std::exception &e) {
        std::stringstream msg;
        msg << "Error: could not read structure file paths from batch file "
            << batch_path << ": " << e.what();
        parser.insert_error("batch", msg.str());
      }
    }
  }

  // read "selection" (if primclex provided)
  if (primclex_ptr != nullptr) {
    // unpack requirements for reading selection
    DB::Database<Configuration> &configuration_db =
        primclex_ptr->db<Configuration>();
    DirectoryStructure const &dir = primclex_ptr->dir();
    std::string calctype = primclex_ptr->settings().default_clex().calctype;

    if (parser.self.contains("selection")) {
      fs::path selection_path;
      parser.optional(selection_path, "selection");
      if (!selection_path.empty()) {
        try {
          DB::Selection<Configuration> selection{configuration_db,
                                                 selection_path};
          auto it = selection.selected().begin();
          auto end = selection.selected().end();
          for (; it != end; ++it) {
            structure_paths.push_back(
                dir.calculated_properties(it.name(), calctype));
          }
        } catch (std::exception &e) {
          std::stringstream msg;
          msg << "Error: could not read selection " << selection_path << ": "
              << e.what();
          parser.insert_error("selection", msg.str());
        }
      }
    }
  }

  if (!structure_paths.size()) {
    std::stringstream msg;
    msg << "Error: no structure files";
    parser.error.insert(msg.str());
  }

  return structure_paths;
}

/// Read xtal::SimpleStructure from structure file paths
///
/// Reads xtal::SimpleStructure from files according to:
/// - If extension == ".json" or ".JSON", read as CASM xtal::SimpleStructure
/// directly from JSON file
/// - Otherwise, assume VASP POSCAR file and read as xtal::BasicStructure and
/// then convert to xtal::SimpleStructure
///
/// This method attempts to catch any reading errors and store error messages
/// inside the validator.
///
/// \param structure_paths Vector of structure file paths to read
/// \param validator Collects failure messages. (This is a base class of
/// InputParser).
///
/// \return std::map of structure file path to xtal::SimpleStructure
///
std::map<fs::path, xtal::SimpleStructure> make_structures(
    std::vector<fs::path> const &structure_paths, Validator &validator) {
  std::map<fs::path, xtal::SimpleStructure> structures;

  for (auto const &p : structure_paths) {
    if (!fs::exists(p)) {
      std::stringstream msg;
      msg << "Error: Failed to read structure file:" << p
          << ": file does not exist";
      validator.error.insert(msg.str());
      continue;
    }

    // read JSON files as xtal::SimpleStructure directly
    if (p.extension() == ".json" || p.extension() == ".JSON") {
      try {
        xtal::SimpleStructure simple_structure;
        jsonParser json{p};
        from_json(simple_structure, json);
        structures.emplace(p, simple_structure);
      } catch (std::exception &e) {
        std::stringstream msg;
        msg << "Error: Failed to read structure file as a CASM JSON structure "
            << "file:" << p << ": " << e.what();
        validator.error.insert(msg.str());
      }
    } else {
      // read others files as VASP POSCAR
      try {
        fs::ifstream sin{p};
        xtal::BasicStructure basic_structure =
            xtal::BasicStructure::from_poscar_stream(sin);
        xtal::SimpleStructure simple_structure =
            make_simple_structure(basic_structure);
        structures.emplace(p, simple_structure);
      } catch (std::exception &e) {
        std::stringstream msg;
        msg << "Error: Failed to read structure file as a VASP POSCAR format "
            << "file:" << p << ": " << e.what();
        validator.error.insert(msg.str());
      }
    }
  }
  return structures;
}

struct StructureImportData;
struct MappingImportData;
struct ConfigurationImportData;

struct StructureImportData {
  /// Location of structure file
  fs::path structure_path;

  /// Structure being mapped
  xtal::SimpleStructure structure;

  /// ConfigMapper results
  ConfigMapperResult mapping_result;
};

/// Used to hold or access configuration <-> structure mapping information
///
/// This will be constructed *after* StructureImportData and
/// ConfigurationImportData?
struct MappingImportData {
  // /// Set to true for mappings that existed prior to this import batch
  // bool preexisting = false;

  /// Properties determined for this mapping, including the name of
  /// configuration mapped to
  MappedProperties properties;

  /// If all PrimClex required properties exist
  bool has_all_required_properties;

  /// If all PrimClex required properties exist, this is the tie-breaking score
  /// as determined by the score_method for the configuration this object
  /// describes the mapping to. Only has value if `has_all_required_properties`
  /// is true.
  std::optional<double> properties_score;

  /// MappingImportData are compared primarily by comparing
  /// `properties_score`. Any Mapping that `has_all_required_properties`
  /// compares less than any that does not.
  bool operator<(MappingImportData const &other) const {
    if (this->properties_score.has_value() &&
        other.properties_score.has_value()) {
      return this->properties_score.value() < other.properties_score.value();
    } else if (this->properties_score.has_value()) {
      return true;
    } else if (other.properties_score.has_value()) {
      return false;
    } else {
      // fallback: compare total mapping cost
      return this->properties.scalar("total_cost") <
             other.properties.scalar("total_cost");
    }
  }
};

struct ConfigurationImportData {
  /// Insert configurations into the configuration database and hold the result
  /// here. Depending on ImportSettings, changes may not be committed.
  ConfigInsertResult insert_result;

  /// Include currently imported and existing imported structures that map to
  /// this configuration
  std::set<MappingImportData> mappings;

  /// Method used to calculate the properties_score of successful mappings. This
  /// determines which successful mappings to this configuration is considered
  /// the best (best is lowest score).
  ScoreMappedProperties mapped_properties_score_method;
};

}  // namespace import_configurations_impl

// COORD_TYPE coordinate_mode

jsonParser &to_json(
    const import_configurations_impl::StructureImportData &structure_data,
    jsonParser &json, COORD_TYPE coordinate_mode) {
  to_json(structure_data.structure_path.string(), json["structure_path"]);
  to_json(structure_data.structure, json["structure"], {},
          coordinate_mode);  // unmapped structure
  to_json(structure_data.mapping_result, json["mapping_result"],
          coordinate_mode);

  return json;
}

jsonParser &to_json(
    const import_configurations_impl::MappingImportData &mapping_data,
    jsonParser &json) {
  to_json(mapping_data.properties, json["properties"]);
  to_json(mapping_data.has_all_required_properties,
          json["has_all_required_properties"]);
  to_json(mapping_data.properties_score, json["properties_score"]);
  // to_json(mapping_data.preexisting, json["preexisting"]);

  return json;
}

jsonParser &to_json(
    const import_configurations_impl::ConfigurationImportData &config_data,
    jsonParser &json, DB::Database<Configuration> const &configuration_db,
    COORD_TYPE coordinate_mode) {
  if (config_data.insert_result.canonical_it != configuration_db.end()) {
    Configuration const &canonical_configuration =
        *config_data.insert_result.canonical_it;
    to_json(canonical_configuration, json["canonical_configuration"]);
    to_json(canonical_configuration.name(),
            json["canonical_configuration_name"]);
    to_json(config_data.insert_result.insert_canonical,
            json["canonical_configuration_is_new"]);
  }

  Configuration const &primitive_configuration =
      *config_data.insert_result.primitive_it;

  to_json(primitive_configuration, json["primitive_configuration"]);
  to_json(primitive_configuration.name(), json["primitive_configuration_name"]);
  to_json(config_data.insert_result.insert_primitive,
          json["primitive_configuration_is_new"]);

  if (config_data.mappings.size()) {
    json["mappings"] = jsonParser::array();
    for (auto const &mapping : config_data.mappings) {
      jsonParser mapping_json;
      to_json(mapping, mapping_json);

      if (config_data.insert_result.canonical_it != configuration_db.end()) {
        Configuration const &canonical_configuration =
            *config_data.insert_result.canonical_it;

        xtal::SimpleStructure simple_structure =
            make_simple_structure(canonical_configuration.supercell(),
                                  canonical_configuration.configdof());

        xtal::SimpleStructure simple_structure_relaxed = make_simple_structure(
            canonical_configuration.supercell(),
            canonical_configuration.configdof(), mapping.properties);

        to_json(simple_structure, mapping_json["structure"], {},
                coordinate_mode);
        to_json(simple_structure_relaxed, mapping_json["relaxed_structure"], {},
                coordinate_mode);
      }
      json["mappings"].push_back(mapping_json);
    }
  }

  to_json(config_data.mapped_properties_score_method,
          json["properties_score_method"]);

  return json;
}

void check_equal(Eigen::MatrixXd const &A, Eigen::MatrixXd const &B,
                 std::string message) {
  if (!almost_equal(A, B)) {
    std::cout << "A: \n" << A << std::endl;
    std::cout << "B: \n" << B << std::endl << std::endl;
    throw std::runtime_error(message);
  }
}

void ImportConfigurations::run(PrimClex &primclex,
                               jsonParser const &json_options,
                               jsonParser const &cli_options_as_json) const {
  std::cout << "run `import --type config`" << std::endl;

  using namespace import_configurations_impl;

  Log &log = CASM::log();

  // Getting data from primclex that will be used below
  // -- using default calctype
  std::shared_ptr<Structure const> shared_prim = primclex.shared_prim();
  DataFormatterDictionary<Supercell> const &supercell_query_dict =
      primclex.settings().query_handler<Supercell>().dict();
  std::string calctype = primclex.settings().default_clex().calctype;
  DB::Database<Configuration> &configuration_db = primclex.db<Configuration>();
  DB::Database<Supercell> &supercell_db = primclex.db<Supercell>();
  DB::PropertiesDatabase &properties_db =
      primclex.db_props<Configuration>(calctype);
  std::vector<std::string> required_properties =
      primclex.settings().required_properties(traits<Configuration>::name,
                                              calctype);

  COORD_TYPE coordinate_mode_for_structure_output = CART;

  ParentInputParser parser =
      make_import_parent_parser(log, json_options, cli_options_as_json);
  std::runtime_error error_if_invalid{
      "Error reading `import --type config` JSON input"};

  log.custom("Checking input");

  // 1) Parse structure file paths and make structures --------------

  std::vector<fs::path> structure_paths =
      parse_structure_paths(parser, &primclex);
  report_and_throw_if_invalid(parser, log, error_if_invalid);

  std::map<fs::path, xtal::SimpleStructure> structures =
      make_structures(structure_paths, parser);
  report_and_throw_if_invalid(parser, log, error_if_invalid);

  // std::cout << "showing structures:\n";
  // for (auto const &value : structures) {
  //   fs::path path = value.first;
  //   xtal::SimpleStructure const &simple_structure = value.second;
  //
  //   std::cout << path << ":" << std::endl;
  //   jsonParser json;
  //   to_json(simple_structure, json);
  //   std::cout << json << std::endl << std::endl;
  // }

  // 2) Parse structure mapping settings ----

  ConfigMapperSettings config_mapper_settings;
  parser.optional(config_mapper_settings, "mapping", shared_prim,
                  supercell_query_dict);
  report_and_throw_if_invalid(parser, log, error_if_invalid);

  log << "input structures: " << structures.size() << std::endl;
  log << "DONE" << std::endl << std::endl;

  // 3) Parse import settings

  DB::ImportSettings import_settings;
  parser.optional(import_settings, "import");

  // 3) Construct ConfigMapper and map structures

  // input:
  // - shared_prim
  // - config_mapper_settings
  // - structures
  // output:
  // - structure_import_data

  ConfigMapper config_mapper{shared_prim, config_mapper_settings};

  log.custom("Mapping structures");
  std::map<fs::path, StructureImportData> structure_import_data;
  for (auto const &value : structures) {
    fs::path path = value.first;
    xtal::SimpleStructure const &simple_structure = value.second;

    log << "Working on " << path << std::endl;

    Configuration const *hint_config_ptr = nullptr;  // TODO: allow this
    ConfigMapperResult mapping_result =
        config_mapper.import_structure(simple_structure, hint_config_ptr);

    StructureImportData _sdata;
    _sdata.structure_path = path;
    _sdata.structure = simple_structure;
    _sdata.mapping_result = mapping_result;

    structure_import_data.emplace(path, _sdata);
  }
  log << "DONE" << std::endl << std::endl;

  for (auto const &value : structure_import_data) {
    fs::path structure_path = value.first;
    StructureImportData const &structure_data = value.second;
    xtal::SimpleStructure const &unmapped_child = structure_data.structure;

    std::cout << "origin: " << value.first.string()
              << "  mappings: " << value.second.mapping_result.maps.size()
              << std::endl;
    jsonParser json;
    to_json(value.second, json, coordinate_mode_for_structure_output);
    std::cout << json << std::endl << std::endl;

    ConfigMapperResult const &mapping_result = structure_data.mapping_result;
    for (auto const &mapping : mapping_result.maps) {
      // collect data
      xtal::SimpleStructure const &unmapped_child = mapping.unmapped_child;
      xtal::MappingNode const &mapping_node = mapping.mapping;
      xtal::LatticeNode const &lattice_node = mapping_node.lattice_node;
      xtal::SimpleStructure const &mapped_child = mapping.mapped_child;
      Configuration const &mapped_configuration = mapping.mapped_configuration;
      Configuration const &final_configuration = mapping.final_configuration;
      // xtal::SimpleStructure const &resolved_struc = indiv.resolved_struc;

      // // LatticeNode:
      // // parent == stretch * isometry * unmapped_child
      // // LatticeMap: unmapped_child == deformation_gradient*parent*N
      // //   parent * N = deformation_gradient.inverse * unmapped_child
      // // paper: S_1*N == L_1*T_1*N == V_N * Q_N * L_2
      //
      // // relationships that should hold:
      // // 1) S1 == lattice_node.stretch * lattice_node.isometry * L2
      //
      // StrucMapper results

      check_equal(mapped_configuration.ideal_lattice().lat_column_mat(),
                  lattice_node.parent.superlattice().lat_column_mat(),
                  "ImportConfigurations check -2a");

      check_equal(mapped_configuration.ideal_lattice().lat_column_mat(),
                  lattice_node.stretch * lattice_node.isometry *
                      unmapped_child.lat_column_mat *
                      lattice_node.child.transformation_matrix_to_super()
                          .cast<double>(),
                  "ImportConfigurations check -2b");

      check_equal(lattice_node.parent.superlattice().lat_column_mat(),
                  lattice_node.stretch * lattice_node.isometry *
                      unmapped_child.lat_column_mat *
                      lattice_node.child.transformation_matrix_to_super()
                          .cast<double>(),
                  "ImportConfigurations check 0");

      // ConfigMapper results
      check_equal(mapped_child.lat_column_mat,
                  lattice_node.isometry * unmapped_child.lat_column_mat,
                  "check 1");
    }
  }

  // input:
  // - structure_import_data
  // - primitive configs?
  // output:
  // - config_import_data

  log.custom("Checking results");

  // collect Configuration data for all configurations that were mapped to
  std::map<Configuration, ConfigurationImportData> config_import_data;

  Index successful_mappings = 0;
  Index failed_mappings = 0;
  Index structures_with_multiple_optimal_mappings = 0;
  for (auto const &value : structure_import_data) {
    fs::path structure_path = value.first;
    StructureImportData const &structure_data = value.second;
    ConfigMapperResult const &mapping_result = structure_data.mapping_result;

    if (mapping_result.success()) {
      successful_mappings++;
      if (mapping_result.n_optimal() > 1) {
        // should this be a failure?
        structures_with_multiple_optimal_mappings++;
      }
    } else {
      failed_mappings++;
      continue;
    }

    for (auto const &mapping : mapping_result.maps) {
      xtal::MappingNode const &mapping_node = mapping.mapping;
      Configuration const &final_configuration = mapping.final_configuration;
      MappedProperties const &final_properties = mapping.final_properties;

      ConfigInsertResult insert_result = make_canonical_and_insert(
          final_configuration, supercell_db, configuration_db,
          import_settings.primitive_only);

      // find existing or create new ConfigurationImportData
      // -- one per unique configuration mapped to by a structure
      auto it = config_import_data.find(*insert_result.canonical_it);
      if (it == config_import_data.end()) {
        it =
            config_import_data
                .emplace(*insert_result.canonical_it, ConfigurationImportData())
                .first;

        // these should only be set once, the first time the configuration is
        // found and ConfigurationImportData constructed
        it->second.insert_result = insert_result;
        it->second.mapped_properties_score_method =
            properties_db.score_method(insert_result.canonical_it->name());

        // get preexisting mappings?
      }
      Configuration const &config = it->first;
      ConfigurationImportData &_cdata = it->second;

      // Construct and store mapping data for each mapping
      MappingImportData _mdata;
      _mdata.properties = final_properties;
      _mdata.properties.to = config.name();
      _mdata.properties.origin = structure_path.string();
      _mdata.properties.file_data = structure_path.string();
      // set 'init_config'?
      _mdata.has_all_required_properties =
          has_all_required_properties(_mdata.properties, required_properties);

      if (_mdata.has_all_required_properties) {
        _mdata.properties_score =
            _cdata.mapped_properties_score_method(_mdata.properties);
      }

      _cdata.mappings.insert(_mdata);
    }

    // if (!map_result.success()) {
    //   res.fail_msg = map_result.fail_msg;
    //   *result++ = std::move(res);
    //   return result;
    // }
    //
    // if (map_result.n_optimal() > 1) {
    //   res.fail_msg = "There were " + std::to_string(map_result.n_optimal()) +
    //                  " optimal mappings, when only one was expected.";
    //   *result++ = std::move(res);
    //   return result;
    // }
  }
  log << "structures successfully mapped: " << successful_mappings << std::endl;
  log << "structures with multiple optimal mappings: "
      << structures_with_multiple_optimal_mappings << std::endl;
  log << "structures that could not be mapped: " << failed_mappings << std::endl
      << std::endl;

  for (auto const &value : config_import_data) {
    std::cout << "config: " << value.first.name()
              << "  mappings: " << value.second.mappings.size() << std::endl;
    jsonParser json;
    to_json(value.second, json, configuration_db,
            coordinate_mode_for_structure_output);
    std::cout << json << std::endl << std::endl;
  }
}

}  // namespace CASM
