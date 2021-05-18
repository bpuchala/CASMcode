#include "casm/app/import/methods/ImportConfigurations.hh"

#include "casm/app/DirectoryStructure.hh"
#include "casm/app/ProjectSettings.hh"
#include "casm/app/QueryHandler.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/dataformatter/DataFormatter_impl.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/ConfigMapping.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/MappedPropertiesTools.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/io/json/ConfigMapping_json_io.hh"
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
    std::copy(array_of_paths.begin(), array_of_paths.end(),
              std::back_inserter(structure_paths));
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
          structure_paths.push_back(tpath);
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
  /// Pointer to an element of ConfigMapperResult::maps with results that mapped
  /// to this configuration - if structure part of current import
  std::map<xtal::MappingNode, ConfigMapperResult::Individual>::value_type const
      *mapping_ptr;

  /// Properties determined for this mapping, including the name of
  /// configuration mapped to
  MappedProperties properties;

  /// If all PrimClex required properties exist
  bool has_all_required_properties;

  /// If all PrimClex required properties exist, this is the tie-breaking score
  /// as determined by PropertiesDatabase score_method for the configuration
  /// this object describes the mapping to
  std::optional<double> mapped_properties_score;

  bool operator<(MappingImportData const &other) const {
    if (this->mapped_properties_score.has_value() &&
        other.mapped_properties_score.has_value()) {
      return this->mapped_properties_score.value() <
             other.mapped_properties_score.value();
    } else if (this->mapped_properties_score.has_value()) {
      return true;
    } else if (other.mapped_properties_score.has_value()) {
      return false;
    } else {
      // compare MappingNode
      return this->mapping_ptr->first < other.mapping_ptr->first;
    }
  }
};

struct ConfigurationImportData {
  /// Insert configurations and hold the result here. Depending on
  /// ImportSettings, changes may not be committed.
  ConfigInsertResult insert_result;

  /// Include currently imported and existing imported structures that map to
  /// this configuration
  std::set<MappingImportData> mappings;
};

}  // namespace import_configurations_impl

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
  // for (auto const &pair : structures) {
  //   fs::path path = pair.first;
  //   xtal::SimpleStructure const &simple_structure = pair.second;
  //
  //   std::cout << path << ":" << std::endl;
  //   jsonParser json;
  //   to_json(simple_structure, json);
  //   std::cout << json << std::endl << std::endl;
  // }

  // 2) Parse structure mapping settings ----

  ConfigMapping::Settings config_mapper_settings;
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

  ConfigMapper config_mapper{shared_prim, config_mapper_settings,
                             shared_prim->lattice().tol()};

  log.custom("Mapping structures");
  std::map<fs::path, StructureImportData> structure_import_data;
  for (auto const &pair : structures) {
    fs::path path = pair.first;
    xtal::SimpleStructure const &simple_structure = pair.second;

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
  for (auto const &pair : structure_import_data) {
    fs::path structure_path = pair.first;
    StructureImportData const &structure_data = pair.second;
    ConfigMapperResult const &mapping_result = structure_data.mapping_result;

    if (mapping_result.success()) {
      successful_mappings++;
      if (mapping_result.n_optimal() > 1) {
        // should this fail?
        structures_with_multiple_optimal_mappings++;
      }
    } else {
      failed_mappings++;
      continue;
    }

    for (auto const &mapping : mapping_result.maps) {
      xtal::MappingNode const &mapping_node = mapping.first;
      ConfigMapperResult::Individual const &indiv = mapping.second;

      ConfigInsertResult insert_result = make_canonical_and_insert(
          indiv.config, supercell_db, configuration_db,
          import_settings.primitive_only);

      // find existing or create new ConfigurationImportData
      // -- one per unique configuration mapped to by a structure
      auto it = config_import_data.find(*insert_result.canonical_it);
      if (it == config_import_data.end()) {
        it =
            config_import_data
                .emplace(*insert_result.canonical_it, ConfigurationImportData())
                .first;
      }
      Configuration const &config = it->first;
      ConfigurationImportData &_cdata = it->second;
      _cdata.insert_result = insert_result;

      // Construct and store mapping data
      MappingImportData _mdata;
      _mdata.mapping_ptr = &mapping;
      _mdata.properties = make_mapped_properties(mapping_node, indiv);
      _mdata.properties.to = config.name();
      _mdata.properties.origin = structure_path.string();
      _mdata.properties.file_data = structure_path.string();
      // set 'init_config'?
      _mdata.has_all_required_properties =
          has_all_required_properties(_mdata.properties, required_properties);

      if (_mdata.has_all_required_properties) {
        _mdata.mapped_properties_score = properties_db.score(_mdata.properties);
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
  }
}

}  // namespace CASM
