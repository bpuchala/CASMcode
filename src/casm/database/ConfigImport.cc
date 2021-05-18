#include "casm/database/ConfigImport.hh"

#include "casm/app/DirectoryStructure.hh"
#include "casm/app/QueryHandler.hh"
#include "casm/casm_io/dataformatter/DataFormatter_impl.hh"
#include "casm/clex/ConfigMapping.hh"
#include "casm/clex/Configuration_impl.hh"
#include "casm/clex/io/json/ConfigMapping_json_io.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/io/SimpleStructureIO.hh"
#include "casm/database/ConfigDatabase.hh"
#include "casm/database/Import_impl.hh"
#include "casm/database/PropertiesDatabase.hh"
#include "casm/database/ScelDatabase.hh"
#include "casm/database/Selection_impl.hh"
#include "casm/database/Update_impl.hh"

namespace CASM {

template class DataFormatter<DB::ConfigIO::Result>;

namespace DB {

// --- Configuration specializations ---------------------------------------

// --- Import<Configuration> ---

ConfigMapping::Settings const &StructureMap<Configuration>::settings() const {
  return m_configmapper->settings();
}

/// Construct with PrimClex and ConfigMapping::Settings (see Import / Update
/// desc)
StructureMap<Configuration>::StructureMap(ConfigMapping::Settings const &_set,
                                          const PrimClex &primclex,
                                          bool _primitive_only)
    : m_primclex_ptr(&primclex), m_primitive_only(_primitive_only) {
  auto const &shared_prim = primclex.shared_prim();

  // -- construct ConfigMapper --
  m_configmapper.reset(
      new ConfigMapper(shared_prim, _set, shared_prim->lattice().tol()));
}

/// \brief Specialized import method for ConfigType
///
/// \param p Path to structure or properties.calc.json file. Not guaranteed to
/// exist or be valid. \param hint std::unique_ptr<Configuration> to 'from'
/// config for 'casm update', or null if unknown as with 'casm import'. \param
/// result Insert iterator of Result objects to output mapping results
///
/// - Should output one or more mapping results from the structure located at
/// specied path
/// - >1 result handles case of non-primitive configurations
/// - responsible for filling in Result data structure
StructureMap<Configuration>::map_result_inserter
StructureMap<Configuration>::map(
    fs::path p, std::vector<std::string> const &req_properties,
    std::unique_ptr<Configuration> const &hint_config,
    map_result_inserter result) const {
  // need to set Result data (w/ defaults):
  // - std::string pos = "";
  // - MappedProperties mapped_props {origin:"", to:"", unmapped:{}, mapped:{}};
  // - bool has_all_required_properties = false;
  // - bool is_new_config = false;
  // - std::string fail_msg = "";
  ConfigIO::Result res;
  res.pos_path = p.string();

  if (!fs::exists(res.pos_path)) {
    res.fail_msg = "Specified file does not exist!";
  }
  // read from structure file or properties.calc.json file (if exists)
  SimpleStructure sstruc = this->_make_structure(res.pos_path);

  // do mapping
  ConfigMapperResult map_result =
      m_configmapper->import_structure(sstruc, hint_config.get());

  if (!map_result.success()) {
    res.fail_msg = map_result.fail_msg;
    *result++ = std::move(res);
    return result;
  }

  if (map_result.n_optimal() > 1) {
    res.fail_msg = "There were " + std::to_string(map_result.n_optimal()) +
                   " optimal mappings, when only one was expected.";
    *result++ = std::move(res);
    return result;
  }

  for (auto const &map : map_result.maps) {
    map.second.config.supercell().set_primclex(m_primclex_ptr);

    // insert in database (note that this also/only inserts primitive)
    ConfigInsertResult insert_result =
        map.second.config.insert(m_primitive_only);

    res.is_new_config = insert_result.insert_canonical;

    res.properties = make_mapped_properties(map.first, map.second);
    res.properties.file_data = p.string();
    res.properties.to = insert_result.canonical_it.name();
    res.properties.origin = p.string();

    if (hint_config) {
      res.properties.init_config = hint_config->name();
    }

    res.has_all_required_properties = true;
    for (std::string const &propname : req_properties) {
      if (!res.properties.global.count(propname) &&
          !res.properties.site.count(propname)) {
        res.has_all_required_properties = false;
        break;
      }
    }

    // at this point, the mapped structure result is complete
    *result++ = res;

    // it may be the structure was not primitive:
    // - in which case we need to create a result indicating that the primitive
    //   was also inserted in the database,
    // - but don't try to scale the data for the primitive
    if (insert_result.canonical_it != insert_result.primitive_it) {
      ConfigIO::Result prim_res;
      prim_res.pos_path = res.pos_path;
      prim_res.properties.file_data = res.properties.file_data;
      prim_res.properties.origin =
          "prim:" +
          res.properties.origin;  // insert_result.primitive_it.name();
      prim_res.properties.to = insert_result.primitive_it.name();
      prim_res.is_new_config = insert_result.insert_primitive;
      // at this point, the mapped structure result is complete
      *result++ = prim_res;
    }
  }
  return result;
}

/// \brief Read BasicStructure to be imported
///
/// If 'p.extension()' == ".json" or ".JSON", read as properties.calc.json
/// Else, read as VASP POSCAR
SimpleStructure StructureMap<Configuration>::_make_structure(
    const fs::path &p) const {
  SimpleStructure sstruc;
  if (p.extension() == ".json" || p.extension() == ".JSON") {
    jsonParser json(p);
    from_json(sstruc, json);
  } else {
    fs::ifstream struc_stream(p);
    BasicStructure struc = BasicStructure::from_poscar_stream(struc_stream);
    sstruc = make_simple_structure(struc);
  }
  return sstruc;
}

// --- Import<Configuration> ---

/// \brief Constructor
Import<Configuration>::Import(const PrimClex &primclex,
                              const StructureMap<Configuration> &mapper,
                              ImportSettings const &_set,
                              std::string const &report_dir)
    :

      ImportT(primclex, mapper, _set, report_dir) {}

const std::string Import<Configuration>::desc =

    "Import Configuration: \n\n"

    "  'casm import' of Configuration proceeds in two steps: \n\n"

    "  1) For each file: \n"
    "     - Read structure from VASP POSCAR type file or CASM \n"
    "       properties.calc.json file. \n"
    "     - Map the structure onto a Configuration of the primitive crystal \n"
    "       structure. \n"
    "     - Record relaxation data (lattice & basis deformation cost). \n\n"

    "  2) If data import is requested, iterate over each import record and\n\n"
    "     do the following:\n"
    "     - If multiple imported structures map onto a configuration for  \n"
    "       which there is no calculation data, import calculation data   \n"
    "       from the structure with the lowest \"conflict_score\".        \n"
    "     - If one or more imported structures map onto a configuration   \n"
    "       for which calculation data already exist, do not import any   \n"
    "       new data unless the \"overwrite\" option is given, in which   \n"
    "       case data and files are imported if the \"conflict_score\"    \n"
    "       will be improved.                                             \n"
    "     - If data is imported, the corresponding properties.calc.json   \n"
    "       file is copied into the directory of the mapped configuration.\n"
    "       Optionally, additional files in the directory of the imported \n"
    "       structure file may also be copied.                            \n"
    "     - Reports are generated detailing the results of the import:    \n"
    "       - import_map_fail: Structures that could not be mapped onto   \n"
    "         the primitive crystal structure.                            \n"
    "       - import_map_success: Configurations that were successfully   \n"
    "         mapped and imported into the Configuration database (or     \n"
    "         already existed).                                           \n"
    "       - import_data_fail: Structures with data that would be        \n"
    "         imported except preexisting data prevents it.               \n"
    "       - import_conflict: Configurations that were mapped to by      \n"
    "         multiple structures.                                        \n\n"

    "Settings: \n\n"

    "  mapping: JSON object (optional)\n"
    "      A JSON object containing the following options controlling the \n"
    "      structure-mapping algorithm:\n\n"

    "    primitive_only: bool (default=false)\n"
    "        By convention, primitive configurations are always imported  \n"
    "        along with non-primitive configurations. If false, only the  \n"
    "        primitive configuration will be imported. Note: data from    \n"
    "        non-primitive configurations is never used for primitive     \n"
    "        configurations.\n\n"

    "    lattice_weight: number in range [0.0, 1.0] (default=0.5)\n"
    "        Candidate configurations are compared using                  \n"
    "        \"deformation_cost\" to determine the best mapping of the    \n"
    "        import structure to a configuration. The \"lattice_weight\"  \n"
    "        determines the relative weight of the                        \n"
    "        \"lattice_deformation_cost\" and \"atomic_deformation_cost\" \n"
    "        when calculating the total \"deformation_cost\".             \n\n"

    "    max_vol_change: number (default=0.3)\n"
    "        Adjusts range of SCEL volumes searched while mapping imported\n"
    "        structure onto ideal crystal (only necessary if the presence \n"
    "        of vacancies makes the volume ambiguous). Default is +/- 30% \n"
    "        of the rounded integer volume (relaxed volume / primitive    \n"
    "        unit cell volume) of the structure. Smaller values yield     \n"
    "        faster execution, larger values may yield more accurate      \n"
    "        mapping.\n\n"

    "    max_va_frac: number (default=0.5)\n"
    "        Places upper bound on the fraction of sites that are allowed \n"
    "        to be vacant after relaxed structure is mapped onto the ideal\n"
    "        crystal. Smaller values yield faster execution, larger values\n"
    "        may yield more accurate mapping. Has no effect if supercell  \n"
    "        volume can be inferred from the number of atoms in the       \n"
    "        structure. Default value allows up to 50% of sites to be     \n"
    "        vacant.\n\n"

    "    min_va_frac: number (default=0.0)\n"
    "        Places lower bound on the fraction of sites that are allowed \n"
    "        to be vacant after relaxed structure is mapped onto the ideal\n"
    "        crystal. Nonzero values may yield faster execution if        \n"
    "        updating configurations that are known to have a large number\n"
    "        of vacancies, at potential sacrifice of mapping accuracy. Has\n"
    "        no effect if supercell volume can be inferred from the number\n"
    "        of atoms in the structure. Default value allows as few as 0% \n"
    "        of sites to be vacant.\n\n"

    "    ideal: bool (optional, default=false)\n"
    "        Assume imported structures are in the setting of the ideal   \n"
    "        crystal. This results in faster mapping, but may not identify\n"
    "        the ideal mapping. If large mapping costs are encountered,   \n"
    "        try re-running with ideal: false\n\n"

    "    robust: bool (optional, default=false)\n"
    "        Perform additional checks to determine if mapping is         \n"
    "        degenerate in cost to other mappings, which can occur if the \n"
    "        imported structure has symmetry that is incompatible with    \n"
    "        prim.json. Results in slower execution.\n\n"

    "    filter: string (optional) \n"
    "        Restricts the import to only consider supercells that match a\n"
    "        provided casm supercell query expression.\n\n"

    "    forced_lattices: array of strings (optional) \n"
    "        Restricts the import to only consider supercells provided via\n"
    "        a list of their conventional names (i.e.,                    \n"
    "        \"SCEL2_2_1_1_1_1_0\").\n\n"

    "    k_best: int (optional, default=1) \n"
    "        Specify the number, k, of k-best mappings to include in the  \n"
    "        solution set (default is 1).\n\n"

    "  data: JSON object (optional)\n"
    "      A JSON object containing the following options controlling when\n"
    "      calculated properties are updated.\n\n"

    "    import: bool (optional, default=false)\n"
    "        If true (default), attempt to import calculation results.    \n"
    "        Results are added from a \"properties.calc.json\" file which \n"
    "        is checked for in the following locations, relative to the   \n"
    "        input file path, 'pos':                                      \n"
    "        1) Is 'pos' a JSON file? If 'pos' ends in \".json\" or       \n"
    "           \".JSON\", then it is assumed to be a                     \n"
    "           'properties.calc.json' file.                              \n"
    "        2) If / path / to / pos, checks for / path / to / "
    "calctype.current / properties.calc.json\n"
    "        3) If / path / to / pos, checks for / path / to / "
    "properties.calc.json \n\n"

    "        If false, only configuration structure is imported.\n\n"

    "    additional_files: bool (optional, default = false)               \n"
    "        If true, attempt to copy all files & directories in the same \n"
    "        directory as the structure file or , if it exists, the       \n"
    "        properties.calc.json file. Files & directories will only by  \n"
    "        copied if there are no existing files or directories in the  \n"
    "        'training_data' directory for the configuration the structure\n"
    "        as been mapped to or \"overwrite\"=true.\n\n"

    "    overwrite: bool (optional, default=false)\n"
    "        If true, data and files will be imported that overwrite      \n"
    "        existing data and files, if the score calculated by the      \n"
    "        \"conflict_score\" for the configuration being mapped to will\n"
    "        be improved.\n\n"

    "Deformation cost:\n"

    "  The \"deformation_cost\" is:\n\n"

    "      deformation_cost = w*lattice_deformation_cost + \n"
    "                         (1.0-w)*atomic_deformation_cost,\n\n"

    "  where \"w\" is the \"lattice_weight\" factor (default=0.5),\n\n"

    "  the \"lattice_deformation_cost\" is the mean-square displacement of a \n"
    "  points on the surface of a sphere of volume equal to the atomic volume\n"
    "  when  it is deformed by the volume preserving deviatoric deformation \n"
    "  matrix, F_deviatoric:\n\n"

    "      V = relaxed_atomic_volume;\n"
    "      F = deformation matrix (lattice_relaxed = F*lattice_ideal);\n"
    "      F_deviatoric = F/pow(F.determinant(), 1./3.);\n"
    "      I = 3x3 identity matrix;\n\n"

    "      lattice_deformation_cost = pow( 3.*V / (4.*pi), 2.0/3.0) / 3.0 * \n"
    "          (0.5 * (F.t * F / pow(std::abs(F.determinant()), 2.0/3.0) - "
    "I)).squaredNorm()\n\n"

    "  and the \"atomic_deformation_cost\" is a cost function for the amount "
    "of\n"
    "  basis site relaxation:\n\n"

    "      D = 3xN matrix of basis site displacements (displacements are "
    "applied \n"
    "          before strain)\n"
    "      Natoms = number of atoms in configuration\n"
    "      atomic_deformation_cost = (F*D * D.transpose() * "
    "F.transpose()).trace() \n"
    "          / (max(Natoms, 1.))\n\n";

int Import<Configuration>::run(const PrimClex &primclex,
                               const jsonParser &kwargs,
                               const Completer::ImportOption &import_opt) {
  // -- collect input settings --

  std::shared_ptr<Structure const> shared_prim = primclex.shared_prim();
  DataFormatterDictionary<Supercell> const &supercell_query_dict =
      primclex.settings().query_handler<Supercell>().dict();

  ConfigMapping::Settings map_settings;
  jsonParser mapping_json;
  kwargs.get_else(mapping_json, "mapping", jsonParser());
  from_json(map_settings, mapping_json, shared_prim, supercell_query_dict);

  bool primitive_only = false;
  if (kwargs.contains("mapping")) {
    kwargs.get_else(primitive_only, "primitive_only", false);
  }

  ImportSettings import_settings;
  {
    jsonParser combined_json = jsonParser::object();
    if (kwargs.contains("data")) {
      combined_json = kwargs["data"];
    }
    auto const &vm = import_opt.vm();
    if (vm.count("data")) {
      combined_json["import_properties"] = true;
    }
    if (vm.count("copy-additional-files")) {
      combined_json["copy_structure_files"] = true;
      combined_json["copy_additional_files"] = true;
    }
    from_json(import_settings, combined_json);
  }

  // get input report_dir, check if exists, and create new report_dir.i if
  // necessary
  std::string report_dir =
      (fs::path(primclex.dir().reports_dir()) / "import_report").string();
  report_dir = create_report_dir(report_dir);

  // 'mapping' subsettings are used to construct ConfigMapper, and also returns
  // the 'used' settings
  StructureMap<Configuration> mapper(map_settings, primclex, primitive_only);

  jsonParser used_settings;
  used_settings["mapping"] = mapping_json;
  used_settings["data"] = import_settings;

  // -- print used settings --
  Log &log = CASM::log();
  log.read("Settings");
  log << used_settings << std::endl << std::endl;

  // -- construct Import --
  Import<Configuration> f(primclex, mapper, import_settings, report_dir);

  // -- read structure file paths --
  std::set<fs::path> pos;
  auto res =
      construct_pos_paths(primclex, import_opt, std::inserter(pos, pos.end()));
  if (res.second) {
    return res.second;
  }

  // -- read structure file paths --
  f.import(pos.begin(), pos.end());

  return 0;
}

const std::string Update<Configuration>::desc =

    "Update Configuration calculation results: \n\n"

    "  'casm update' of Configuration calculation results proceeds as follows: "
    "\n\n"

    "  For each Configuration in the input selection: \n"
    "   - Read properties.calc.json file from training_data directory.        "
    "\n"
    "   - Map the relaxed structure onto a Configuration of the primitive "
    "crystal\n"
    "     structure. \n"
    "   - Record relaxation data: \n"
    "     - Lattice & basis deformation cost \n"
    "     - Initial configuration and relaxed configuration \n\n"
    "   - If multiple configurations relax onto a configuration for which "
    "there \n"
    "     is no calculation data, the calculation data from the with the "
    "lowest \n"
    "     conflict resolution score is used for the relaxed configuration.\n\n"

    "Settings: \n\n"

    "  force: bool (optional, default=false) \n"
    "    Force update all specified Configuration, else use timestamps to      "
    " \n"
    "    determine which to update. \n"

    "Settings: \n\n"

    "  mapping: JSON object (optional)\n"
    "      A JSON object containing the following options controlling the \n"
    "      structure-mapping algorithm:\n\n"

    "    primitive_only: bool (default=false)\n"
    "        By convention, primitive configurations are always imported  \n"
    "        along with non-primitive configurations. If false, only the  \n"
    "        primitive configuration will be imported. Note: data from    \n"
    "        non-primitive configurations is never used for primitive     \n"
    "        configurations.\n\n"

    "    lattice_weight: number in range [0.0, 1.0] (default=0.5)\n"
    "        Candidate configurations are compared using                  \n"
    "        \"deformation_cost\" to determine the best mapping of the    \n"
    "        import structure to a configuration. The \"lattice_weight\"  \n"
    "        determines the relative weight of the                        \n"
    "        \"lattice_deformation_cost\" and \"atomic_deformation_cost\" \n"
    "        when calculating the total \"deformation_cost\".             \n\n"

    "    max_vol_change: number (default=0.3)\n"
    "        Adjusts range of SCEL volumes searched while mapping imported\n"
    "        structure onto ideal crystal (only necessary if the presence \n"
    "        of vacancies makes the volume ambiguous). Default is +/- 30% \n"
    "        of the rounded integer volume (relaxed volume / primitive    \n"
    "        unit cell volume) of the structure. Smaller values yield     \n"
    "        faster execution, larger values may yield more accurate      \n"
    "        mapping.\n\n"

    "    max_va_frac: number (default=0.5)\n"
    "        Places upper bound on the fraction of sites that are allowed \n"
    "        to be vacant after relaxed structure is mapped onto the ideal\n"
    "        crystal. Smaller values yield faster execution, larger values\n"
    "        may yield more accurate mapping. Has no effect if supercell  \n"
    "        volume can be inferred from the number of atoms in the       \n"
    "        structure. Default value allows up to 50% of sites to be     \n"
    "        vacant.\n\n"

    "    min_va_frac: number (default=0.0)\n"
    "        Places lower bound on the fraction of sites that are allowed \n"
    "        to be vacant after relaxed structure is mapped onto the ideal\n"
    "        crystal. Nonzero values may yield faster execution if        \n"
    "        updating configurations that are known to have a large number\n"
    "        of vacancies, at potential sacrifice of mapping accuracy. Has\n"
    "        no effect if supercell volume can be inferred from the number\n"
    "        of atoms in the structure. Default value allows as few as 0% \n"
    "        of sites to be vacant.\n\n"

    "    ideal: bool (optional, default=false)\n"
    "        Assume imported structures are in the setting of the ideal   \n"
    "        crystal. This results in faster mapping, but may not identify\n"
    "        the ideal mapping. If large mapping costs are encountered,   \n"
    "        try re-running with ideal: false\n\n"

    "    robust: bool (optional, default=false)\n"
    "        Perform additional checks to determine if mapping is         \n"
    "        degenerate in cost to other mappings, which can occur if the \n"
    "        imported structure has symmetry that is incompatible with    \n"
    "        prim.json. Results in slower execution.\n\n"

    "    filter: string (optional) \n"
    "        Restricts the import to only consider supercells that match a\n"
    "        provided casm supercell query expression.\n\n"

    "    forced_lattices: array of strings (optional) \n"
    "        Restricts the import to only consider supercells provided via\n"
    "        a list of their conventional names (i.e.,                    \n"
    "        \"SCEL2_2_1_1_1_1_0\").\n\n"

    "    k_best: int (optional, default=1) \n"
    "        Specify the number, k, of k-best mappings to include in the  \n"
    "        solution set (default is 1).\n\n"

    // "  data: JSON object (optional)\n"
    // "      A JSON object containing the following options controlling how \n"
    // "      calculated properties are updated. After structural relaxation,\n"
    // "      a structure that began as structure 'A' may map more closely   \n"
    // "      onto a different structure, 'B' In some cases, multiple        \n"
    // "      structures may map onto the same 'B', and a scoring metric is  \n"
    // "      used to specify which set of calculation data is associated    \n"
    // "      with configuration 'B'. The following values determine the     \n"
    // "      scoring metric:\n\n"
    //
    // "    \"deformation_cost\":\n"
    // "       \"lattice_weight\": number, in range [0, 1.0]\n"
    //
    // "       Uses a weighted sum of cost functions for lattice and basis \n"
    // "       deformation. See below for complete definition. Ex: \n"
    // "         {\"method\":\"deformation_cost\", \"lattice_weight\":0.5} \n\n"
    //
    // "    \"minimum\":\n"
    // "       \"property\": property name (ex: \"energy\")\n"
    //
    // "       Reads the specified property from the mapped properties and\n"
    // "       selects the minimum to be the best mapping. Ex: \n"
    // "         {\"method\":\"minimum\", \"property\": \"energy\"} \n\n"
    //
    // "    \"maximum\":\n"
    // "       \"property\": property name (ex: \"energy\")\n"
    //
    // "       Reads the specified property from the mapped properties and\n"
    // "       selects the maximum to be the best mapping. Ex: \n"
    // "         {\"method\":\"maximum\", \"property\": \"energy\"} \n\n"
    //
    // "    \"direct_selection\":\n"
    // "       \"name\": configname to force as 'best' (ex: \n"
    // "       \"SCEL3_1_1_3_0_0_0/4\")\n"
    //
    // "       Directly specify which configuration's properties should be used.
    // " "Ex: \n" "         {                                     \n" "
    // \"method\":\"direct_selection\",    \n" "           \"name\":
    // \"SCEL3_1_1_3_0_0_0/4\"   \n" "         } \n\n"
    //
    // "  The default value used by CASM is:                                 \n"
    // "    {\"method\": \"minimum\", \"property\": \"energy\"} \n\n"

    "Deformation cost:\n"

    "  The \"deformation_cost\" is:\n\n"

    "      deformation_cost = w*lattice_deformation_cost + \n"
    "                         (1.0-w)*atomic_deformation_cost,\n\n"

    "  where \"w\" is the \"lattice_weight\" factor (default=0.5),\n\n"

    "  the \"lattice_deformation_cost\" is the mean-square displacement of a \n"
    "  points on the surface of a sphere of volume equal to the atomic volume\n"
    "  when  it is deformed by the volume preserving deviatoric deformation \n"
    "  matrix, F_deviatoric:\n\n"

    "      V = relaxed_atomic_volume;\n"
    "      F = deformation matrix (lattice_relaxed = F*lattice_ideal);\n"
    "      F_deviatoric = F/pow(F.determinant(), 1./3.);\n"
    "      I = 3x3 identity matrix;\n\n"

    "      lattice_deformation_cost = pow( 3.*V / (4.*pi), 2.0/3.0) / 3.0 * \n"
    "          (0.5 * (F.t * F / pow(std::abs(F.determinant()), 2.0/3.0) - "
    "I)).squaredNorm()\n\n"

    "  and the \"atomic_deformation_cost\" is a cost function for the amount "
    "of\n"
    "  basis site relaxation:\n\n"

    "      D = 3xN matrix of basis site displacements (displacements are "
    "applied \n"
    "          before strain)\n"
    "      Natoms = number of atoms in configuration\n"
    "      atomic_deformation_cost = (F*D * D.transpose() * "
    "F.transpose()).trace() \n"
    "          / (max(Natoms, 1.))\n\n";

/// Allow ConfigType to specialize the report formatting for 'import'
DataFormatter<ConfigIO::Result> Import<Configuration>::_import_formatter()
    const {
  DataFormatterDictionary<ConfigIO::Result> dict;
  ConfigIO::default_import_formatters(dict, db_props());

  std::vector<std::string> col = {"initial_path",
                                  "selected",
                                  "to_configname",
                                  "final_path",
                                  "is_new_config",
                                  "has_all_required_properties",
                                  "preexisting_properties",
                                  "preexisting_files",
                                  "did_insert_properties",
                                  "did_copy_structure_file",
                                  "did_copy_additional_files",
                                  "score",
                                  "best_score",
                                  "has_best_scoring_mapped_properties",
                                  "lattice_deformation_cost",
                                  "atomic_deformation_cost",
                                  "energy"};

  return dict.parse(col);
}

// --- Update<Configuration> ---

/// \brief Constructor
Update<Configuration>::Update(const PrimClex &primclex,
                              const StructureMap<Configuration> &mapper,
                              UpdateSettings const &_set,
                              std::string const &report_dir)
    : UpdateT(primclex, mapper, _set, report_dir) {}

int Update<Configuration>::run(const PrimClex &primclex,
                               const jsonParser &kwargs,
                               const Completer::UpdateOption &update_opt) {
  // -- collect input settings --

  const po::variables_map &vm = update_opt.vm();
  jsonParser used;

  std::shared_ptr<Structure const> shared_prim = primclex.shared_prim();
  DataFormatterDictionary<Supercell> const &supercell_query_dict =
      primclex.settings().query_handler<Supercell>().dict();

  // general settings
  bool force;
  if (vm.count("force")) {
    force = true;
  } else {
    kwargs.get_else(force, "force", false);
  }
  used["force"] = force;

  // TODO: this could take more settings, for now output_as_json fixed true
  UpdateSettings update_settings;

  // 'mapping' subsettings are used to construct ConfigMapper and return 'used'
  // settings values still need to figure out how to specify this in general
  ConfigMapping::Settings map_settings;
  jsonParser mapping_json;
  kwargs.get_else(mapping_json, "mapping", jsonParser());
  from_json(map_settings, mapping_json, shared_prim, supercell_query_dict);

  bool primitive_only = false;
  if (kwargs.contains("mapping")) {
    kwargs.get_else(primitive_only, "primitive_only", false);
  }

  StructureMap<Configuration> mapper(map_settings, primclex, primitive_only);
  used["mapping"] = mapping_json;

  // 'data' subsettings
  // bool import_data = true;
  // bool import_additional_files = false;
  // bool overwrite = false;

  // -- print used settings --
  Log &log = CASM::log();
  log.read("Settings");
  log << used << std::endl << std::endl;

  // get input report_dir, check if exists, and create new report_dir.i if
  // necessary
  std::string report_dir =
      (fs::path(primclex.dir().reports_dir()) / "update_report").string();
  report_dir = create_report_dir(report_dir);

  // -- construct Update --
  Update<Configuration> f(primclex, mapper, update_settings, report_dir);

  // -- read selection --
  DB::Selection<Configuration> sel(primclex, update_opt.selection_path());

  // -- update --
  f.update(sel, force);

  return 0;
}

// Allow ConfigType to specialize the report formatting for 'update'
DataFormatter<ConfigIO::Result> Update<Configuration>::_update_formatter()
    const {
  DataFormatterDictionary<ConfigIO::Result> dict;
  ConfigIO::default_update_formatters(dict, db_props());

  std::vector<std::string> col = {"data_origin",
                                  "selected",
                                  "to_configname",
                                  "has_all_required_properties",
                                  "score",
                                  "best_score",
                                  "has_best_scoring_mapped_properties",
                                  "lattice_deformation_cost",
                                  "atomic_deformation_cost",
                                  "energy"};

  return dict.parse(col);
}

}  // namespace DB
}  // namespace CASM
