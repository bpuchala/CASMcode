#include "casm/clex/io/json/ConfigMapping_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/casm_io/dataformatter/DataFormatter_impl.hh"
#include "casm/casm_io/enum/io_traits.hh"
#include "casm/clex/ConfigMapping.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/io/json/Configuration_json_io.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/io/SimpleStructureIO.hh"
#include "casm/crystallography/io/json/StrucMapping_json_io.hh"
#include "casm/global/definitions.hh"
#include "casm/symmetry/SupercellSymInfo.hh"

namespace CASM {

//--ConfigMapping::Settings------------------------------------------------

// jsonParser &to_json(ConfigMapping::Settings const &_set, jsonParser &_json) {
//   _json["lattice_weight"] = _set.lattice_weight;
//   _json["ideal"] = _set.ideal;
//   _json["strict"] = _set.strict;
//   // _json["primitive_only"] = _set.primitive_only;
//   _json["robust"] = _set.robust;
//   _json["fix_volume"] = _set.fix_volume;
//   _json["fix_lattice"] = _set.fix_lattice;
//   _json["k_best"] = _set.k_best;
//   if (!_set.forced_lattices.empty())
//     _json["forced_lattices"] = _set.forced_lattices;
//   if (!_set.filter.empty()) _json["filter"] = _set.filter;
//   _json["cost_tol"] = _set.cost_tol;
//   _json["min_va_frac"] = _set.min_va_frac;
//   _json["max_va_frac"] = _set.max_va_frac;
//   _json["max_vol_change"] = _set.max_vol_change;
//
//   return _json;
// }

/// Read ConfigMapping::Settings from JSON
///
/// \param _set ConfigMapping::Settings to be assigned from JSON
/// \param _json jsonParser JSON Input
/// \param shared_prim Prim used for parent structure
/// \param supercell_query_dict DataFormatterDictionary<Supercell> for parsing
///     filter expressions
jsonParser const &from_json(
    ConfigMapping::Settings &_set, jsonParser const &_json,
    std::shared_ptr<Structure const> const &shared_prim,
    DataFormatterDictionary<Supercell> const &supercell_query_dict) {
  _set.set_default();

  if (_json.contains("lattice_weight"))
    _set.lattice_weight = _json["lattice_weight"].get<double>();

  if (_json.contains("ideal")) _set.ideal = _json["ideal"].get<bool>();

  if (_json.contains("strict")) _set.strict = _json["strict"].get<bool>();

  if (_json.contains("robust")) _set.robust = _json["robust"].get<bool>();

  // if (_json.contains("primitive_only"))
  //   _set.primitive_only = _json["primitive_only"].get<bool>();

  if (_json.contains("fix_volume"))
    _set.fix_volume = _json["fix_volume"].get<bool>();

  if (_json.contains("fix_lattice"))
    _set.fix_lattice = _json["fix_lattice"].get<bool>();

  if (_json.contains("k_best")) {
    _set.k_best = _json["k_best"].get<Index>();
  }

  if (_json.contains("forced_lattices")) {
    // _set.forced_lattices =
    //     _json["forced_lattices"].get<std::vector<std::string> >();

    _set.forced_lattices.clear();

    std::vector<std::string> scelnames =
        _json["forced_lattices"].get<std::vector<std::string> >();

    for (std::string scelname : scelnames) {
      xtal::Superlattice super_lat = make_superlattice_from_supercell_name(
          shared_prim->factor_group(), shared_prim->lattice(), scelname);
      _set.forced_lattices.push_back(super_lat.superlattice());
    }
  }

  if (_json.contains("filter")) {
    /// If a filter string is specified, construct the Supercell query filter
    /// for restricting potential supercells for mapping

    // _set.filter = _json["filter"].get<std::string>();
    // dict = primclex.settings().query_handler<Supercell>().dict();

    std::string filter_expression = _json["filter"].get<std::string>();

    DataFormatter<Supercell> formatter =
        supercell_query_dict.parse(filter_expression);
    auto filter = [formatter, shared_prim](Lattice const &parent,
                                           Lattice const &child) -> bool {
      ValueDataStream<bool> check_stream;
      check_stream << formatter(Supercell(shared_prim, parent));
      return check_stream.value();
    };
  }

  if (_json.contains("cost_tol"))
    _set.cost_tol = _json["cost_tol"].get<double>();

  if (_json.contains("min_va_frac"))
    _set.min_va_frac = _json["min_va_frac"].get<double>();

  if (_json.contains("max_va_frac"))
    _set.max_va_frac = _json["max_va_frac"].get<double>();

  if (_json.contains("max_vol_change"))
    _set.max_vol_change = _json["max_vol_change"].get<double>();

  return _json;
}

ENUM_TRAITS(ConfigMapperResult::HintStatus)

namespace {
jsonParser &to_json(ConfigMapperResult::Individual const &individual_result,
                    jsonParser &json, COORD_TYPE coordinate_mode) {
  json["configuration"] = individual_result.config;
  to_json(individual_result.resolved_struc, json["structure"], {},
          coordinate_mode);
  if (individual_result.hint_status != ConfigMapperResult::HintStatus::None) {
    json["hint_status"] = to_string(individual_result.hint_status);
    json["hint_cost"] = individual_result.hint_cost;
  }
  return json;
}
}  // anonymous namespace

jsonParser &to_json(ConfigMapperResult const &config_mapping_result,
                    jsonParser &json, COORD_TYPE coordinate_mode) {
  json["success"] = config_mapping_result.success();
  json["n_optimal_mappings"] = config_mapping_result.n_optimal();
  if (!config_mapping_result.success()) {
    json["fail_msg"] = config_mapping_result.fail_msg;
  } else {
    json["maps"] = jsonParser::array();
    for (auto const &mapping : config_mapping_result.maps) {
      jsonParser mapping_json;

      // xtal::MappingNode
      to_json(mapping.first, mapping_json["mapping"]);

      // ConfigMapperResult::Individual
      to_json(mapping.second, mapping_json["result"], coordinate_mode);

      json["maps"].push_back(mapping_json);
    }
  }

  return json;
}

const std::string traits<ConfigMapperResult::HintStatus>::name = "hint_status";

const std::multimap<ConfigMapperResult::HintStatus, std::vector<std::string> >
    traits<ConfigMapperResult::HintStatus>::strval = {
        {ConfigMapperResult::HintStatus::None, {"None"}},
        {ConfigMapperResult::HintStatus::Derivative, {"Derivative"}},
        {ConfigMapperResult::HintStatus::Equivalent, {"Equivalent"}},
        {ConfigMapperResult::HintStatus::Identical, {"Identical"}},
        {ConfigMapperResult::HintStatus::Derivative, {"NewOcc"}},
        {ConfigMapperResult::HintStatus::Derivative, {"NewScel"}}};

}  // namespace CASM
