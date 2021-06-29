#include "casm/clex/io/json/ConfigMapping_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/casm_io/dataformatter/DataFormatter_impl.hh"
#include "casm/casm_io/enum/io_traits.hh"
#include "casm/casm_io/json/InputParser_impl.hh"
#include "casm/clex/ConfigMapping.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/io/json/Configuration_json_io.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/io/SimpleStructureIO.hh"
#include "casm/crystallography/io/json/StrucMapping_json_io.hh"
#include "casm/global/definitions.hh"
#include "casm/symmetry/SupercellSymInfo.hh"
#include "casm/symmetry/io/json/SuperlatticeIO_json_io.hh"

namespace CASM {

/// Parse ConfigMapperSettings from JSON
///
/// \param _set ConfigMapperSettings to be assigned from JSON
/// \param _json jsonParser JSON Input
/// \param shared_prim Prim used for parent structure
/// \param supercell_query_dict DataFormatterDictionary<Supercell> for parsing
///     filter expressions
///
/// Lattice mapping method options (choose one, default is
/// auto_lattice_volume_range):
///
/// fix_ideal: bool (optional)
///     Force lattice mapping solutions of the form `L1 * T = Q * L2`, where
///         L1 = the prim lattice
///         L2 = child.lat_column_mat
///          T = an integer transformation matrix, to be determined
///          Q = isometry matrix, constrained to be an elemenet of the prim
///              factor group, to be determined
///
/// fix_lattice_mapping: bool (optional)
///     Force lattice mapping solutions of the form `S1 = V * Q * L2`, where
///         S1 = this->configuration_lattice.value()
///         L2 = child.lat_column_mat
///       V, Q = stretch and isometry matrices, to be determined
///
/// fix_lattice: bool (optional)
///     Force lattice mapping solutions of the form `S1 * N = V * Q * L2`, where
///         S1 = this->configuration_lattice.value()
///         L2 = child.lat_column_mat
///          N = unimodular matrix, generates non-identical but equivalent
///              parent superlattices, to be determined
///       V, Q = stretch and isometry matrices, to be determined
///
/// fix_lattice_volume_range: bool (optional)
///     Force lattice mapping solutions of the form `S1 * N = V * Q * L2`, where
///         this->min_lattice_volume <=
///             S1.determinant() <= this->max_lattice_volume,
///         S1 = lattice equivalent to the configuration lattice, to be
///              determined
///         L2 = child.lat_column_mat
///          N = unimodular matrix, generates non-identical but equivalent
///              parent superlattices, to be determined
///       V, Q = stretch and isometry matrices, to be determined
///
/// auto_lattice_volume_range: bool (optional)
///     Equivalent to `fix_lattice_volume_range`, with a volume range that is
///     determined by constaints on the Va fraction and the change in volume
///     between the input structure and mapped structure.
///
///
/// Structure mapping parameters:
///
/// use_symmetry_breaking_strain_cost: bool (default=false)
///     If true, use the symmetry-breaking strain cost.
///
/// use_symmetry_breaking_displacement_cost: bool (default=false)
///     If true, use the symmetry-breaking displacement cost.
///
/// lattice_weight: number (default=0.5)
///     Specifies the cost function in terms of lattice deformation
///     cost and atomic deformation cost:
///         total_cost = lattice_weight*lattice_deformation_cost +
///                      (1-lattice_weight)*atomic_deformation_cost
///
/// cost_tol: number (default=1e-5)
///     Tolerance used to determine if two mappings have identical cost.
///
/// robust: bool (default=false)
///     True invokes additional checks which determine whether there are any
///     other mappings that are distinct from the best mapping, but have the
///     same cost.
///
/// k_best: int (default=1)
///     Specify the number, k, of k-best mappings to include in solution set
///
/// configuration_lattice: Superlattice (optional)
///     Superlattice of the prim used by `fix_lattice` and
///     `fix_lattice_mapping` methods. This is a required value if one
///     of those methods is chosen. It must be a superlattice of the prim
///     lattice. It does not have to be canonical. The Superlattice format is
///     described below.
///
/// allowed_lattices: array of Superlattice (optional)
///     List of superlattices of the prim to consider when searching for
///     mappings. If any provided, only these lattices will be considered, and
///     only if they fall within the volume range being searched. The
///     Superlattice format is described below.
///
/// filter: str (optional)
///     If a filter string is specified, it is interpreted as a `casm query`-
///     style expression used to query supercells. If provided, used to filter
///     list of potential supercells of the parent structure to search over.
///
///     Applies to the `fixed_lattice_volume_range` and
///     `auto_lattice_volume_range` lattice mapping methods. The filter is also
///     applied to lattices specified by `allowed_lattices`.
///
/// min_lattice_volume: int (default=1)
///     The minimum configuration lattice volume considered. Applies to the
///     `fixed_lattice_volume_range` lattice mapping method only.
///
/// max_lattice_volume: int (default=1)
///     The maximum configuration lattice volume considered. Applies to the
///     `fixed_lattice_volume_range` lattice mapping method only.
///
/// min_va_frac: number (default=0.)
///     The minimum fraction of vacant sites. Below this fraction a mapping will
///     not be considered. Applies to the `auto_lattice_volume_range` lattice
///     mapping method only.
///
/// max_va_frac: number (default=1.)
///     The maximum fraction of vacant sites. Above this fraction a mapping will
///     not be considered. Applies to the `auto_lattice_volume_range` lattice
///     mapping method only.
///
/// max_volume_change: number (default=0.3)
///     Constrains the search space by setting a limit on allowed volume change.
///     Applies to the `auto_lattice_volume_range` lattice mapping method only.
///
///
/// Equivalent configuration selection parameters:
///
/// finalize_strict: bool (default=false)
///     If true, as a post-processing step, the "mapped configuration" is
///     transformed to a "final configuration", in the canonical supercell,
///     which preserves the orientation of the original unmapped structure as
///     much as possible. Otherwise, the "final configuration" is the canonical
///     equivalent configuration in the canonical supercell.
///
///
/// Reading superlattices for `configuration_lattice` or `allowed_lattices`:
///
/// Superlattices are expected to be JSON objects, from which a Superlattice
/// can be specified using any of the following attributes. When searching in
/// the order listed below, the first found will be used. If none found, the
/// prim lattice is used.
///
///   transformation_matrix_to_super: 3x3 array of integer (optional)
///     Transformation matrix T, defining the supercell lattice vectors
///     S, in terms of the prim lattice vectors, P: `S = P * T`, where
///     S and P are column vector matrices.
///
///   supercell_lattice_row_vectors: 3x3 array of integer (optional)
///     Supercell lattice vectors, as a row vector matrix.
///
///   supercell_lattice_column_matrix: 3x3 array of integer (optional)
///     Supercell lattice vectors, as a column vector matrix.
///
///   supercell_name: string (optional)
///     A name given to all equivalent super lattices of the prim
///     lattice. For the canonical super lattice, the name is
///     constructed from the hermite normal form of
///     `transformation_matrix_to_super`. For a non-canonical super
///     lattice, the name is the constructed from the name of the
///     canonical super lattice and the index of the prim factor group
///     operation that (excluding the shift) transforms the canonical
///     super lattice into this super lattice.
///
///     Example 1: Canonical supercell name
///
///         \"SCEL8_4_2_1_1_3_2\"
///
///     Example 2: Non-canonical supercell name representing a
///     re-orientation by application of prim factor group operation
///     with index 2 (indexing starting at 0).
///
///         \"SCEL8_4_2_1_1_3_2.2\"
///
///   make_canonical: bool (optional, default=false)
///     If \"true\", the canonical equivalent supercell is used.
///
void parse(InputParser<ConfigMapperSettings> &parser,
           std::shared_ptr<Structure const> const &shared_prim,
           DataFormatterDictionary<Supercell> const &supercell_query_dict) {
  ConfigMapperSettings settings;

  // parse lattice mapping method:
  parser.optional(settings.fix_ideal, "fix_ideal");
  parser.optional(settings.fix_lattice_mapping, "fix_lattice_mapping");
  parser.optional(settings.fix_lattice, "fix_lattice");
  parser.optional(settings.fix_lattice_volume_range,
                  "fix_lattice_volume_range");
  parser.optional(settings.auto_lattice_volume_range,
                  "auto_lattice_volume_range");

  int count = 0;
  count += (settings.fix_ideal == true);
  count += (settings.fix_lattice_mapping == true);
  count += (settings.fix_lattice == true);
  count += (settings.fix_lattice_volume_range == true);
  count += (settings.auto_lattice_volume_range == true);
  if (count == 0) {
    settings.auto_lattice_volume_range = true;
  } else if (count > 1) {
    std::stringstream msg;
    msg << "Error: One or zero lattice mapping method may be used. Options are "
           "\"fix_ideal\", \"fix_lattice_mapping\", \"fix_lattice\", "
           "\"fix_lattice_volume_range\", \"auto_lattice_volume_range\" "
           "(default).";
    parser.error.insert(msg.str());
  }

  // parse structure mapping parameters:
  parser.optional(settings.use_symmetry_breaking_strain_cost,
                  "use_symmetry_breaking_strain_cost");
  parser.optional(settings.use_symmetry_breaking_displacement_cost,
                  "use_symmetry_breaking_displacement_cost");
  parser.optional(settings.lattice_weight, "lattice_weight");
  parser.optional(settings.cost_tol, "cost_tol");
  parser.optional(settings.robust, "robust");
  parser.optional(settings.k_best, "k_best");
  parser.optional(settings.min_lattice_volume, "min_lattice_volume");
  parser.optional(settings.max_lattice_volume, "max_lattice_volume");
  parser.optional(settings.min_va_frac, "min_va_frac");
  parser.optional(settings.max_va_frac, "max_va_frac");
  parser.optional(settings.max_volume_change, "max_volume_change");

  auto superlattice_subparser =
      parser.subparse_if<SuperlatticeIO>("configuration_lattice", shared_prim);
  if (superlattice_subparser->value != nullptr) {
    settings.configuration_lattice =
        superlattice_subparser->value->superlattice.superlattice();
  }

  if (parser.self.contains("allowed_lattices")) {
    settings.allowed_lattices.clear();

    if (!parser.self["allowed_lattices"].is_array()) {
      std::stringstream msg;
      msg << "Error: \"allowed_lattices\" must be an array.";
      parser.insert_error("allowed_lattices", msg.str());
    }

    for (Index i = 0; i < parser.self["allowed_lattices"].size(); ++i) {
      fs::path option = "allowed_lattices";
      option = option / std::to_string(i);
      auto superlattice_subparser =
          parser.subparse_if<SuperlatticeIO>(option, shared_prim);
      if (superlattice_subparser->value != nullptr) {
        settings.allowed_lattices.push_back(
            superlattice_subparser->value->superlattice.superlattice());
      }
    }
  }

  if (parser.self.contains("filter")) {
    /// If a filter string is specified, construct the Supercell query filter
    /// for restricting potential supercells for mapping

    // _set.filter = _json["filter"].get<std::string>();
    // dict = primclex.settings().query_handler<Supercell>().dict();
    std::string filter_expression;
    parser.optional(filter_expression, "filter");

    try {
      DataFormatter<Supercell> formatter =
          supercell_query_dict.parse(filter_expression);
      auto filter = [formatter, shared_prim](Lattice const &parent,
                                             Lattice const &child) -> bool {
        ValueDataStream<bool> check_stream;
        check_stream << formatter(Supercell(shared_prim, parent));
        return check_stream.value();
      };
      settings.filter = filter;
    } catch (std::exception &e) {
      std::stringstream msg;
      msg << "Error: could not construct filter from expression: '"
          << filter_expression << "'.";
      parser.insert_error("filter", msg.str());
    }
  }

  // parse equivalent configuration selection parameters:
  parser.optional(settings.finalize_strict, "finalize_strict");

  parser.value = notstd::clone(settings);
}

ENUM_TRAITS(ConfigComparison)

namespace {

jsonParser &to_json(SymOp const &op, jsonParser &json) {
  to_json(get_matrix(op), json["matrix"]);
  to_json_array(get_translation(op), json["translation"]);
  to_json(get_time_reversal(op), json["time_reversal"]);
  return json;
}

jsonParser &to_json(ConfigurationMapping const &mapping, jsonParser &json,
                    COORD_TYPE coordinate_mode) {
  to_json(mapping.mapping, json["mapping"]);  // xtal::MappingNode
  to_json(mapping.mapped_child, json["mapped"]["structure"], {},
          coordinate_mode);
  json["mapped"]["configuration"] = mapping.mapped_configuration;
  json["mapped"]["properties"] = mapping.mapped_properties;

  to_json(mapping.symop_to_canon_scel, json["symop_to_canon_scel"]);
  to_json(mapping.permutation_to_canon_scel.perm_array(),
          json["permutation_to_canon_scel"]);
  json["transformation_matrix_to_canon_scel"] =
      mapping.transformation_matrix_to_canon_scel;
  to_json(mapping.structure_in_canon_scel,
          json["in_canonical_supercell"]["structure"], {}, coordinate_mode);
  json["in_canonical_supercell"]["configuration"] =
      mapping.configuration_in_canon_scel;
  json["in_canonical_supercell"]["properties"] =
      mapping.properties_in_canon_scel;

  to_json(mapping.symop_to_final, json["symop_to_final"]);
  to_json(mapping.permutation_to_final.perm_array(),
          json["permutation_to_final"]);
  to_json(mapping.final_structure, json["final"]["structure"], {},
          coordinate_mode);
  json["final"]["configuration"] = mapping.final_configuration;
  json["final"]["properties"] = mapping.final_properties;

  if (mapping.hint_status != ConfigComparison::None) {
    json["hint_status"] = to_string(mapping.hint_status);
    json["hint_cost"] = mapping.hint_cost;
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
      to_json(mapping, mapping_json, coordinate_mode);
      json["maps"].push_back(mapping_json);
    }
  }

  return json;
}

const std::string traits<ConfigComparison>::name = "hint_status";

const std::multimap<ConfigComparison, std::vector<std::string> >
    traits<ConfigComparison>::strval = {
        {ConfigComparison::None, {"None"}},
        {ConfigComparison::Derivative, {"Derivative"}},
        {ConfigComparison::Equivalent, {"Equivalent"}},
        {ConfigComparison::Identical, {"Identical"}},
        {ConfigComparison::Derivative, {"NewOcc"}},
        {ConfigComparison::Derivative, {"NewScel"}}};

}  // namespace CASM
