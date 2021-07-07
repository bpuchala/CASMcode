#ifndef CASM_ConfigMapping
#define CASM_ConfigMapping

#include "casm/clex/Configuration.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/SimpleStrucMapCalculator.hh"
#include "casm/crystallography/StrucMapping.hh"
#include "casm/global/definitions.hh"
#include "casm/misc/CASM_math.hh"

namespace CASM {
namespace xtal {
class Lattice;
class Site;
class BasicStructure;
struct MappingNode;
class SimpleStructure;
class SimpleStrucMapCalculator;
class StrucMapper;
}  // namespace xtal

class Supercell;
class PermuteIterator;
class PrimClex;
class Configuration;
class ConfigDoF;
class SupercellSymInfo;

/// \brief Settings for ConfigMapper
struct ConfigMapperSettings {
  ConfigMapperSettings()
      : lattice_weight(0.5),
        cost_tol(TOL),
        k_best(1),
        min_lattice_volume(1),
        max_lattice_volume(1),
        min_va_frac(0.),
        max_va_frac(1.),
        max_volume_change(0.3) {}

  // --- Lattice mapping method choices: ---

  /// \brief Force lattice mapping assuming no deformation
  ///
  /// Force lattice mapping solutions of the form \f$L_1 * T_1 = Q * L_2\f$,
  /// where
  ///
  /// - \f$L_1\f$ = the prim lattice, as a column vector matrix
  /// - \f$L_2\f$ = `child.lat_column_mat`, the unmapped child structure's
  ///   lattice vectors, as a column vector matrix
  /// - \f$T_1\f$ = an integer transformation matrix, to be determined
  /// - \f$Q\f$ = isometry matrix, constrained to be an elemenet of the prim
  ///   factor group, to be determined
  bool fix_ideal = false;

  /// \brief Force lattice mapping to particular lattice vectors
  ///
  /// Force lattice mapping solutions of the form \f$S_1 = V * Q * L_2\f$, where
  ///
  /// - \f$S_1\f$ = `this->configuration_lattice.value()`, lattice vectors, as
  ///   a column vector matrix, of an exact supercell of the prim lattice
  ///   vectors
  /// - \f$L_2\f$ = `child.lat_column_mat`, the unmapped child structure's
  ///   lattice vectors, as a column vector matrix
  /// - \f$V\f$, \f$Q\f$ = stretch and isometry matrices, to be determined
  bool fix_lattice_mapping = false;

  /// \brief Force lattice mapping to an equivalent of a particular lattice
  ///
  /// Force lattice mapping solutions of the form \f$S_1 * N = V * Q * L_2\f$,
  /// where
  ///
  /// - \f$S_1\f$ = `this->configuration_lattice.value()`, lattice vectors, as
  ///   a column vector matrix, of an exact supercell of the prim lattice
  ///   vectors
  /// - \f$L_2\f$ = `child.lat_column_mat`, the unmapped child structure's
  ///   lattice vectors, as a column vector matrix
  /// - \f$N\f$ = unimodular matrix, generates non-identical but equivalent
  ///   parent superlattices, to be determined
  /// - \f$V\f$, \f$Q\f$ = stretch and isometry matrices, to be determined
  bool fix_lattice = false;

  /// \brief Force lattice mapping to a range of superlattice volumes
  ///
  /// Force lattice mapping solutions of the form \f$L_1 * T_1 * N = V * Q *
  /// L_2\f$, where
  ///
  ///     this->min_lattice_volume <=
  ///         (L_1 * T_1 * N).determinant() <= this->max_lattice_volume,
  ///
  /// - \f$L_1\f$ = the prim lattice, as a column vector matrix
  /// - \f$L_2\f$ = `child.lat_column_mat`, the unmapped child structure's
  ///   lattice vectors, as a column vector matrix
  /// - \f$T_1\f$ = an integer transformation matrix, to be determined
  /// - \f$N\f$ = unimodular matrix, generates non-identical but equivalent
  ///   parent superlattices, to be determined
  /// - \f$V\f$, \f$Q\f$ = stretch and isometry matrices, to be determined
  bool fix_lattice_volume_range = false;

  /// \brief Force lattice mapping to a range of superlattice volumes
  /// determined by specifying an allowed Va concentration range and a maximum
  /// allowed volume change
  ///
  /// Equivalent to `fix_lattice_volume_range`, with a volume range that is
  /// determined by constaints on the Va fraction and the change in volume
  /// between the input structure and mapped structure.
  bool auto_lattice_volume_range = false;

  // --- Structure mapping parameters: ---

  /// If true, use the symmetry-breaking strain cost.
  bool use_symmetry_breaking_strain_cost = false;

  /// If true, use the symmetry-breaking displacement cost.
  bool use_symmetry_breaking_displacement_cost = false;

  /// \brief Set the lattice deformation cost weight (= 1 - atomic deformation
  /// cost weight) in the total cost
  ///
  /// Specifies the cost function in terms of lattice deformation
  /// cost and atomic deformation cost:
  ///
  ///     total_cost = lattice_weight*lattice_deformation_cost +
  ///                  (1-lattice_weight)*atomic_deformation_cost
  ///
  /// Default: 0.5
  double lattice_weight;

  /// Tolerance used to determine if two mappings have identical cost
  ///
  /// Default: TOL
  double cost_tol;

  /// \brief True invokes additional checks which determine whether there are
  /// any other mappings that are distinct from the best mapping, but have the
  /// same cost
  ///
  /// Default: false
  bool robust = false;

  /// Specify the target number of mappings to include in solution set
  ///
  /// Note that the actual number of mappings returned may vary from k, if
  /// valid mappings are not found or mapping costs are tied.
  ///
  /// Default: 1
  Index k_best;

  /// Lattice used by `fix_lattice` and `fix_lattice_mapping` methods.
  ///
  /// This is required to have a value of one of those methods is chosen. It
  /// must be a superlattice of the prim lattice. It does not have to be
  /// canonical.
  /// Default: empty
  std::optional<Lattice> configuration_lattice;

  /// \brief List of particular superlattices of the parent structure to
  /// consider when searching for mappings within a range of volumes.
  ///
  /// If any provided, only these lattices will be considered, and
  /// only if they fall within the volume range being searched.
  ///
  /// Applies to the `fixed_lattice_volume_range` and
  /// `auto_lattice_volume_range` lattice mapping methods.
  /// Default: empty
  std::vector<Lattice> allowed_lattices;

  /// \brief If provided, used to filter list of potential supercells of the
  /// parent structure to search over when searching for mappings within a
  /// range of volumes.
  ///
  /// The filter function is of the form
  /// `bool filter(parent_prim_lattice, proposed_parent_superlattice)`, where
  /// `parent_prim_lattice` is the lattice of the primitive parent structure,
  /// and `proposed_parent_superlattice` is a proposed superlattice of the
  /// parent structure.
  ///
  /// Applies to the `fixed_lattice_volume_range` and
  /// `auto_lattice_volume_range` lattice mapping methods. The filter is also
  /// applied to lattices specified by `allowed_lattices`.
  /// Default: empty
  xtal::StrucMapping::LatticeFilterFunction filter;

  /// The minimum configuration lattice volume considered.
  ///
  /// Applies to the `fixed_lattice_volume_range` lattice mapping method only.
  /// Default: 1
  Index min_lattice_volume;

  /// The maximum configuration lattice volume considered.
  ///
  /// Applies to the `fixed_lattice_volume_range` lattice mapping method only.
  /// Default: 1
  Index max_lattice_volume;

  /// The minimum fraction of vacant sites.
  ///
  /// Below this fraction a lattice mapping will not be considered.
  ///
  /// Applies to the `auto_lattice_volume_range` lattice mapping method only.
  /// Default: 0.
  double min_va_frac;

  /// The maximum fraction of vacant sites.
  ///
  /// Above this fraction a lattice mapping will not be considered
  ///
  /// Applies to the `auto_lattice_volume_range` lattice mapping method only.
  /// Default: 1.
  double max_va_frac;

  /// Constrains the search space by setting a limit on allowed volume change.
  ///
  /// Applies to the `auto_lattice_volume_range` lattice mapping method only.
  /// Default: 0.3
  double max_volume_change;

  /// If true, ensures that if no supercell volume satisfies vacancy
  /// constraints, the smallest possible volume is used. Otherwise, the result
  /// will be no valid mapping.
  ///
  /// Applies to the `auto_lattice_volume_range` lattice mapping method only.
  /// Default: false
  bool soft_va_limit = false;

  // TODO: specify a particular translation
  // std::optional<Eigen::Vector3d> fix_translation;

  // --- Equivalent configuration selection parameters: ---

  /// If true, as a post-processing step, the `mapped_configuration` is
  /// transformed to a `final_configuration` in the canonical supercell which
  /// preserves the orientation of the original unmapped structure as much as
  /// possible. Otherwise, the `final_configuration` is the canonical
  /// equivalent configuration in the canonical supercell
  bool finalize_strict = false;
};

/// Specify degree to which the hinted configuration matches the final mapped
/// configuration:
/// - None : unspecified/unknown
/// - Derivate : same occupation, but other DoFs are different
/// - Equivalent : same occupation and DoFs, but related to Hint by a SymOp
///   of the ideal crystal
/// - Identical : Exact same occupation and DoFs
/// - NewOcc : Occupation has no relation to Hint (no statement about other
///   DoFs, but they should be presumed different)
/// - NewScel : Mapped configuration corresponds to different supercell than
///   Hint
enum class ConfigComparison {
  None,
  Derivative,
  Equivalent,
  Identical,
  NewOcc,
  NewScel
};

/// Data structure storing an individual ConfigMapper solution
///
/// See member comments for definitions. See ConfigMapper for more detailed
/// context.
struct ConfigurationMapping {
  ConfigurationMapping(
      xtal::SimpleStructure const &_unmapped_child,
      xtal::MappingNode const &_mapping,
      xtal::SimpleStructure const &_mapped_child,
      Configuration const &_mapped_configuration,
      MappedProperties const &_mapped_properties,
      SymOp const &_symop_to_canon_scel,
      Permutation const &_permutation_to_canon_scel,
      Eigen::Matrix3l const &_transformation_matrix_to_canon_scel,
      xtal::SimpleStructure const &_structure_in_canon_scel,
      Configuration const &_configuration_in_canon_scel,
      MappedProperties const &_properties_in_canon_scel,
      SymOp const &_symop_to_final, Permutation const &_permutation_to_final,
      Eigen::Matrix3l const &_transformation_matrix_to_final,
      xtal::SimpleStructure const &_final_structure,
      Configuration const &_final_configuration,
      MappedProperties const &_final_properties, ConfigComparison _hint_status,
      double _hint_cost)
      : unmapped_child(_unmapped_child),
        mapping(_mapping),
        mapped_child(_mapped_child),
        mapped_configuration(_mapped_configuration),
        mapped_properties(_mapped_properties),
        symop_to_canon_scel(_symop_to_canon_scel),
        permutation_to_canon_scel(_permutation_to_canon_scel),
        transformation_matrix_to_canon_scel(
            _transformation_matrix_to_canon_scel),
        structure_in_canon_scel(_structure_in_canon_scel),
        configuration_in_canon_scel(_configuration_in_canon_scel),
        properties_in_canon_scel(_properties_in_canon_scel),
        symop_to_final(_symop_to_final),
        permutation_to_final(_permutation_to_final),
        transformation_matrix_to_final(_transformation_matrix_to_final),
        final_structure(_final_structure),
        final_configuration(_final_configuration),
        final_properties(_final_properties),
        hint_status(_hint_status),
        hint_cost(_hint_cost) {}

  /// \brief The original child structure, without any mapping or modifications
  xtal::SimpleStructure unmapped_child;

  /// \brief The StrucMapper solution mapping from `unmapped_child` to
  /// `mapped_child`
  xtal::MappingNode mapping;

  /// \brief The mapped child structure
  ///
  /// The mapped child structure is equivalent to `unmapped_child` but mapped by
  /// StrucMapper and transformed according to the results stored in
  /// `mapping` to the setting of the `parent` superstructure.
  ///
  /// Additionally, the following global properties are added:
  ///
  ///     mapped_child.properties["Ustrain"] =
  ///         StrainConverter("Ustrain").unroll_E(
  ///             mapping.lattice_node.stretch.inverse());
  ///
  ///     mapped_child.mol_info.properties["disp"] =
  ///         // displacements associated with molecules (mol_displacement in
  ///         //  MappingNode documentation)
  ///
  ///     mapped_child.properties["isometry"] =
  ///         // unrolled 9-element vector from mapping.lattice_node.isometry
  ///
  xtal::SimpleStructure mapped_child;

  /// \brief The configuration directly constructed from `mapped_child`
  ///
  /// Relationships that hold:
  ///
  ///     mapped_configuration.ideal_lattice().lat_column_mat() =
  ///         mapped_child.lat_column_mat
  ///
  ///     mapped_configuration.occ(i) = mapping.mol_labels[i].second
  ///
  /// Continuous DoF are read directly from `mapped_child.properties` and
  /// `mapped_child.mol_info.properties`.
  Configuration mapped_configuration;

  /// \brief Properties, as determined from `mapping`, and transformed to match
  /// the `mapped_child` and `mapped_configuration`
  ///
  /// Includes:
  /// - Global properties from `mapped_child.properties` in
  /// `mapped_properties.global` (reminder: stored as name:vector)
  /// - Local properties from `mapped_child.mol_info.properties` in
  /// `mapped_properties.site` (reminder: stored as name:matrix, with one
  /// column per site)
  ///
  /// Additionally:
  /// - If "disp" is not a DoF, then it will also store the molecular
  /// coordinates using:
  ///     mapped_properties.site["coordinate"] = mapped_child.mol_info.coords
  ///
  /// - If "*strain" (any flavor of strain) is not a DoF, then it will also
  /// store the lattice vectors:
  ///     mapped_properties.global["latvec"] = mapped_child.lat_column_mat;
  ///
  MappedProperties mapped_properties;

  /// \brief The transformation from the mapped_configuration to
  /// configuration_in_canon_scel
  ///
  /// The transformation from mapped configuration and structure to the
  /// configuration in the canonical supercell and corresponding structure.
  /// An element of the prim factor group. Lattice vectors are related
  /// according to:
  ///
  ///     configuration_in_canon_scel.ideal_lattice().lat_column_mat() =
  ///         symop_to_canon_scel.matrix() *
  ///         mapped_configuration.ideal_lattice().lat_column_mat() *
  ///         transformation_matrix_to_canon_scel
  ///
  /// DoFs, represented in the prim basis, are transformed according to
  /// representations stored for each prim factor group opreation. For
  /// continuous DoFs:
  ///
  ///     v_in_canon_scel = matrix * v_mapped,
  ///
  /// For global DoF, matrix is stored as:
  ///
  ///     symop_to_canon_scel.representation(
  ///         global_dof_info.symrep_ID()).MatrixXd()
  ///
  /// For local DoF, matrix is stored as:
  ///
  ///     symop_to_canon_scel.representation(
  ///          site_dof_info[b].symrep_ID()).MatrixXd()
  ///
  /// where `b` is the prim sublattice the site DoF begins on.
  ///
  /// Occupation values for anisotropic occupants are permuted according to
  /// the permutation representation, stored as:
  ///
  ///     occ[site_index] = occ_perm_b[occ[site_index]]
  ///     occ_perm_b = symop_to_canon_scel.get_permutation_rep(
  ///         m_occupation.symrep_IDs()[b])).permtuation(),
  ///
  /// where `occ_perm_b` is the permutation representation for occupants that
  /// begin on sublattice `b` of the prim.
  ///
  /// Properties, represented in the \link xtal::DoFSet "standard
  /// basis"\endlink, are transformed according to:
  ///
  ///     v_in_canon_scel = matrix * v_mapped,
  ///
  /// where
  ///
  ///     matrix = AnisoValTraits(type).symop_to_matrix(
  ///                  symop_to_canon_scel.matrix,
  ///                  symop_to_canon_scel.translation,
  ///                  symop_to_canon_scel.time_reversal),
  ///
  /// Site mapping is described by `permutation_to_canon_scel`.
  ///
  SymOp symop_to_canon_scel;

  /// \brief The permutation that describes site mapping that occurs for the
  /// transformation from mapped_configuration configuration_in_canon_scel
  ///
  /// The value `j = permutation_to_canon_scel[i]` indicates that values on
  /// site `j` in `mapped_configuration` are transformed according to the
  /// appropriate symmetry representation for `symop_to_canon_scel`, and then
  /// mapped to site `i` in `configuration_in_canon_scel`. The same
  /// permutation is used to map `mol_info` properties from `mapped_properties`
  /// to `properties_in_canon_scel`.
  Permutation permutation_to_canon_scel;

  /// \brief The transformation of the supercell lattice vectors going from
  /// mapped_configuration to configuration_in_canon_scel.
  ///
  /// The transformation of the supercell lattice vectors going from
  /// mapped_configuration to configuration_in_canon_scel satisfies the
  /// relation:
  ///
  ///     configuration_in_canon_scel.ideal_lattice().lat_column_mat() ==
  ///         symop_to_canon_scel.matrix() *
  ///         mapped_configuration.ideal_lattice().lat_column_mat() *
  ///         transformation_matrix_to_canon_scel
  ///
  Eigen::Matrix3l transformation_matrix_to_canon_scel;

  /// \brief The structure corresponding configuration_in_canon_scel
  xtal::SimpleStructure structure_in_canon_scel;

  /// \brief The configuration after mapping to the canonical supercell
  Configuration configuration_in_canon_scel;

  /// \brief Properties after mapping to the canonical supercell;
  MappedProperties properties_in_canon_scel;

  /// \brief The transformation from configuration_in_canon_scel to
  /// final_configuration
  ///
  /// The transformation from configuration DoF and properties in the
  /// canonical supercell to final configuration DoF and properties is made
  /// according to `symop_to_final` and `permutation_to_final`, using the
  /// same relationships described for `symop_to_canon_scel` and
  /// `permutation_to_canon_scel`.
  SymOp symop_to_final;

  /// \brief The permutation that describes site mapping from
  /// configuration_in_canon_scel to final_configuration.
  ///
  /// Note:
  /// - This is a different type of permutation than
  /// `permutation_to_canon_scel` because it acts on ConfigDoF instead of
  /// SimpleStructure coordinates and properties
  ///
  /// After being transformed, site properties are mapped using:
  ///
  ///     result.site[property_name].col(i) =
  ///         (transformed) input.site[property_name].col(permutation[i]).
  ///
  Permutation permutation_to_final;

  /// \brief The transformation of the supercell lattice vectors going from
  /// configuration_in_canon_scel to final_configuration.
  ///
  /// While the supercell of final_configuration is still the canonical
  /// supercell, to satisfy the following relationship a non-identity
  /// transformation matrix may be necessary:
  ///
  ///     final_configuration.ideal_lattice().lat_column_mat() ==
  ///         symop_to_final.matrix() *
  ///         configuration_in_canon_scel.ideal_lattice().lat_column_mat() *
  ///         transformation_matrix_to_final
  ///
  Eigen::Matrix3l transformation_matrix_to_final;

  /// \brief The structure corresponding to final_configuration
  ///
  /// The final_structure is symmetrically equivalent to mapped_child,
  /// but possibly with a different choice of supercell and permutation.
  xtal::SimpleStructure final_structure;

  /// \brief The final configuration
  ///
  /// The final_configuration is symmetrically equivalent to
  /// mapped_configuration, but possibly with a different choice of supercell
  /// and permutation.
  Configuration final_configuration;

  /// \brief Properties after mapping to the final configuration;
  MappedProperties final_properties;

  /// \brief The degree to which the hinted configuration matches
  /// final_configuration
  ConfigComparison hint_status;

  /// \brief The structure mapping score for the hinted configuration
  double hint_cost;

  /// \brief Compares ConfigurationMapping using MappingNode (this->mapping)
  bool operator<(ConfigurationMapping const &other) const {
    return this->mapping < other.mapping;
  }
};

/// Data structure holding results of the ConfigMapper algorithm
///
/// Consists of:
/// - std::set<ConfigurationMapping> maps: The set of structure-
///   to-configuration mapping solutions found by ConfigMapper. See
///   ConfigurationMapping for details.
/// - std::string fail_msg: Failure message if ConfigMapper could not map input
///   structure to the prim.
struct ConfigMapperResult {
  ConfigMapperResult() {}

  /// Returns true if 1 or more solutions was found
  bool success() const { return !maps.empty(); }

  /// The number of solutions with score equal to the best scoring solution
  Index n_optimal(double tol = TOL) const;

  /// The set of mapping solutions from input structure to configuration
  std::set<ConfigurationMapping> maps;

  /// Failure message if could not map input structure to prim
  std::string fail_msg;
};

/// \brief Implements a method for mapping structures to configurations
class ConfigMapper {
 public:
  /// ConfigMapper constructor
  ///
  /// \param _shared_prim The reference structure that input structures should
  ///   be mapped to
  /// \param _settings ConfigMapping::Settings, all collected parameters
  ///   controlling the configuration mapping
  ///
  ConfigMapper(std::shared_ptr<Structure const> const &_shared_prim,
               ConfigMapperSettings const &_settings);

  std::shared_ptr<Structure const> const &shared_prim() const {
    return m_shared_prim;
  }

  ConfigMapperSettings const &settings() const { return m_settings; }

  xtal::StrucMapper const &struc_mapper() const { return m_struc_mapper; }

  /// \brief Map structure to prim
  ///
  /// \param _struc Structure to map to prim
  /// \param hint_ptr[in]
  /// \parblock
  ///                Provides a suggestion for which Configuration _struc should
  ///                map onto. The hint is used to reduce search times, score
  ///                the mapping to a particular configuration, and also used
  ///                in combination with the ConfigMapping::Settings
  ///                parameter 'strict' to try to map onto a particular
  ///                orientation of a configuration.
  /// \endparblock
  /// \param hint_dofs DoFs to include when structure-mapping `_struc` to a
  ///   xtal::SimpleStructure constructed from `hint_ptr`.
  ///
  ConfigMapperResult import_structure(xtal::SimpleStructure const &_struc,
                                      Configuration const *hint_ptr = nullptr,
                                      std::vector<DoFKey> const &_hint_dofs = {
                                          "occ"}) const;

 private:
  /// Compare configurations
  ConfigComparison _make_hint_status(
      Configuration const *hint_ptr,
      Configuration const &final_configuration) const;

  /// Prim structure being mapped to (parent)
  std::shared_ptr<Structure const> m_shared_prim;

  /// Implements the structure mapping
  xtal::StrucMapper m_struc_mapper;

  /// See ConfigMapperSettings members for parameter descriptions
  ConfigMapperSettings m_settings;
};

class PrimStrucMapCalculator : public xtal::SimpleStrucMapCalculator {
 public:
  PrimStrucMapCalculator(xtal::BasicStructure const &_prim,
                         std::vector<xtal::SymOp> const &symgroup = {},
                         xtal::SimpleStructure::SpeciesMode _species_mode =
                             xtal::SimpleStructure::SpeciesMode::ATOM);

 private:
  /// \brief Make an exact copy of the calculator (including any initialized
  /// members)
  xtal::StrucMapCalculatorInterface *_clone() const override {
    return new PrimStrucMapCalculator(*this);
  }

  xtal::BasicStructure m_prim;
};

}  // namespace CASM

#endif
