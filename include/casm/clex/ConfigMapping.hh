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

/// \brief Struct with optional parameters for Config Mapping
/// Specifies default parameters for all values, in order to simplify
/// parsing from JSON
///
/// There are 5 choices for the lattice-mapping portion of the mapping method.
/// The choice of one of these methods determines which other parameters have
/// an effect. One and only one must be chosen. The methods are:
///
/// 1. fix_ideal = true: Force lattice mapping solutions of the form
///     L1 * T = Q * L2, where
///     L1 = the prim lattice
///     L2 = child.lat_column_mat
///      T = an integer transformation matrix, to be determined
///      Q = isometry matrix, constrained to be an elemenet of the prim
///          factor group, to be determined
///
/// 2. fix_lattice_mapping = true: Force lattice mapping solutions of the form
///     S1 = V * Q * L2, where
///     S1 = this->configuration_lattice
///     L2 = child.lat_column_mat
///   V, Q = stretch and isometry matrices, to be determined
///
/// 3. fix_lattice = true: Force lattice mapping solutions of the form
///     S1 * N = V * Q * L2, where
///     S1 = this->configuration_lattice
///     L2 = child.lat_column_mat
///      N = unimodular matrix, generates non-identical but equivalent
///          parent superlattices, to be determined
///   V, Q = stretch and isometry matrices, to be determined
///
/// 4. fix_lattice_volume_range = true: Force lattice mapping solutions of the
/// form
///     S1 * N = V * Q * L2, where
///     this->min_lattice_volume <=
///         S1.determinant() <= this->max_lattice_volume,
///     S1 = lattice equivalent to the configuration lattice, to be determined
///     L2 = child.lat_column_mat
///      N = unimodular matrix, generates non-identical but equivalent
///          parent superlattices, to be determined
///   V, Q = stretch and isometry matrices, to be determined
///
/// 5. auto_lattice_volume_range = true: The most general mapping approach, it
///   is equivalent `fix_lattice_volume_range` with a volume range that is
///   determined by constaints on the Va frac and change in volume between the
///   input structure and mapped structure.
///
/// If the mapping is succesful, one or more mapping solutions will be found
/// and the unmapped input structure will be mapped to the prim, creating a
/// `mapped_child` structure, `mapped_configuration`, and `mapped_properties`,
/// where the `mapped_configuration` and `mapped_properties` have the same
/// lattice and orientation as `mapped_child`.
///
/// As a post-processing step, the `mapped_configuration` is transformed to the
/// `final_configuration` which is (a) always the canonical supercells and (b)
/// either:
/// 1. finalize_strict = true: (the 'strict' method) the configuration which
/// preserves the orientation of the original unmapped structure as much as
/// possible.
/// 2. otherwise, the finalized configuration is the canonical equivalent
/// configuration in the canonical supercell
///
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

  int options() const {
    int opt = 0;
    if (robust) opt |= xtal::StrucMapper::robust;
    return opt;
  }

  // --- Lattice mapping method choices: ---

  /// Force lattice mapping solutions of the form `L1 * T = Q * L2`, where
  ///     L1 = the prim lattice
  ///     L2 = child.lat_column_mat
  ///      T = an integer transformation matrix, to be determined
  ///      Q = isometry matrix, constrained to be an elemenet of the prim
  ///          factor group, to be determined
  bool fix_ideal = false;

  /// Force lattice mapping solutions of the form `S1 = V * Q * L2`, where
  ///     S1 = this->configuration_lattice.value()
  ///     L2 = child.lat_column_mat
  ///   V, Q = stretch and isometry matrices, to be determined
  bool fix_lattice_mapping = false;

  /// Force lattice mapping solutions of the form `S1 * N = V * Q * L2`, where
  ///     S1 = this->configuration_lattice.value()
  ///     L2 = child.lat_column_mat
  ///      N = unimodular matrix, generates non-identical but equivalent
  ///          parent superlattices, to be determined
  ///   V, Q = stretch and isometry matrices, to be determined
  bool fix_lattice = false;

  /// Force lattice mapping solutions of the form `S1 * N = V * Q * L2`, where
  ///     this->min_lattice_volume <=
  ///         S1.determinant() <= this->max_lattice_volume,
  ///     S1 = lattice equivalent to the configuration lattice, to be determined
  ///     L2 = child.lat_column_mat
  ///      N = unimodular matrix, generates non-identical but equivalent
  ///          parent superlattices, to be determined
  ///   V, Q = stretch and isometry matrices, to be determined
  bool fix_lattice_volume_range = false;

  /// Equivalent to `fix_lattice_volume_range`, with a volume range that is
  /// determined by constaints on the Va fraction and the change in volume
  /// between the input structure and mapped structure.
  bool auto_lattice_volume_range = false;

  // --- Structure mapping parameters: ---

  /// If true, use the symmetry-breaking strain cost.
  bool use_symmetry_breaking_strain_cost = false;

  /// If true, use the symmetry-breaking displacement cost.
  bool use_symmetry_breaking_displacement_cost = false;

  /// Specifies the cost function in terms of lattice deformation
  /// cost and atomic deformation cost:
  ///     total_cost = lattice_weight*lattice_deformation_cost +
  ///                  (1-lattice_weight)*atomic_deformation_cost
  /// Default: 0.5
  double lattice_weight;

  /// Tolerance used to determine if two mappings have identical cost
  /// Default: TOL
  double cost_tol;

  /// True invokes additional checks which determine whether there are any other
  /// mappings that are distinct from the best mapping, but have the same cost
  /// Default: false
  bool robust = false;

  /// Specify the number, k, of k-best mappings to include in solution set
  /// Default: 1
  Index k_best;

  /// Lattice used by `fix_lattice` and `fix_lattice_mapping` methods.
  /// This is required to have a value of one of those methods is chosen. It
  /// must be a superlattice of the prim lattice. It does not have to be
  /// canonical.
  /// Default: empty
  std::optional<Lattice> configuration_lattice;

  /// List of superlattices of parent structure to consider when searching for
  /// mappings. If any provided, only these lattices will be considered, and
  /// only if they fall within the volume range being searched.
  ///
  /// Applies to the `fixed_lattice_volume_range` and
  /// `auto_lattice_volume_range` lattice mapping methods.
  /// Default: empty
  std::vector<Lattice> allowed_lattices;

  /// If provided, used to filter list of potential supercells of the parent
  /// structure to search over.
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

  /// The minimum fraction of vacant sites. Below this fraction a mapping will
  /// not be considered.
  ///
  /// Applies to the `auto_lattice_volume_range` lattice mapping method only.
  /// Default: 0.
  double min_va_frac;

  /// The maximum fraction of vacant sites. Above this fraction a mapping will
  /// not be considered
  ///
  /// Applies to the `auto_lattice_volume_range` lattice mapping method only.
  /// Default: 1.
  double max_va_frac;

  /// Constrains the search space by setting a limit on allowed volume change.
  ///
  /// Applies to the `auto_lattice_volume_range` lattice mapping method only.
  /// Default: 0.3
  double max_volume_change;

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

/// Data structure holding results of ConfigMapper algorithm
struct ConfigMapperResult {
  typedef ConfigComparison HintStatus;

  /// Data structure storing a structure -> configuration mapping
  ///
  /// See member comments for definitions
  struct ConfigurationMapping {
    ConfigurationMapping(xtal::SimpleStructure const &_unmapped_child,
                         xtal::MappingNode const &_mapping,
                         xtal::SimpleStructure const &_mapped_child,
                         Configuration const &_mapped_configuration,
                         MappedProperties const &_mapped_properties,
                         xtal::SymOp const &_symop_to_final,
                         Eigen::Matrix3l const &_transformation_matrix_to_final,
                         Configuration const &_final_configuration,
                         MappedProperties const &_final_properties,
                         HintStatus _hint_status, double _hint_cost)
        : unmapped_child(_unmapped_child),
          mapping(_mapping),
          mapped_child(_mapped_child),
          mapped_configuration(_mapped_configuration),
          mapped_properties(_mapped_properties),
          symop_to_final(_symop_to_final),
          transformation_matrix_to_final(_transformation_matrix_to_final),
          final_configuration(_final_configuration),
          final_properties(_final_properties),
          hint_status(_hint_status),
          hint_cost(_hint_cost) {}

    /// Original child structure, without any mapping or modifications
    xtal::SimpleStructure unmapped_child;

    /// StrucMapper output, describes mapping from `unmapped_child` structure to
    /// `mapped_child` structure via a lattice mapping and an atomic assignment
    /// mapping
    xtal::MappingNode mapping;

    /// Mapped child structure, equivalent to `unmapped_child`, but mapped by
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
    ///         // `MappingNode` documentation)
    ///
    ///     mapped_child.properties["isometry"] =
    ///         // unrolled 9-element vector from mapping.lattice_node.isometry
    ///
    xtal::SimpleStructure mapped_child;

    /// Configuration, as directly constructed from `mapped_child`
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

    /// Properties, as determined from mapping, and transformed to match the
    /// `mapped_child` and `mapped_configuration`
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

    /// The transformation from mapped configuration DoF and properties to
    /// final configuration DoF and properties is made according to:
    ///
    ///     v_final = matrix * v_mapped,
    ///     where matrix = AnisoValTraits(type).symop_to_matrix(
    ///                        symop_to_final.matrix,
    ///                        symop_to_final.translation,
    ///                        symop_to_final.time_reversal)
    ///
    /// The site DoF and property values are similarly transformed, and also
    /// moved from site with coordinate `r_mapped` to site with coordinate
    /// `r_final` (and then placed within the periodic boundaries of the
    /// supercell) according to:
    ///
    ///     r_final = symop_to_final.matrix * r_mapped +
    ///               symop_to_final.translation
    ///     TODO: update relationship to show coordinate->index conversions
    ///
    xtal::SymOp symop_to_final;

    /// The transformation from mapped to final supercell lattice vectors is
    /// made according to:
    ///
    ///     final_configuration.ideal_lattice().lat_column_mat() =
    ///         symop_to_final.matrix *
    ///         final_configuration.ideal_lattice().lat_column_mat() *
    ///         transformation_matrix_to_final
    ///
    Eigen::Matrix3l transformation_matrix_to_final;

    /// Final configuration, symmetrically equivalent to `mapped_configuration`,
    /// but possibly with a different choice of supercell and permutation.
    Configuration final_configuration;

    /// Properties, symmetrically equivalent to `mapped_properties`,
    /// transformed to match `final_configuration`
    MappedProperties final_properties;

    /// Degree to which the hinted configuration matches final_configuration
    HintStatus hint_status;

    /// Structure mapping score for hinted configuration matches
    /// final_configuration
    double hint_cost;

    /// Compare using MappingNode (this->mapping)
    bool operator<(ConfigurationMapping const &other) const {
      return this->mapping < other.mapping;
    }
  };

  ConfigMapperResult() {}

  bool success() const { return !maps.empty(); }

  Index n_optimal(double tol = TOL) const;

  /// The set of mapping solutions from input structure to configuration
  std::set<ConfigurationMapping> maps;

  /// Failure message if could not map input structure to prim
  std::string fail_msg;
};

/// A class for mapping an arbitrary crystal structure as a configuration of a
/// crystal template as described by a PrimClex.  ConfigMapper manages options
/// for the mapping algorithm and mapping cost function It also caches some
/// information about supercell lattices so that batch imports are more
/// efficient
///
/// \ingroup Configuration
class ConfigMapper {
 public:
  using HintStatus = ConfigMapperResult::HintStatus;

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
  HintStatus _make_hint_status(Configuration const *hint_ptr,
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
