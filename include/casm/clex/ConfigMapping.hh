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

namespace ConfigMapping {

/// \brief Struct with optional parameters for Config Mapping
/// Specifies default parameters for all values, in order to simplify
/// parsing from JSON
///
/// There are essentially 5 modes for structure mapping:
///
/// Method 1: hint && 'ideal' option (exactly known lattice mapping)
/// Method 2: hint && 'fix_lattice' option (exactly known lattice)
/// Method 3: hint && 'fix_volume' option (exactly known volume)
/// Method 4: no hint && 'ideal' option (ideal integer supercell of prim)
/// Method 5: otherwise, most general mapping

struct Settings {
  Settings(double _lattice_weight = 0.5, bool _ideal = false,
           bool _strict = false, bool _robust = false,
           // bool _primitive_only = false,
           bool _fix_volume = false, bool _fix_lattice = false,
           Index _k_best = 1, std::vector<Lattice> _forced_lattices = {},
           xtal::StrucMapping::LatticeFilterFunction _filter =
               xtal::StrucMapping::LatticeFilterFunction(),
           double _cost_tol = CASM::TOL, double _min_va_frac = 0.,
           double _max_va_frac = 0.5, double _max_vol_change = 0.3)
      : lattice_weight(_lattice_weight),
        ideal(_ideal),
        strict(_strict),
        robust(_robust),
        // primitive_only(_primitive_only),
        fix_volume(_fix_volume),
        fix_lattice(_fix_lattice),
        k_best(_k_best),
        forced_lattices(_forced_lattices),
        filter(_filter),
        cost_tol(_cost_tol),
        min_va_frac(_min_va_frac),
        max_va_frac(_max_va_frac),
        max_vol_change(_max_vol_change) {}

  int options() const {
    int opt = 0;
    if (robust) opt |= xtal::StrucMapper::robust;
    return opt;
  }

  void set_default() { *this = Settings(); }

  /// lattice_weight specifies the cost function in terms of lattice deformation
  /// cost and atomic deformation cost (i.e., atomic displacement) cost =
  /// lattice_weight*lattice_cost + (1l-lattice_weight)*atomic_displacement_cost
  double lattice_weight;

  /// True if child structure's lattice should be assumed to be ideal integer
  /// supercell of parent structure
  bool ideal;

  /// True invokes post-processing step to find a symmetry operation of the
  /// parent structure that preserves setting of child structure as much as
  /// possible after mapping
  bool strict;

  /// True invokes additional checks which determine whether there are any other
  /// mappings that are distinct from the best mapping, but have the same cost
  bool robust;

  /// If true, search for potential mappings will be constrained to the
  /// supercell volume of the hint config
  bool fix_volume;

  /// If true, search for potential mappings will be constrained to the exact
  /// supercell of the hint config
  bool fix_lattice;

  /// Specify the number, k, of k-best mappings to include in solution set
  /// (default is 1)
  Index k_best;

  /// List of superlattices of parent structure to consider when searching for
  /// mappings
  std::vector<Lattice> forced_lattices;

  /// If provided, used to filter list of potential supercells of the parent
  /// structure to search over
  xtal::StrucMapping::LatticeFilterFunction filter;

  /// Tolerance used to determine if two mappings have identical cost
  double cost_tol;

  /// minimum fraction of vacant sites, below this fraction a mapping will not
  /// be considered
  double min_va_frac;

  /// maximum fraction of vacant sites, above this fraction a mapping will not
  /// be considered
  double max_va_frac;

  /// constrains the search space by assuming a limit on allowed volume change
  /// only taken into account when non-interstitial vacancies are allowed in
  /// parent structure
  double max_vol_change;
};
}  // namespace ConfigMapping

/// Data structure holding results of ConfigMapper algorithm
struct ConfigMapperResult {
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
  enum class HintStatus {
    None,
    Derivative,
    Equivalent,
    Identical,
    NewOcc,
    NewScel
  };

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
  /// \param _tol tolerance for mapping comparisons (typically,
  ///   _shared_prim->crystallography_tol())
  ///
  ConfigMapper(std::shared_ptr<Structure const> const &_shared_prim,
               ConfigMapping::Settings const &_settings, double _tol);

  std::shared_ptr<Structure const> const &shared_prim() const {
    return m_shared_prim;
  }

  ConfigMapping::Settings const &settings() const { return m_settings; }

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

  /// See ConfigMapping::Settings members for parameter descriptions
  ConfigMapping::Settings m_settings;
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
