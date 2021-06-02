#ifndef CASM_SimpleStrucMapCalculator
#define CASM_SimpleStrucMapCalculator

#include "casm/crystallography/Adapter.hh"
#include "casm/crystallography/StrucMapCalculatorInterface.hh"

namespace CASM {
namespace xtal {
class SimpleStructure;
struct SymOp;
typedef std::vector<SymOp> SymOpVector;

// In this file:
struct MappingNode;
class StrucMapCalculatorInterface;
class SimpleStrucMapCalculator;

// TODO:
// - Explain _factor_group default / non-default behavior
// - Does species_mode SpeciesMode::MOL work? how does it behave?

class SimpleStrucMapCalculator : public StrucMapCalculatorInterface {
 public:
  /// StrucMapCalculatorInterface constructor
  ///
  /// \param _parent Reference structure to be mapped to
  /// \param _factor_group Factor group of the parent structure.
  /// \param _species_mode Specifies whether to map to parent atoms or
  ///     molecules. Use `SimpleStructure::SpeciesMode::ATOM`.
  /// \param allowed_species Names of allowed species on each parent structure
  ///     site. If empty, `allowed_species[site_index]` is set to the 1-element
  ///     vector with the name of the current species on the parent structure
  ///     site. Use `allowed_molecule_names` to use the names of
  ///     `xtal::Molecule` allowed on `xtal::BasicStructure` sites.
  ///
  SimpleStrucMapCalculator(
      SimpleStructure _parent,
      SymOpVector const &_factor_group = {SymOp::identity()},
      SimpleStructure::SpeciesMode species_mode =
          SimpleStructure::SpeciesMode::ATOM,
      StrucMapping::AllowedSpecies allowed_species = {})
      : StrucMapCalculatorInterface(std::move(_parent), _factor_group,
                                    species_mode, std::move(allowed_species)) {}

  /// StrucMapCalculatorInterface constructor
  ///
  /// \param _parent Reference structure to be mapped to
  /// \param _factor_group Factor group of the parent structure.
  /// \param _species_mode Specifies whether to map to parent atoms or
  ///     molecules. Use `SimpleStructure::SpeciesMode::ATOM`.
  /// \param allowed_species Names of allowed species on each parent structure
  ///     site. If empty, `allowed_species[site_index]` is set to the 1-element
  ///     vector with the name of the current species on the parent structure
  ///     site. Use `allowed_molecule_names` to use the names of
  ///     `xtal::Molecule` allowed on `xtal::BasicStructure` sites.
  template <typename ExternSymOpVector>
  SimpleStrucMapCalculator(
      SimpleStructure _parent,
      ExternSymOpVector const &_factor_group = {SymOp::identity()},
      SimpleStructure::SpeciesMode species_mode =
          SimpleStructure::SpeciesMode::ATOM,
      StrucMapping::AllowedSpecies allowed_species = {})
      : SimpleStrucMapCalculator(
            _parent,
            adapter::Adapter<SymOpVector, ExternSymOpVector>()(_factor_group),
            species_mode, allowed_species) {}

  virtual ~SimpleStrucMapCalculator() {}

  /// Constructs a list of prospective mapping translations
  std::vector<Eigen::Vector3d> translations(
      MappingNode const &_node,
      SimpleStructure const &child_struc) const override;

  /// Creates a copy of the child structure and applies mapping
  virtual SimpleStructure resolve_setting(
      MappingNode const &_node,
      SimpleStructure const &_child_struc) const override;

  /// Sets MappingNode data based on lattice and atomic mapping results
  void finalize(MappingNode &_node, SimpleStructure const &child_struc,
                bool const &symmetrize_atomic_cost = false) const override;

  /// Populates the cost matrix for the atomic assignment problem
  bool populate_cost_mat(MappingNode &_node,
                         SimpleStructure const &child_struc) const override;

 private:
  /// \brief Make an exact copy of the calculator (including any initialized
  /// members)
  virtual StrucMapCalculatorInterface *_clone() const override {
    return new SimpleStrucMapCalculator(*this);
  }

  /// \brief Make an exact copy of the calculator (including any initialized
  /// members)
  virtual StrucMapCalculatorInterface *_quasi_clone(
      SimpleStructure _parent,
      SymOpVector const &_factor_group = {SymOp::identity()},
      SimpleStructure::SpeciesMode _species_mode =
          SimpleStructure::SpeciesMode::ATOM,
      StrucMapping::AllowedSpecies _allowed_species = {}) const override {
    return new SimpleStrucMapCalculator(std::move(_parent), _factor_group,
                                        _species_mode,
                                        std::move(_allowed_species));
  }
};
}  // namespace xtal
}  // namespace CASM
#endif
