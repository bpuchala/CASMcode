#include "casm/clex/SimpleStructureTools.hh"

#include <string>

#include "casm/basis_set/DoFTraits.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/Supercell.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {

namespace clex_SimpleStructureTools_impl {

/// \brief Imposes DoF values from ConfigDoF _config onto _sstruc, using using
/// any necessary information contained in _reference
void _apply_dofs(xtal::SimpleStructure &_sstruc, ConfigDoF const &_config,
                 xtal::BasicStructure const &_reference,
                 std::vector<DoFKey> which_dofs);

/// \brief Imposes DoF and property values from `dofs_and_properties` onto
/// simplestructure.
void _apply_dofs_and_properties(
    xtal::SimpleStructure &simplestructure, ConfigDoF const &configdof,
    xtal::BasicStructure const &prim,
    std::map<std::string, Eigen::MatrixXd> dofs_and_properties);

}  // namespace clex_SimpleStructureTools_impl

/// Construct xtal::SimpleStructure from ConfigDoF
///
/// \param supercell Supercell of the configuration
/// \param configdof ConfigDoF to make the structure from
/// \param which_dofs Names of DoFs to include in the resulting structure. If
///   empty, all DoFs are included. To exclude all DoFs from the result, use
///  `{"none"}`.
xtal::SimpleStructure make_simple_structure(
    Supercell const &supercell, ConfigDoF const &configdof,
    std::vector<DoFKey> const &which_dofs) {
  return make_simple_structure(supercell, configdof, MappedProperties(),
                               which_dofs);
}

/// Construct xtal::SimpleStructure from Configuration (DoF only)
///
/// Equivalent to:
/// \code
/// make_simple_structure(configuration.supercell(),
///     configuration.configdof(), which_dofs);
/// \endcode
xtal::SimpleStructure make_simple_structure(
    Configuration const &configuration, std::vector<DoFKey> const &which_dofs) {
  return make_simple_structure(configuration.supercell(),
                               configuration.configdof(), which_dofs);
}

/// Construct xtal::SimpleStructure from ConfigDoF and MappedProperties
///
/// \param supercell Supercell of the configuration
/// \param configdof ConfigDoF to make the structure from
/// \param properties MappedProperties to include in the resulting structure
/// \param which_dofs Names of DoFs to include in the resulting structure. If
///   empty, all DoFs are included. To exclude all DoFs from the result, use
///  `{"none"}`.
///
/// Note:
/// - First, a reference structure is construced from the ideal supercell
/// lattice, site coordinates, and the occupants listed at each site by
/// configdof.
/// - Then, properties and DoF are:
///
///   1) copied to structure properties (global) and mol_info.properties
///        (local)
///   2) and applied, if DoFTraits exist, using the `apply_standard_values`
///        method of DoFTraits according to the `must_apply_before` /
///        `must_apply_after` order specified by AnisoValTraits.
xtal::SimpleStructure make_simple_structure(
    Supercell const &supercell, ConfigDoF const &configdof,
    MappedProperties const &properties, std::vector<DoFKey> const &which_dofs) {
  auto const &prim = supercell.prim();

  // create reference structure
  xtal::SimpleStructure result;
  result.lat_column_mat = supercell.lattice().lat_column_mat();
  result.mol_info.resize(configdof.size());
  for (Index b = 0, l = 0; b < configdof.n_sublat(); ++b) {
    for (Index v = 0; v < configdof.n_vol(); ++v, ++l) {
      result.mol_info.cart_coord(l) = supercell.coord(l).const_cart();
    }
  }

  // copy molecule attributes, including names
  for (Index b = 0, l = 0; b < configdof.n_sublat(); ++b) {
    for (Index v = 0; v < configdof.n_vol(); ++v, ++l) {
      Molecule const &mol = prim.basis()[b].occupant_dof()[configdof.occ(l)];

      // Fill up the molecule's SpeciesAttributes
      for (auto const &attr : mol.attributes()) {
        // Has this attribute been encountered yet??
        auto it = result.mol_info.properties.find(attr.first);
        /// If not, initialize it
        if (it == result.mol_info.properties.end()) {
          // Iterator now points to initialized matrix
          it = result.mol_info.properties
                   .emplace(attr.first,
                            Eigen::MatrixXd::Zero(attr.second.traits().dim(),
                                                  configdof.size()))
                   .first;
        }
        it->second.col(l) = attr.second.value();
      }

      // Record name
      result.mol_info.names[l] = mol.name();
    }
  }

  std::map<std::string, Eigen::MatrixXd> dofs_and_properties;
  std::set<std::string> which_dof_set{which_dofs.begin(), which_dofs.end()};

  // add local dofs
  for (std::string const &dof : xtal::continuous_local_dof_types(prim)) {
    if (which_dof_set.empty() || which_dof_set.count(dof)) {
      dofs_and_properties.emplace(dof,
                                  configdof.local_dof(dof).standard_values());
    }
  }

  // add global dofs
  for (std::string const &dof : xtal::global_dof_types(prim)) {
    if (which_dof_set.empty() || which_dof_set.count(dof)) {
      dofs_and_properties.emplace(dof,
                                  configdof.global_dof(dof).standard_values());
    }
  }

  // add local properties
  for (auto const &site_property : properties.site) {
    auto insert_result = dofs_and_properties.insert(site_property);
    if (!insert_result.second) {
      std::stringstream msg;
      msg << "Error in make_simple_structure: site property '"
          << site_property.first << "' conflicts with a DoF.";
      throw std::runtime_error(msg.str());
    }
  }

  // add global properties
  for (auto const &global_property : properties.global) {
    auto insert_result = dofs_and_properties.insert(global_property);
    if (!insert_result.second) {
      std::stringstream msg;
      msg << "Error in make_simple_structure: global property '"
          << global_property.first << "' conflicts with a DoF.";
      throw std::runtime_error(msg.str());
    }
  }

  clex_SimpleStructureTools_impl::_apply_dofs_and_properties(
      result, configdof, prim, dofs_and_properties);
  return result;
}

/// Construct xtal::SimpleStructure from Configuration and MappedProperties
///
/// Equivalent to:
/// \code
/// make_simple_structure(configuration.supercell(),
///     configuration.configdof(), properties, which_dofs);
/// \endcode

xtal::SimpleStructure make_simple_structure(
    Configuration const &configuration, MappedProperties const &properties,
    std::vector<DoFKey> const &which_dofs) {
  return make_simple_structure(configuration.supercell(),
                               configuration.configdof(), properties,
                               which_dofs);
}

//***************************************************************************

std::vector<std::set<Index> > mol_site_compatibility(
    xtal::SimpleStructure const &sstruc, Configuration const &_config) {
  std::vector<std::set<Index> > result;
  result.reserve(sstruc.mol_info.names.size());
  for (std::string const &sp : sstruc.mol_info.names) {
    result.push_back({});
    for (Index l = 0; l < _config.size(); ++l) {
      if (_config.mol(l).name() == sp) {
        result.back().insert(l);
      }
    }
  }
  return result;
}

std::vector<std::set<Index> > atom_site_compatibility(
    xtal::SimpleStructure const &sstruc, Configuration const &_config) {
  std::vector<std::set<Index> > result;
  result.reserve(sstruc.atom_info.names.size());
  for (std::string const &sp : sstruc.atom_info.names) {
    result.push_back({});
    for (Index l = 0; l < _config.size(); ++l) {
      if (_config.mol(l).contains(sp)) {
        result.back().insert(l);
      }
    }
  }
  return result;
}

namespace clex_SimpleStructureTools_impl {

/// Class that helps manage which order DoF are applied to a SimpleStructure
/// when building it from a Configuration. Each TransformDirective applies one
/// DoF type or "atomizes" molecules (populating SimpleStructure::atom_info from
/// SimpleStructure::mol_info). It is meant to be stored in a
/// std::set<TransformDirective> and has a comparison operator and uses
/// information from AnisoValTraits to appropriately order TransformDirective in
/// the set. Then they can be applied sequentially (using `transform`) to build
/// the SimpleStructure.
class TransformDirective {
 public:
  /// \brief consturct from transformation or DoF type name, ConfigDoF, and prim
  TransformDirective(std::string const &_name,
                     Eigen::MatrixXd const &_standard_values,
                     ConfigDoF const &_configdof,
                     xtal::BasicStructure const &_prim);

  /// \brief Name of DoFType or transformation
  std::string const &name() const { return m_name; }

  /// \brief Compare with _other TransformDirective. Returns true if this
  /// TransformDirective has precedence
  bool operator<(TransformDirective const &_other) const;

  /// \brief Applies transformation to _struc using information contained in
  /// _config
  void transform(ConfigDoF const &_config,
                 xtal::BasicStructure const &_reference,
                 xtal::SimpleStructure &_struc) const;

  /// \brief Applies transformation to _struc using information contained in
  /// _config.
  ///
  /// This version enables applying properties and dofs on the same footing,
  /// using `m_standard_values`.
  void transform(xtal::SimpleStructure &_struc) const;

 private:
  /// \brief Build m_before object by recursively traversing DoF dependencies
  void _accumulate_before(std::set<std::string> const &_queue,
                          std::set<std::string> &_result) const;

  /// \brief Build m_after object by recursively traversing DoF dependencies
  void _accumulate_after(std::set<std::string> const &_queue,
                         std::set<std::string> &_result) const;

  std::string m_name;
  Eigen::MatrixXd m_standard_values;
  ConfigDoF const &m_configdof;
  xtal::BasicStructure const &m_prim;
  std::set<std::string> m_before;
  std::set<std::string> m_after;

  DoFType::Traits const *m_traits_ptr;
};

void _apply_dofs(xtal::SimpleStructure &_sstruc, ConfigDoF const &_config,
                 xtal::BasicStructure const &_reference,
                 std::vector<DoFKey> which_dofs) {
  std::set<TransformDirective> tformers;
  tformers.emplace("atomize", Eigen::MatrixXd(), _config, _reference);
  if (which_dofs.empty()) {
    for (std::string const &dof : continuous_local_dof_types(_reference))
      which_dofs.push_back(dof);
    for (std::string const &dof : global_dof_types(_reference))
      which_dofs.push_back(dof);
  }

  // this version gets dof values from _config in the transform call
  for (DoFKey const &dof : which_dofs) {
    if (dof != "none" && dof != "occ")
      tformers.emplace(dof, Eigen::MatrixXd(), _config, _reference);
  }

  for (TransformDirective const &tformer : tformers) {
    tformer.transform(_config, _reference, _sstruc);
  }
}

/// \brief Imposes DoF and property values from `dofs_and_properties` onto
/// simplestructure.
///
/// \param simplestructure Structure being transformed
/// \param configdof Configuration DoF, used to set "atomize" and set
///     occupation, and provided to DoFTraits implementations along with DoF
///     and property values.
/// \param prim Prim structure, used to set "atomize" and set
///     occupation, and provided to DoFTraits implementations along with DoF
///     and property values.
/// \param dofs_and_properties {name, standard value}, this is what is applied
///
void _apply_dofs_and_properties(
    xtal::SimpleStructure &simplestructure, ConfigDoF const &configdof,
    xtal::BasicStructure const &prim,
    std::map<std::string, Eigen::MatrixXd> dofs_and_properties) {
  // Create set of transformers, initialize with "atomize"
  std::set<TransformDirective> transformers;
  transformers.emplace("atomize", Eigen::MatrixXd(), configdof, prim);

  // Create all other transformers to apply DoFs and properties
  // this version gets dof values from dofs_and_properties
  for (auto const &dof_or_property : dofs_and_properties) {
    std::string dof = dof_or_property.first;
    if (dof != "none" && dof != "occ") {
      Eigen::MatrixXd const &standard_values = dof_or_property.second;
      transformers.emplace(dof, standard_values, configdof, prim);
    }
  }

  // Apply DoFs and properties
  for (TransformDirective const &transformer : transformers) {
    transformer.transform(simplestructure);
  }
}

TransformDirective::TransformDirective(std::string const &_name,
                                       Eigen::MatrixXd const &_standard_values,
                                       ConfigDoF const &_configdof,
                                       xtal::BasicStructure const &_prim)
    : m_name(_name),
      m_standard_values(_standard_values),
      m_configdof(_configdof),
      m_prim(_prim),
      m_traits_ptr(nullptr) {
  if (name() != "atomize") {
    if (DoFType::traits_dict().contains(name())) {
      m_traits_ptr = &DoFType::traits(name());
    }
    _accumulate_before({_name}, m_before);
    _accumulate_after({_name}, m_after);
    if (m_after.count("atomize") == 0) m_before.insert("atomize");
  }
}

bool TransformDirective::operator<(TransformDirective const &_other) const {
  if (m_before.count(_other.name()) || _other.m_after.count(name())) {
    return false;
  }
  if (m_after.count(_other.name()) || _other.m_before.count(name())) {
    return true;
  }
  return name() < _other.name();
}

void TransformDirective::_accumulate_before(
    std::set<std::string> const &_queue, std::set<std::string> &_result) const {
  for (std::string const &el : _queue) {
    if (el != name()) _result.insert(el);
    if (el != "atomize")
      _accumulate_before(AnisoValTraits(el).must_apply_before(), _result);
  }
}

void TransformDirective::_accumulate_after(
    std::set<std::string> const &_queue, std::set<std::string> &_result) const {
  for (std::string const &el : _queue) {
    if (el != name()) _result.insert(el);
    if (el != "atomize")
      _accumulate_after(AnisoValTraits(el).must_apply_after(), _result);
  }
}

// This is the original version, uses DoF values from _dof
void TransformDirective::transform(ConfigDoF const &_dof,
                                   xtal::BasicStructure const &_reference,
                                   xtal::SimpleStructure &_struc) const {
  if (name() == "atomize") {
    _atomize(_struc, _dof.occupation(), _reference);
  } else if (m_traits_ptr) {
    if (m_traits_ptr->val_traits().global())
      _struc.properties[m_traits_ptr->name()] =
          _dof.global_dof(m_traits_ptr->name()).standard_values();
    else {
      _struc.mol_info.properties[m_traits_ptr->name()] =
          _dof.local_dof(m_traits_ptr->name()).standard_values();
    }
    m_traits_ptr->apply_dof(_dof, _reference, _struc);
  }
}

// This version enables applying properties and dofs on the same footing
void TransformDirective::transform(xtal::SimpleStructure &_struc) const {
  if (name() == "atomize") {
    _atomize(_struc, m_configdof.occupation(), m_prim);
  } else {
    if (AnisoValTraits(name()).global())
      _struc.properties[name()] = m_standard_values;
    else {
      _struc.mol_info.properties[name()] = m_standard_values;
    }
    if (m_traits_ptr) {
      m_traits_ptr->apply_standard_values(m_standard_values, _struc);
    }
  }
}

}  // namespace clex_SimpleStructureTools_impl
}  // namespace CASM
