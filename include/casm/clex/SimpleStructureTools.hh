#ifndef CLEX_SIMPLESTRUCTURETOOLS_HH
#define CLEX_SIMPLESTRUCTURETOOLS_HH

#include <set>
#include <vector>

#include "casm/crystallography/DoFDecl.hh"
#include "casm/global/definitions.hh"

namespace CASM {

namespace xtal {
class SimpleStructure;
}

class ConfigDoF;
class Configuration;
struct MappedProperties;
class Supercell;

/// Construct xtal::SimpleStructure from ConfigDoF (DoF only)
xtal::SimpleStructure make_simple_structure(
    Supercell const &supercell, ConfigDoF const &configdof,
    std::vector<DoFKey> const &which_dofs = {});

/// Construct xtal::SimpleStructure from Configuration (DoF only)
xtal::SimpleStructure make_simple_structure(
    Configuration const &configuration,
    std::vector<DoFKey> const &which_dofs = {});

/// Construct xtal::SimpleStructure from ConfigDoF and MappedProperties
xtal::SimpleStructure make_simple_structure(
    Supercell const &supercell, ConfigDoF const &configdof,
    MappedProperties const &properties,
    std::vector<DoFKey> const &which_dofs = {});

/// Construct xtal::SimpleStructure from Configuration and MappedProperties
xtal::SimpleStructure make_simple_structure(
    Configuration const &configuration, MappedProperties const &properties,
    std::vector<DoFKey> const &which_dofs = {});

/// \brief Determine which sites of a Configuration can host each atom of a
/// SimpleStructure result[i] is set of site indices in @param _config that can
/// host atom 'i' of @param sstruc
std::vector<std::set<Index> > atom_site_compatibility(
    xtal::SimpleStructure const &sstruc, Configuration const &_config);

/// \brief Determine which sites of a Configuration can host each molecule of a
/// SimpleStructure result[i] is set of site indices in @param _config that can
/// host molecule 'i' of @param sstruc
std::vector<std::set<Index> > mol_site_compatibility(
    xtal::SimpleStructure const &sstruc, Configuration const &_config);

}  // namespace CASM

#endif
