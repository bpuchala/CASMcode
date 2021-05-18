#include "casm/clex/MappedPropertiesTools.hh"

#include "casm/clex/MappedProperties.hh"
#include "casm/crystallography/AnisoValTraits.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/symmetry/PermuteIterator.hh"

namespace CASM {
MappedProperties copy_apply(PermuteIterator const &op,
                            MappedProperties const &props) {
  SymOp symop = op.sym_op();

  MappedProperties result;
  for (auto it = props.global.begin(); it != props.global.end(); ++it) {
    AnisoValTraits traits(AnisoValTraits::name_suffix(it->first));
    result.global[it->first] =
        traits.symop_to_matrix(symop.matrix(), symop.tau(),
                               symop.time_reversal()) *
        it->second;
  }

  Permutation tperm(op.combined_permute());
  for (auto it = props.site.begin(); it != props.site.end(); ++it) {
    AnisoValTraits traits(AnisoValTraits::name_suffix(it->first));
    Eigen::MatrixXd new_matrix =
        traits.symop_to_matrix(symop.matrix(), symop.tau(),
                               symop.time_reversal()) *
        it->second;
    auto it2 = result.site
                   .emplace(std::make_pair(
                       it->first,
                       Eigen::MatrixXd(new_matrix.rows(), new_matrix.cols())))
                   .first;
    for (Index i = 0; i < (it->second).cols(); i++) {
      (it2->second).col(i) = new_matrix.col(tperm[i]);
    }
  }
  return result;
}

/// Construct MappedProperties from a mapped SimpleStructure
///
/// \param simple_structure A SimpleStructure solution of the mapping algorithm.
/// This represents a way the configuration mapping algorithm found to map the
/// SimpleStructure to a configuratio of the prim.
/// \param dof_managed_properties A list of the SimpleStructure properties, both
/// site and global, that correspond to prim degrees of freedom and thus should
/// not be included in MappedProperties
///
/// Notes:
/// - Property names in "simple_structure" and "dof_managed_properties" must
/// follow CASM property naming conventions as documented for AnisoValTraits.
/// - If "disp" is a property of simple_structure (i.e. it is in
/// `simple_structure.mol_info.properties`), rather than a DoF, then store
/// store the mol_info coordinates using:
///     result.site["coordinate"] = simple_structure.mol_info.coords;
/// - If "*strain" is a global property (i.e. "*strain" is in
/// `simple_structure.properties`), rather than a DoF, we will also store the
/// lattice vectors using:
///     result.global["latvec"] = simple_structure.lat_column_mat;
/// - This only sets the `site` and `global` members of MappedProperties
///
MappedProperties make_mapped_properties(
    xtal::SimpleStructure const &simple_structure,
    std::set<std::string> const &dof_managed_properties) {
  MappedProperties result;

  for (auto const &prop : simple_structure.properties) {
    if (!dof_managed_properties.count(prop.first)) {
      result.global[prop.first] = prop.second;
      // If "*strain" is a property, rather than a DoF, we will also store the
      // lattice
      if (prop.first.find("strain") != std::string::npos) {
        result.global["latvec"] = simple_structure.lat_column_mat;
      }
    }
  }

  for (auto const &prop : simple_structure.mol_info.properties) {
    if (!dof_managed_properties.count(prop.first)) {
      result.site[prop.first] = prop.second;
      // If "disp" is a property, rather than a DoF, we will also store the
      // coordinates
      if (prop.first == "disp") {
        result.site["coordinate"] = simple_structure.mol_info.coords;
      }
    }
  }

  return result;
}

/// Check if all required properties are included in MappedProperties
bool has_all_required_properties(
    MappedProperties const &properties,
    std::vector<std::string> const &required_properties) {
  for (std::string const &propname : required_properties) {
    if (!properties.global.count(propname) &&
        !properties.site.count(propname)) {
      return false;
    }
  }
  return true;
}

}  // namespace CASM
