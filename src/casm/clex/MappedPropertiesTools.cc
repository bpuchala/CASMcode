#include "casm/clex/MappedPropertiesTools.hh"

#include "casm/basis_set/DoFTraits.hh"
#include "casm/clex/MappedProperties.hh"
#include "casm/crystallography/AnisoValTraits.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/symmetry/PermuteIterator.hh"

namespace CASM {

namespace sym {

/// \brief Apply symmetry to MappedProperties, using permutation for mapping
/// site properties
///
/// \param symop Symmetry operation used to transform site and global properties
/// \param permutation Permutation used for mapping site properties.
///
/// Transformed site properties are mapped using:
///     result.site[property_name].col(i) =
///         (transformed) input.site[property_name].col(permutation[i]).
///
MappedProperties copy_apply(SymOp const &symop, Permutation const &permutation,
                            MappedProperties const &props) {
  MappedProperties result;
  for (auto it = props.global.begin(); it != props.global.end(); ++it) {
    AnisoValTraits traits(AnisoValTraits::name_suffix(it->first));
    Eigen::MatrixXd matrix = traits.symop_to_matrix(symop.matrix(), symop.tau(),
                                                    symop.time_reversal());
    result.global[it->first] = matrix * it->second;
  }

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
      (it2->second).col(i) = new_matrix.col(permutation[i]);
    }
  }
  return result;
}

/// \brief Apply symmetry to MappedProperties, using combined permutation for
/// mapping site properties
MappedProperties copy_apply(PermuteIterator const &op,
                            MappedProperties const &props) {
  return copy_apply(op.sym_op(), op.combined_permute(), props);
}

}  // namespace sym

namespace {

/// Return set of xtal::SimpleStructure properties that are used for DoF
///
/// \throws if DoF types are not found in structure properties
std::set<std::string> make_dof_managed_properties(
    xtal::SimpleStructure const &child, xtal::BasicStructure const &prim) {
  std::set<std::string> dof_managed_properties;

  for (auto const &dof_key : xtal::global_dof_types(prim)) {
    auto val = DoFType::traits(dof_key).find_values(child.properties);
    dof_managed_properties.insert(val.second.begin(), val.second.end());
  }

  for (auto const &dof_key : xtal::continuous_local_dof_types(prim)) {
    auto val = DoFType::traits(dof_key).find_values(child.mol_info.properties);
    dof_managed_properties.insert(val.second.begin(), val.second.end());
  }

  return dof_managed_properties;
}

Eigen::VectorXd as_column_major_vector(Eigen::MatrixXd const &matrix) {
  Eigen::MatrixXd M = matrix;
  return Eigen::Map<Eigen::VectorXd>(M.data(), M.size());
}

Eigen::MatrixXd as_column_major_matrix(Eigen::VectorXd vector, int rows,
                                       int cols) {
  return Eigen::Map<Eigen::MatrixXd>(vector.data(), rows, cols);
}

}  // namespace

/// Construct MappedProperties from a mapped SimpleStructure
///
/// \param simple_structure A SimpleStructure solution of the mapping algorithm.
/// This represents a way the configuration mapping algorithm found to map the
/// SimpleStructure to a configuratio of the prim.
/// \param shared_prim The prim the structure was mapped to. Determines which
/// properties are DoF and which are MappedProperties
///
/// Notes:
/// - DoFTraits, specifically the DoFType::traits(dof_key).find_values methods,
/// are used to determine which global and mol_info "mapped_structure"
/// properties are represented by DoF, and the rest are copied into the
/// `global` and `site` members of the resulting MappedProperties object.
/// - Property names in "mapped_structure" must follow CASM property naming
/// conventions as documented for AnisoValTraits.
///
MappedProperties make_mapped_properties(
    xtal::SimpleStructure const &mapped_structure,
    xtal::BasicStructure const &prim) {
  MappedProperties result;
  auto dof_managed_properties =
      make_dof_managed_properties(mapped_structure, prim);

  for (auto const &prop : mapped_structure.properties) {
    if (!dof_managed_properties.count(prop.first)) {
      result.global[prop.first] = prop.second;
    }
  }

  for (auto const &prop : mapped_structure.mol_info.properties) {
    if (!dof_managed_properties.count(prop.first)) {
      result.site[prop.first] = prop.second;
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
