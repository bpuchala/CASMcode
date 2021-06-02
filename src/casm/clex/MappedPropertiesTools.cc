#include "casm/clex/MappedPropertiesTools.hh"

#include "casm/basis_set/DoFTraits.hh"
#include "casm/clex/MappedProperties.hh"
#include "casm/crystallography/AnisoValTraits.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/symmetry/PermuteIterator.hh"

namespace CASM {
MappedProperties copy_apply(PermuteIterator const &op,
                            MappedProperties const &props) {
  SymOp symop = op.sym_op();

  MappedProperties result;
  for (auto it = props.global.begin(); it != props.global.end(); ++it) {
    AnisoValTraits traits(AnisoValTraits::name_suffix(it->first));
    Eigen::MatrixXd matrix = traits.symop_to_matrix(symop.matrix(), symop.tau(),
                                                    symop.time_reversal());
    result.global[it->first] = matrix * it->second;
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

namespace {

/// Return set of xtal::SimpleStructure properties that are used for DoF
///
/// \throws if DoF types are not found in structure properties
std::set<std::string> make_dof_managed_properties(
    xtal::SimpleStructure const &child, xtal::BasicStructure const &prim) {
  std::set<std::string> dof_managed_properties;

  for (auto const &dof_key : xtal::global_dof_types(prim)) {
    auto val = DoFType::traits(dof_key).find_values(child.properties);

    // TODO: what case would require multiple property names? currently
    // only one is ever returned
    dof_managed_properties.insert(val.second.begin(), val.second.end());
  }

  for (auto const &dof_key : xtal::continuous_local_dof_types(prim)) {
    auto val = DoFType::traits(dof_key).find_values(child.mol_info.properties);

    // TODO: what case would require multiple property names? currently
    // only one is ever returned
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
/// - Property names in "mapped_structure" must follow CASM property naming
/// conventions as documented for AnisoValTraits.
/// - If "disp" is a property of mapped_structure (i.e. it is in
/// `simple_structure.mol_info.properties`), rather than a DoF, then store
/// store the mol_info coordinates using:
///     result.site["coordinate"] = mapped_structure.mol_info.coords;
/// - If "*strain" is a global property (i.e. "*strain" is in
/// `mapped_structure.properties`), rather than a DoF, we will also store the
/// lattice vectors using:
///     result.global["latvec"] = mapped_structure.lat_column_mat;
/// - This only sets the `site` and `global` members of MappedProperties
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
      // If "*strain" is a property, rather than a DoF, we will also store the
      // lattice
      if (prop.first.find("strain") != std::string::npos) {
        result.global["latvec"] =
            as_column_major_vector(mapped_structure.lat_column_mat);
      }
    }
  }

  for (auto const &prop : mapped_structure.mol_info.properties) {
    if (!dof_managed_properties.count(prop.first)) {
      result.site[prop.first] = prop.second;
      // If "disp" is a property, rather than a DoF, we will also store the
      // coordinates
      if (prop.first == "disp") {
        result.site["coordinate"] = mapped_structure.mol_info.coords;
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
