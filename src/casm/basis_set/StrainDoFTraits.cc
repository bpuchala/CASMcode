#include "casm/basis_set/StrainDoFTraits.hh"

#include "casm/basis_set/Adapter.hh"
#include "casm/basis_set/BasisSet.hh"
#include "casm/basis_set/DoFSet.hh"
#include "casm/basis_set/FunctionVisitor.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/clex/ClexBasis.hh"
#include "casm/clex/ConfigDoF.hh"
#include "casm/clex/NeighborList.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/clusterography/IntegralCluster.hh"
#include "casm/crystallography/Adapter.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/Site.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/strain/StrainConverter.hh"

namespace CASM {
namespace DoF_impl {

/// \brief Retrieve the standard values for a DoF from dictionary of properties
/// from a SimpleStructure or MappedProperties object
///  Returns matrix with standard values, and names of properties that were used
///  to construct the matrix
std::pair<Eigen::MatrixXd, std::set<std::string> > StrainDoFTraits::find_values(
    std::map<std::string, Eigen::MatrixXd> const &values) const {
  std::pair<Eigen::MatrixXd, std::set<std::string> > result;
  for (auto const &val : values) {
    std::string valname = AnisoValTraits::name_suffix(val.first);
    auto pos = valname.find("strain");
    if (pos != std::string::npos) {
      StrainConverter c(valname.substr(0, pos));
      Eigen::Matrix3d F = c.unrolled_strain_metric_to_F(val.second);
      c.set_mode(m_metric);

      result.first = c.unrolled_strain_metric(F);
      result.second.insert(val.first);
      return result;
    }
  }
  throw std::runtime_error(
      "Could not identify DoF values for DoF '" + (this->name()) +
      "' from provided list of tabulated structure properties.");
  return result;
}

/// \brief Construct the site basis (if DOF_MODE is LOCAL) for a DoF, given its
/// site
std::vector<BasisSet> StrainDoFTraits::construct_site_bases(
    Structure const &_prim,
    std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster> > > &_asym_unit,
    BasisFunctionSpecs const &_basis_function_specs) const {
  // std::cout << "Using " << func_type << " site basis functions." << std::endl
  // << std::endl;
  if (_prim.structure().global_dofs().find(name()) ==
      _prim.structure().global_dofs().end())
    return std::vector<BasisSet>();

  std::vector<BasisSet> result(1);
  result[0].set_variable_basis(adapter::Adapter<CASM::DoFSet, xtal::DoFSet>()(
      _prim.structure().global_dof(name()),
      _prim.global_dof_symrep_ID(name())));

  return result;
}

/// \brief Apply DoF values for this DoF to _struc
void StrainDoFTraits::apply_dof(ConfigDoF const &_dof,
                                BasicStructure const &_reference,
                                SimpleStructure &_struc) const {
  this->apply_standard_values(_dof.global_dof(name()).standard_values(),
                              _struc);
}

/// \brief Apply DoF or property values to _struc
void StrainDoFTraits::apply_standard_values(
    Eigen::MatrixXd const &standard_values, SimpleStructure &_struc) const {
  Eigen::VectorXd const &unrolled_metric = standard_values;
  StrainConverter c(m_metric);
  Eigen::Matrix3d F = c.unrolled_strain_metric_to_F(unrolled_metric);

  _struc.lat_column_mat = F * _struc.lat_column_mat;
  _struc.mol_info.coords = F * _struc.mol_info.coords;
  _struc.atom_info.coords = F * _struc.atom_info.coords;
}

jsonParser StrainDoFTraits::dof_to_json(
    ConfigDoF const &_dof, BasicStructure const &_reference) const {
  Eigen::VectorXd unrolled_metric = _dof.global_dof(name()).standard_values();
  StrainConverter c(m_metric);
  Eigen::Matrix3d F = c.unrolled_strain_metric_to_F(unrolled_metric);

  jsonParser json;
  json["deformation"] = F;
  to_json_array(unrolled_metric, json["value"]);
  return json;
}

}  // namespace DoF_impl

namespace DoFType {
DoF_impl::StrainDoFTraits GLstrain() { return DoF_impl::StrainDoFTraits("GL"); }

DoF_impl::StrainDoFTraits EAstrain() { return DoF_impl::StrainDoFTraits("EA"); }

DoF_impl::StrainDoFTraits Hstrain() { return DoF_impl::StrainDoFTraits("H"); }

}  // namespace DoFType
}  // namespace CASM
