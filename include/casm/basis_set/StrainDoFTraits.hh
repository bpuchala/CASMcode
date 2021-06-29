#ifndef CASM_StrainDoFTraits
#define CASM_StrainDoFTraits
#include "casm/basis_set/DoFTraits.hh"

namespace CASM {
namespace DoF_impl {
class StrainDoFTraits : public DoFType::Traits {
 public:
  StrainDoFTraits(std::string const &_metric)
      : DoFType::Traits(AnisoValTraits::strain(_metric)), m_metric(_metric) {}

  /// \brief Retrieve the standard values for a DoF from dictionary of
  /// properties from properties.calc.json
  ///  Returns matrix with standard values, and names of properties that were
  ///  used to construct the matrix
  std::pair<Eigen::MatrixXd, std::set<std::string> > find_values(
      std::map<std::string, Eigen::MatrixXd> const &values) const override;

  /// \brief Construct the site basis (if DOF_MODE is LOCAL) for a DoF, given
  /// its site
  std::vector<BasisSet> construct_site_bases(
      Structure const &_prim,
      std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster> > > &_asym_unit,
      BasisFunctionSpecs const &_basis_function_specs) const override;

  /// \brief Serialize strain DoF values from ConfigDoF
  jsonParser dof_to_json(ConfigDoF const &_dof,
                         BasicStructure const &_reference) const override;

  /// \brief Apply DoF values for this DoF to _struc
  void apply_dof(ConfigDoF const &_dof, BasicStructure const &_reference,
                 SimpleStructure &_struc) const override;

  /// \brief Apply DoF or property values to _struc
  void apply_standard_values(Eigen::MatrixXd const &standard_values,
                             SimpleStructure &_struc) const override;

 protected:
  DoFType::Traits *_clone() const override {
    return new StrainDoFTraits(*this);
  }

  std::string const m_metric;
};
}  // namespace DoF_impl

namespace DoFType {
DoF_impl::StrainDoFTraits GLstrain();

DoF_impl::StrainDoFTraits EAstrain();

DoF_impl::StrainDoFTraits Hstrain();

}  // namespace DoFType
}  // namespace CASM
#endif
