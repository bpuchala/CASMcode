#ifndef CASM_DisplacementDoFTraits
#define CASM_DisplacementDoFTraits

#include "casm/basis_set/DoFTraits.hh"

namespace CASM {
namespace DoF_impl {
class DisplacementDoFTraits : public DoFType::Traits {
 public:
  DisplacementDoFTraits() : DoFType::Traits(AnisoValTraits::disp()) {}

  /// \brief Apply DoF values for this DoF to _struc
  void apply_dof(ConfigDoF const &_dof, BasicStructure const &_reference,
                 SimpleStructure &_struc) const override;

  /// \brief Apply DoF or property values to _struc
  void apply_standard_values(Eigen::MatrixXd const &standard_values,
                             SimpleStructure &_struc) const override;

  std::vector<BasisSet> construct_site_bases(
      Structure const &_prim,
      std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster> > > &_asym_unit,
      BasisFunctionSpecs const &_basis_function_specs) const override;

 protected:
  DoFType::Traits *_clone() const override;
};
}  // namespace DoF_impl

namespace DoFType {
DoF_impl::DisplacementDoFTraits displacement();
}

}  // namespace CASM
#endif
