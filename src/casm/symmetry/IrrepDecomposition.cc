#include "casm/symmetry/IrrepDecomposition.hh"

#include <iostream>

#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/symmetry/IrrepDecompositionImpl.hh"

namespace CASM {

namespace SymRepTools_v2 {

IrrepInfo::IrrepInfo(Eigen::MatrixXcd _trans_mat, Eigen::VectorXcd _characters)
    : trans_mat(std::move(_trans_mat)),
      characters(std::move(_characters)),
      complex(!almost_zero(trans_mat.imag())),
      pseudo_irrep(false),
      index(0) {}

/// Construct a "dummy" IrrepInfo with user specified transformtion matrix
///
/// The "dummy" IrrepInfo is constructed with specified transformtion matrix
/// and character vector of [(dim,0)] where 'dim' is the dimension of irrep
/// (number of rows of `trans_mat`)
IrrepInfo make_dummy_irrep_info(Eigen::MatrixXcd const &trans_mat) {
  Eigen::VectorXcd tchar(1);
  tchar(0) = std::complex<double>(double(trans_mat.rows()), 0.);
  return IrrepInfo(trans_mat, tchar);
}

/// \brief Assumes that irreps are real, and concatenates their individual
/// trans_mats to form larger trans_mat
Eigen::MatrixXd full_trans_mat(std::vector<IrrepInfo> const &irreps) {
  Index row = 0;
  Index col = 0;
  for (auto const &irrep : irreps) {
    col = irrep.vector_dim();
    row += irrep.irrep_dim();
  }
  Eigen::MatrixXd trans_mat(row, col);
  row = 0;
  for (auto const &irrep : irreps) {
    trans_mat.block(row, 0, irrep.irrep_dim(), irrep.vector_dim()) =
        irrep.trans_mat.real();
    row += irrep.irrep_dim();
  }
  return trans_mat;
}

/// IrrepDecomposition constructor
///
/// \param rep Full space matrix representation (rep[0].rows() ==
///     init_subspace.rows())
/// \param head_group Group for which the irreps are to be found
/// \param _subspace Input subspace in which irreps are to be found. Will be
///     expanded by application of rep and orthogonalization to form an
///     invariant subspace (i.e. column space dimension is not increased by
///     application of elements in head_group)
///
IrrepDecomposition::IrrepDecomposition(
    MatrixRep const &_fullspace_rep, GroupIndices const &_head_group,
    Eigen::MatrixXd const &init_subspace,
    GroupIndicesOrbitVector const &_cyclic_subgroups,
    GroupIndicesOrbitVector const &_all_subgroups, bool allow_complex)
    : fullspace_rep(_fullspace_rep),
      head_group(_head_group),
      cyclic_subgroups(_cyclic_subgroups),
      all_subgroups(_all_subgroups) {
  using namespace IrrepDecompositionImpl;

  // 1) Expand subspace by application of group, and orthonormalization
  subspace = make_invariant_space(fullspace_rep, head_group, init_subspace);

  // 2) Create `subspace_rep`, a transformed copy of `fullspace_rep` that acts
  //    on coordinates with `subspace` columns as their basis. Matrices in
  //    `subspace_rep` are shape (subspace.cols() x subspace.cols())
  subspace_rep = make_subspace_rep(fullspace_rep, subspace);

  std::cout << "fullspace_rep dim: " << fullspace_rep[0].rows() << std::endl;
  std::cout << "subspace_rep dim: " << subspace_rep[0].rows() << std::endl;

  // 3) Perform irrep_decomposition
  std::vector<IrrepInfo> subspace_irreps =
      irrep_decomposition(subspace_rep, head_group, allow_complex);

  // 4) Symmetrize subspace irreps
  std::vector<IrrepInfo> symmetrized_subspace_irreps =
      symmetrize_irreps(subspace_rep, head_group, subspace_irreps,
                        cyclic_subgroups, all_subgroups);

  for (auto const &irrep : symmetrized_subspace_irreps) {
    std::cout << "---" << std::endl;
    std::cout << "symmetrized irrep: " << std::endl;
    std::cout << "index: " << irrep.index << std::endl;
    std::cout << "characters: " << prettyc(irrep.characters.transpose())
              << std::endl;
    std::cout << "subspace: \n" << prettyc(irrep.trans_mat) << std::endl;
    std::cout << "directions.size() (number of orbits): "
              << irrep.directions.size() << std::endl;
    Index orbit_index = 0;
    for (auto const &orbit : irrep.directions) {
      std::cout << "-" << std::endl;
      std::cout << "orbit: " << orbit_index << std::endl;
      for (auto const &direction : orbit) {
        std::cout << prettyc(direction.transpose()) << std::endl;
      }
      ++orbit_index;
    }
  }

  // 5) Transform to fullspace irreps
  irreps = make_fullspace_irreps(symmetrized_subspace_irreps, subspace);

  // 6) Combine to form symmetry adapted subspace
  symmetry_adapted_subspace = full_trans_mat(irreps).adjoint();
}

}  // namespace SymRepTools_v2

}  // namespace CASM