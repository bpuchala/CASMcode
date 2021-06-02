#ifndef LATTICEMAP_HH
#define LATTICEMAP_HH

#include "casm/container/Counter.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/global/definitions.hh"

namespace CASM {
namespace xtal {

struct SymOp;
typedef std::vector<SymOp> SymOpVector;

/** \ingroup Lattice
 *  @{
 */

class StrainCostCalculator {
 public:
  StrainCostCalculator(Eigen::Ref<const Eigen::MatrixXd> const
                           &strain_gram_mat = Eigen::MatrixXd::Identity(9, 9));

  //\brief Isotropic strain cost, without gram matrix
  static double isotropic_strain_cost(
      Eigen::Matrix3d const &_deformation_gradient);

  //\brief Isotropic strain cost, without gram matrix
  static double isotropic_strain_cost(
      Eigen::Matrix3d const &_deformation_gradient, double _vol_factor);

  // \brief Volumetric factor :
  // pow(abs(_deformation_gradient.determinant()),1./3.), used to normalize the
  // strain cost to make it volume-independent
  static double vol_factor(Eigen::Matrix3d const &_deformation_gradient) {
    return pow(std::abs(_deformation_gradient.determinant()), 1. / 3.);
  }

  //\brief Anisotropic strain cost; utilizes stored gram matrix to compute
  // strain cost
  double strain_cost(Eigen::Matrix3d const &_deformation_gradient) const;

  //\brief Anisotropic strain cost; utilizes stored gram matrix to compute
  // strain cost
  double strain_cost(Eigen::Matrix3d const &_deformation_gradient,
                     double _vol_factor) const;

  //\brief Symmetrized strain cost; Utilizes the parent point group symmetry to
  // calculate only the symmetry breaking lattice cost
  double strain_cost(Eigen::Matrix3d const &_deformation_gradient,
                     SymOpVector const &parent_point_group) const;

 private:
  Eigen::MatrixXd m_gram_mat;
  bool m_sym_cost;

  mutable Eigen::Matrix3d m_cache;
  mutable Eigen::Matrix3d m_cache_inv;
};

/// \brief LatticeMap finds optimal mappings between lattices
class LatticeMap {
 public:
  typedef Eigen::Matrix<double, 3, 3, Eigen::DontAlign> DMatType;
  typedef Eigen::Matrix<int, 3, 3, Eigen::DontAlign> IMatType;

  /// LatticeMap constructor
  LatticeMap(Lattice const &_parent, Lattice const &_child, int _range,
             SymOpVector const &_parent_point_group,
             SymOpVector const &_child_point_group,
             double _init_better_than = 1e20,
             bool _symmetrize_strain_cost = false, double _cost_tol = TOL);

  /// LatticeMap constructor
  LatticeMap(Eigen::Ref<const DMatType> const &_parent,
             Eigen::Ref<const DMatType> const &_child, int _range,
             SymOpVector const &_parent_point_group,
             SymOpVector const &_child_point_group,
             double _init_better_than = 1e20,
             bool _symmetrize_strain_cost = false, double _cost_tol = TOL);

  /// Iterate until all possible solutions have been considered
  LatticeMap const &best_strain_mapping() const;

  /// \brief Iterate until the next solution (N, F^{N}) with lattice mapping
  /// score less than `max_cost` is found.
  LatticeMap const &next_mapping_better_than(double max_cost) const;

  /// The cost of the current solution
  double strain_cost() const { return m_cost; }

  /// The name of the method used to calculate the lattice deformation cost
  std::string cost_method() const;

  /// This is N for the current solution
  const DMatType &matrixN() const { return m_N; }

  /// \brief This is F_reverse^{N} (parent to child deformation) for the
  /// current solution
  ///
  ///     F_reverse^{N} * (L1 * T1) * N = L2
  const DMatType &deformation_gradient() const {
    return m_deformation_gradient;
  }

  /// This is the parent lattice column matrix (L1 * T1)
  const DMatType &parent_matrix() const { return m_parent; }

  /// This is the child lattice column matrix (L2)
  const DMatType &child_matrix() const { return m_child; }

  /// This is the reduced cell of the parent lattice column matrix
  const DMatType &reduced_parent_matrix() const { return m_reduced_parent; }

  /// This is the reduced cell of the child lattice column matrix
  const DMatType &reduced_child_matrix() const { return m_reduced_child; }

  /// This is the tolerance used for cost comparisons
  double cost_tol() const { return m_cost_tol; }

  /// \brief If true, use the symmetry breaking strain cost; else use the
  /// isotropic strain cost
  bool symmetrize_strain_cost() const { return m_symmetrize_strain_cost; }

 private:
  /// These are the original, not reduced, parent and child lattice column
  /// matrices
  DMatType m_parent, m_child;

  /// These are reduced cell parent and child lattice column matrices
  DMatType m_reduced_parent, m_reduced_child;

  // Conversion matrices:
  // m_reduced_parent = m_parent * m_transformation_matrix_to_reduced_parent
  // m_reduced_child * m_transformation_matrix_to_reduced_child_inv = m_child
  // N = m_transformation_matrix_to_reduced_parent * N_reduced *
  //     m_transformation_matrix_to_reduced_child.inverse(),
  // where N_reduced == inv_mat().inverse()
  DMatType m_transformation_matrix_to_reduced_parent;
  DMatType m_transformation_matrix_to_reduced_child_inv;

  // The absolute value of an element of N is not allowed be be greater than
  // m_range.
  int m_range;

  // pointer to static list of unimodular matrices (used as N.inverse())
  std::vector<Eigen::Matrix3i> const *m_mvec_ptr;

  // parent point group matrices, in fractional coordinates
  std::vector<Eigen::Matrix3i> m_parent_fsym_mats;

  // parent point group in cartesian coordinates
  SymOpVector m_parent_point_group;

  // child point group matrices, in fractional coordinates
  std::vector<Eigen::Matrix3i> m_child_fsym_mats;

  // flag indicating if the symmetrized strain cost should be used while
  // searching for the best lattice maps
  bool m_symmetrize_strain_cost;
  double m_cost_tol;

  mutable double m_cost;
  mutable Index m_currmat;
  mutable DMatType m_deformation_gradient, m_N, m_dcache;
  mutable IMatType m_icache;

  void _reset(double _better_than = 1e20);

  /// Calculate the strain cost of a deformation
  double _calc_strain_cost(const Eigen::Matrix3d &deformation_gradient) const;

  ///\brief Returns the inverse of the current transformation matrix under
  /// consideration
  // We treat the unimodular matrices as the inverse of the transformation
  // matrices that we are actively considering, allowing fewer matrix inversions
  Eigen::Matrix3i const &inv_mat() const { return (*m_mvec_ptr)[m_currmat]; }

  /// Number of unimodular matrices
  Index n_mat() const { return m_mvec_ptr->size(); }

  /// Returns true if current N matrix is the canonical equivalent
  bool _check_canonical() const;

  /// \brief Iterate until the next solution (N, F^{N}) with lattice mapping
  /// score less than `max_cost` is found.
  LatticeMap const &_next_mapping_better_than(double max_cost) const;
};

/// \brief Returns the volume-normalized strain cost, calculated to be
/// invariant to which structure is the child/parent.
double isotropic_strain_cost(Eigen::Matrix3d const &deformation_gradient);

/// \brief Returns the symmetrized stretch
Eigen::Matrix3d symmetrized_stretch(Eigen::Matrix3d const &deformation_gradient,
                                    SymOpVector const &parent_point_group);

/// \brief Returns the symmetry-breaking component of the stretch
Eigen::Matrix3d symmetry_breaking_stretch(
    Eigen::Matrix3d const &deformation_gradient,
    SymOpVector const &parent_point_group);

/// \brief Returns the symmetry-breaking strain cost, calculated to be
/// invariant to which structure is the child/parent.
double symmetry_breaking_strain_cost(
    Eigen::Matrix3d const &deformation_gradient,
    SymOpVector const &parent_point_group);

/** @} */
}  // namespace xtal
}  // namespace CASM
#endif
