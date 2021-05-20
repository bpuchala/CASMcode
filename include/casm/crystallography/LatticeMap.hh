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

  /// \brief Volumetric factor, used to normalize the strain cost to make it
  /// volume-independent
  ///
  /// BP: I think there might be a inconsistency here...
  ///
  /// using deformation_gradient of child -> parent,
  ///     parent * N == deformation_gradient * child
  ///
  /// vol_factor == pow(abs(deformation_gradient.determinant()), 1./3.)
  ///            ?= pow(std::abs(volume(child) / volume(parent)), 1./3.)
  static double vol_factor(Eigen::Matrix3d const &deformation_gradient) {
    return pow(std::abs(deformation_gradient.determinant()), 1. / 3.);
  }

  /// \brief Volumetric factor, used to normalize the strain cost to make it
  /// volume-independent
  ///
  /// BP: I think there might be a inconsistency here...
  ///
  /// using deformation_gradient of child -> parent,
  ///     parent * N == deformation_gradient * child
  ///
  /// vol_factor == pow(abs(deformation_gradient.determinant()), 1./3.)
  ///            == pow(std::abs(volume(child) / volume(parent)), 1./3.)
  static double vol_factor(Lattice const &parent, Lattice const &child) {
    return pow(std::abs(volume(child) / volume(parent)), 1. / 3.);
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

  /// Return name of method used by `strain_cost(Eigen::Matrix3d const &)`
  std::string cost_method() const;

 private:
  Eigen::MatrixXd m_gram_mat;
  bool m_sym_cost;

  mutable Eigen::Matrix3d m_cache;
  mutable Eigen::Matrix3d m_cache_inv;
};

/// Find the parent mapping of Lattice _parent onto Lattice _child
/// Denoting _parent.lat_column_mat() as 'parent' and _child.lat_column_mat() as
/// 'child', we want a mapping
///            child = deformation_gradient*parent*N
/// where deformation_gradient is an arbitrary 3x3 matrix and 'N' is a 3x3
/// unimodular (i.e., determinant=+/-1) integer matrix such that
/// cost(deformation_gradient) is minimized with respect to the matrix 'N' For
/// the cost function we use:
///         cost(deformation_gradient) = (_child.volume()/num_atoms)^(2/3) *
///         trace(D.transpose()*D) / g
/// where 'D' is the isovolumetric strain
///         D = deformation_gradient/det(deformation_gradient)^(1/3)-Identity
/// and 'g' is a const geometric factor.  The cost function corresponds to the
/// mean-square displacement of a point on the surface of a sphere having
/// V=_child.volume()/num_atoms (i.e., the atomic volume of the child crystal)
/// when the sphere is deformed at constant volume by
/// deformation_gradient/det(deformation_gradient)^(1/3)

/// LatticeMap finds mappings between lattices that minimize a strain cost.
///
/// Note: this documentation is attempting to use the notation and conventions
/// from the paper `Comparing crystal structures with symmetry and geometry`,
/// by John C Thomas, Anirudh Raju Natarajan, Anton Van der Ven
///
/// Goal:
/// - Find mappings (F^{N}, N) between lattices (L1 * T1) (parent) and L2
///   (child), of the form
///       (L1 * T1) * N = F^{N} * L2,
///   that minimize a strain cost, strain_cost(F^{N}).
///
/// Variables:
/// - parent (L1 * T1), child (L2) (Eigen::Matrix3d): the lattices to be
///   mapped, represented as 3x3 matrices whose columns are the lattice vectors
/// - N (Eigen::Matrix3d): a unimodular matrix (integer valued, though
///   represented with a floating point matrix for multiplication, with
///   N.determinant()==1) that generates lattice vectors (L1 * T1 * N) for
///   lattices that are equivalent to the parent lattice (L1 * T1).
/// - F^{N} (Eigen::Matrixd): a deformation gradient,
///        (L1 * T1) * N == F^{N} * L2,
///   from child lattice vectors, L2, to lattice vectors (L1 * T1 * N).
///
/// Some relations:
///   The deformation gradient can be decomposed into a stretch tensor and
///   isometry,
///       F^{N} == V^{N} * Q^{N}.
///   It can also be defined in the reverse direction (parent -> child),
///       F_rev^{N} * (L1 * T1) * N == L2,
///   and decomposed into the reverse stretch tensor and isometry
///       F_rev^{N} = V_rev^{N} * Q_rev^{N},
///   Which are related via:
///       V_rev^{N} == (V^{N}).inverse()
///       Q_rev^{N} == (Q^{N}).inverse() == (Q^{N}).transpose()
///
/// Strain cost:
/// - If `_symmetrize_strain_cost==false`, use a volume-normalized strain cost,
///
///       \tilde{V}^N = \frac{1}{\det{V}^{1/3}} V^{N}
///       \tilde{B}^{N} = \tilde{V}^N - I
///       \tilde{B_rev}^{N} = \tilde{V_rev}^N - I
///
/// The search for minimum cost mappings is best done using the reduced cells
/// of the parent and child lattices:
///
///     reduced_parent = parent.reduced_cell(),
///     reduced_parent == parent * transformation_matrix_to_reduced_parent,
///
///     reduced_child = child.reduced_cell(),
///     reduced_child == child * transformation_matrix_to_reduced_child,
///
/// Then,
///     reduced_parent * N_reduced == F_reduced^{N} * reduced_child
///     parent * transformation_matrix_to_reduced_parent * N_reduced ==
///         F_reduced^{N} * child * transformation_matrix_to_reduced_child
///
/// Thus,
///     N == transformation_matrix_to_reduced_parent * N_reduced *
///         transformation_matrix_to_reduced_child.inverse(),
///     F^{N} == F_reduced^{N}
///
///
class LatticeMap {
 public:
  typedef Eigen::Matrix<double, 3, 3, Eigen::DontAlign> DMatType;
  typedef Eigen::Matrix<int, 3, 3, Eigen::DontAlign> IMatType;

  LatticeMap(Lattice const &_parent, Lattice const &_child, Index _num_atoms,
             int _range, SymOpVector const &_parent_point_group,
             SymOpVector const &_child_point_group,
             Eigen::Ref<const Eigen::MatrixXd> const &strain_gram_mat =
                 Eigen::MatrixXd::Identity(9, 9),
             double _init_better_than = 1e20,
             bool _symmetrize_strain_cost = false, double _xtal_tol = TOL);

  LatticeMap(Eigen::Ref<const DMatType> const &_parent,
             Eigen::Ref<const DMatType> const &_child, Index _num_atoms,
             int _range, SymOpVector const &_parent_point_group,
             SymOpVector const &_child_point_group,
             Eigen::Ref<const Eigen::MatrixXd> const &strain_gram_mat =
                 Eigen::MatrixXd::Identity(9, 9),
             double _init_better_than = 1e20,
             bool _symmetrize_strain_cost = false, double _xtal_tol = TOL);

  void reset(double _better_than = 1e20);

  // Finds the smallest strain tensor (in terms of Frobenius norm) that deforms
  // (*this) into a lattice symmetrically equivalent to 'child_lattice'
  LatticeMap const &best_strain_mapping() const;

  LatticeMap const &next_mapping_better_than(double max_cost) const;

  double strain_cost() const { return m_cost; }

  double calc_strain_cost(const Eigen::Matrix3d &deformation_gradient) const;

  /// The name of the method used to calculate the lattice deformation cost
  std::string cost_method() const;

  /// This is N
  const DMatType &matrixN() const { return m_N; }

  /// This is F_rev^{N}
  ///
  ///     (L1 * T1) * N = F^{N} * L2
  ///     F_rev^{N} * (L1 * T1) * N = L2
  const DMatType &deformation_gradient() const {
    return m_deformation_gradient;
  }

  /// This is reduced_parent
  const DMatType &parent_matrix() const { return m_parent; }

  /// This is reduced_child
  const DMatType &child_matrix() const { return m_child; }

  double xtal_tol() const { return m_xtal_tol; }

  /// \brief If true, use cost of symmetry breaking strain; else use isotropic
  /// strain cost or anisotropic strain cost, depending on strain_gram_mat
  bool symmetrize_strain_cost() const { return m_symmetrize_strain_cost; }

 private:
  /// These are reduced_parent, and reduced_child
  DMatType m_parent, m_child;

  // Conversion matrices:
  // N = m_transformation_matrix_to_reduced_parent * N_reduced *
  //     m_transformation_matrix_to_reduced_child.inverse(),
  // where N_reduced == inv_mat().inverse()
  DMatType m_transformation_matrix_to_reduced_parent;
  DMatType m_transformation_matrix_to_reduced_child_inv;

  StrainCostCalculator m_calc;

  // m_vol_factor == pow( std::abs(volume(child) / volume(parent)), 1. / 3.)
  //              == pow( (V^{N}).determinant(), 1./3.)
  double m_vol_factor;

  int m_range;

  // pointer to static list of unimodular matrices
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
  double m_xtal_tol;

  mutable double m_cost;
  mutable Index m_currmat;
  mutable DMatType m_deformation_gradient, m_N, m_dcache;
  mutable IMatType m_icache;

  ///\brief Returns the inverse of the current transformation matrix under
  /// consideration
  // We treat the unimodular matrices as the inverse of the transformation
  // matrices that we are actively considering, allowing fewer matrix inversions
  Eigen::Matrix3i const &inv_mat() const { return (*m_mvec_ptr)[m_currmat]; }

  ///\brief Number of unimodular matrices
  Index n_mat() const { return m_mvec_ptr->size(); }

  /// \brief Returns true if current transformation is the canonical equivalent
  bool _check_canonical() const;

  LatticeMap const &_next_mapping_better_than(double max_cost) const;

  // use m_deformation_gradient to calculate strain cost
  double _calc_strain_cost() const;
};

/** @} */
}  // namespace xtal
}  // namespace CASM
#endif
