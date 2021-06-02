#include "casm/crystallography/LatticeMap.hh"

#include <iterator>

#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/LatticeIsEquivalent.hh"
#include "casm/crystallography/Strain.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/external/Eigen/src/Core/Matrix.h"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {
namespace xtal {

StrainCostCalculator::StrainCostCalculator(
    Eigen::Ref<const Eigen::MatrixXd> const
        &strain_gram_mat /*= Eigen::MatrixXd::Identity(9,9)*/) {
  if (strain_gram_mat.size() == 0 || strain_gram_mat.isIdentity(1e-9)) {
    m_sym_cost = false;
  } else {
    m_sym_cost = true;
    m_gram_mat.resize(6, 6);
    if (strain_gram_mat.rows() == 6 && strain_gram_mat.cols() == 6) {
      std::vector<Index> map({0, 5, 4, 1, 3, 2});
      for (Index i = 0; i < 6; ++i) {
        for (Index j = 0; j < 6; ++j) {
          m_gram_mat(i, j) = strain_gram_mat(map[i], map[j]);
          if (i > 2) m_gram_mat(i, j) *= sqrt(2.);
          if (j > 2) m_gram_mat(i, j) *= sqrt(2.);
        }
      }
    }
    if (strain_gram_mat.rows() == 9 && strain_gram_mat.cols() == 9) {
      Index m = 0;
      for (Index i = 0; i < 3; ++i) {
        for (Index j = i; j < 3; ++j, ++m) {
          Index n = 0;
          for (Index k = 0; k < 3; ++k) {
            for (Index l = k; l < 3; ++l, ++n) {
              m_gram_mat(m, n) = strain_gram_mat(i * 3 + j, k * 3 + l);
              if (m > 2) m_gram_mat(m, n) *= sqrt(2.);
              if (n > 2) m_gram_mat(m, n) *= sqrt(2.);
            }
          }
        }
      }
    }
  }
}

//*******************************************************************************************
// static function
double StrainCostCalculator::isotropic_strain_cost(
    const Eigen::Matrix3d &_deformation_gradient, double _vol_factor) {
  Eigen::Matrix3d tmat =
      polar_decomposition(_deformation_gradient / _vol_factor);

  // -> epsilon=(_deformation_gradient_deviatoric-identity)
  return ((tmat - Eigen::Matrix3d::Identity(3, 3)).squaredNorm() +
          (tmat.inverse() - Eigen::Matrix3d::Identity(3, 3)).squaredNorm()) /
         6.;
}

//*******************************************************************************************
// static function
double StrainCostCalculator::isotropic_strain_cost(
    const Eigen::Matrix3d &_deformation_gradient) {
  return isotropic_strain_cost(_deformation_gradient,
                               vol_factor(_deformation_gradient));
}

//*******************************************************************************************

// strain_cost is the mean-square displacement of a point on the surface of a
// unit sphere when it is deformed by the volume-preserving deformation
// _deformation_gradient_deviatoric =
// _deformation_gradient/det(_deformation_gradient)^(1/3)
double StrainCostCalculator::strain_cost(
    const Eigen::Matrix3d &_deformation_gradient, double _vol_factor) const {
  if (m_sym_cost) {
    double cost = 0;
    m_cache = polar_decomposition(_deformation_gradient / _vol_factor);
    m_cache_inv = m_cache.inverse() - Eigen::Matrix3d::Identity(3, 3);
    m_cache -= Eigen::Matrix3d::Identity(3, 3);
    Index m = 0;
    for (Index i = 0; i < 3; ++i) {
      for (Index j = i; j < 3; ++j, ++m) {
        Index n = 0;
        for (Index k = 0; k < 3; ++k) {
          for (Index l = k; l < 3; ++l, ++n) {
            cost += m_gram_mat(m, n) *
                    (m_cache(i, j) * m_cache(j, k) +
                     m_cache_inv(i, j) * m_cache_inv(j, k)) /
                    6.;
          }
        }
      }
    }
    // geometric factor: (3*V/(4*pi))^(2/3)/3 = V^(2/3)/7.795554179
    return cost;
  }

  return isotropic_strain_cost(_deformation_gradient, _vol_factor);
}

//*******************************************************************************************

double StrainCostCalculator::strain_cost(
    const Eigen::Matrix3d &_deformation_gradient) const {
  return strain_cost(_deformation_gradient, vol_factor(_deformation_gradient));
}

//*******************************************************************************************
double StrainCostCalculator::strain_cost(
    Eigen::Matrix3d const &_deformation_gradient,
    SymOpVector const &parent_point_group) const {
  // Apply the sym op to the deformation gradient and average
  Eigen::Matrix3d stretch_aggregate = Eigen::Matrix3d::Zero();
  Eigen::Matrix3d stretch = strain::right_stretch_tensor(_deformation_gradient);
  for (auto const &op : parent_point_group) {
    stretch_aggregate += op.matrix * stretch * op.matrix.inverse();
  }
  stretch_aggregate = stretch_aggregate / double(parent_point_group.size());
  return strain_cost(stretch - stretch_aggregate + Eigen::Matrix3d::Identity(),
                     1.0);
}

/// \class LatticeMap
/// \brief LatticeMap finds optimal mappings between lattices
///
/// Note: this documentation uses the notation and conventions from the paper
/// `Comparing crystal structures with symmetry and geometry`,
/// by John C Thomas, Anirudh Raju Natarajan, Anton Van der Ven
///
/// Goal:
/// - Find mapping solutions, (N, F^{N}), between lattices (L1 * T1) (parent)
///   and L2 (child), of the form
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
///       F_reverse^{N} * (L1 * T1) * N == L2,
///   and decomposed into the reverse stretch tensor and isometry
///       F_reverse^{N} = V_reverse^{N} * Q_reverse^{N},
///   Which are related via:
///       V_reverse^{N} == (V^{N}).inverse()
///       Q_reverse^{N} == (Q^{N}).inverse() == (Q^{N}).transpose()
///   The Biot strain tensor, B^{N}, and reverse:
///       B^{N} = V^{N} - I
///       B_reverse^{N} = V_reverse^{N} - I
///
/// Strain cost:
/// - 'isotropic_strain_cost': (default, used if
///  `_symmetrize_strain_cost==false`) a volume-normalized strain cost,
///   calculated to be invariant to which structure is the child/parent. The
///  `\tilde` indicates that the value is normalized to be volume invariant.
///
///       \tilde{V} = \frac{1}{\det{V}^{1/3}} V
///       \tilde{B} = \tilde{V} - I
///
///   For reverse cost, use V.inverse() = U_reverse
///       \tilde{U_reverse} = \frac{1}{\det{U_reverse}^{1/3}} U_reverse
///       \tilde{B_reverse} = \tilde{U_reverse} - I
///
///       strain_cost = (1./2.)*(
///           (1./3.)*\tr{\tilde{B}^{2}} +
///           (1./3.)*\tr{\tilde{B_reverse}^{2}})
///
/// - 'symmetry_breaking_strain_cost': (used if
///   `_symmetrize_strain_cost==true`) a strain cost, calculated with only the
///   the components of the strain that break the symmetry of the parent
///   structure.
///
///       symmetry_breaking_strain_cost =
///           (1./2.)*(
///               (1./3.)*trace( (B_sym_break)^2 ) +
///               (1./3.)*trace( (B_sym_break_reverse)^2 ) )
///
///       B_sym_break = B - B_sym, the symmetry-breaking Biot strain
///       B_sym = (1./N_G1) * sum_i ( G1_i * B * G1_i.transpose() ),
///       where G1_i are point group operations of parent structure, N_G1 is
///       the total number of operations,
///       and similarly for B_sym_break_reverse from B_reverse.
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

/// LatticeMap constructor
///
/// \param _parent Reference lattice (L1 * T1)
/// \param _child Lattice to be mapped to _parent (L2)
/// \param num_atoms Unused. (TODO: remove)
/// \param _range Determines range of N matrices to be searched when optimizing
///     the lattice mapping. The absolute value of an element of N is not
///     allowed be be greater than `_range`. Typically 1 is a good enough
///     choice.
/// \param _parent_point_group Point group of the parent (i.e. crystal point
///     group), used to identify canonical N matrices and reduce the number of
///     operations. Also used to calculate the symmetry-breaking strain cost if
///     `_symmetrize_strain_cost==true`.
/// \param _child_point_group Point group of the child (structure), used to
///     identify canonical N matrices and reduce the number of operations.
///     (Typically, just use {SymOp::identity()})
/// \param _init_better_than Initializes the LatticeMap to the first mapping
///     with a cost less than `_init_better_than + _cost_tol`. Use the default
///     value (1e20) to initialize to the first valid mapping.
/// \param _symmetrize_strain_cost Boolean flag (default=false), which if true
///     indicates that the `symmetry_breaking_strain_cost` should be used to
///     score lattice mappings.
/// \param _cost_tol Tolerance used for cost comparisons
LatticeMap::LatticeMap(const Lattice &_parent, const Lattice &_child,
                       int _range, SymOpVector const &_parent_point_group,
                       SymOpVector const &_child_point_group,
                       double _init_better_than /* = 1e20 */,
                       bool _symmetrize_strain_cost, double _cost_tol)
    : m_parent(_parent.lat_column_mat()),
      m_child(_child.lat_column_mat()),
      m_range(_range),
      m_symmetrize_strain_cost(_symmetrize_strain_cost),
      m_cost_tol(_cost_tol),
      m_cost(1e20),
      m_currmat(0) {
  Lattice reduced_parent = _parent.reduced_cell();
  m_reduced_parent = reduced_parent.lat_column_mat();

  Lattice reduced_child = _child.reduced_cell();
  m_reduced_child = reduced_child.lat_column_mat();

  // m_reduced_parent = m_parent * m_transformation_matrix_to_reduced_parent
  m_transformation_matrix_to_reduced_parent =
      _parent.inv_lat_column_mat() * m_reduced_parent;
  // m_reduced_child * m_transformation_matrix_to_reduced_child_inv = m_child
  m_transformation_matrix_to_reduced_child_inv =
      m_reduced_child.inverse() * _child.lat_column_mat();

  if (_range == 1)
    m_mvec_ptr = &unimodular_matrices<1>();
  else if (_range == 2)
    m_mvec_ptr = &unimodular_matrices<2>();
  else if (_range == 3)
    m_mvec_ptr = &unimodular_matrices<3>();
  else if (_range == 4)
    m_mvec_ptr = &unimodular_matrices<4>();
  else
    throw std::runtime_error(
        "LatticeMap cannot currently be invoked for range>4");

  // Construct inverse fractional symops for parent
  {
    IsPointGroupOp symcheck(reduced_parent);
    m_parent_fsym_mats.reserve(_parent_point_group.size());
    for (auto const &op : _parent_point_group) {
      if (!symcheck(op)) continue;
      m_parent_fsym_mats.push_back(
          iround(reduced_parent.inv_lat_column_mat() * op.matrix.transpose() *
                 reduced_parent.lat_column_mat()));
      for (Index i = 0; i < (m_parent_fsym_mats.size() - 1); ++i) {
        if (m_parent_fsym_mats[i] == m_parent_fsym_mats.back()) {
          m_parent_fsym_mats.pop_back();
          break;
        }
      }
    }
  }

  // Store the parent symmetry operations
  m_parent_point_group = _parent_point_group;

  // Construct fractional symops for child
  {
    IsPointGroupOp symcheck(reduced_child);
    m_child_fsym_mats.reserve(_child_point_group.size());
    for (auto const &op : _child_point_group) {
      if (!symcheck(op)) continue;
      m_child_fsym_mats.push_back(
          iround(reduced_child.inv_lat_column_mat() * op.matrix *
                 reduced_child.lat_column_mat()));
      for (Index i = 0; i < (m_child_fsym_mats.size() - 1); ++i) {
        if (m_child_fsym_mats[i] == m_child_fsym_mats.back()) {
          m_child_fsym_mats.pop_back();
          break;
        }
      }
    }
  }

  _reset(_init_better_than);
}

LatticeMap::LatticeMap(Eigen::Ref<const LatticeMap::DMatType> const &_parent,
                       Eigen::Ref<const LatticeMap::DMatType> const &_child,
                       int _range, SymOpVector const &_parent_point_group,
                       SymOpVector const &_child_point_group,
                       double _init_better_than /* = 1e20 */,
                       bool _symmetrize_strain_cost, double _cost_tol)
    : LatticeMap(Lattice(_parent), Lattice(_child), _range, _parent_point_group,
                 _child_point_group, _init_better_than, _symmetrize_strain_cost,
                 _cost_tol) {}

void LatticeMap::_reset(double _better_than) {
  m_currmat = 0;

  // From relation F * parent * inv_mat.inverse() = child
  m_deformation_gradient =
      m_reduced_child * inv_mat().cast<double>() *
      m_reduced_parent.inverse();  // -> _deformation_gradient

  double tcost = _calc_strain_cost(m_deformation_gradient);

  // Initialize to first valid mapping
  if (tcost <= _better_than && _check_canonical()) {
    m_cost = tcost;
    // reconstruct correct N for unreduced lattice
    m_N = m_transformation_matrix_to_reduced_parent *
          inv_mat().cast<double>().inverse() *
          m_transformation_matrix_to_reduced_child_inv;
  } else
    next_mapping_better_than(_better_than);
}

const LatticeMap &LatticeMap::best_strain_mapping() const {
  m_currmat = 0;

  // Get an upper bound on the best mapping by starting with no lattice
  // equivalence
  m_N = DMatType::Identity(3, 3);
  // m_dcache -> value of inv_mat() that gives m_N = identity;
  m_dcache = m_transformation_matrix_to_reduced_child_inv *
             m_transformation_matrix_to_reduced_parent;
  m_deformation_gradient =
      m_reduced_child * m_dcache * m_reduced_parent.inverse();

  double best_cost = _calc_strain_cost(m_deformation_gradient);

  while (next_mapping_better_than(best_cost).strain_cost() < best_cost) {
    best_cost = strain_cost();
  }

  m_cost = best_cost;
  return *this;
}

/// Calculate the strain cost of a deformation
///
/// \param deformation_gradient Parent to child deformation gradient,
///    F_reverse^{N}, where F_reverse^{N} * (L1 * T1) * N = L2.
///
double LatticeMap::_calc_strain_cost(
    const Eigen::Matrix3d &deformation_gradient) const {
  if (symmetrize_strain_cost())
    return symmetry_breaking_strain_cost(deformation_gradient,
                                         m_parent_point_group);
  else
    return isotropic_strain_cost(deformation_gradient);
}

/// The name of the method used to calculate the lattice deformation cost
///
/// \returns "isotropic_strain_cost", "anisotropic_strain_cost", or
///   "symmetry_breaking_strain_cost", as determined by constructor arguments
///
std::string LatticeMap::cost_method() const {
  if (symmetrize_strain_cost()) {
    return "symmetry_breaking_strain_cost";
  } else {
    return "isotropic_strain_cost";
  }
}

/// \brief Iterate until the next solution (N, F^{N}) with lattice mapping
/// score less than `max_cost` is found.
///
/// The current solution is:
/// - F_reverse^{N} = this->deformation_cost()
/// - N = this->matrixN()
///
/// More solutions may be found as long as
///     `this->strain_cost() < max_cost + this->cost_tol()`.
/// Otherwise, iteration is complete.
///
const LatticeMap &LatticeMap::next_mapping_better_than(double max_cost) const {
  m_cost = 1e20;
  return _next_mapping_better_than(max_cost);
}

/// \brief Iterate until the next solution (N, F^{N}) with lattice mapping
/// score less than `max_cost` is found.
const LatticeMap &LatticeMap::_next_mapping_better_than(double max_cost) const {
  DMatType init_deformation_gradient(m_deformation_gradient);
  // tcost initial value shouldn't matter unles m_inv_count is invalid
  double tcost = max_cost;

  while (++m_currmat < n_mat()) {
    if (!_check_canonical()) {
      continue;
    }

    // From relation _deformation_gradient * parent * inv_mat.inverse() = child
    m_deformation_gradient =
        m_reduced_child * inv_mat().cast<double>() *
        m_reduced_parent.inverse();  // -> _deformation_gradient
    tcost = _calc_strain_cost(m_deformation_gradient);
    if (std::abs(tcost) < (std::abs(max_cost) + std::abs(cost_tol()))) {
      m_cost = tcost;

      // need to undo the effect of transformation to reduced cell on 'N'
      // Maybe better to get m_N from m_deformation_gradient instead?
      // m_transformation_matrix_to_reduced_parent and
      // m_transformation_matrix_to_reduced_child_inv depend on the lattice
      // reduction that was performed in the constructor, so we would need to
      // store "non-reduced" parent and child
      m_N = m_transformation_matrix_to_reduced_parent *
            inv_mat().cast<double>().inverse() *
            m_transformation_matrix_to_reduced_child_inv;
      // std::cout << "N:\n" << m_N << "\n";
      //  We already have:
      //        m_deformation_gradient = m_reduced_child *
      //        inv_mat().cast<double>() * m_reduced_parent.inverse();
      break;
    }
  }
  if (!(std::abs(tcost) < (std::abs(max_cost) + std::abs(cost_tol())))) {
    // If no good mappings were found, uncache the starting value of
    // m_deformation_gradient
    m_deformation_gradient = init_deformation_gradient;
    // m_N hasn't changed if tcost>max_cost
    // m_cost hasn't changed either
  }
  // m_N, m_deformation_gradient, and m_cost will describe the best mapping
  // encountered, even if nothing better than max_cost was encountered

  return *this;
}

/// Returns true if current N matrix is the canonical equivalent
bool LatticeMap::_check_canonical() const {
  // Purpose of jmin is to exclude (i,j)=(0,0) element
  // jmin is set to 0 at end of i=0 pass;
  Index jmin = 1;
  for (Index i = 0; i < m_parent_fsym_mats.size(); ++i, jmin = 0) {
    auto const &inv_parent_op = m_parent_fsym_mats[i];
    for (Index j = jmin; j < m_child_fsym_mats.size(); ++j) {
      auto const &child_op = m_child_fsym_mats[j];
      m_icache = child_op * inv_mat() * inv_parent_op;
      // Skip ops that transform matrix out of range; they won't be enumerated
      if (std::abs(m_icache(0, 0)) > m_range ||
          std::abs(m_icache(0, 1)) > m_range ||
          std::abs(m_icache(0, 2)) > m_range ||
          std::abs(m_icache(1, 0)) > m_range ||
          std::abs(m_icache(1, 1)) > m_range ||
          std::abs(m_icache(1, 2)) > m_range ||
          std::abs(m_icache(2, 0)) > m_range ||
          std::abs(m_icache(2, 1)) > m_range ||
          std::abs(m_icache(2, 2)) > m_range)
        continue;
      if (std::lexicographical_compare(m_icache.data(), m_icache.data() + 9,
                                       inv_mat().data(), inv_mat().data() + 9))
        return false;
    }
  }
  return true;
}

/// \brief Returns the volume-normalized strain cost, calculated to be
/// invariant to which structure is the child/parent.
///
/// Given a lattice mapping:
///     (L1 * T1) * N = V^{N} * Q^{N} * L2 * T2
///
/// Or equivalently,
///     V_reverse^{N} * Q_reverse^{N} * (L1 * T1) * N = L2 * T2
///
/// This calculates the cost as:
///
/// Child to parent deformation cost:
///       \tilde{V} = \frac{1}{\det{V}^{1/3}} V
///       \tilde{B} = \tilde{V} - I
///       cost = (1./3.)*\tr{\tilde{B}^{2}}
///
/// Parent to child deformation cost, using V.inverse() = U_reverse:
///       \tilde{U_reverse} = \frac{1}{\det{U_reverse}^{1/3}} U_reverse
///       \tilde{B_reverse} = \tilde{U_reverse} - I
///       cost = (1./3.)*\tr{\tilde{B_reverse}^{2}}
///
/// Direction invariant cost:
///       strain_cost = (1./2.)*(
///           (1./3.)*\tr{\tilde{B}^{2}} +
///           (1./3.)*\tr{\tilde{B_reverse}^{2}})
///
/// In the above, the `\tilde` indicates that the value is normalized to be
/// volume invariant.
///
/// \param deformation_gradient The deformation gradient, F or F_reverse. The
///     result is equivalent whether this is the parent to child deformation or
///     child to parent deformation.
///
double isotropic_strain_cost(Eigen::Matrix3d const &deformation_gradient) {
  // written using convention B = V - I of the mapping paper:
  Eigen::Matrix3d const &F_reverse = deformation_gradient;
  double vol_factor = pow(std::abs(F_reverse.determinant()), 1. / 3.);

  Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
  Eigen::Matrix3d U_reverse_normalized =
      strain::right_stretch_tensor(F_reverse / vol_factor);
  Eigen::Matrix3d V_normalized = U_reverse_normalized.inverse();

  // M.squaredNorm() = squared Frobenius norm of M = tr(M*M.transpose())
  return ((U_reverse_normalized - I).squaredNorm() +
          (V_normalized - I).squaredNorm()) /
         6.;
}

/// \brief Returns the symmetrized stretch
///
/// U_symmetrized = (1/NG) * sum_i (G_i * U * G_i.inverse()),
/// where NG: parent_point_group.size(),
/// G_i: elements of parent_point_group,
/// U = right stretch tensor of deformation_gradient
///
/// \param deformation_gradient The deformation gradient.
/// \param parent_point_group Parent point group. Use the point group of the
///     parent structure if mapping structures.
///
/// \return U_symmetrized
///
Eigen::Matrix3d symmetrized_stretch(Eigen::Matrix3d const &deformation_gradient,
                                    SymOpVector const &parent_point_group) {
  Eigen::Matrix3d const &F = deformation_gradient;
  Eigen::Matrix3d U = strain::right_stretch_tensor(F);

  Eigen::Matrix3d U_aggregate = Eigen::Matrix3d::Zero();
  for (auto const &op : parent_point_group) {
    U_aggregate += op.matrix * U * op.matrix.inverse();
  }
  return U_aggregate / double(parent_point_group.size());
}

/// \brief Returns the symmetry-breaking component of the stretch
///
/// U_symmetrized = (1/NG) * sum_i (G_i * U * G_i.inverse()),
/// where NG: parent_point_group.size(),
/// G_i: elements of parent_point_group,
/// U = right stretch tensor of deformation_gradient,
/// I = identity matrix
///
/// \param deformation_gradient The deformation gradient.
/// \param parent_point_group Parent point group. Use the point group of the
///     parent structure if mapping structures.
///
/// \return U - U_symmetrized + I
///
Eigen::Matrix3d symmetry_breaking_stretch(
    Eigen::Matrix3d const &deformation_gradient,
    SymOpVector const &parent_point_group) {
  Eigen::Matrix3d const &F = deformation_gradient;
  Eigen::Matrix3d U = strain::right_stretch_tensor(F);

  Eigen::Matrix3d U_aggregate = Eigen::Matrix3d::Zero();
  for (auto const &op : parent_point_group) {
    U_aggregate += op.matrix * U * op.matrix.inverse();
  }

  Eigen::Matrix3d U_sym = U_aggregate / double(parent_point_group.size());

  Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d U_sym_breaking = U - U_sym + I;
  return U_sym_breaking;
}

/// \brief Returns the symmetry-breaking strain cost, calculated to be
/// invariant to which structure is the child/parent.
///
/// The symmetry-breaking strain cost is calculated:
///
///     symmetry_breaking_strain_cost =
///         (1./2.)*(
///             (1./3.)*trace( (B_sym_break)^2 ) +
///             (1./3.)*trace( (B_reverse_sym_break)^2 ) )
///
/// Where:
///     B_sym_break = B - B_sym, the symmetry-breaking Biot strain, with
///         B = V - I as defined in mapping paper
///     B_sym = (1./N_G1) * sum_i ( G1_i * B * G1_i.transpose() ),
///     G1_i are point group operations of parent structure,
///     N_G1 is the total number of operations,
///     B_reverse_sym_break = B_reverse - B_reverse_sym is calculated similarly
///         to B_sym_break, but using B_reverse = U_reverse - I in place of B.
///
/// \param deformation_gradient The deformation gradient, F or F_reverse. The
///     result is equivalent whether this is the parent to child deformation or
///     child to parent deformation.
/// \param parent_point_group Parent point group. Use the point group of the
///     parent structure if mapping structures.
///
double symmetry_breaking_strain_cost(
    Eigen::Matrix3d const &deformation_gradient,
    SymOpVector const &parent_point_group) {
  // written using convention B = V - I of the mapping paper:

  // deformation_gradient = F_reverse (parent to child deformation)
  Eigen::Matrix3d const &F_reverse = deformation_gradient;
  Eigen::Matrix3d U_reverse_sym_breaking =
      symmetry_breaking_stretch(F_reverse, parent_point_group);
  Eigen::Matrix3d V_sym_breaking = U_reverse_sym_breaking.inverse();

  // M.squaredNorm() = squared Frobenius norm of M = tr(M*M.transpose())
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
  return ((V_sym_breaking - I).squaredNorm() +
          (U_reverse_sym_breaking - I).squaredNorm()) /
         6.;
}

}  // namespace xtal
}  // namespace CASM
