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

/// Return name of method used by `strain_cost(Eigen::Matrix3d const &)`
///
/// \returns One of "isotropic_strain_cost" or "anisotropic_strain_cost",
/// depending on `strain_gram_mat` constructor argument.
///
/// Note:
/// - Currently to use the symmetrized strain cost the method
///   `strain_cost(Eigen::Matrix3d const &, SymOpVector const &)` must be used
///   explicitly.
///
/// TODO:
/// - make StrainCostCalculator a unary functor with method chosen at
///   construction time, and then this could return "symmetrized_strain_cost"
///   also.
std::string StrainCostCalculator::cost_method() const {
  if (m_sym_cost) {
    return "anisotropic_strain_cost";
  } else {
    return "strain_cost";
  }
}

//*******************************************************************************************

/// LatticeMap constructor
///
/// \param _parent Reference lattice (L1 * T1)
/// \param _child Lattice to be mapped to _parent (L2)
/// \param num_atoms Unused. (TODO: remove)
/// \param _range Determines range of N matrices to be searched when optimizing
/// the lattice mapping. The absolute value of an element of N is not allowed be
/// be greater than `_range`. Typically 1 is a good enough choice. \param
/// _parent_point_group Point group of the parent (i.e. crystal point group),
/// used to identify canonical N matrices and reduce the number of operations.
/// Also used to calculate the symmetry-breaking strain cost if
/// `_symmetrize_strain_cost==true`. \param _child_point_group Point group of
/// the child (structure), used to identify canonical N matrices and reduce the
/// number of operations. (TODO: seems to typically just be
/// {SymOp::identity()}?) \param strain_gram_mat Gram matrix, may be used to
/// calculate an anisotropic strain cost. Use default value
/// Eigen::MatrixXd::Identity(9, 9) otherwise. \param _init_better_than
/// Initializes the LatticeMap to the first mapping with a cost less than
/// `_init_better_than + _xtal_tol`. Use the default value (1e20) to initialize
/// to the first valid mapping. \param _symmetrize_strain_cost Boolean flag
/// (default=false), which if true indicates that the `strain-breaking strain
/// cost` should be used to score lattice mappings. \param _xtal_tol Tolerance
/// used for crystallography comparisons
LatticeMap::LatticeMap(const Lattice &_parent, const Lattice &_child,
                       Index num_atoms, int _range,
                       SymOpVector const &_parent_point_group,
                       SymOpVector const &_child_point_group,
                       Eigen::Ref<const Eigen::MatrixXd> const &strain_gram_mat,
                       double _init_better_than /* = 1e20 */,
                       bool _symmetrize_strain_cost, double _xtal_tol)
    : m_parent(_parent.lat_column_mat()),
      m_child(_child.lat_column_mat()),
      m_calc(strain_gram_mat),
      m_vol_factor(pow(std::abs(volume(_child) / volume(_parent)), 1. / 3.)),
      m_range(_range),
      m_cost(1e20),
      m_currmat(0),
      m_symmetrize_strain_cost(_symmetrize_strain_cost),
      m_xtal_tol(_xtal_tol) {
  Lattice reduced_parent = _parent.reduced_cell();
  m_reduced_parent = reduced_parent.lat_column_mat();

  Lattice reduced_child = _child.reduced_cell();
  m_reduced_child = reduced_child.lat_column_mat();

  m_transformation_matrix_to_reduced_parent =
      _parent.inv_lat_column_mat() * m_reduced_parent;
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
                       Index _num_atoms, int _range,
                       SymOpVector const &_parent_point_group,
                       SymOpVector const &_child_point_group,
                       Eigen::Ref<const Eigen::MatrixXd> const &strain_gram_mat,
                       double _init_better_than /* = 1e20 */,
                       bool _symmetrize_strain_cost, double _xtal_tol)
    : LatticeMap(Lattice(_parent), Lattice(_child), _num_atoms, _range,
                 _parent_point_group, _child_point_group, strain_gram_mat,
                 _init_better_than, _symmetrize_strain_cost, _xtal_tol) {}

void LatticeMap::_reset(double _better_than) {
  m_currmat = 0;

  // From relation F * parent * inv_mat.inverse() = child
  m_deformation_gradient =
      m_reduced_child * inv_mat().cast<double>() *
      m_reduced_parent.inverse();  // -> _deformation_gradient

  double tcost = calc_strain_cost(m_deformation_gradient);

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

  double best_cost = calc_strain_cost(m_deformation_gradient);

  while (next_mapping_better_than(best_cost).strain_cost() < best_cost) {
    best_cost = strain_cost();
  }

  m_cost = best_cost;
  return *this;
}

/// \brief Use the current strain cost method (determined by constructor
/// arguments) to calculate the strain cost of a deformation
///
/// \param deformation_gradient Parent to child deformation gradient,
///    F_rev^{N}, where F_rev^{N} * (L1 * T1) * N == L2.
///
double LatticeMap::calc_strain_cost(
    const Eigen::Matrix3d &deformation_gradient) const {
  if (symmetrize_strain_cost())
    return m_calc.strain_cost(deformation_gradient, m_parent_point_group);
  else
    return m_calc.strain_cost(deformation_gradient, m_vol_factor);
}

/// The name of the method used to calculate the lattice deformation cost
///
/// \returns "strain_cost", "anisotropic_strain_cost", or
///   "symmetry_breaking_strain_cost", as determined by constructor arguments
///
std::string LatticeMap::cost_method() const {
  if (symmetrize_strain_cost()) {
    return "symmetry_breaking_strain_cost";
  } else {
    return m_calc.cost_method();
  }
}

/// \brief Iterate until the next solution (N, F^{N}) with lattice mapping
/// score less than `max_cost` is found.
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
    tcost = calc_strain_cost(m_deformation_gradient);
    if (std::abs(tcost) < (std::abs(max_cost) + std::abs(xtal_tol()))) {
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
  if (!(std::abs(tcost) < (std::abs(max_cost) + std::abs(xtal_tol())))) {
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
}  // namespace xtal
}  // namespace CASM
