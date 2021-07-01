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
/// Note: This documentation uses the notation and conventions from the paper
/// "Comparing crystal structures with symmetry and geometry",
/// by John C. Thomas, Anirudh Raju Natarajan, Anton Van der Ven.
///
/// Goal:
/// - Find lattice mapping solutions, \f$(N, F^{N})\f$, that minimize a
///   cost, \f$\mathrm{cost}(F^{N})\f$, of straining a "child" lattice
///   (\f$L_2\f$) to a reference "parent" lattice (\f$L_1 * T_1\f$)
///   according to:
///   \f[
///       L_1 * T_1 * N = F^{N} * L_2.
///   \f]
///
/// \note Within the context of LatticeMap, \f$L_1\f$ and \f$T_1\f$ are
///   never separated and the LatticeMap documentation could be written using
///   only \f$L_1\f$ for the parent lattice. In practice, LatticeMap is used by
///   the structure mapping algorithms implemented by StrucMapper which choose
///   or search for \f$T_1\f$ in order to map a child lattice to a superlattice
///   of the parent lattice. The LatticeMap documentation is written using
///   \f$(L_1 * T_1)\f$ for the parent lattice to be clearer within the
///   context of the StrucMapper documentation.
///
/// Variables:
/// - \f$(L_1 * T_1)\f$: the reference parent lattice to be
///   mapped to, represented as a 3x3 matrix whose columns are the lattice
///   vectors.
/// - \f$L_2\f$: the child lattice to be mapped, represented
///   as a 3x3 matrix whose columns are the lattice vectors
/// - \f$N\f$: a unimodular matrix (integer valued, with \f$\det{N}=1\f$,
///   though represented with a floating point matrix for multiplication
///   purposes) that generates lattice vectors \f$(L_1 * T_1 * N)\f$ of
///   lattices that are equivalent to the input parent lattice \f$(L_1 *
///   T_1)\f$. Elsewhere, the superscript \f$N\f$ is used to indicate values
///   that depend on the choice of \f$N\f$.
/// - \f$F^{N}\f$: the deformation gradient of mapping from the child lattice
///   vectors to the parent lattice vectors.
///
///
/// ### Example: Construction, mapping to a BasicStructure ###
///
/// When mapping to a BasicStructure, the crystal point group should be used:
/// \code
/// // parent, use crystal point group
/// xtal::BasicStructure parent = ...;
/// double tol = TOL;
/// std::vector<xtal::SymOp> parent_factor_group =
///     xtal::make_factor_group(parent, tol);
/// std::vector<xtal::SymOp> parent_crystal_point_group =
///     xtal::make_crystal_point_group(parent_factor_group, tol);
///
/// // child, use identity_group
/// xtal::SimpleStructure unmapped_child = ...;
/// std::vector<xtal::SymOp> child_point_group = xtal::identity_group();
///
/// // range of elements of N matrix (controls number of potential mappings to
/// // be considered... larger is more, typically 1 is sufficient)
/// int unimodular_element_range = 1;
///
/// double max_lattice_cost = 1e20;
///
/// // Use isotropic_strain_cost
/// bool use_symmetry_breaking_strain_cost = false;
///
/// xtal::LatticeMap lattice_map{parent.lattice(),
///                              Lattice(unmapped_child.lat_column_mat)
///                              unimodular_element_range,
///                              parent_crystal_point_group,
///                              child_point_group,
///                              max_lattice_cost,
///                              use_symmetry_breaking_strain_cost};
/// \endcode
///
/// ### Example: Construction, mapping lattices without a basis ###
///
/// When mapping lattices without a basis, the lattice point groups should be
/// used:
/// \code
/// // parent, use lattice point group
/// Eigen::Matrix3d L1;
/// L1 << 3.233986860000, -1.616993430000, 0.000000000000,
///     0.000000000000, 2.800714770000, 0.000000000000,
///     0.000000000000, 0.000000000000, 10.337356680000;
///
/// xtal::Lattice L1_lattice(L1);
/// auto L1_point_group = xtal::make_point_group(L1_lattice);
///
/// // child, use lattice point group
/// Eigen::Matrix3d L2;
/// L2 << 3.269930775653, 0.000000000000, 0.000000000000,  //
///     -1.634965387827, 2.831843113861, 0.000000000000,   //
///     0.000000000000, 0.000000000000, 10.464806115486;   //
///
/// xtal::Lattice L2_lattice(L2);
/// auto L2_point_group = xtal::make_point_group(L2_lattice);
///
/// // range of elements of N matrix (controls number of potential mappings to
/// // be considered... larger is more, typically 1 is sufficient)
/// int unimodular_element_range = 1;
///
/// double max_lattice_cost = 1e20;
///
/// // Use isotropic_strain_cost
/// bool use_symmetry_breaking_strain_cost = false;
///
/// xtal::LatticeMap lattice_map{L1_lattice,
///                              L2_lattice,
///                              unimodular_element_range,
///                              L1_point_group,
///                              L2_point_group,
///                              max_lattice_cost,
///                              use_symmetry_breaking_strain_cost};
/// \endcode
///
/// ### Example: Find the optimal mapping ###
///
/// Use LatticeMap::best_strain_mapping to find the optimal mapping according
/// to the chosen cost method.
/// \code
/// // find optimal mapping
/// lattice_map.best_strain_mapping();
///
/// // print solution
/// std::cout << "N: \n" << lattice_map.matrixN() << std::endl;
/// std::cout << "F^{N}: \n" << lattice_map.deformation_gradient() << std::endl;
/// std::cout << "cost method: \n" << lattice_map.cost_method() << std::endl;
/// std::cout << "cost: \n" << lattice_map.strain_cost() << std::endl;
/// \endcode
///
/// ### Example: Find sub-optimal mappings ###
///
/// Find all mappings with cost better than some `max_cost` by using
/// LatticeMap::next_mapping_better_than:
/// \code
/// double max_cost = 0.4;
/// std::vector<Eigen::Matrix3d> N;
/// std::vector<Eigen::Matrix3d> F;
/// std::vector<double> cost;
/// while (lattice_map.next_mapping_better_than(max_cost)) {
///   N.push_back(lattice_map.matrixN());
///   F.push_back(lattice_map.deformation_gradient());
///   cost.push_back(lattice_map.strain_cost());
/// }
/// \endcode
///
/// ### Some useful relations and definitions ###
///
/// The deformation gradient can be decomposed into a stretch tensor (\f$V\f$,
/// the "left stretch tensor", or \f$U\f$ the "right stretch tensor") and an
/// isometry (distance-preserving transformation) matrix (\f$Q\f$),
/// \f[
///     F^{N} = V^{N} * Q^{N} = Q^{N} * U^{N}.
/// \f]
/// It can also be defined in the reverse direction (parent to child
/// deformation),
/// \f[
///     F_{reverse}^{N} * (L_1 * T_1) * N = L_2,
/// \f]
/// and decomposed into the reverse stretch tensor and isometry matrices
/// \f[
///     F_{reverse}^{N} = V_{reverse}^{N} * Q_{reverse}^{N} = Q_{reverse}^N *
///         U_{reverse}^N.
/// \f]
/// There are related via:
/// \f[
///     U_{reverse}^{N} = (V^{N})^{-1} \\
///     U^{N} = (V_{reverse}^{N})^{-1} \\
///     Q_{reverse}^{N} = (Q^{N})^{-1} = (Q^{N})^{\top}
/// \f]
/// The Biot strain tensor, \f$B^{N}\f$, and reverse are defined:
/// \f[
///     B^{N} = V^{N} - I \\
///     B_{reverse}^{N} = V_{reverse}^{N} - I
/// \f]
///
///
/// ### Strain cost ###
///
/// The strain cost may be calculated using one of:
/// - `"isotropic_strain_cost"`: (default, used if constructed with
///  `_symmetrize_strain_cost==false`, calculated by isotropic_strain_cost) a
///  volume-normalized strain cost, calculated to be invariant to which
///  structure is the child/parent. The use of \f$\tilde{}\f$ indicates that
///  the value is normalized to be volume invariant.
///  \f[
///       \mathrm{isotropic\ strain\ cost} = \frac{1}{2}*\left(
///           \frac{1}{3}*\mathrm{tr}\left(\tilde{B}^{2} \right) +
///           \frac{1}{3}*\mathrm{tr}\left(\tilde{B}_{reverse}^{2} \right)
///       \right) \\
///       \tilde{V} = \frac{1}{\det{V}^{1/3}} V \\
///       \tilde{B} = \tilde{V} - I
///  \f]
///   For reverse cost, use \f$V^{-1} = U_{reverse}\f$:
///  \f[
///       \tilde{U}_{reverse} = \frac{1}{\det{U_{reverse}}^{1/3}} U_{reverse} \\
///       \tilde{B}_{reverse} = \tilde{U}_{reverse} - I
///  \f]
/// - `"symmetry_breaking_strain_cost"`: (used if constructed with
///   `_symmetrize_strain_cost==true`, calculated by
///   symmetry_breaking_strain_cost) a strain cost, calculated with only the
///   the components of the strain that break the symmetry of the parent
///   structure.
///  \f[
///       \mathrm{symmetry\ breaking\ strain\ cost} =
///           \frac{1}{2}*\left(
///               \frac{1}{3}*\mathrm{tr}\left( B_{sym-break}^2 \right) +
///               \frac{1}{2}*\mathrm{tr}\left( B_{sym-break-reverse}^2 \right)
///           \right) \\
///       B_{sym-break} = B - B_{sym} \\
///       B_{sym} = \frac{1}{N_{G_1}} * \sum_i ( G_1(i) * B * G_1^{\top}(i) )
///  \f]
///   where \f$B_{sym-break}\f$, is the symmetry-breaking Biot strain,
///   \f$G_1(i)\f$ are point group operations of the parent structure, and
///   \f$N_{G_1}\f$ is the total number of operations. Similar relations hold
///   for calculating \f$B_{sym-break-reverse}\f$ from \f$B_{reverse}\f$.
///
/// ### Implementation note: Use of reduced cells ###
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
///
///     reduced_parent * N_reduced == F_reduced * reduced_child
///     parent * transformation_matrix_to_reduced_parent * N_reduced ==
///         F_reduced * child * transformation_matrix_to_reduced_child
///
/// Finally, this gives:
///
///     N == transformation_matrix_to_reduced_parent * N_reduced *
///         transformation_matrix_to_reduced_child.inverse(),
///     F == F_reduced

/// \brief LatticeMap constructor
///
/// \param _parent Reference lattice (\f$L_1 * T_1\f$)
/// \param _child Lattice to be mapped to _parent (\f$L_2\f$)
/// \param _range Determines range of \f$N\f$ matrices to be searched when
///     optimizing the lattice mapping. The absolute value of an element of N
///     is not allowed be be greater than `_range`. Typically 1 is a good enough
///     choice.
/// \param _parent_point_group Point group of the parent (i.e. crystal point
///     group), used to identify canonical N matrices and reduce the number of
///     operations. Also used to calculate the symmetry-breaking strain cost if
///     `_symmetrize_strain_cost==true`.
/// \param _child_point_group Point group of the child (structure), used to
///     identify canonical \f$N\f$ matrices and reduce the number of operations.
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
      m_has_current_solution(false),
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

void LatticeMap::_reset(double _better_than) {
  m_currmat = 0;

  // From relation F * parent * inv_mat.inverse() = child
  m_deformation_gradient =
      m_reduced_child * inv_mat().cast<double>() *
      m_reduced_parent.inverse();  // -> _deformation_gradient

  double tcost = _calc_strain_cost(m_deformation_gradient);

  // Initialize to first valid mapping
  if (tcost <= _better_than && _check_canonical()) {
    m_has_current_solution = true;
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

  m_has_current_solution = true;
  m_cost = best_cost;
  return *this;
}

/// Calculate the strain cost of a deformation
///
/// \param deformation_gradient Parent to child deformation gradient,
///    \f$F_{reverse}^{N}\f$, where \f$F_{reverse}^{N} * (L_1 * T_1) * N =
///    L_2\f$.
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

/// \brief Iterate until the next solution \f$(N, F^{N})\f$ with lattice
/// mapping score less than `max_cost` is found.
///
/// The current solution is:
/// - \f$F_{reverse}^{N}\f$ = this->deformation_cost()
/// - \f$N\f$ = this->matrixN()
///
/// The current solution exists if
/// \code
///     static_cast<bool>(*this) == true.
/// \endcode
/// Otherwise, iteration is complete.
///
const LatticeMap &LatticeMap::next_mapping_better_than(double max_cost) const {
  m_has_current_solution = false;
  m_cost = 1e20;
  return _next_mapping_better_than(max_cost);
}

/// \brief Iterate until the next solution \f$(N, F^{N})\f$ with lattice mapping
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
      m_has_current_solution = true;
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
/// \f[
///     (L_1 * T_1) * N = V^{N} * Q^{N} * L_2 * T_2
/// \f]
///
/// Or equivalently,
/// \f[
///     V_{reverse}^{N} * Q_{reverse}^{N} * (L_1 * T_1) * N = L_2 * T_2
/// \f]
///
/// This calculates the cost as:
///
/// Child to parent deformation cost:
/// \f[
///       \tilde{V} = \frac{1}{\det{V}^{1/3}} V
///       \tilde{B} = \tilde{V} - I
///       cost = (1./3.)*\mathmr{Tr}{\tilde{B}^{2}}
/// \f]
///
/// Parent to child deformation cost, using \f$V^{-1} = U_{reverse}\f$:
/// \f[
///       \tilde{U}_{reverse} = \frac{1}{\det{U_{reverse}}^{1/3}} U_{reverse}
///       \tilde{B}_{reverse} = \tilde{U}_{reverse} - I
///       cost = (1./3.)*\mathrm{tr}(\tilde{B}_{reverse}^{2})
/// \f]
///
/// Direction invariant cost:
/// \f[
///       strain_cost = (1./2.)*(
///           (1./3.)*\mathrm{tr}(\tilde{B}^{2}) +
///           (1./3.)*\mathrm{tr}(\tilde{B}_{reverse}^{2}))
/// \f]
///
/// In the above, the \f$\tilde\f$ indicates that the value is normalized to be
/// volume invariant.
///
/// \param deformation_gradient The deformation gradient, \f$F or F_reverse\f$.
///     The result is equivalent whether this is the parent to child
///     deformation or child to parent deformation.
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
/// \param deformation_gradient The deformation gradient.
/// \param parent_point_group Parent point group. Use the point group of the
///     parent structure if mapping structures.
///
/// \return \f$U_{symmetrized}\f$, where:
///   - \f$U_{symmetrized} = (1/N_{G_1}) * \sum_i (G_1(i) * U * G_1(i)^{-1})\f$,
///   - \f$U\f$: right stretch tensor of the deformation gradient
///   - \f$N_{G_1}\f$: Parent point group size
///   - \f$G_1(i)\f$: Element \f$i\f$ of the parent point group
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
/// \param deformation_gradient The deformation gradient.
/// \param parent_point_group Parent point group. Use the point group of the
///     parent structure if mapping structures.
///
/// \return \f$U - U_{symmetrized} + I\f$, where:
///   - \f$U\f$: Right stretch tensor of the deformation gradient
///   - \f$U_{symmetrized} = (1/N_{G_1}) * \sum_i (G_1(i) * U * G_1(i)^{-1})\f$,
///   - \f$N_{G_1}\f$: Parent point group size
///   - \f$G_1(i)\f$: Element \f$i\f$ of the parent point group
///   - \f$I\f$: Identity matrix
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
/// \f[
///     \frac{1}{2}*\left(
///         \frac{1}{3}*\mathrm{tr}\left( B_{sym-break}^2 \right) +
///         \frac{1}{3}*\mathrm{tr}\left( B_{reverse-sym-break}^2 \right)
///     \right)
/// \f]
///
/// Where:
/// - \f$B = V - I\f$ Biot strain, as defined in the mapping paper
/// - \f$B_{sym} = (1./N_{G_1}) * \sum_i ( G_1(i) * B * G_1(i)^{-1} )\f$,
/// - \f$B_{sym-break} = B - B_{sym}\f$, the symmetry-breaking Biot strain
/// - \f$G_1(i)\f$: Element \f$i\f$ of the parent point group
/// - \f$N_{G_1}\f$: Parent point group size
/// - \f$B_{reverse-sym-break} = B_{reverse} - B_{reverse-sym}\f$ is calculated
///   similarly to \f$B_{sym-break}\f$, but using \f$B_{reverse} = U_{reverse}
///   - I\f$ in place of \f$B\f$.
///
/// \param deformation_gradient The deformation gradient, \f$F\f$ or
///     \f$F_{reverse}\f$. The result is equivalent whether this is the parent
///     to child deformation or child to parent deformation.
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
