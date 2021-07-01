#include "casm/crystallography/StrucMapping.hh"

#include "casm/container/algorithm.hh"
#include "casm/crystallography/CanonicalForm.hh"
#include "casm/crystallography/Coordinate.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/LatticeMap.hh"
#include "casm/crystallography/Strain.hh"
#include "casm/crystallography/StrucMapCalculatorInterface.hh"
#include "casm/crystallography/SuperlatticeEnumerator.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/crystallography/io/VaspIO.hh"
#include "casm/external/Eigen/src/Core/Map.h"
#include "casm/external/Eigen/src/Core/PermutationMatrix.h"
#include "casm/external/Eigen/src/Core/util/Constants.h"
#include "casm/external/Eigen/src/Core/util/Meta.h"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {
namespace xtal {

namespace Local {

static bool lex_lt(Eigen::Matrix<long, 3, 3> const &A,
                   Eigen::Matrix<long, 3, 3> const &B) {
  return std::lexicographical_compare(A.data(), A.data() + 9, B.data(),
                                      B.data() + 9);
}
}  // namespace Local
//*******************************************************************************************

namespace StrucMapping {
double atomic_cost_child(const MappingNode &mapped_result, Index Nsites) {
  Nsites = max(Nsites, Index(1));
  // mean square displacement distance in deformed coordinate system
  double atomic_vol =
      mapped_result.lattice_node.parent.superlattice().volume() /
      double(Nsites) / mapped_result.lattice_node.stretch.determinant();
  return pow(3. * abs(atomic_vol) / (4. * M_PI), -2. / 3.) *
         (mapped_result.lattice_node.stretch.inverse() *
          mapped_result.atom_displacement)
             .squaredNorm() /
         double(Nsites);
}
//*******************************************************************************************

double atomic_cost_parent(const MappingNode &mapped_result, Index Nsites) {
  Nsites = max(Nsites, Index(1));
  // mean square displacement distance in deformed coordinate system
  double atomic_vol =
      mapped_result.lattice_node.parent.superlattice().volume() /
      double(Nsites);
  return pow(3. * abs(atomic_vol) / (4. * M_PI), -2. / 3.) *
         (mapped_result.atom_displacement).squaredNorm() / double(Nsites);
}

//*******************************************************************************************

double atomic_cost(const MappingNode &mapped_result, Index Nsites) {
  // mean square displacement distance in deformed coordinate system
  return (atomic_cost_child(mapped_result, Nsites) +
          atomic_cost_parent(mapped_result, Nsites)) /
         2.;
}

//*******************************************************************************************
double atomic_cost(
    const MappingNode &basic_mapping_node, SymOpVector &factor_group,
    const std::vector<Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic,
                                               Index>> &permutation_group,
    Index Nsites) {
  const auto disp_matrix = basic_mapping_node.atom_displacement;
  Eigen::MatrixXd symmetry_preserving_displacement =
      Eigen::MatrixXd::Zero(disp_matrix.rows(), disp_matrix.cols());
  for (int i = 0; i < factor_group.size(); ++i) {
    auto transformed_disp = factor_group[i].matrix * disp_matrix;
    Eigen::MatrixXd transformed_and_permuted_disp =
        transformed_disp * permutation_group[i];
    symmetry_preserving_displacement += transformed_and_permuted_disp;
  }
  symmetry_preserving_displacement =
      symmetry_preserving_displacement / factor_group.size();
  auto new_report = basic_mapping_node;
  new_report.atom_displacement = disp_matrix - symmetry_preserving_displacement;
  return atomic_cost_parent(new_report, Nsites);
};

}  // namespace StrucMapping

//*******************************************************************************************
namespace Local {
// Local helper function for StrucMapper::k_best_maps_better_than
template <typename OutputIterator>
static bool initial_atomic_maps(SimpleStructure child_struc,
                                MappingNode const &seed,
                                StrucMapCalculatorInterface const &calculator,
                                double max_cost,
                                bool const &symmetrize_atomic_cost,
                                OutputIterator it) {
  // derotate first
  child_struc.rotate_coords(seed.isometry());

  // Then undeform by inverse of right stretch
  child_struc.deform_coords(seed.stretch());

  // We want to get rid of translations.
  // define translation such that:
  //    IDEAL = RELAXED + translation
  // and use it when calculating cost matrix

  for (Eigen::Vector3d const &translation :
       calculator.translations(seed, child_struc)) {
    MappingNode node = seed;
    node.atomic_node.translation = translation;
    if (!calculator.populate_cost_mat(node, child_struc)) {
      // Indicates that structure is incompatible with supercell, regardless of
      // translation so return false
      return false;
    }

    // The mapping routine is called here
    node.calc();

    // if assignment is smaller than child_struc.basis().size(), then
    // child_struc is incompattible with supercell (assignment.size()==0 if the
    // hungarian routine detects an incompatibility, regardless of translation)
    if (!node.is_viable) {
      return false;
    }
    // Now we are filling up displacements
    calculator.finalize(node, child_struc, symmetrize_atomic_cost);

    if (node.cost < max_cost) {
      *it = node;
    }
  }

  return true;
}

//*******************************************************************************************
// Local helper function for StrucMapper::k_best_maps_better_than
template <typename OutputIterator>
static void partition_node(MappingNode const &_node,
                           StrucMapCalculatorInterface const &_calculator,
                           SimpleStructure child_struc,
                           bool const &symmetrize_atomic_cost,
                           OutputIterator it) {
  // derotate first
  child_struc.rotate_coords(_node.isometry());

  // Then undeform by inverse of right stretch
  child_struc.deform_coords(_node.stretch());

  Index cN = _node.lattice_node.child.size() *
             _calculator.struc_info(child_struc).size();

  // We will increment from i=0 to number of rows in cost matrix of '_node'
  // For each (i,j) assignment of '_node', we will will spawn a new node, where
  // the (i,j) assignment is forced OFF and all (k,l) assignments with k<i will
  // be forced ON. This involves
  //   (1) setting the cost of element (i,j) to infinity (forccing it off)
  //   (2) striking the first 'k' rows from the staring cost matrix, and all
  //   columns corresponding to
  //       their atomic assignments.
  //   (3) recording the (i,j) pairs that have been turned on
  //   (4) resolving the reduced assignment problem with the new cost matrix
  Index old_j, new_j, deleted_j;
  Index n = _node.atomic_node.assignment.size();
  MappingNode t1(_node), t2(_node);

  _node.is_partitioned = true;
  t1.is_partitioned = false;
  t2.is_partitioned = false;

  MappingNode *p1 = &t1;
  MappingNode *p2 = &t2;
  for (Index m = n; 1 < m && (p1->is_viable || p2->is_viable); --m) {
    AssignmentNode &n1(p1->atomic_node);
    AssignmentNode &n2(p2->atomic_node);
    // clear assignment and cost_mat
    n2.assignment.clear();
    n2.irow.clear();
    n2.icol.clear();
    n2.cost_mat.resize(m - 1, m - 1);
    n2.forced_on = n1.forced_on;

    // We are forcing on the first site assignment in t1: i.e.,  [0,deleted_j]
    // this involves striking row 0 and col deleted_j from t2's cost_mat
    deleted_j = n1.assignment[0];

    // [0,deleted_j] have local context only. We store the forced assignment in
    // forced_on using the original indexing from n1
    n2.forced_on.emplace(n1.irow[0], n1.icol[deleted_j]);

    // Strike row 0 and col deleted_j to form new AssignmentNode for t2
    // (i,j) indexes the starting cost_mat, (i-1, new_j) indexes the resulting
    // cost_mat
    n2.irow = std::vector<Index>(++n1.irow.begin(), n1.irow.end());
    // We will also store an updated assignment vector in t2, which will be
    // used to construct next node of partition
    n2.assignment =
        std::vector<Index>(++n1.assignment.begin(), n1.assignment.end());
    for (old_j = 0, new_j = 0; old_j < m; ++old_j, ++new_j) {
      if (old_j == deleted_j) {
        --new_j;
        continue;
      }
      n2.icol.push_back(n1.icol[old_j]);

      // We will also store an updated assignment vector in t2, which will be
      // used to construct next node of partition
      if (n2.assignment[new_j] > deleted_j) n2.assignment[new_j]--;

      // Fill col new_j of t2's cost mat
      for (Index i = 1; i < m; ++i)
        n2.cost_mat(i - 1, new_j) = n1.cost_mat(i, old_j);
    }
    // t2 properly initialized; we can now force OFF [0,deleted_j] in t1, and
    // add it to node list
    n1.cost_mat(0, deleted_j) = StrucMapping::big_inf();
    // IMPORTANT: If n1.icol[deleted_j]=cN, it is a virtual vacancy.
    // If we exclude a single virtual vacancy from occupying this parent site,
    // the marginal cost of assigning a different virtual vacancy to the same
    // site is zero, and the mapping is equivalent. So, we need to exclude ALL
    // virtual vacancies from occupying this parent site:
    if (n1.icol[deleted_j] >= cN) {
      for (old_j = n1.icol.size() - 1; old_j >= 0 && n1.icol[old_j] >= cN;
           --old_j) {
        n1.cost_mat(0, old_j) = StrucMapping::big_inf();
      }
    }
    n1.assignment.clear();
    p1->is_viable = true;
    p1->calc();
    if (p1->is_viable) {
      // even if p1 is unviable, p2 may still be viable, so we continue
      _calculator.finalize(*p1, child_struc, symmetrize_atomic_cost);
      it = *p1;
    }
    std::swap(p1, p2);
  }
}
}  // namespace Local

namespace {
void check_equal(Eigen::MatrixXd const &A, Eigen::MatrixXd const &B,
                 std::string message) {
  if (!almost_equal(A, B)) {
    std::cout << "A: \n" << A << std::endl;
    std::cout << "B: \n" << B << std::endl << std::endl;
    throw std::runtime_error(message);
  }
}
}  // namespace

/// \struct LatticeNode
/// \brief Data structure describing a lattice mapping relationship
///
/// A general map for a child structure onto a parent structure may require
/// forming a supercell of the parent structure (most commonly) and/or of the
/// child structure. As such, the LatticeNode is specified in terms of
/// superlattices of both the parent and the child, as well as deformation and
/// rotation information sufficient to fully define the lattice mapping
/// transformation.
///
/// ### LatticeNode definitions
///
/// LatticeNode encodes the general lattice mapping relationship:
///
/// \f[
///     L_1 * T_1 * N = V^{N} * Q^{N} * L_2 * T_2
/// \f]
///
/// - \f$L_1\f$: The parent (reference structure) lattice vectors, as a column
///   vector matrix. It is equal to
///
///       parent.prim_lattice().lat_column_mat()
///
/// - \f$T_1\f$: An integer transformation matrix to a parent superlattice.
/// - \f$N\f$: A unimodular matrix (integer matrix, with determinant 1)
///   transforms parent superlattice vectors to create an equivalent
///   superlattice. Elsewhere, the superscript \f$N\f$ is used to indicate
///   values that depend on the choice of N.
/// - \f$L_2\f$: The child lattice (lattice of structure to be mapped), as a
///   column vector matrix.
/// - \f$T_2\f$: An integer transformation matrix to a child superlattice,
///   \f$(L_2 * T_2)\f$. When the parent lattice is the lattice of the
///   primitive structure, then it will be the case that \f$T_2 = I\f$.
/// - \f$V^{N}\f$: The stretch matrix, a symmetric matrix that describes the
///   deformation that maps a de-rotated child superlattice to a parent
///   superlattice. It is equal to:
///
///       stretch
///
/// - \f$Q^{N}\f$: The isometry (distance-preserving transformation) matrix,
///   which describes the transformation that de-rotates (de-reflects, etc.) a
///   superlattice of child. A property of \f$Q^{N}\f$ is that \f$(Q^{N})^{-1}
///   == (Q^{N})^{T}\f$. It is equal to:
///
///       isometry
///
///
/// Thus:
/// - \f$L_1 * T_1 * N\f$: A parent supercell lattice, as a column vector
///   matrix. It is equal to:
///
///       parent.superlattice().lat_column_mat()
///
/// - \f$V^{N} * Q^{N} * L_2\f$: The mapped and un-deformed child prim lattice,
///   as a column vector matrix. It is equal to:
///
///       child.prim_lattice().lat_column_mat()
///
/// - \f$V^{N} * Q^{N} * L_2 * T_2\f$: The mapped and un-deformed child
///   supercell lattice, as a column vector matrix. It is equal to:
///
///       child.superlattice().lat_column_mat()
///
///
/// ### Deformation gradient
///
/// The deformation gradient describes the transformation that maps two sets of
/// lattice vectors. The deformation gradient for transformation of a
/// superlattice of the child lattice to a superlattice of the parent lattice
/// can be written:
/// \f[
///     F_{child \to parent} = F = V * Q = Q * U
/// \f]
/// The reverse transformation (parent to child transformation) can be written:
/// \f[
///     F_{parent \to child} = F_{reverse} = V_{reverse} *
///     Q_{reverse} = Q_{reverse} * U_{reverse}
/// \f]
/// where \f$Q^{-1} = Q^{\top}\f$, using \f$^{\top}\f$ to indicate matrix
/// transpose, and \f$V\f$ and \f$U\f$ are symmetric stretch tensors (left and
/// right, respectively). Similar relations hold for the "reverse"
/// transformations. Then
/// \f[
///     F^{\top} * F = U^{\top} * Q^{-1} * Q * U = U^{2}
/// \f]
/// Thus the \link polar_decomposition polar decomposition\endlink of \f$F\f$
/// can be used to obtain \f$U\f$:
///
///     U = polar_decomposition(F)
///     U_reverse = polar_decomposition(F_reverse)
///
/// Relations between "forward" and the "reverse" definitions:
/// \f[
///     F_{reverse} = F^{-1} \\
///     Q_{reverse} * U_{reverse} = Q^{-1} * V^{-1}  \\
///     Q_{reverse} = Q^{-1} = Q^{\top}  \\
///     V = U_{reverse}^{-1}  \\
///     Q = (F_{reverse} * V)^{\top}
/// \f]
///
///
/// \note Elsewhere in CASM the strain associated with a configuration, whether
/// as a DoF or in MappedProperties, is defined using the "reverse" sense, as a
/// transformation of parent (supercell of prim) to child:
/// \f[
///     F_{parent \to child} * L_1 * T_1 * N = L_2 * T_2
/// \f]
/// For example, the \f$U\f$ of `"Ustrain"`, is \f$U_{reverse} = $U_{parent \to
/// child}\f$.
///
///
/// ### Relating LatticeNode and LatticeMap definitions:
///
/// The LatticeMap class performs searches over possible lattice mappings in
/// order to find low deformation cost mappings. It generates solutions \f$(N,
/// F_{reverse}^{N})\f$ that satisfy:
/// \f[
///     F_{reverse}^{N} * L_1 * T_1 * N = L_2 * T_2
/// \f]
/// - \f$F_{reverse}^{N}\f$: \code lattice_map.deformation_gradient() \endcode
/// - \f$L_1 * T_1\f$: \code lattice_map.parent_matrix() \endcode
/// - \f$N\f$: \code lattice_map.matrixN() \endcode
/// - \f$L_2 * T_2\f$: \code lattice_map.child_matrix() \endcode
///
/// So relating xtal::LatticeMap definitions and xtal::LatticeNode defintions:
///
///     lattice_node.stretch = polar_decomposition(
///                                lattice_map.deformation_gradient()).inverse();
///     lattice_node.isometry = (lattice_map.deformation_gradient() *
///                              lattice_node.stretch).transpose();
///

/// Construct LatticeNode, setting all members directly
LatticeNode::LatticeNode(Superlattice _parent, Superlattice _child,
                         Eigen::Matrix3d _stretch, Eigen::Matrix3d _isometry,
                         double _cost, std::string _cost_method)
    : stretch(_stretch),
      isometry(_isometry),
      parent(_parent),
      child(_child),
      cost(_cost),
      cost_method(_cost_method) {
  check_equal(
      parent.superlattice().lat_column_mat(),
      child.superlattice().lat_column_mat(),
      "LatticeNode constructor error: _parent.superlattice().lat_column_mat() "
      "!= _child.superlattice().lat_column_mat()");
}

/// \brief Construct a LatticeNode by calculating the deformation tensor that
/// maps a particular child superlattice to a particular parent superlattice
/// [deprecated]
///
/// \param parent_prim primitive lattice being mapped to (\f$L_1\f$)
/// \param parent_scel exact integral multiple of parent_prim (\f$L_1 * T_1 *
///     N\f$)
/// \param unmapped_child_prim primitive lattice being mapped (\f$L_2\f$)
/// \param unmapped_child_scel exact integral multiple of child_prim (\f$L_2 *
///     T_2\f$)
/// \param child_N_atom is number of sites in the child (Not used)
/// \param _cost is used to specify mapping cost (in default case -- big_inf()
/// -- cost will be calculated from scratch)
///
/// Note: This method is deprecated. Prefer using \ref make_lattice_node_1.
LatticeNode::LatticeNode(Lattice const &parent_prim, Lattice const &parent_scel,
                         Lattice const &unmapped_child_prim,
                         Lattice const &unmapped_child_scel, Index child_N_atom,
                         double _cost /*=StrucMapping::big_inf()*/)
    : parent(parent_prim, parent_scel),
      // Transform child_prim lattice to its idealized state using same
      // F.inverse as below, but inline:
      child(Lattice((parent_scel.lat_column_mat() *
                     unmapped_child_scel.inv_lat_column_mat()) *
                        unmapped_child_prim.lat_column_mat(),
                    parent_prim.tol()),
            parent_scel),
      cost(_cost) {
  // see LatticeNode class documentation for more on relations

  // parent_prim = L1
  // parent_scel = L1 * T1 * N
  // child_prim = L2
  // child_scel = L2 * T2
  // F_reverse * L1 * T1 * N = L2 * T2
  Eigen::Matrix3d F_reverse =
      unmapped_child_scel.lat_column_mat() * parent_scel.inv_lat_column_mat();

  // V = U_reverse.inverse()
  stretch = strain::right_stretch_tensor(F_reverse).inverse();

  // Q = (F_reverse * V).transpose()
  isometry = (F_reverse * stretch).transpose();

  if (StrucMapping::is_inf(cost)) {
    cost = isotropic_strain_cost(stretch);
    cost_method = "isotropic_strain_cost";
  } else {
    cost_method = "unknown";
  }

  check_equal(parent.superlattice().lat_column_mat(),
              stretch * isometry * unmapped_child_scel.lat_column_mat(),
              "LatticeNode constructor error: "
              "parent.superlattice().lat_column_mat() != "
              "stretch * isometry * unmapped_child_scel.lat_column_mat()");
}

/// \brief Construct a LatticeNode using the mapping calculated by LatticeMap
/// [deprecated]
///
/// \param lattice_map The lattice mapping is used to specify
/// the parent superlattice and the current solution specifies the deformation
/// gradient and choice of lattice vectors that map a supercell of the unmapped
/// child to a supercell of the parent. Specifically:
/// - \f$ L_1 * T_1 \f$ = `lattice_map.parent_matrix()`
/// - \f$F_{parent \to child}\f$ = `lattice_map.deformation_gradient()`
/// - \f$N\f$ = `lattice_map.matrixN`
/// \param parent_prim primitive lattice being mapped to (\f$L_1\f$)
/// \param unmapped_child_prim primitive lattice being mapped (\f$L_2\f$)
///
/// Note:
/// - The lattice deformation cost is calculated using the method specified by
///   `lattice_map`.
///
/// Note: This method is deprecated. Prefer using \ref make_lattice_node_2.
LatticeNode::LatticeNode(LatticeMap const &lattice_map,
                         Lattice const &parent_prim,
                         Lattice const &unmapped_child_prim)
    :  // see LatticeNode class documentation for more on relations
       // V = U_reverse.inverse()
      stretch(
          polar_decomposition(lattice_map.deformation_gradient()).inverse()),
      // Q = (F_reverse * V).transpose()
      isometry((lattice_map.deformation_gradient() * stretch).transpose()),
      // parent.prim_lattice() = L1
      // lattice_map.parent_matrix() = L1 * T1
      // parent.superlattice() = L1 * T1 * N
      parent(parent_prim,
             Lattice(lattice_map.parent_matrix() * lattice_map.matrixN(),
                     parent_prim.tol())),
      // (mapped) child.prim_lattice() = F * L2
      // (mapped) child.superlattice() = parent.superlattice()
      child(Lattice(lattice_map.deformation_gradient().inverse() *
                        unmapped_child_prim.lat_column_mat(),
                    parent_prim.tol()),
            parent.superlattice()),
      cost(lattice_map.strain_cost()),
      cost_method(lattice_map.cost_method()) {
  check_equal(
      parent.superlattice().lat_column_mat(),
      child.superlattice().lat_column_mat(),
      "LatticeNode constructor error: parent.superlattice().lat_column_mat() "
      "!= child.superlattice().lat_column_mat()");
  check_equal(
      lattice_map.deformation_gradient().inverse(), stretch * isometry,
      "LatticeNode constructor error: "
      "lattice_map.deformation_gradient().inverse() != stretch * isometry");
}

/// \brief Construct a LatticeNode by calculating the deformation tensor that
/// maps a particular child superlattice to a particular parent superlattice
///
/// \param parent_prim primitive lattice being mapped to (\f$L_1\f$)
/// \param parent_scel exact integral multiple of parent_prim (\f$L_1 * T_1 *
///     N\f$)
/// \param unmapped_child_prim primitive lattice being mapped (\f$L_2\f$)
/// \param unmapped_child_scel exact integral multiple of child_prim (\f$L_2 *
///     T_2\f$)
///
/// Note:
/// - In result: `parent_scel = stretch * isometry * unmapped_child_scel'
/// - The lattice deformation cost is calculated using
///   `isotropic_strain_cost(stretch)`
///
/// \anchor make_lattice_node_1
/// \relates LatticeNode
///
LatticeNode make_lattice_node(Lattice const &parent_prim,
                              Lattice const &parent_scel,
                              Lattice const &unmapped_child_prim,
                              Lattice const &unmapped_child_scel) {
  // see LatticeNode class documentation for more on relations

  // parent_prim = L1
  // parent_scel = L1 * T1 * N
  // child_prim = L2
  // child_scel = L2 * T2

  // parent.prim_lattice() = L1
  // parent.superlattice() = L1 * T1 * N
  Superlattice parent{parent_prim, parent_scel};

  // F_reverse * L1 * T1 * N = L2 * T2
  Eigen::Matrix3d F_reverse =
      unmapped_child_scel.lat_column_mat() * parent_scel.inv_lat_column_mat();

  // L1 * T1 * N = F * L2 * T2, F = F_reverse.inverse()
  Eigen::Matrix3d F =
      parent_scel.lat_column_mat() * unmapped_child_scel.inv_lat_column_mat();

  // mapped_child_prim_lattice = F * L2
  Lattice mapped_child_prim_lattice{F * unmapped_child_prim.lat_column_mat(),
                                    parent_prim.tol()};

  // (mapped) child.prim_lattice() = F * L2
  // (mapped) child.superlattice() = parent.superlattice()
  Superlattice mapped_child{mapped_child_prim_lattice, parent_scel};

  // V = U_reverse.inverse()
  Eigen::Matrix3d stretch = strain::right_stretch_tensor(F_reverse).inverse();

  // Q = (F_reverse * V).transpose()
  Eigen::Matrix3d isometry = (F_reverse * stretch).transpose();

  double cost = isotropic_strain_cost(stretch);

  check_equal(
      parent.superlattice().lat_column_mat(),
      stretch * isometry * unmapped_child_scel.lat_column_mat(),
      "Error in make_lattice_node: parent.superlattice().lat_column_mat() != "
      "stretch * isometry * unmapped_child_scel.lat_column_mat()");

  return LatticeNode(parent, mapped_child, stretch, isometry, cost,
                     "isotropic_strain_cost");
}

/// \brief Construct a LatticeNode using the mapping calculated by LatticeMap
///
/// \param lattice_map The lattice mapping is used to specify
/// the parent superlattice and the current solution specifies the deformation
/// gradient and choice of lattice vectors that map a supercell of the unmapped
/// child to a supercell of the parent. Specifically, it specifies:
/// - \f$ L_1 * T_1 \f$ = `lattice_map.parent_matrix()`
/// - \f$F_{parent \to child}\f$ = `lattice_map.deformation_gradient()`
/// - \f$N\f$ = `lattice_map.matrixN`
/// \param parent_prim primitive lattice being mapped to (\f$L_1\f$)
/// \param unmapped_child_prim primitive lattice being mapped (\f$L_2\f$)
///
/// Note:
/// - The lattice deformation cost is calculated using the method specified by
///   `lattice_map`.
///
/// \anchor make_lattice_node_2
/// \relates LatticeNode
///
LatticeNode make_lattice_node(LatticeMap const &lattice_map,
                              Lattice const &parent_prim,
                              Lattice const &unmapped_child_prim) {
  // see LatticeNode class documentation for more on relations

  // F_reverse * L1 * T1 * N = L2 * T2
  Eigen::Matrix3d F_reverse = lattice_map.deformation_gradient();

  // V = U_reverse.inverse()
  Eigen::Matrix3d stretch = polar_decomposition(F_reverse).inverse();

  // Q = (F_reverse * V).transpose()
  Eigen::Matrix3d isometry = (F_reverse * stretch).transpose();

  // parent.prim_lattice() = L1
  // parent.superlattice() = L1 * T1 * N
  Lattice parent_scel{lattice_map.parent_matrix() * lattice_map.matrixN(),
                      parent_prim.tol()};
  Superlattice parent{parent_prim, parent_scel};

  // mapped_child_prim = F * L2 = F_reverse.inverse() * L2
  Lattice mapped_child_prim{lattice_map.deformation_gradient().inverse() *
                            unmapped_child_prim.lat_column_mat()};

  // (mapped) child.prim_lattice() = F * L2
  // (mapped) child.superlattice() = parent.superlattice()
  Superlattice mapped_child{mapped_child_prim, parent_scel};

  double cost = lattice_map.strain_cost();
  std::string cost_method = lattice_map.cost_method();

  return LatticeNode(parent, mapped_child, stretch, isometry, cost,
                     cost_method);
}

//*******************************************************************************************

/// \brief Compare two LatticeMap objects, based on their mapping cost
/// \relates LatticeNode
bool less(LatticeNode const &A, LatticeNode const &B, double cost_tol) {
  if (!almost_equal(A.cost, B.cost, cost_tol)) return A.cost < B.cost;
  if (A.child.transformation_matrix_to_super() !=
      B.child.transformation_matrix_to_super())
    return Local::lex_lt(A.child.transformation_matrix_to_super(),
                         B.child.transformation_matrix_to_super());
  if (A.parent.transformation_matrix_to_super() !=
      B.parent.transformation_matrix_to_super())
    return Local::lex_lt(A.parent.transformation_matrix_to_super(),
                         B.parent.transformation_matrix_to_super());
  return false;
}

//*******************************************************************************************

/// \relates LatticeNode
bool identical(LatticeNode const &A, LatticeNode const &B, double cost_tol) {
  if (!almost_equal(A.cost, B.cost, cost_tol)) return false;
  if (A.parent.transformation_matrix_to_super() !=
      B.parent.transformation_matrix_to_super())
    return false;
  if (A.child.transformation_matrix_to_super() !=
      B.child.transformation_matrix_to_super())
    return false;
  return true;
}

//*******************************************************************************************

bool AssignmentNode::operator<(AssignmentNode const &B) const {
  if (empty() != B.empty()) return empty();
  if (time_reversal != B.time_reversal) return B.time_reversal;
  if (!almost_equal(translation, B.translation, 1e-6))
    return float_lexicographical_compare(translation, B.translation, 1e-6);
  return false;
}

//*******************************************************************************************

/// \relates AssignmentNode
bool identical(AssignmentNode const &A, AssignmentNode const &B) {
  if (A.empty() != B.empty()) return false;
  if (A.time_reversal != B.time_reversal) return false;
  if (!almost_equal(A.translation, B.translation, 1e-6)) return false;
  return true;
}

/// \class MappingNode
/// \brief Data structure holding a mapping between structures
///
/// Contains structure mapping results of the methods implemented by
/// StrucMapper.
///
/// ### Interpreting the mapping results ###
///
/// Mapping inputs:
/// - SimpleStructure parent: Reference structure begin mapped to. It is
///   specified via StrucMapper constructor arguments.
/// - SimpleStructure unmapped_child: Structure being mapped. It is input to a
///   StrucMapper `map_X_struc[_Y]` method.
///
/// Mapping instance:
/// - StrucMapper struc_mapper: StrucMapper instance, for mapping
///   unmapped_child to parent.
///
/// Results:
/// - MappingNode mapping: An individual solution returned by a StrucMapper
///   `map_X_struc[_Y]` method. Solutions are returned as std::set<MappingNode>
///   because, depending on the parameters used, multiple solutions may be
///   returned.
/// - SimpleStructure mapped_child: Resulting mapped structure. It is not
///   returned directly as part of the mapping results, but it can be
///   constructed using:
///
///       struc_mapper.calculator().resolve_setting(mapping_node,
///                                                 unmapped_child)
///
/// - SimpleStructure parent_superstructure: The ideal super structure of the
///   reference parent structure the child was mapped to. It is not
///   returned directly as part of the mapping results, but it can be
///   constructed using:
///
///       make_superstructure(mapping.parent.transformation_matrix_to_super(),
///                           parent)
///
///
/// ### Lattice mapping results ###
///
/// Lattice mapping results satisfy the relation:
///
///     mapped_child.lat_column_mat == mapping.lattice_node.stretch *
///         mapping.lattice_node.isometry *
///         unmapped_child.lat_column_mat *
///         mapping.lattice_node.child.transformation_matrix_to_super()
///             .cast<double>();
///
/// ### Atomic assignment results ###
///
/// Atomic assignment mapping results satisfy the relations:
///
///     mapped_child.atom_info.names[i] ==
///         unmapped_child.atom_info.names[perm[i]],
///
///     mapped_child.atom_info.coords.col(i) ==
///         parent_superstructure.coords(i) +
///         mapping.atom_displacement.col(i) ==
///         mapping.lattice_node.stretch * mapping.lattice_node.isometry *
///         unmapped_child.atom_info.coords.col(perm[i]) +
///             mapping.atomic_node.translation,
///
/// where
///
///     perm = mapping.atom_permutation,
///
/// and due to mapping within periodic boundaries, coordinate comparisons
/// must be made checking for equality up to a lattice translation.
///
///
/// ### Atom to molecule assignment results ###
///
/// \note Currently only structures with single atom molecules can be mapped.
/// These values are structured with future developments in mind and currently
/// used even though molecules only contain a single atom.
///
/// Additionally, \link MappingNode::mol_map mol_map\endlink, and \link
/// MappingNode::mol_labels mol_labels\endlink hold information used for
/// molecule mapping:
///
/// - \link MappingNode::mol_map mol_map\endlink[i] (std::set<Index>): the set
///   of atom indices of the parent superstructure that comprise the molecule
///   at its i-th site.
///
/// - \link MappingNode::mol_labels mol_labels\endlink[i]
///   (std::pair<std::string, Index>): the name and occupant index of the
///   molecule on the the i-th site of the parent superstructure
///
/// So the following hold:
///
///     mapped_child.mol_info.names[i] = mapping.mol_labels[i].first
///
///     mapped_child.mol_info.coords.col(i) =
///         mapping.lattice_node.stretch.inverse() *
///         (parent_superstructure.coords(i) + mol_displacement.col(i))
///
/// where mol_displacement is a 3xN matrix in which the i-th column gives the
/// mean displacement of the atoms in the molecule on that site:
///
///     mol_displacement.col(i) =
///         (mean of mapping.atom_displacement.col(j),
///         for j in mapping.mol_map[i]),
///
/// N is number of sites in parent_superstructure, and due to mapping within
/// periodic boundaries, coordinate comparisons must be made checking for
/// equality up to a lattice translation.
///
///
/// ### Degree of freedom and property mapping ###
///
/// DoF / properties are mapped using a transformation matrix constructed by
/// AnisoValTraits from the mapping results.
///
/// For global DoF/properties, the value is mapped using:
///
///     v_mapped = matrix * v_unmapped,
///
/// For site DoF/properties, the value on site i is mapped using:
///
///     v_mapped[i] = matrix * v_unmapped[perm[i]],
///
/// where
///
///     matrix = AnisoValTraits(type).symop_to_matrix(
///                        mapping.lattice_node.isometry,
///                        mapping.lattice_node.stretch.inverse() *
///                            mapping.atomic_node.translation,
///                        mapping.atomic_node.time_reversal)
/// and
///
///     perm = mapping.atom_permutation
///
///
/// ### Conversion to equivalent SymOp ###
///
/// If the same SimpleStructure is used for both the parent and child, and the
/// lattice of the solutions is constrained to be identical to the input
/// lattice, then the resulting mapping solutions give the factor group of the
/// structure. The set of MappingNode solutions can be converted to vector of
/// SymOp using:
///
/// \code
/// std::set<MappingNode> mappings = ... result from StrucMapper ...
/// // converts std::set<MappingNode> to
/// //     xtal::SymOpVector (= std::vector<xtal::SymOp>)
/// adapter::Adapter<xtal::SymOpVector, std::set<MappingNode>> adapter();
/// xtal::SymOpVector factor_group = adapter(mappings);
/// \endcode
///
/// See StrucMapper for a complete example getting the factor group of a
/// SimpleStructure.
///

//*******************************************************************************************

MappingNode MappingNode::invalid() {
  static MappingNode result(LatticeNode(Lattice::cubic(), Lattice::cubic(),
                                        Lattice::cubic(), Lattice::cubic(), 1),
                            0.5);
  result.is_viable = false;
  result.is_valid = false;
  result.is_partitioned = false;
  return result;
}

//*******************************************************************************************

void MappingNode::calc() {
  if (is_viable) {
    if (atomic_node.irow.empty())
      atomic_node.irow = sequence<Index>(0, atomic_node.cost_mat.rows() - 1);
    if (atomic_node.icol.empty())
      atomic_node.icol = sequence<Index>(0, atomic_node.cost_mat.cols() - 1);
    double tcost =
        hungarian_method(atomic_node.cost_mat, atomic_node.assignment,
                         cost_tol());  // + atomic_node.cost_offset;
    if (StrucMapping::is_inf(tcost)) {
      is_viable = false;
      cost = StrucMapping::big_inf();
    }
  } else
    cost = StrucMapping::big_inf();
}

//*******************************************************************************************

bool MappingNode::operator<(MappingNode const &B) const {
  double _cost_tol = max(this->cost_tol(), B.cost_tol());
  if (!almost_equal(this->cost, B.cost, _cost_tol)) {
    return this->cost < B.cost;
  }
  if (!almost_equal(this->lattice_node.cost, B.lattice_node.cost, _cost_tol)) {
    return this->lattice_node.cost < B.lattice_node.cost;
  }
  if (this->atomic_node.empty() != B.atomic_node.empty()) {
    return this->atomic_node.empty();
  }
  if (!identical(this->lattice_node, B.lattice_node, _cost_tol)) {
    return less(this->lattice_node, B.lattice_node, _cost_tol);
  }
  if (!identical(this->atomic_node, B.atomic_node)) {
    return this->atomic_node < B.atomic_node;
  }
  if (atom_permutation != B.atom_permutation)
    return std::lexicographical_compare(
        this->atom_permutation.begin(), this->atom_permutation.end(),
        B.atom_permutation.begin(), B.atom_permutation.end());

  return false;
}

/// \class StrucMapper
/// \brief Implements a method for mapping a "child" crystal structure as a
/// deformation of a reference "parent" crystal structure.
///
/// Note: This documentation uses the notation and conventions from the paper
/// "Comparing crystal structures with symmetry and geometry",
/// by John C. Thomas, Anirudh Raju Natarajan, Anton Van der Ven.
///
/// ### Example 1: Basic usage ###
///
/// In many cases, the default structure mapping options will work fine, but
/// perhaps be less efficient than if some known constraints can be applied.
/// \code
/// // parent: (reference structure)
/// xtal::BasicStructure basicstructure = ...;
/// xtal::SimpleStructure parent = make_simple_structure(basicstructure);
/// xtal::SymOpVector parent_factor_group =
///     xtal::make_factor_group(basicstructure);
///
/// // child: (structure being mapped to parent)
/// xtal::SimpleStructure child = ...;
/// xtal::SymOpVector child_factor_group = identity_group();
///
/// xtal::SimpleStrucMapCalculator calculator{
///     parent, parent_factor_group,
///     CASM::xtal::SimpleStructure::SpeciesMode::ATOM,
///     allowed_molecule_names(basicstructure)};
/// xtal::StrucMapper mapper(calculator);
///
/// std::set<xtal::MappingNode> mappings = mapper.map_deformed_struc(child);
/// \endcode
///
/// ### Example 2: All options, most general method ###
///
/// The following lays out all possible options for the most general mapping
/// approach.
/// \code
/// // --- Parent and child structures ---
///
/// // parent: (reference structure)
/// xtal::BasicStructure basicstructure = ...;
/// xtal::SimpleStructure parent = make_simple_structure(basicstructure);
/// xtal::SymOpVector parent_factor_group =
///     xtal::make_factor_group(basicstructure);
///
/// // child: (structure being mapped to parent)
/// xtal::SimpleStructure child = ...;
/// xtal::SymOpVector child_factor_group = identity_group();
///
/// // --- StrucMapper constructor parameters ---
///
/// // Specifies parent structure, factor group, allowed species on each
/// // structure site (default: SimpleStrucMapCalculator(parent));
/// xtal::SimpleStrucMapCalculator calculator{
///     parent, parent_factor_group,
///     CASM::xtal::SimpleStructure::SpeciesMode::ATOM,
///     allowed_molecule_names(basicstructure)};
///
/// // Total cost weighting, lattice deformation component (default: 0.5)
/// // total_cost = lattice_weight * lattice_deformation_cost +
/// //              (1-lattice_weight) * atomic_deformation_cost
/// double lattice_weight = 0.5;
///
/// // Potential parent superlattice volume constraint (default: 0.5)
/// double max_volume_change = 0.5;
///
/// // Perform additional checks to determine if mapping is degenerate in cost
/// // to other mappings, which can occur if the imported structure has
/// // symmetry that is incompatible with the parent structure. Results in
/// // slower execution. (default: false)
/// bool robust = false;
///
/// // If true, ensures that if no supercell volume satisfies vacancy
/// // constraints, the smallest possible volume is used.
/// // Default behavior results in no valid mapping. (default: false)
/// bool soft_va_limit = false;
///
/// // Mapping cost comparison tolerance (default: TOL)
/// double cost_tol = TOL;
///
/// // Potential parent superlattice volume constraint (default: 0.)
/// double min_va_frac = 0.0;
///
/// // Potential parent superlattice volume constraint (default: 1.)
/// double max_va_frac = 1.0;
///
/// xtal::StrucMapper mapper{calculator, lattice_weight, max_volume_change,
///                          robust, soft_va_limit, cost_tol, min_va_frac,
///                          max_va_frac};
///
/// // --- Standard mapping method parameters ---
///
/// // How many solutions to return (approximately, see detailed description
/// // for StrucMapper::map_deformed_struc) (default: 1)
/// Index k_best = 1;
///
/// // Max cost solution to keep (default: xtal::StrucMapping::big_inf())
/// double max_cost = xtal::StrucMapping::big_inf();
///
/// // All mappings better than min_cost will be reported, without contributing
/// // to 'k'. (default: -TOL)
/// double min_cost = -TOL;
///
/// // Keep invalid solutions due to failed assignment? (default: false)
/// bool keep_invalid = false;
///
/// std::set<xtal::MappingNode> mappings = mapper.map_deformed_struc(
///     child, k_best, max_cost, min_cost, keep_invalid, child_factor_group);
/// \endcode
///
/// ### Example: Using StrucMapper to find the factor group ###
///
/// If the same SimpleStructure is used for both the parent and child, and the
/// lattice of the solutions is constrained to be identical to the input
/// lattice, the resulting mapping solutions give the factor group of the
/// structure.
/// \code
/// xtal::SimpleStructure structure = ...;
///
/// // note: no "parent_factor_group" provided here
/// xtal::StrucMapper mapper{xtal::SimpleStrucMapCalculator(structure)};
///
/// double max_cost = xtal::StrucMapping::big_inf();
/// double min_cost = 1e-3; // all mappings better than min_cost are kept
/// int k = 0; // don't keep any mappings that are worse than 'min_cost'
///
/// std::set<MappingNode> mappings = mapper.map_deformed_struc_impose_lattice(
///     structure, xtal::Lattice(structure.lat_column_mat), k, max_cost,
///     min_cost);
///
/// // converts std::set<MappingNode> to
/// //     xtal::SymOpVector (= std::vector<xtal::SymOp>)
/// adapter::Adapter<xtal::SymOpVector, std::set<MappingNode>> adapter();
/// xtal::SymOpVector factor_group = adapter(mappings);
/// \endcode
///
/// ### Structure mapping summary ###
///
/// In general, structure mapping may require mapping a superstructure of the
/// child structure to a superstructure of the parent structure. In this
/// method, lattice mapping is done first, then atomic assignment for each
/// potential lattice mapping solution. A variety of settings control which
/// lattice mappings are considered, how many solutions are stored, and how
/// they are scored. For those specifics, see the StrucMapper constructor
/// documentation and the `map_X_struc[_Y]` methods which run the algorithm
/// in general or under different simplifying assumptions:
///
/// - `StrucMapper::map_deformed_struc`: (Most general) Find the k-best
///   mappings of an arbitrary child structure onto the parent structure,
///   without simplifying assumptions
/// - `StrucMapper::map_deformed_struc_impose_lattice_vols`: Find the k-best
///   mappings of an arbitrary child structure onto the parent structure,
///   specifying the range of parent superlattice volumes considered
/// - `StrucMapper::map_deformed_struc_impose_lattice`: Find the k-best
///   mappings of an arbitrary child structure onto the parent structure,
///   specifying the parent superlattice exactly
/// - `StrucMapper::map_deformed_struc_impose_lattice_node`: Find the k-best
///   mappings of an arbitrary child structure onto the parent structure,
///   specifying the lattice mapping exactly
/// - `StrucMapper::map_ideal_struc`: Find the k-best mappings of a child
///   structure onto the parent structure, assuming that the child lattice and
///   parent lattice are related by an integer transformation and a parent
///   structure point group operation
///
/// ### Mapping results ###
///
/// The `map_X_struc[_Y]` method results are returned as a set of MappingNode,
/// which may contain more than one potential solution depending on the method
/// arguments. The MappingNode data structure contains:
/// - LatticeNode \link MappingNode::lattice_node lattice_node\endlink: Holds
///   lattice mapping results
/// - AssignmentNode \link MappingNode::atomic_node atomic_node\endlink: Holds
///   atomic assignment results
/// - Additional data on the displacements, permutation, and cost determined
///   from the lattice and assignment mapping solutions.
///
/// ### Lattice mapping ###
///
/// \note In the following, `mapping` is a MappingNode result, and
/// `unmapped_child` is a `map_X_struc[_Y]` method argument.
///
/// The lattice mapping portion of the structure mapping algorithm finds
/// solutions (\f$N, V^{N}, Q^{N}\f$) of:
/// \f[
///     L_1 * T_1 * N = V^{N} * Q^{N} * L_2 * T_2,
/// \f]
/// where:
/// - \f$L_1\f$: Parent lattice (lattice of reference structure). It is equal to
///
///       this->parent().lat_column_mat()
///
/// - \f$T_1\f$: Integer transformation matrix to parent superlattice.
/// - \f$N\f$: Unimodular matrix (integer matrix, with \f$\det{N}=1\f$)
///   that transforms parent superlattice vectors to create an equivalent
///   lattice.
/// - \f$L_1 * T_1 * N\f$: The combination is the parent superlattice the child
///   has been mapped to. It is equal to
///
///       mapping.lattice_node.parent.super_lattice().lat_column_mat()
///
/// - \f$L_2\f$: Child lattice (lattice of structure to be mapped). It is equal
///   to
///
///       unmapped_child.lat_column_mat
///
/// - \f$T_2\f$: Integer transformation matrix to child superlattice (\f$L_2 *
///   T_2\f$). When the parent lattice is the lattice of the primitive
///   structure, then \f$T_2 = I\f$.
/// - \f$V^{N}\f$: Stretch, a symmetric matrix that describes the deformation
///   that maps a de-rotated child superlattice to a parent superlattice. It is
///   equal to
///
///       mapping.lattice_node.stretch
///
/// - \f$Q^{N}\f$: Isometry (distance-preserving transformation) matrix,
///   that describes the transformation that de-rotates (de-reflects, etc.) a
///   superlattice of child to align it with a parent superlattice. A property
///   of \f$Q^{N}\f$ is that \f$(Q^{N})^{-1}  == (Q^{N})^{\top}\f$. It is
///   equal to
///
///       mapping.lattice_node.isometry
///
///
/// Spelling it all out:
/// - \f$L_1\f$: The (unmapped) child lattice
/// - \f$L_1 * T_1\f$: A parent superlattice
/// - \f$L_1 * T_1 * N\f$: lattice vectors for an equivalent lattice to a parent
///   superlattice
/// - \f$L_2 * T_2\f$: A child superlattice. When the parent lattice is the
///   lattice of the primitive structure, then \f$T_2 = I\f$.
/// - \f$Q^{N} * L_2 * T_2\f$: a de-rotated child superlattice
/// - \f$V^{N} * Q^{N} * L_2 * T_2\f$: a de-rotated and undeformed child
///   superlattice
/// - \f$F^{N} = V^{N} * Q^{N}\f$: the deformation gradient for the
///   transformation of a superlattice of the child lattice to a superlattice
///   of the parent lattice (\f$F^{N} = F_{child \to parent}\f$). Note that the
///   strain associated with a configuration, whether as a DoF or in
///   MappedProperties is defined in the reverse sense, as a transformation of
///   parent (supercell of prim) to child:
///   \f[
///       F_{parent \to child} * L_1 * T_1 * N = L_2 * T_2, \\
///       F_{parent \to child} = F_{child \to parent}^{-1}
///   \f]
///   For example, the \f$U\f$ of `"Ustrain"` is defined: \f$F_{parent \to
///   child} = Q * U\f$. See LatticeNode for a more detailed description.
///
/// In general, lattice mapping proceeds using the following steps:
///
/// 1) Propose pairs of possible parent and child superlattices, \f$S_1 = L_1 *
/// T_1\f$, and \f$S_2 = L_2 * T_2\f$, respectively. For mapping to a prim,
/// \f$S_2\f$ is fixed to \f$S_2 = L_2\f$. The possible values of \f$S_1\f$ may
/// be constrained based on whether or not vacancies are allowed, and if so how
/// many (through the fixed_volume, min_va_frac, and max_va_frac parameters).
///
/// 2) For each pair \f$(S_1, S_2)\f$, find the \f$(N, F^{N})\f$, where \f$S_1
/// * N = F^{N} * S_2\f$, that minimize the strain cost, by looping over
/// unimodular matrices \f$N\f$ with elements in the range -1/+1, -2/+2, etc.
/// and solving for \f$F^{N}\f$. If requested, the k-best solutions \f$(N,
/// F^{N})\f$ for a given \f$(S_1, S_2)\f$ can be saved. (Note that the
/// solution is performed by LatticeMap, and that
/// LatticeMap::deformation_gradient() is F_{parent \to child}. See
/// isotropic_strain_cost and symmetry_breaking_strain_cost for details on
/// the strain cost calculation, which can be performed including the entire
/// lattice deformation, or just the part of the lattice deformation that
/// breaks the symmetry of the parent structure.
///
/// Sub-optimal solutions: Sub-optimal lattice mapping solutions are of interest
/// for several reasons: they might still be optimal total mapping soluations,
/// they might be important low energy transformation pathways, they might be
/// optimal mapping solutions when considering only symmetry-breaking lattice
/// deformations. Therefore, CASM also allows for keeping the k-best scoring
/// lattice mapping solutions. To avoid keeping any solution, the
/// "map_X_struc[_Y]" methods also allow for specifying a "max_cost" above
/// which mapping solutions are not kept.
///
///
/// ### Atomic assignment ###
///
/// The assignment portion of the structure mapping algorithm finds
/// solutions \f$(p_i, \vec{t}, \vec{d}(i))\f$ of:
/// \f[
///     \vec{r_1}(i) + \vec{d}(i) = V * Q * \vec{r_2}(p_i) + \vec{t}
/// \f]
///
/// where:
/// - \f$\vec{r_1}(i)\f$: Vector of coordinates of atoms in the parent
///   superstructure. The value \f$\vec{r_1}(i)\f$ represents the Cartesian
///   coordinate of the \f$i\f$-th atom in the parent structure. The parent
///   superstructure is not returned directly as part of the mapping results,
///   but it can be constructed using:
///
///       xtal::Simpletructure parent_superstructure = make_superstructure(
///           mapping.parent.transformation_matrix_to_super(), this->parent());
///
///   Then the \f$i\f$-th atom coordinate, \f$\vec{r_1}(i)\f$, is equal to:
///
///       parent_superstructure.atom_info.coords.col(i)
///
/// - \f$\vec{r_2}(i)\f$: Vector of coordinates of atoms in the unmapped child
///   superstructure. The value \f$\vec{r_1}(i)\f$ represents the Cartesian
///   coordinate of the \f$i\f$-th atom in the unmapped child structure.
/// - \f$V * Q\f$: Lattice transformation, from the unmapped child superlattice
///   to the parent superlattice, as determined by a solution to the lattice
///   mapping problem.
/// - \f$p_i\f$: A permutation vector, describes which atom in the unmapped
///   structure (\f$p_i\f$) is mapped to the i-th site of the mapped structure.
///   Values of \f$p_i\f$ greater than the number of atoms in the unmapped
///   structure indicate inferred vacancies.
/// - \f$\vec{t}\f$: A translation vector, in Cartesian coordinates, of the de-
///   rotated and undeformed (mapped) child superstructure that minimizes the
///   atomic deformation cost.
/// - \f$\vec{d}(i)\f$: The displacement associated with the atom at the i-th
///   site in parent superstructure.
///
/// In general, given a lattice mapping solution, the assignment problem is
/// solved using the following steps:
///
/// 1) Considering, in turn, possible translations of minority atoms in the de-
/// rotated and undeformed child superstructure to sites in which that atom
/// type is allowed in the parent superstructure, translate the child
/// superstrucuture into registry with (at least one site of) the parent
/// superstructure.
///
/// 2) Use the Hungarian Method to find the optimal assignment permutation,
/// specifying the cost of assignment as infinite if atom type is not allowed
/// on the parent superstructure site, otherwise with a cost based on either
/// geometric distance or symmetry-breaking distance. See atomic_cost for more
/// details on the atomic deformation cost calculations.
///
/// 3) After the optimal assignment is determined, subtract the average
/// displacement from the initial proposed  translation (checking that the
/// average total displacement is zero and performing corrections for periodic
/// boundary effects if this is not the case), and re-calculate the atomic
/// deformation cost using the same assignment. This results in an optimal
/// assignment solution \f$(p_i, \vec{t}, \vec{d}(i))\f$.
///
/// Sub-optimal solutions: As with lattice mapping, sub-optimal assingment
/// solutions may be required. In this case, Murty's algorithm is applied (which
/// systematically forces some assignments) to enumerate the k-best assignment
/// solutions \f$(p_i, \vec{t}, \vec{d}(i))\f$.
///
///
/// ### Total mapping cost ###
///
/// The optimal total mapping solution, a combination of a lattice mapping, and
/// an assignment mapping calculated given that lattice mapping, is determined
/// by minimizing a total cost function which combines the lattice deformation
/// cost and atomic deformation costs:
///
///     total_cost = lattice_weight * lattice_deformation_cost +
///                  atomic_weight * atomic_deformation_cost
///
///

//*******************************************************************************************

StrucMapper::StrucMapper(StrucMapCalculatorInterface const &calculator,
                         double _lattice_weight /*= 0.5*/,
                         double _max_volume_change /*= 0.5*/,
                         bool _robust /*= false*/,
                         bool _soft_va_limit /*= false*/,
                         double _cost_tol /*= TOL*/,
                         double _min_va_frac /*= 0.*/,
                         double _max_va_frac /*= 1.*/)
    : m_calc_ptr(calculator.clone()),
      m_max_volume_change(_max_volume_change),
      m_robust(_robust),
      m_soft_va_limit(_soft_va_limit),
      m_cost_tol(max(1e-10, _cost_tol)),
      m_xtal_tol(TOL),
      m_lattice_transformation_range(1),
      m_symmetrize_lattice_cost(false),
      m_symmetrize_atomic_cost(false),
      m_filtered(false) {
  set_min_va_frac(_min_va_frac);
  set_max_va_frac(_max_va_frac);

  // squeeze lattice_weight into (0,1] if necessary
  set_lattice_weight(_lattice_weight);

  // Make sure that max_volume_change is positive
  m_max_volume_change = max(3 * xtal_tol(), _max_volume_change);
}

//*******************************************************************************************

Index StrucMapper::_n_species(SimpleStructure const &sstruc) const {
  return calculator().struc_info(sstruc).size();
}

//*******************************************************************************************

/// \brief Calculate min_vol, max_vol range from min_va_frac, max_va_frac, and
/// properties of the child and parent structures and calculator.
///
/// Calculation inputs are:
/// - min_va_frac
/// - max_va_frac
/// - parent structure
/// - child structure
/// - soft_va_limit
///
std::pair<Index, Index> StrucMapper::_vol_range(
    const SimpleStructure &unmapped_child) const {
  Index min_vol(0), max_vol(0);
  // mapped_result.clear();
  SimpleStructure::Info const &c_info(calculator().struc_info(unmapped_child));

  // Using _n_species() instead of parent().n_mol() here.
  // For mapping molecules, may need to push this entire routine into the
  // StrucMapCalculator?
  double N_sites_p = double(_n_species(parent()));

  double N_sites_c = double(_n_species(unmapped_child));

  if (calculator().fixed_species().size() > 0) {
    std::string tcompon = calculator().fixed_species().begin()->first;
    int ncompon(0);
    for (std::string const &sp : c_info.names) {
      if (sp == tcompon) ncompon++;
    }
    min_vol = ncompon / int(calculator().fixed_species().begin()->second);
    max_vol = min_vol;
  } else {
    // Try to narrow the range of supercell volumes -- the best bounds are
    // obtained from the convex hull of the end-members, but we need to wait for
    // improvements to convex hull routines

    int max_n_va = calculator().max_n_va();

    // Absolute largest Va fraction
    double max_va_frac_limit = double(max_n_va) / N_sites_p;

    double t_min_va_frac = min(min_va_frac(), max_va_frac_limit);
    double t_max_va_frac = min(max_va_frac(), max_va_frac_limit);

    // Use vacancy range to narrow volume range.

    // the number of sites implied in order to realize va_frac >= min_va_frac
    // (TOL is not super important):
    double min_n_sites = floor(N_sites_c / (1. - t_min_va_frac) - TOL);

    // the number of sites implied to realize va_frac <= max_va_frac  (TOL is
    // not super important):
    double max_n_sites = ceil(N_sites_c / (1. - t_max_va_frac) + TOL);

    // min_vol assumes min number vacancies -- best case scenario (buffer by
    // half an atom)
    min_vol = ceil((min_n_sites - 0.5) / N_sites_p);

    // max_vol assumes min number vacancies -- best case scenario (buffer by
    // half an atom)
    max_vol = floor((max_n_sites + 0.5) / N_sites_p);

    // If volume is not fully determined, use volume ratio of child to parent to
    // narrow search
    if (min_vol != max_vol) {
      // Nvol is approximate integer volume-- assume that actual integer is
      // within m_max_volume_change*100% of this volume, and use it to tighten
      // our bounds
      double Nvol = (std::abs(unmapped_child.lat_column_mat.determinant() /
                              parent().lat_column_mat.determinant()));
      min_vol = min(
          max_vol, max<Index>(round((1.0 - m_max_volume_change) * double(Nvol)),
                              min_vol));
      max_vol = max(
          min_vol, min<Index>(round((1.0 + m_max_volume_change) * double(Nvol)),
                              max_vol));
    }
  }

  // make sure volume range is above zero
  if (m_soft_va_limit) {
    Index smallest_possible_vol = ceil((N_sites_c - 0.5) / N_sites_p);
    min_vol = max<Index>(min_vol, smallest_possible_vol);
    max_vol = max<Index>(max_vol, smallest_possible_vol);
  }

  return std::pair<Index, Index>(min_vol, max_vol);
}

//*******************************************************************************************

SimpleStructure const &StrucMapper::parent() const {
  return calculator().parent();
}

//*******************************************************************************************

/// \brief Calculate min_vol, max_vol range from min_va_frac, max_va_frac, and
/// properties of the child and parent structures and calculator.
std::set<MappingNode> StrucMapper::_seed_from_vol_range(
    SimpleStructure const &unmapped_child, Index k, Index min_vol,
    Index max_vol, double max_lattice_cost, double min_lattice_cost,
    SymOpVector const &child_factor_group) const {
  if (!valid_index(min_vol) || !valid_index(max_vol) || max_vol < min_vol) {
    std::stringstream msg;
    msg << "Error in StrucMapper: invalid volume range [" << min_vol << ","
        << max_vol << "]";
    throw std::runtime_error(msg.str());
  }

  // There is nothing to enumerate, don't even bother.
  if (max_vol < 1) {
    return {};
  }
  // Ensure that you don't try to enumerate size zero supercells
  min_vol = std::max(min_vol, Index{1});

  Lattice child_lat(unmapped_child.lat_column_mat, xtal_tol());
  std::set<MappingNode> mapping_seed;
  for (Index i_vol = min_vol; i_vol <= max_vol; i_vol++) {
    std::vector<Lattice> lat_vec;
    for (Lattice const &lat : _lattices_of_vol(i_vol)) {
      if (m_filtered && !_filter_lat(lat, child_lat)) {
        continue;
      }
      lat_vec.push_back(lat);
    }

    std::set<MappingNode> t_seed = _seed_k_best_from_super_lats(
        unmapped_child, lat_vec,
        {Lattice(unmapped_child.lat_column_mat, xtal_tol())}, k,
        max_lattice_cost, max(min_lattice_cost, cost_tol()),
        child_factor_group);

    mapping_seed.insert(std::make_move_iterator(t_seed.begin()),
                        std::make_move_iterator(t_seed.end()));
  }
  return mapping_seed;
}

//*******************************************************************************************
/*
 * Given a structure and a mapping node, find a perfect supercell of the prim
 * that is equivalent to structure's lattice and then try to map the structure's
 * basis onto that supercell
 *
 * Returns false if no mapping is possible, or if the lattice is not ideal
 *
 * What this does NOT do:
 *    -Check if the imported Structure is the same as one in a smaller Supercell
 *
 */
//*******************************************************************************************

/// \brief Find the k-best mappings of a child structure onto the parent
/// structure, assuming that the child lattice and parent lattice are related
/// by an integer transformation and a parent structure point group operation
///
/// Similar to `map_deformed_struc`, but with the following modification:
/// - A single lattice mapping, `imposed_node`, is calculated from the
///   assumption that the child lattice and parent lattice are related
///   by an integer transformation, and the structures are related by a parent
///   structure point group operation. Thus, the mapping will satisfy L1 * T1 =
///   Q * L2 * T2, where Q is restricted to parent structure point group
///   operations.
/// - Parameters `min_va_frac`, `max_va_frac`, `max_volumne_change`,
///   `soft_va_limit`, `add_allowed_lattice`, and `set_filter` have no effect.
/// - All assignment is the same as in `map_deformed_struc`.
///
/// \param unmapped_child Input structure to be mapped onto parent structure
/// \param k Number of k-best mapping relations to return
/// \param max_cost Search will terminate once no mappings better than
///     max_cost are found
/// \param min_cost All mappings better than min_cost will be reported,
///     without contributing to 'k'
/// \param keep_invalid If true, invalid mappings (when an atomic assignment
///     maps an atom of one type to a different type) are retained in output;
///     otherwise they are discarded
///
std::set<MappingNode> StrucMapper::map_ideal_struc(
    const SimpleStructure &unmapped_child, Index k,
    double max_cost /*=StrucMapping::big_inf()*/, double min_cost /*=-TOL*/,
    bool keep_invalid /*=false*/) const {
  // xtal::is_superlattice isn't very smart right now, and will return
  // false if the two lattices differ by a rigid rotation
  // In the future this may not be the case, so we will assume that
  // unmapped_child may be rigidly rotated relative to prim
  Eigen::Matrix3d trans_mat;

  // c_lat must be an ideal supercell of the parent lattice, but it need not be
  // canonical We will account for the difference in orientation between c_lat
  // and the canonical supercell, which must be related by a point group
  // operation
  Lattice c_lat(unmapped_child.lat_column_mat, xtal_tol());

  bool c_lat_is_supercell_of_parent;
  std::tie(c_lat_is_supercell_of_parent, trans_mat) = xtal::is_superlattice(
      c_lat, Lattice(parent().lat_column_mat), xtal_tol());
  if (!c_lat_is_supercell_of_parent) {
    return {};
  }

  // We know unmapped_child.lattice() is a supercell of the prim, now we have to
  // reorient 'unmapped_child' by a point-group operation of the parent to match
  // canonical lattice vectors This may not be a rotation in the child
  // structure's point group
  Lattice derot_c_lat(canonical::equivalent(
      Lattice(parent().lat_column_mat * trans_mat, xtal_tol()),
      calculator().point_group()));

  // We now find a transformation matrix of c_lat so that, after transformation,
  // it is related to derot_c_lat by rigid rotation only. Following line finds R
  // and T such that derot_c_lat = R*c_lat*T
  auto res = xtal::is_equivalent_superlattice(
      derot_c_lat, c_lat, calculator().point_group().begin(),
      calculator().point_group().end(), xtal_tol());
  LatticeNode lattice_node(
      Lattice(parent().lat_column_mat, xtal_tol()), derot_c_lat, c_lat,
      Lattice(unmapped_child.lat_column_mat * res.second.cast<double>(),
              xtal_tol()),
      _n_species(unmapped_child), 0. /*strain_cost is zero in ideal case*/);

  return map_deformed_struc_impose_lattice_node(
      unmapped_child, lattice_node, k, max_cost, min_cost, keep_invalid);
}

//*******************************************************************************************

/// \brief Find the k-best mappings of an arbitrary child structure onto the
/// parent structure, without simplifying assumptions
///
/// Besides parameters passed via the StrucMapper consturctor and
/// `map_deformed_struc` method call, this method is also controlled by
/// additional parameters that may be set on StrucMapper. The StrucMapper
/// modifiers that control these parameters are listed here, and the effects
/// are described in the notes below.
/// - `StrucMapper::add_allowed_lattices`
/// - `StrucMapper::set_filter`
/// - `StrucMapper::set_symmetrize_lattice_cost`
/// - `StrucMapper::set_symmetrize_atomic_cost`
///
/// Description:
/// - The value `k` is a target, but more or fewer mapping solutions may be
/// returned. Solutions with identical cost will not be excluded just for
/// hitting the `k` limit. For example, if k and k+1, k+2, etc mappings have
/// identical cost, k will be increased until k+n mapping has cost greater than
/// mapping k. Solutions with cost great than (worse than) `max_cost` will never
/// be included. Solutions with cost less than (better than) `min_cost` will
/// always be included and do not count against `k`.
/// - Candidate initial parent superlattices (\f$L_1 * T_1\f$) are determined
///   from:
///   - The range of possible volumes (integer volume multiple of the parent
///     lattice volume) is calculated to be consistent with `min_va_frac`,
///     `max_va_frac`, `max_volume_change`, the parent structure and allowed
///     species, the child structure, and `soft_va_limit`.
///   - If lattices have been added via "add_allowed_lattice", then for each
///     volume in the range, only "allowed_lattices" are considered; else
///     at each volume enumerate canonical parent superlattices (using the point
///     group of the provided parent structure factor group).
///   - If a lattice filter has been set, then enumerated lattices must pass
///     the lattice filter.
/// - Given an initial parent superlattice, the k-best scoring lattice mapping
///   solutions \f$(N, F^{N})\f$, where \f$(L_1 * T_1) * N = F^{N} * L_2\f$ are
///   kept and used as the starting point for atomic assignment. The symmetry-
///   breaking lattice deformation cost is used if
///   `symmetrize_lattice_cost==true`.
/// - Given a lattice mapping, candidate atomic assignments are determined by:
///   - First, finding a site in child whose occupant has the lowest number of
///     compatible sites in the parent. This will minimize translations that
///     need to be considered.
///   - Then, generating the translations from that site in the derotated and
///     undeformed child structure (according to the current lattice mapping
///     solution) to all chemically compatible sites in the parent (removing
///     symmetrically equivalent translations).
///   - If partitioning is allowed ('robust' option or k > 1), Murty's
///     algorithm is also used to systematically force some assignments in
///     order to find the k-best assignments.
///   - For each translation/assignment considered, assignments are calculated
///     by the Hungarian Method, and then the translation is adjusted so that
///     the mean displacement is zero. Then the atomic deformation is scored.
///     The symmetry-breaking atomic deformation cost is used if
///     `symmetrize_atomic_cost==true`.
/// - Total mapping solutions are scored (mapping_node.cost) using the weighted
///   average of `lattice_node.cost` (with weight `lattice_weight`) and
///   `atomic_node.cost` (with weight `1-lattice_weight`).
///
/// Sub-optimal solutions: As with lattice mapping, sub-optimal assignment
/// solutions may be required. In this case, Murty's algorithm is applied (which
/// systematically forces some assignments) to enumerate the k-best assignment
/// solutions (perm, trans).
/// - Suboptimal candidates are generated based on the initial atomic mapping
///
/// \param unmapped_child Input structure to be mapped onto parent structure
/// \param k Number of k-best mapping relations to return
/// \param max_cost Search will terminate once no mappings better than
/// max_cost are found
/// \param min_cost All mappings better than min_cost will be reported,
///     without contributing to 'k'
/// \param keep_invalid If true, invalid mappings (when an atomic assignment
///     maps an atom of one type to a different type) are retained in output;
///     otherwise they are discarded
/// \param child_factor_group TODO: Document this. In practice, just use the
///     default.
///
/// \result std::set<MappingNode> a list of valid mappings, sorted first by
/// cost, and then other attributes
std::set<MappingNode> StrucMapper::map_deformed_struc(
    const SimpleStructure &unmapped_child, Index k /*=1*/,
    double max_cost /*=StrucMapping::big_inf()*/, double min_cost /*=-TOL*/,
    bool keep_invalid /*=false*/, SymOpVector const &child_factor_group) const {
  auto vols = _vol_range(unmapped_child);
  return map_deformed_struc_impose_lattice_vols(
      unmapped_child, vols.first, vols.second, k, max_cost, min_cost,
      keep_invalid, child_factor_group);
}

//*******************************************************************************************

/// \brief Find the k-best mappings of an arbitrary child structure onto the
/// parent structure, specifying the range of parent superlattice volumes
/// considered
///
/// Similar to `map_deformed_struc`, but with the following modification:
/// - The inclusive range [min_vol, max_vol] of candidate parent superlattice
///   volumes can be specified directly by `min_vol` and `max_vol`.
/// - If invalid values of `min_vol` or `max_vol` are provided (negative values
///   or max_vol < min_vol), then this method throws.
/// - Parameters `min_va_frac`, `max_va_frac`, `max_volume_change`,
///   and `soft_va_limit` have no effect.
/// - Constraints set by `add_allowed_lattice` and `set_filter` are still in
///   effect, and all assignment is the same as in `map_deformed_struc`.
///
/// \param unmapped_child Input structure to be mapped onto parent structure
/// \param min_vol,max_vol The inclusive range [min_vol, max_vol] of candidate
///      parent superlattice volumes considered.
/// \param k Number of k-best mapping relations to return
/// \param max_cost Search will terminate once no mappings better than
/// max_cost are found
/// \param min_cost All mappings better than min_cost will be reported,
///     without contributing to 'k'
/// \param keep_invalid If true, invalid mappings (when an atomic assignment
///     maps an atom of one type to a different type) are retained in output;
///     otherwise they are discarded
/// \param child_factor_group TODO: Document this. In practice, just use the
///     default.
///
std::set<MappingNode> StrucMapper::map_deformed_struc_impose_lattice_vols(
    const SimpleStructure &unmapped_child, Index min_vol, Index max_vol,
    Index k /*=1*/, double max_cost /*=StrucMapping::big_inf()*/,
    double min_cost /*=-TOL*/, bool keep_invalid /*=false*/,
    SymOpVector const &child_factor_group) const {
  int seed_k = 10 + 5 * k;
  std::set<MappingNode> mapping_seed = _seed_from_vol_range(
      unmapped_child, seed_k, min_vol, max_vol,
      max_cost / (this->lattice_weight()),
      max(min_cost / (this->lattice_weight()), cost_tol()), child_factor_group);

  bool no_partition = !m_robust && k <= 1;
  k_best_maps_better_than(unmapped_child, mapping_seed, k, max_cost, min_cost,
                          keep_invalid, false, no_partition);

  return mapping_seed;
}

/// \brief Find the k-best mappings of an arbitrary child structure onto the
/// parent structure, specifying the parent superlattice exactly, but not the
/// way the child lattice maps to the parent superlattice
///
/// Similar to `map_deformed_struc`, but with the following modification:
/// - Only 1 parent superlattice, `imposed_lat` is considered.
/// - Parameters `min_va_frac`, `max_va_frac`, `max_volumne_change`,
///   `soft_va_limit`, `add_allowed_lattice`, and `set_filter` have no effect.
/// - Parent superlattices that are symmetrically equivalent to `imposed_lat`
///   via parent structure point group operations, but not identical are not
///   considered.
/// - The child structure may be mapped to a parent superlattice that is
///   identical to `imposed_lat`, but with differing lattice vectors (i.e. the
///   lattice mapping may be \f$L_1 * T_1 * N = F^{N} * L_2\f$, where \f$L_1 *
///   T_1\f$ equals `imposed_lat`, and \f$N\f$ is a unimodular matrix that may
///   not be \f$I\f$.)
/// - All assignment is the same as in `map_deformed_struc`.
///
/// \param unmapped_child Input structure to be mapped onto parent structure
/// \param imposed_lat The imposed parent superlattice, `imposed_lat` = \f$L1 *
///     T1\f$.
/// \param k Number of k-best mapping relations to return
/// \param max_cost Search will terminate once no mappings better than
///     max_cost are found
/// \param min_cost All mappings better than min_cost will be reported,
///     without contributing to 'k'
/// \param keep_invalid If true, invalid mappings (when an atomic assignment
///     maps an atom of one type to a different type) are retained in output;
///     otherwise they are discarded
/// \param child_factor_group TODO: Document this. In practice, just use the
///     default.
///
std::set<MappingNode> StrucMapper::map_deformed_struc_impose_lattice(
    const SimpleStructure &unmapped_child, const Lattice &imposed_lat,
    Index k /*=1*/, double max_cost /*=StrucMapping::big_inf()*/,
    double min_cost /*=-TOL*/, bool keep_invalid /*=false*/,
    SymOpVector const &child_factor_group) const {
  std::set<MappingNode> mapping_seed = _seed_k_best_from_super_lats(
      unmapped_child, {imposed_lat},
      {Lattice(unmapped_child.lat_column_mat, xtal_tol())}, k,
      max_cost / (this->lattice_weight()),
      max(min_cost / (this->lattice_weight()), cost_tol()), child_factor_group);

  bool no_partition = !m_robust && k <= 1;
  k_best_maps_better_than(unmapped_child, mapping_seed, k, max_cost, min_cost,
                          keep_invalid, false, no_partition);
  return mapping_seed;
}

/// \brief Find the k-best mappings of an arbitrary child structure onto the
/// parent structure, specifying the lattice mapping exactly
///
/// Similar to `map_deformed_struc`, but with the following modification:
/// - A single lattice mapping, `imposed_node`, is specified exactly.
/// - Parameters `min_va_frac`, `max_va_frac`, `max_volumne_change`,
///   `soft_va_limit`, `add_allowed_lattice`, and `set_filter` have no effect.
/// - All assignment is the same as in `map_deformed_struc`.
///
/// \param unmapped_child Input structure to be mapped onto parent structure
/// \param imposed_node The imposed lattice mapping, encodes \f$T_1, N, F^{N\f$
///     in the lattice mapping relationship \f$L_1 * T_1 * N = F^{N} * L_2\f$,
///     along with the lattice deformation cost.
/// \param k Number of k-best mapping relations to return
/// \param max_cost Search will terminate once no mappings better than
/// max_cost are found
/// \param min_cost All mappings better than min_cost will be reported,
///     without contributing to 'k'
/// \param keep_invalid If true, invalid mappings (when an atomic assignment
///     maps an atom of one type to a different type) are retained in output;
///     otherwise they are discarded
/// \param child_factor_group TODO: Document this. In practice, just use the
///     default.
///
std::set<MappingNode> StrucMapper::map_deformed_struc_impose_lattice_node(
    const SimpleStructure &unmapped_child, const LatticeNode &imposed_node,
    Index k /*=1*/, double max_cost /*=StrucMapping::big_inf()*/,
    double min_cost /*=-TOL*/, bool keep_invalid /*=false*/) const {
  std::set<MappingNode> mapping_seed;
  mapping_seed.emplace(imposed_node, m_lattice_weight);
  bool no_partition = !m_robust && k <= 1;
  k_best_maps_better_than(unmapped_child, mapping_seed, k, max_cost, min_cost,
                          keep_invalid, false, no_partition);
  return mapping_seed;
}

//*******************************************************************************************

std::vector<Lattice> StrucMapper::_lattices_of_vol(Index prim_vol) const {
  if (!valid_index(prim_vol)) {
    throw std::runtime_error("Cannot enumerate lattice of volume " +
                             std::to_string(prim_vol) +
                             ", which is out of bounds.\n");
  }

  // If you specified that you wanted certain lattices, return those, otherwise
  // do the usual enumeration
  if (this->lattices_constrained()) {
    // This may very well return an empty vector, saving painful time
    // enumerating things
    return m_allowed_superlat_map[prim_vol];
  }

  // If we already have candidate lattices for the given volume, return those
  auto it = m_superlat_map.find(prim_vol);
  if (it != m_superlat_map.end()) return it->second;

  // We don't have any lattices for the provided volume, enumerate them all!!!
  std::vector<Lattice> &lat_vec = m_superlat_map[prim_vol];

  auto pg = calculator().point_group();
  SuperlatticeEnumerator enumerator(
      pg.begin(), pg.end(), Lattice(parent().lat_column_mat, xtal_tol()),
      ScelEnumProps(prim_vol, prim_vol + 1));

  for (auto it = enumerator.begin(); it != enumerator.end(); ++it) {
    Lattice canon_lat = *it;
    if (canonical::check(canon_lat, calculator().point_group())) {
      canon_lat = canonical::equivalent(canon_lat, calculator().point_group());
    }
    lat_vec.push_back(canon_lat);
  }

  return lat_vec;
}

/// Find k-best mappings
///
/// This function gets called by all of the `map_X_struc[_impose_Y]`
/// functions to perform assignment mappings and finalize the total mappings.
///
/// Starting from a child structure and a queue of (possibly unfinalized)
/// mapping nodes, visit each node in the queue, test its fitness, and add any
/// of its viable derivative nodes to the queue.
///
/// The queue is ordered by the lower bound of a node's mapping cost. This lower
/// cost is less than or equal to its finalized mapping cost or the finalized
/// mapping costs of any of its derivatives node. As such, the queue only needs
/// to be iterated over once.
///
/// At the end of the process, the queue may be pruned depending on the choice
/// of input parameters.
///
///\param unmapped_child Input structure to be mapped onto parent structure
///\param queue
///\parblock
///  List of partial mappings (i.e., MappingNode with LatticeNode only) to
///  seed search. Partial mappings are removed from list as they are searched
///  finalized mappings are inserted into list as they are found.
///\endparblock
///
/// \param k Number of k-best mapping relations to return
/// \param max_cost Search will terminate once no mappings better than
///     max_cost are found
/// \param min_cost All mappings better than min_cost will be reported,
///     without contributing to 'k'
/// \param keep_invalid If true, invalid mappings are retained in output;
///     otherwise they are discarded
/// \param keep_tail If true, any partial or unresolved mappings at back of
///     queue will be retained (they are deleted otherwise)
/// \param no_partition
/// \parblock
///   If true, search over suboptimal basis permutations will be skipped (The
///   optimal mapping will be found, but suboptimal mappings may be excluded
///   in the k-best case. Degenerate permutations of atoms will not be
///   identified.)
/// \endparblock
///
///\result Index specifying the number of complete, valid mappings in the
/// solution set
///
Index StrucMapper::k_best_maps_better_than(
    SimpleStructure const &unmapped_child, std::set<MappingNode> &queue,
    Index k, double max_cost /*=StrucMapping::big_inf()*/,
    double min_cost /*=-TOL*/, bool keep_invalid /*=false*/,
    bool keep_tail /*= false*/, bool no_partition /*= false*/) const {
  int nfound = 0;
  // Track pairs of supercell volumes that are chemically incompatible
  std::set<std::pair<Index, Index>> vol_mismatch;

  if (k == 0) {
    // If k==0, then we only keep values less than min_cost
    // However, max_cost controls search loop, so we set max_cost to min_cost
    max_cost = min_cost;
  }

  auto it = queue.begin();
  while (it != queue.end()) {
    bool erase = true;
    auto current = it;

    if (it->cost <= (max_cost + this->cost_tol())) {
      // If supercell volumes have already been determined incompatible, we do
      // nothing; current node is deleted
      if (!vol_mismatch.count(current->vol_pair())) {
        // Consider two exlusive cases (and base case, in which current node
        // isn't even viable)
        if (current->atomic_node.empty()) {
          // Case 1: Current node only describes a lattice mapping with no basis
          // mapping
          //         Perform basis mapping, insert result into queue, and erase
          //         current (new node must have cost greather than the current
          //         node, so will
          //          appear later in the queue)
          if (!Local::initial_atomic_maps(
                  unmapped_child, *current, calculator(), max_cost,
                  symmetrize_atomic_cost(), std::inserter(queue, current))) {
            // If no basis maps are viable, it indicates volume mismatch; add to
            // vol_mismatch
            vol_mismatch.insert(current->vol_pair());
          }
        } else if (current->is_viable) {
          // Case 2: Current node is a complete mapping and is viable
          //         Either it is a valid node, and thus part of the solution
          //         set, or it is invalid, but we must add its partition to the
          //         queue because it may have a suboptimal derivative mapping
          //         that is valid and is part of the solution set

          if (current->is_valid && current->cost > min_cost && nfound < k) {
            // current node is part of solution set
            ++nfound;

            // Need to account for case where mapping k+1 is as good as mapping
            // k A few ways to make this happen, but we will do it by increasing
            // k when nfound==k, and shrink max_cost to be current cost
            if (nfound == k) {
              ++k;
              max_cost = current->cost;  // + tol();
            }
          }

          // Regardless of validity, we partition current node and add the
          // derivative nodes to queue (but only if we haven't reached stopping
          // condition) Skip partitioning if node is already partitioned, or if
          // caller has asked not to
          if (nfound < k || current->cost <= min_cost) {
            if (!(no_partition || current->is_partitioned)) {
              Local::partition_node(*current, calculator(), unmapped_child,
                                    symmetrize_atomic_cost(),
                                    std::inserter(queue, current));
            }

            // Keep current node if it is in the solution set if we have been
            // asked to keep invalids
            if (current->is_valid || keep_invalid) erase = false;
          }
        }
        // Never keep unviable nodes or incomplete nodes
      }
    } else {
      erase = !keep_tail;
    }
    // Safe to increment here:
    //  1) No continue/break statements
    //  2) Nothing has been deleted yet
    //  3) Everything that needs to be inserted has been inserted
    ++it;

    // Erase current if no longer needed
    if (erase) queue.erase(current);
  }

  return nfound;
}

//****************************************************************************************************************

/// \brief Generate a set of mapping seeds (i.e., MappingNode with
/// LatticeNode only) from a list of supercells of the parent structure and
/// a list of supercells of the child structure
///
/// \param k Number of k-best mapping relations to return
/// \param max_lattice_cost Search will terminate once no lattice mappings
/// better than max_lattice_cost are found
/// \param min_lattice_cost All lattice mappings better than min_lattice_cost
/// will be returned, without contributing to 'k'
///
/// Note:
/// - This method ignores min_va_frac, max_va_frac,
std::set<MappingNode> StrucMapper::_seed_k_best_from_super_lats(
    SimpleStructure const &unmapped_child,
    std::vector<Lattice> const &_parent_scels,
    std::vector<Lattice> const &_child_scels, Index k,
    double max_lattice_cost /*=StrucMapping::small_inf()*/,
    double min_lattice_cost /*=1e-6*/,
    SymOpVector const &child_factor_group) const {
  Lattice p_prim_lat(parent().lat_column_mat, xtal_tol());
  Lattice c_prim_lat(unmapped_child.lat_column_mat, xtal_tol());
  std::set<MappingNode> result;

  if (k == 0 || !valid_index(k)) {
    max_lattice_cost = min_lattice_cost;
  }
  for (Lattice const &c_lat : _child_scels) {
    auto pg_indices = invariant_subgroup_indices(c_lat, child_factor_group);
    SymOpVector c_lat_factor_group;
    c_lat_factor_group.reserve(pg_indices.size());
    for (Index i : pg_indices) {
      c_lat_factor_group.push_back(child_factor_group[i]);
    }

    for (Lattice const &p_lat : _parent_scels) {
      LatticeMap lattice_map(p_lat, c_lat, this->lattice_transformation_range(),
                             calculator().point_group(), c_lat_factor_group,
                             max_lattice_cost, symmetrize_lattice_cost(),
                             cost_tol());

      // lattice_map is initialized to first mapping better than
      // 'max_lattice_cost', if such a mapping exists We will continue checking
      // possibilities until all such mappings are exhausted
      while (lattice_map.strain_cost() < (max_lattice_cost + cost_tol())) {
        // Mappings worse than min_lattice_cost count against k
        if (lattice_map.strain_cost() > min_lattice_cost) {
          if (k == 0 || !valid_index(k)) {
            // If k is already depleted, we will still add this mapping, but
            // adjust max_lattice_cost to avoid adding anything worse
            max_lattice_cost = max(min_lattice_cost, lattice_map.strain_cost());
          } else {
            // If k is not depleted, decrement k
            k--;
          }
        }
        LatticeNode lattice_node{lattice_map, p_prim_lat, c_prim_lat};
        result.emplace(lattice_node, this->lattice_weight());
        lattice_map.next_mapping_better_than(max_lattice_cost);
      }
    }
  }
  return result;
}

}  // namespace xtal
}  // namespace CASM
