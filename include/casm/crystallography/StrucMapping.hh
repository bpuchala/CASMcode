#ifndef CASM_StrucMapping
#define CASM_StrucMapping

#include <unordered_set>
#include <vector>

#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/StrucMapCalculatorInterface.hh"
#include "casm/crystallography/Superlattice.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/external/Eigen/Core"
#include "casm/global/definitions.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/misc/CASM_math.hh"
#include "casm/misc/cloneable_ptr.hh"

namespace CASM {
namespace xtal {
class Lattice;
class LatticeMap;
class SimpleStructure;
class StrucMapCalculatorInterface;

// In this file:
struct LatticeNode;
struct AssignmentNode;
struct MappingNode;
class StrucMapper;

namespace StrucMapping {

// TODO:
// - Pass all StrucMapper parameters via constructor or via
// `map_X` function call and remove StrucMapper parameter modifiers such as
// `StrucMapper::set_filter` or `StrucMapper::set_max_va_frac`.
// - Find an alternative to the use of StrucMapCalculatorInterface, or reduce
// it to only the virtual members and make it optional.
// - Document more explicitly the conversion of std::set<MappingNode> results
// to SymOpVector
// - Document use of k=0 mapping
// - Consider 'isometry' vs 'rigid_rotation' vs.
//   'distance_preserving_transformation' etc.

/// Lattice filter function for structure mapping
///
/// The filter function is of the form
/// `bool filter(parent_prim_lattice, proposed_parent_superlattice)`, where
/// `parent_prim_lattice` is the lattice of the primitive parent structure, and
/// `proposed_parent_superlattice` is a proposed superlattice of the parent
/// structure
typedef std::function<bool(Lattice const &, Lattice const &)>
    LatticeFilterFunction;

typedef std::vector<std::vector<Index>> PermuteOpVector;
/// \brief Very large value used to denote invalid or impossible mapping
inline double big_inf() { return 10E20; }

/// \brief use as default value to initialize mapping costs. Does not indicate
/// ivalidity
inline double small_inf() { return 10E10; }

inline bool is_inf(double _val) { return _val > small_inf() / 2.; }

/// Calculate the atomic deformation cost (using child's coordinates)
///
/// The atomic deformation cost of a MappingNode is calculated as the normalized
/// mean-square displacement of its atoms. The displacement vectors are
/// deformed to the child structure's coordinate system before calculating.
///
/// \param mapped_result a proposed mapping; both lattice_node and atomic_node
///     of mapped_result must be initialized
/// \param Nsites number of atoms (excluding vacancies) in the relaxed
///     structure, for proper normalization
///
/// Result is dimensionless, having been normalized by the squared radius of a
/// sphere having the same atomic volume of child structure
double atomic_cost_child(const MappingNode &mapped_result, Index Nsites);

/// Calculate the atomic deformation cost (using parent's coordinates)
///
/// The atomic deformation cost of a MappingNode is calculated as the normalized
/// mean-square displacement of its atoms. The displacement vectors are
/// deformed to the parent structure's coordinate system before calculating.
///
/// \param mapped_result a proposed mapping; both lattice_node and atomic_node
///     of mapped_result must be initialized
/// \param Nsites number of atoms (excluding vacancies) in the relaxed
///     structure, for proper normalization.
///
/// Result is dimensionless, having been normalized by the squared radius of a
/// sphere having the same atomic volume of parent structure
double atomic_cost_parent(const MappingNode &mapped_result, Index Nsites);

/// \brief Calculate the atomic deformation cost as the average of the
/// `atomic_cost_child` and `atomic_cost_parent`
double atomic_cost(const MappingNode &mapped_config, Index Nsites);

/// \brief Calculate the symmetry breaking atomic cost of a MappingNode
double atomic_cost(
    const MappingNode &basic_mapping_node, SymOpVector &factor_group,
    const std::vector<Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic,
                                               Index>> &permutation_group,
    Index Nsites);

}  // namespace StrucMapping

/// \brief Data structure describing a lattice mapping relationship
struct LatticeNode {  // Note: See full description in StrucMapping.cc

  /// \brief The stretch matrix
  ///
  /// This is \f$V^{N}\f$, the stretch matrix, a symmetric matrix that
  /// describes the deformation that maps a de-rotated child superlattice to a
  /// parent superlattice.
  Eigen::Matrix3d stretch;

  /// \brief The isometry (distance-preserving transformation) matrix
  ///
  /// This is \f$Q^{N}\f$, the isometry matrix, which describes a rigid
  /// transformation that de-rotates (de-reflects, etc.) a superlattice of
  /// child. A property of \f$Q^{N}\f$ is that \f$(Q^{N})^{-1} ==
  /// (Q^{N})^{T}\f$.
  Eigen::Matrix3d isometry;

  /// \brief A superlattice of the parent (reference) structure
  ///
  /// This encodes the parent lattice (\f$L_1\f$) and superlattice of parent
  /// that the child is mapped to (\f$L_1 * T_1 * N\f$). Specifically:
  /// - \f$L_1\f$ = `parent.prim_lattice().lat_column_mat()`
  /// - \f$L_1 * T_1 * N\f$ = `parent.superlattice().lat_column_mat()`
  ///
  Superlattice parent;

  /// \brief A superlattice of the mapped and un-deformed child structure
  ///
  /// This encodes the mapped and un-deformed child lattice and the mapped and
  /// un-deformed child superlattice. Specifically:
  /// - \f$V^{N} * Q^{N} * L_2\f$ = `child.prim_lattice().lat_column_mat()`
  /// - \f$V^{N} * Q^{N} * L_2 * T_2\f$ =
  ///   `child.superlattice().lat_column_mat()`
  ///
  /// Note:
  /// - parent.superlattice() == child.superlattice()
  /// - When the parent lattice is the lattice of the primitive structure, then
  ///   it will be the case that \f$T_2 = I\f$.
  Superlattice child;

  /// \brief The lattice deformation cost for this lattice mapping
  double cost;

  /// \brief The name of the method used to calculate the lattice deformation
  /// cost
  std::string cost_method;

  /// Construct LatticeNode, setting all members directly
  LatticeNode(Superlattice parent, Superlattice child, Eigen::Matrix3d stretch,
              Eigen::Matrix3d isometry, double cost, std::string cost_method);

  // TODO: We should prefer putting the logic of how LatticeNode members are
  // calculated in standalone methods rather than constructors...

  /// \brief Construct a LatticeNode by calculating the deformation tensor that
  /// maps a particular child superlattice to a particular parent superlattice
  /// [deprecated]
  LatticeNode(Lattice const &parent_prim, Lattice const &parent_scel,
              Lattice const &unmapped_child_prim,
              Lattice const &unmapped_child_scel, Index child_N_atom,
              double _cost = StrucMapping::big_inf());

  /// \brief Construct a LatticeNode using the mapping calculated by LatticeMap
  /// [deprecated]
  LatticeNode(LatticeMap const &lattice_map, Lattice const &parent_prim,
              Lattice const &unmapped_child_prim);
};

/// \brief Construct a LatticeNode by calculating the deformation tensor that
/// maps a particular child superlattice to a particular parent superlattice
LatticeNode make_lattice_node(Lattice const &parent_prim,
                              Lattice const &parent_scel,
                              Lattice const &unmapped_child_prim,
                              Lattice const &unmapped_child_scel);

/// \brief Construct a LatticeNode using the mapping calculated by LatticeMap
LatticeNode make_lattice_node(LatticeMap const &_lat_map,
                              Lattice const &parent_prim,
                              Lattice const &unmapped_child_prim);

/// \brief Compare two LatticeMap objects, based on their mapping cost first,
/// followed by PrimGrid transformation matrices
bool less(LatticeNode const &A, LatticeNode const &B, double cost_tol);

/// \brief returns true if cost values and parent/child supercell
/// transformations are same for A and B
bool identical(LatticeNode const &A, LatticeNode const &B, double cost_tol);

/// \brief Atomic assignment solution data structure
///
/// AssignmentNode holds the solution of a constrained atomic assignmnet
/// problem. This describes the permutation, translation, and time-reversal of
/// the atoms of a child structure to bring them into registration with the
/// atoms of a parent structure (assuming periodic boundary conditions). Also
/// records the constrained and unconstrained assignment costs
struct AssignmentNode {
  /// Constructor
  ///
  /// \param _cost_tol The cost comparison tolerance
  ///
  AssignmentNode(double _cost_tol = 1e-6)
      : time_reversal(false), cost(0), m_cost_tol(_cost_tol) {}

  /// \brief Mapping translation from child to parent
  ///
  /// Defined such that
  ///
  ///     translation = parent_coord.col(i)-
  ///         child_coord.col(permutation[i]) + displacement.col(i)
  ///
  /// where coordinate comparisons must be made checking for equality up to a
  /// lattice translation.
  ///
  /// This definition assumes that 'child_coord' has been de-rotated and
  /// un-deformed according to a particular lattice mapping (as defined by a
  /// LatticeNode object)
  Eigen::Vector3d translation;

  /// \brief Time reversal relationship between child and parent
  bool time_reversal;

  /// \brief Parent-to-child site assignments that have been forced on.
  ///
  /// Defined such that forced_on.count(el)==1 denotes that
  /// child_coord.col(el.second) maps onto parent_coord.col(el.first).
  std::set<std::pair<Index, Index>> forced_on;

  /// \brief 'Real' indices of rows in the reduced 'cost_mat'
  ///
  /// When a site assignment {i,j} is added to forced_on, row 'i' and col 'j'
  /// are removed from cost_mat. An element cost_mat(k,l) in the cost_mat
  /// corresponds to the original element at (irow[k],icol[l]) in the original
  /// cost_mat.
  std::vector<Index> irow;

  /// \brief 'Real' indices of columns in the reduced 'cost_mat'
  ///
  /// When a site assignment {i,j} is added to forced_on, row 'i' and col 'j'
  /// are removed from cost_mat. An element cost_mat(k,l) in the cost_mat
  /// corresponds to the original element at (irow[k],icol[l]) in the original
  /// cost_mat.
  std::vector<Index> icol;

  /// \brief Solution of the assignment problem for the reduced 'cost_mat'
  ///
  /// An assignment {k,l} in the reduced problem occurs when assignment[k]==l.
  /// In the unreduced problem, this assignment corresponds to
  /// {irow[k],icol[l]}.
  std::vector<Index> assignment;

  /// \brief Cost matrix for an assignment problem
  ///
  /// - If the parent structure allows vacancies, cost_mat may have more
  /// columns than child sites. These correspond to 'virtual vacancies', which
  /// have zero cost of mapping onto a parent site that allows vacancies and
  /// an infinite cost of mapping onto a parent site that does not allow
  /// vacancies.
  /// - This may be a reduced assignment problem. If forced_on.size()>0, then
  /// cost_mat(i,j) is the cost of mapping child site 'j' onto parent atom 'i'.
  /// - In the case of a reduced assignment problem,
  ///   cost_mat.cols()=icol.size() and cost_mat.rows()=irow.size()
  Eigen::MatrixXd cost_mat;

  /// \brief Total cost of best solution to the constrained assignment problem,
  /// possibly having some forced_on assignments
  double cost;

  /// \brief Name of the method used to calculate the atomic displacement cost
  std::string cost_method;

  /// \brief Return cost comparison tolerance
  double cost_tol() const { return m_cost_tol; }

  /// \brief True if cost matrix and assignment vector are uninitialized
  bool empty() const { return cost_mat.size() == 0 && assignment.empty(); }

  /// \brief Compares time_reversal and translation
  ///
  /// - For time_reversal, false is less than true.
  /// - For translation, elements are compared lexicographically
  bool operator<(AssignmentNode const &other) const;

  /// \brief Combines constrained vector HungarianNode::assignment and
  /// HungarianNode::forced_on to obtain total permutation vector.
  ///
  /// The value of child_after.site[i] after assignment is given by
  /// child_before.site[result[i]] before assignment.
  std::vector<Index> permutation() const {
    std::vector<Index> result(assignment.size() + forced_on.size(), 0);
    for (auto const &pair : forced_on) {
      result[pair.first] = pair.second;
    }
    for (Index i = 0; i < assignment.size(); ++i) {
      result[irow[i]] = icol[assignment[i]];
    }
    return result;
  }

 private:
  double m_cost_tol;
};

/// \brief true if time_reversal and translation are identical
bool identical(AssignmentNode const &A, AssignmentNode const &B);

/// \brief Data structure holding a mapping between structures
struct MappingNode {  // Note: See full description in StrucMapping.cc

  // typedefs to provide flexibility if we eventually change to a
  // Eigen::Matrix<3,Eigen::Dynamic>
  typedef Eigen::MatrixXd DisplacementMatrix;

  // Can treat as a Eigen::VectorXd
  using Displacement = DisplacementMatrix::ColXpr;
  using ConstDisplacement = DisplacementMatrix::ConstColXpr;

  /// Label molecules as name and occupant index
  using MoleculeLabel = std::pair<std::string, Index>;

  using AtomIndexSet = std::set<Index>;
  using MoleculeMap = std::vector<AtomIndexSet>;

  /// \brief Contains lattice mapping results
  LatticeNode lattice_node;

  /// \brief Contains atomic assignment results
  AssignmentNode atomic_node;

  /// \brief Component of total cost due to lattice mapping cost
  double lattice_weight;

  /// \brief Component of total cost due to atomic assignment cost
  double atomic_weight;

  /// \brief True if assignment problem is not yet known to be insoluable
  ///
  /// - Used during structure mapping algorithm
  /// - Default: true
  bool is_viable;

  /// \brief True if assignment is valid
  ///
  /// - Specifies that all atoms are assigned to sites where that atom type is
  /// allowed
  /// - Default before checking is false
  mutable bool is_valid;

  /// \brief True if node has been partitioned
  ///
  /// - Used during structure mapping algorithm
  /// - Partitioning into sub-nodes to find sub-optimal assignment solutions is
  /// performed as part of the generalized k-best assignment problem
  /// - Default: false
  mutable bool is_partitioned;

  /// \brief Total solution cost, including lattice and atomic costs
  ///
  /// - Populated by a StrucMapCalculator
  /// - Not guaranteed to be a linear function of lattice_node.cost and
  /// atomic_node.cost
  double cost;

  /// \brief Method used to calculate the total finalized cost
  std::string cost_method;

  /// \brief 3xN matrix of displacements for all sites in parent supercell (Va
  /// are included, but set to Zero)
  Eigen::MatrixXd atom_displacement;

  /// \brief Atom assignment solution as a permutation
  ///
  /// `atom_permutation` lists indices of sites in input structure, as-read, so
  /// that they constitute particular mapping onto parent structure.
  /// If we define index `j = atom_permutation[i]`, this indicates that atom `j`
  /// of the child superstructure maps onto site `i` of parent superstructure.
  /// If the parent structure has `N` sites and child has `M<N` atoms,
  /// vacancies  are designated by values `j>=M`.
  std::vector<Index> atom_permutation;

  /// \brief Mapping between atomic and molecular representation
  ///
  /// mol_map[j] lists atom indices of parent superstructure that comprise the
  /// molecule at its j'th molecular site
  ///
  /// \note Currently only single atom molecules can be mapped, and the general
  /// implementation is for future support.
  MoleculeMap mol_map;

  /// \brief List of assigned molecule names
  ///
  /// Assigned molecule names, as the pair {species_name, occupant_index},
  /// where occupant_index is the index into allowed_species for the parent
  /// structure basis site the molecule is mapped to (index used for ConfigDoF
  /// occupation).
  std::vector<MoleculeLabel> mol_labels;

  /// \brief Static constructor to build an invalid MappingNode
  ///
  /// - Can be used as return value when no valid mapping exists
  static MappingNode invalid();

  /// \brief Construct with lattice node and lattice_weight.
  ///
  /// - Cost is initialized assuming zero atomic_node cost
  MappingNode(LatticeNode _lattice_node, double _lattice_weight)
      : lattice_node(std::move(_lattice_node)),
        is_viable(true),
        is_valid(false),
        is_partitioned(false) {
    set_lattice_weight(_lattice_weight);
    cost = lattice_weight * lattice_node.cost;
  }

  /// \brief Tolerance used for cost comparison
  double cost_tol() const { return atomic_node.cost_tol(); }

  /// \brief Set the lattice_weight
  ///
  /// Cost is calculated as:
  ///
  ///     cost = lattice_weight*lattice_node.cost +
  ///            atomic_weight*atomic_node.cost
  ///
  /// where `lattice_weight` must be on interval `(0.,1.]`, and `atomic_weight`
  /// is `1.-lattice_weight`.
  void set_lattice_weight(double _lw) {
    lattice_weight = max(min(_lw, 1.0), 1e-9);
    atomic_weight = 1. - lattice_weight;
  }

  /// \brief Return pair of integer volumes {Vp, Vc}, where Vp is parent
  /// supercell volume and Vc is child supercell volume
  std::pair<Index, Index> vol_pair() const {
    return std::pair<Index, Index>(lattice_node.parent.size(),
                                   lattice_node.child.size());
  }

  /// \brief Solves the assignment problem
  ///
  /// - Solves the assignment problem using hungarian_method
  /// - Sets \link MappingNode::is_viable is_viable\endling to false if no
  ///   solution
  void calc();

  /// \brief Convenience method to access MappingNode::lattice_node.isometry
  Eigen::Matrix3d const &isometry() const { return lattice_node.isometry; }

  /// \brief Convenience method to access MappingNode::lattice_node.stretch
  Eigen::Matrix3d const &stretch() const { return lattice_node.stretch; }

  /// \brief Convenience method to access MappingNode::atomic_node.translation
  Eigen::Vector3d const &translation() const { return atomic_node.translation; }

  /// \brief Convenience method to access MappingNode::atomic_node.time_reveral
  bool time_reversal() const { return atomic_node.time_reversal; }

  /// \brief Access the i'th atomic displacement of mapped structure
  Displacement disp(Index i) { return atom_displacement.col(i); }

  /// \brief Access the i'th atomic displacement of mapped structure
  ConstDisplacement disp(Index i) const { return atom_displacement.col(i); }

  /// \brief Compares cost, lattice_node, atomic_node, and permutation
  ///
  /// If costs are tied, compares lattice_node.cost first, then:
  /// - if tied, uninitialized atomic_node comes before initialized atomic_node
  /// - if tied, compares lattice_node, using defined comparator
  /// - if tied, compares atomic_node, using defined comparator
  /// - if tied, does lexicographic comparison of permutation
  ///
  /// \note This order is essential for proper behavior of mapping algorithm
  ///
  bool operator<(MappingNode const &other) const;
};

/// \brief External accessor for isometry, to provide xtal::SymOp adaptability
inline Eigen::Matrix3d const &get_matrix(MappingNode const &_node) {
  return _node.isometry();
}

/// \brief External accessor for translation, to provide xtal::SymOp
/// adaptability
inline Eigen::Vector3d const &get_translation(MappingNode const &_node) {
  return _node.translation();
}

/// \brief External accessor for time_reversal, to provide xtal::SymOp
/// adaptability
inline bool get_time_reversal(MappingNode const &_node) {
  return _node.time_reversal();
}

/// \brief Implements a method for mapping a "child" crystal structure as a
/// deformation of a reference "parent" crystal structure.
class StrucMapper {
 public:
  using LatMapType = std::map<Index, std::vector<Lattice>>;

  ///\brief Construct and initialize a StrucMapper
  ///
  ///\param _calculator
  ///\parblock
  ///          specialization of base class StrucMapCalculatorInterface, which
  ///          controls the way that the cost_matrix and costs are calculated,
  ///          determines validity of proposed mappings, and finalizes the
  ///          representation of proposed mappings
  ///\endparblock
  ///
  ///\param _lattice_weight
  ///\parblock
  ///          free parameter 'w' in the cost function: total_cost =
  ///          w*lattice_deformation+(1-w)*atomic_deformation can vary between 0
  ///          (completely basis-focused) and 1 (completely lattice-focused)
  ///\endparblock
  ///
  ///\param _max_volume_change
  ///\parblock
  ///          constrains the search space by assuming a limit on allowed volume
  ///          change only taken into account when non-interstitial vacancies
  ///          are allowed in parent structure
  ///\endparblock
  ///
  ///\param _robust Perform additional checks to determine if mapping is
  /// degenerate in cost to other mappings, which can occur if the imported
  /// structure has symmetry that is incompatible with parent structure.
  /// Results in slower execution.
  ///
  ///\param _soft_va_limit If true, ensures that if no supercell volume
  /// satisfies vacancy constraints, the smallest possible volume is used.
  /// Default behavior results in no valid mapping.
  ///
  ///\param _cost_tol tolerance for mapping comparisons
  ///
  ///\param _min_va_frac minimum fraction of vacant sites, below this fraction a
  /// mapping will not be considered
  ///
  ///\param _max_va_frac maximum fraction of vacant sites, above this fraction a
  /// mapping will not be considered
  StrucMapper(StrucMapCalculatorInterface const &_calculator,
              double _lattice_weight = 0.5, double _max_volume_change = 0.5,
              bool _robust = false, bool _soft_va_limit = false,
              double _cost_tol = TOL, double _min_va_frac = 0.,
              double _max_va_frac = 1.);

  ///\brief Tolerance for determining if two mapping-cost values are identical
  double cost_tol() const { return m_cost_tol; }

  ///\brief Tolerance for initializing lattices. For now it is initialized to
  /// CASM::TOL
  double xtal_tol() const { return m_xtal_tol; }

  double lattice_weight() const { return m_lattice_weight; }

  double atomic_weight() const { return 1. - m_lattice_weight; }

  void set_lattice_weight(double _lw) {
    m_lattice_weight = max(min(_lw, 1.0), 1e-9);
  }

  /// \brief Max element considered for integer unimodular matrix
  /// transformations (which define orientation relationship of mapping)
  Index lattice_transformation_range() const {
    return m_lattice_transformation_range;
  }

  /// \brief Max element considered for integer unimodular matrix
  /// transformations (which define orientation relationship of mapping)
  void set_lattice_transformation_range(Index _new_range) {
    m_lattice_transformation_range = _new_range;
  }

  /// \brief Flag that enables the calculation of a symmetrized lattice cost
  /// when performing the lattice maps. This cost only accounts for that part of
  /// the deformation gradient that breaks the symmetry of the parent crystal
  /// structure
  void set_symmetrize_lattice_cost(bool _sym_lat_cost) {
    m_symmetrize_lattice_cost = _sym_lat_cost;
  }

  /// \brief Flag that enables the calculation of a symmetrized atomic cost
  /// when performing the atomic maps. This cost only accounts for that part of
  /// the displacement field that breaks the symmetry of the parent crystal
  /// structure
  void set_symmetrize_atomic_cost(
      bool _sym_atomic_cost, const SymOpVector &factor_group,
      const std::vector<Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic,
                                                 Index>> &permutation_group) {
    m_symmetrize_atomic_cost = _sym_atomic_cost;
    m_calc_ptr->set_sym_invariant_displacement_modes(
        generate_invariant_shuffle_modes(factor_group, permutation_group));
  }

  bool symmetrize_lattice_cost() const { return m_symmetrize_lattice_cost; }
  bool symmetrize_atomic_cost() const { return m_symmetrize_atomic_cost; }

  /// \brief Returns the minimum fraction of sites allowed to be vacant in the
  /// mapping relation Vacancy fraction is used to constrain the mapping
  /// supercell search, but is only used when the supercell volume cannot is not
  /// uniquely determined by the number of each species in the child structure
  double min_va_frac() const { return m_min_va_frac; }

  /// \brief Sets the minimum fraction of sites allowed to be vacant in the
  /// mapping relation Vacancy fraction is used to constrain the mapping
  /// supercell search, but is only used when the supercell volume cannot is not
  /// uniquely determined by the number of each species in the child structure
  void set_min_va_frac(double _min_va) { m_min_va_frac = max(_min_va, 0.); }

  /// \brief Returns the maximum fraction of sites allowed to be vacant in the
  /// mapping relation Vacancy fraction is used to constrain the mapping
  /// supercell search, but is only used when the supercell volume cannot is not
  /// uniquely determined by the number of each species in the child structure
  double max_va_frac() const { return m_max_va_frac; }

  /// \brief Sets the maximum fraction of sites allowed to be vacant in the
  /// mapping relation Vacancy fraction is used to constrain the mapping
  /// supercell search, but is only used when the supercell volume cannot is not
  /// uniquely determined by the number of each species in the child structure
  void set_max_va_frac(double _max_va) { m_max_va_frac = min(_max_va, 0.99); }

  /// \brief 'robust' option: see constructor for details
  bool robust() const { return m_robust; }

  /// \brief 'soft_va_limit' option: see constructor for details
  bool soft_va_limit() const { return m_soft_va_limit; }

  /// \brief returns reference to parent structure
  SimpleStructure const &parent() const;

  /// \brief specify a superlattice of the parent to be searched during mapping
  ///
  /// Used if:
  /// -
  void add_allowed_lattice(Lattice const &_lat) {
    m_allowed_superlat_map
        [std::abs(round(volume(_lat) / parent().lat_column_mat.determinant()))]
            .push_back(_lat);
  }

  ///\brief clear the list of allowed parent superlattices;
  /// all superlattices will be generated automatically, as needed (default)
  void clear_allowed_lattices() const { m_allowed_superlat_map.clear(); }

  ///\brief returns true if the search of parent superlattices is constrained to
  /// a pre-specified list
  bool lattices_constrained() const { return m_allowed_superlat_map.size(); }

  ///\brief specify to use filtered lattices for mapping. The filter function is
  /// of the form
  ///  bool filter(parent_lattice, proposed_lattice)
  /// where parent_lattice is the primitive lattice of the parent structure, and
  /// proposed_lattice is a proposed superlattice of the parent structure
  void set_filter(StrucMapping::LatticeFilterFunction _filter_f) {
    m_filtered = true;
    m_filter_f = _filter_f;
    m_superlat_map.clear();
  }

  ///\brief specify not to use filtered lattice for mapping
  void unset_filter() {
    m_filtered = false;
    m_superlat_map.clear();
  }

  // --- The `map_X_struc[_impose_Y]` methods are what run the algorithm ---

  /// \brief Find the k-best mappings of a child structure onto the parent
  /// structure, assuming that the child lattice and parent lattice are related
  /// by an integer transformation and a parent structure point group operation
  std::set<MappingNode> map_ideal_struc(
      const SimpleStructure &unmapped_child, Index k = 1,
      double max_cost = StrucMapping::big_inf(), double min_cost = -TOL,
      bool keep_invalid = false) const;

  /// \brief Find the k-best mappings of an arbitrary child structure onto the
  /// parent structure, without simplifying assumptions
  std::set<MappingNode> map_deformed_struc(
      const SimpleStructure &unmapped_child, Index k = 1,
      double max_cost = StrucMapping::big_inf(), double min_cost = -TOL,
      bool keep_invalid = false,
      SymOpVector const &child_factor_group = {SymOp::identity()}) const;

  /// \brief Find the k-best mappings of an arbitrary child structure onto the
  /// parent structure, specifying the range of parent superlattice volumes
  /// considered
  std::set<MappingNode> map_deformed_struc_impose_lattice_vols(
      const SimpleStructure &unmapped_child, Index min_vol, Index max_vol,
      Index k = 1, double max_cost = StrucMapping::big_inf(),
      double min_cost = -TOL, bool keep_invalid = false,
      SymOpVector const &child_factor_group = {SymOp::identity()}) const;

  /// \brief Find the k-best mappings of an arbitrary child structure onto the
  /// parent structure, specifying the parent superlattice exactly, but not the
  /// way the child lattice maps to the parent superlattice
  std::set<MappingNode> map_deformed_struc_impose_lattice(
      const SimpleStructure &unmapped_child, const Lattice &imposed_lat,
      Index k = 1, double max_cost = StrucMapping::big_inf(),
      double min_cost = -TOL, bool keep_invalid = false,
      SymOpVector const &child_factor_group = {SymOp::identity()}) const;

  /// \brief Find the k-best mappings of an arbitrary child structure onto the
  /// parent structure, specifying the lattice mapping exactly
  std::set<MappingNode> map_deformed_struc_impose_lattice_node(
      const SimpleStructure &unmapped_child, const LatticeNode &imposed_node,
      Index k = 1, double max_cost = StrucMapping::big_inf(),
      double min_cost = -TOL, bool keep_invalid = false) const;

  /// Find k-best mappings
  Index k_best_maps_better_than(SimpleStructure const &unmapped_child,
                                std::set<MappingNode> &queue, Index k = 1,
                                double max_cost = StrucMapping::big_inf(),
                                double min_cost = -TOL,
                                bool keep_invalid = false,
                                bool keep_tail = false,
                                bool no_partiton = false) const;

  StrucMapCalculatorInterface const &calculator() const { return *m_calc_ptr; }

 private:
  /// \brief Generate a set of mapping seeds (i.e., MappingNode with
  /// LatticeNode only) from a list of supercells of the parent structure and
  /// a list of supercells of the child structure
  std::set<MappingNode> _seed_k_best_from_super_lats(
      SimpleStructure const &unmapped_child,
      std::vector<Lattice> const &_parent_scels,
      std::vector<Lattice> const &_child_scels, Index k, double max_strain_cost,
      double min_strain_cost,
      SymOpVector const &child_factor_group = {SymOp::identity()}) const;

  /// \brief construct partial mapping nodes (with uninitialized atomic_node)
  /// based on current settings considers supercells with integer volume between
  /// min_vol and max_vol
  ///
  /// Note:
  std::set<MappingNode> _seed_from_vol_range(
      SimpleStructure const &unmapped_child, Index k, Index min_vol,
      Index max_vol, double max_strain_cost, double min_strain_cost,
      SymOpVector const &child_factor_group = {SymOp::identity()}) const;

  ///\brief returns number of species in a SimpleStructure given the current
  /// calculator settings.
  ///       Use instead of sstruc.n_atom() for consistency
  Index _n_species(SimpleStructure const &sstruc) const;

  // Check lattice filter function
  bool _filter_lat(Lattice const &_parent_lat,
                   Lattice const &_child_lat) const {
    return m_filter_f(_parent_lat, _child_lat);
  }

  /// \brief Calculate min_vol, max_vol range from min_va_frac, max_va_frac, and
  /// properties of the child and parent structures and calculator.
  std::pair<Index, Index> _vol_range(
      const SimpleStructure &unmapped_child) const;

  notstd::cloneable_ptr<StrucMapCalculatorInterface> m_calc_ptr;

  double m_lattice_weight;
  double m_max_volume_change;
  bool m_robust;
  bool m_soft_va_limit;
  double m_cost_tol;
  double m_xtal_tol;
  double m_min_va_frac;
  double m_max_va_frac;

  Index m_lattice_transformation_range;

  bool m_symmetrize_lattice_cost;
  bool m_symmetrize_atomic_cost;

  bool m_filtered;
  StrucMapping::LatticeFilterFunction m_filter_f;

  /// Maps the supercell volume to a vector of Lattices with that volume
  mutable LatMapType m_superlat_map;
  mutable LatMapType m_allowed_superlat_map;

  std::vector<Lattice> _lattices_of_vol(Index prim_vol) const;
};

}  // namespace xtal
}  // namespace CASM
#endif
