#include "casm/clex/ConfigMapping.hh"

#include <vector>

#include "casm/app/ProjectSettings.hh"
#include "casm/app/QueryHandler_impl.hh"
#include "casm/basis_set/DoFTraits.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/dataformatter/DataStream.hh"
#include "casm/clex/ConfigDoFTools.hh"
#include "casm/clex/ConfigIsEquivalent.hh"
#include "casm/clex/Configuration.hh"
#include "casm/clex/FillSupercell.hh"
#include "casm/clex/MappedProperties.hh"
#include "casm/clex/MappedPropertiesTools.hh"
#include "casm/clex/ParamComposition.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/SimpleStructureTools.hh"
#include "casm/clex/Supercell.hh"
#include "casm/completer/Handlers.hh"
#include "casm/crystallography/Adapter.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/CanonicalForm.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/LatticeMap.hh"
#include "casm/crystallography/Niggli.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/database/ScelDatabase.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/strain/StrainConverter.hh"
#include "casm/symmetry/PermuteIterator.hh"

namespace CASM {

/// \class ConfigMapper
/// \brief Implements a method for mapping structures to configurations
/// \ingroup Configuration
///
/// The ConfigMapper class implements a method for mapping an arbitrary crystal
/// structure (xtal::SimpleStructure, called the "child" structure), to a
/// configuration (Configuration) of a crystal template (xtal::BasicStructure,
/// which is used to create the "parent" reference structure) and a set of
/// properties (MappedProperties).
///
///
/// Examples
/// --------
///
/// ConfigMapper maps a "child" xtal::SimpleStructure to a "parent" prim
/// xtal::BasicStructure:
///
///     // structure being mapped
///     xtal::SimpleStructure unmapped_child = ...
///
///     // structure being mapped to
///     std::shared_ptr<xtal::BasicStructure const> shared_prim = ...
///
///
/// ### Example 1: Mapping with "auto" lattice volume range
///
/// When mapping to prim structures that allow vacancies, the volume of the
/// mapped configuration is not trivially determinable from the number of atoms
/// in the input structure. In this example, the \link
/// ConfigMapperSettings::auto_lattice_volume_range
/// "auto_lattice_volume_range"\endlink lattice mapping method is used to
/// constrain the configuration mapping solutions by specifying the allowed Va
/// concentation range in the solution configurations, and overall volume
/// change for the supercell relative the input structure. This is the most
/// general lattice mapping approach.
///
///     ConfigMapperSettings config_mapper_settings;
///     config_mapper_settings.auto_lattice_volume_range = true;
///     config_mapper_settings.min_va_frac = 0.0;       // default = 0.0;
///     config_mapper_settings.max_va_frac = 0.05;      // default = 1.0;
///     config_mapper_settings.max_volume_change = 0.3; // default = 0.3;
///
///     ConfigMapper mapper{shared_prim, config_mapper_settings};
///     ConfigMapperResult result = mapper.import_structure(unmapped_child);
///
///
/// ### Example 2: Mapping to a configuration with a given lattice
///
/// In this example, the \link ConfigMapperSettings::fix_lattice
/// "fix_lattice"\endlink lattice mapping method is used to constrain the
/// configuration mapping solutions to supercells with a lattice that is
/// equivalent to a specified superlattice of the prim structure. The resulting
/// solutions may not have identical lattice vectors as the specified lattice,
/// but they will have equivalent lattices (same lattice points).
///
///     ConfigMapperSettings config_mapper_settings;
///     config_mapper_settings.fix_lattice = true;
///
///     Eigen::Vector3d a, b, c;
///     std::tie(a, b, c) = shared_prim->lattice().vectors();
///     Eigen::Matrix3d S; // supercell lattice column vector matrix
///     S.col(0) = 3*a;
///     S.col(1) = 3*b + a;
///     S.col(2) = 1*c;
///     config_mapper_settings.configuration_lattice = Lattice{S};
///
///     ConfigMapper mapper{shared_prim, config_mapper_settings};
///     ConfigMapperResult result = mapper.import_structure(unmapped_child);
///
/// ### Example 3: Printing results as JSON.
///
/// The ConfigMapperResult data structure holds mapping solutions. It can be
/// printed as a JSON object.
///
///     jsonParser json;
///     json.put_obj();
///     to_json(unmapped_child, json["unmapped_child"], {}, CART);
///     std::cout << json << std::endl;
///
///     json.put_obj();
///     to_json(result, json, CART);
///     std::cout << json << std::endl;
///
///
/// Configuration mapping method summary
/// ------------------------------------
///
/// Configuration mapping steps:
///
/// 1) Generate and score potential mappings of the child structure to the
/// reference parent structure. For configuration mapping, the parent structure
/// (xtal::SimpleStructure) is constructed from a prim (xtal::BasicStructure)
/// which also provides the list of allowed occupants on each site. Lattice
/// mapping is done first, then atomic assignment for each potential lattice
/// mapping solution. A variety of settings control which lattice mappings are
/// considered, how many solutions are stored, and how they are scored. See
/// ConfigMapperSettings for configuration mapping options and
/// xtal::StrucMapper for structure mapping details.
///
/// 2) The best scoring mappings are used to create "mapped structures",
/// which are de-rotated according to the lattice mapping solution, expressed
/// in terms of molecules, ordered according to the atomic assignment solution
/// and to match the ordering of sites in the appropriate supercell, and have
/// displacements, strains, and other DoF and properties appropriately de-
/// rotated, permuted, scaled (if extensive), and expressed in terms of the
/// reference structure, so that mapped configurations and properties may be
/// constructed.
///
///
/// Configuration mapping results
/// -----------------------------
///
/// The ConfigMapperResult data structure holds the configuration mapping
/// results. It contains:
///
/// - \link ConfigurationMapping std::set<ConfigurationMapping>\endlink \link
/// ConfigMapperResult::maps maps\endlink: The set of mapping solutions, sorted
/// by increasing mapping cost. The ConfigurationMapping data structure,
/// summarized below, holds the details of each mapping solution.
/// - std::string \link ConfigMapperResult::fail_msg fail_msg\endlink: Holds
/// messages describing why the mapping algorithm failed if no mapping solutions
/// were found.
///
/// The ConfigurationMapping data structure contains:
/// - xtal::SimpleStructure \link ConfigurationMapping::mapped_child
/// mapped_child\endlink: Resulting mapped structure, equivalent to the input
/// structure, but mapped to the parent by means of rigid transformation,
/// permutation, and translation. It is *not* deformed to match the lattice of
/// the reference structure it is mapped to. For the mapped child structure:
///   - Atoms and inferred vacancies are resolved as molecules
///   - Molecule displacements are set to the average of the coordinates of
///   atoms in the molecule, and then molecules are translated such that average
///   displacement of molecules relative to the ideal sites is 0
///   - Properties are transformed to match the orientation of the mapped
///   structure
///   - Coordinates and properties are permuted so that they match the order of
///   sites in the supercell with lattice vectors equal to those of the mapped
///   structure
/// - xtal::MappingNode \link ConfigurationMapping::mapping mapping\endlink:
/// Holds the transformation (rigid transformation, deformation, permutation,
/// translation) that maps the unmapped child to the mapped child.
/// - Configuration \link ConfigurationMapping::mapped_configuration
/// mapped_configuration\endlink: Constructed with the properties of mapped
/// child that correspond to prim DoF.
/// - Configuration \link ConfigurationMapping::mapped_properties
/// mapped_properties\endlink: Constructed with properties of mapped child that
/// do not correspond to prim DoF.
/// - Additional data describing the configuration and properties mapped to the
/// canonical equivalent supercell, and in an equivalent configuration as
/// selected according to the ConfigMapperSettings.
///
///
/// Relations between the unmapped and mapped structures
/// ----------------------------------------------------
///
/// \note In the following, objects referred to are either members of
/// ConfigMapper (`struc_mapper()`, `shared_prim()`), the input structure being
/// mapped (`unmapped_child`), or members of a ConfigurationMapping solution
/// (`mapping`, `mapped_child`, `mapped_configuration`, `mapped_properties`).
///
/// \note The current implementation is written with planned future support for
/// resolution of atomic properties to molecular properties in mind (for
/// example, resolving atoms in the input structure to be molecules in the
/// mapped structure, and averaging displacements of atoms that make up an
/// inferred molecule to calculate molecular displacements), but currently only
/// structures with single-atom molecules can be mapped.
///
/// ### Lattice vector mapping relations
///
/// When a valid mapping solution is found, the following relations hold between
/// the mapped and unmapped lattice vectors:
///
/// \f[
///     L_1 * T_1 * N = V * Q * L_2,
/// \f]
///
/// - \f$L_1\f$: The prim lattice, as a column vector matrix. It is equal to
///
///       shared_prim()->lattice().lat_column_mat()
///
/// - \f$T_1\f$: An integer transformation matrix
/// - \f$N\f$: A unimodular matrix, generates non-identical but equivalent
///   parent superlattices, to be determined
/// - \f$L_2\f$: The unmapped child structure's lattice vectors, as a column
///   vector matrix. It is equal to
///
///       unmapped_child.lat_column_mat \endcode
///
/// - \f$V\f$, \f$Q\f$: The stretch (deformation) and isometry (rigid
///   transformation) matrices, stored in the mapping results as:
///
///   - \code mapping.lattice_node.stretch \endcode
///   - \code mapping.lattice_node.isometry \endcode
///
/// Thus:
/// - \f$L_1 * T_1 * N\f$: The supercell lattice, as a column vector matrix, of
///   the mapped configuration. The following are all equal to
///   \f$L_1 * T_1 * N\f$:
///
///   - \code
///     mapping.lattice_node.parent.superlattice().lat_column_mat()
///     \endcode
///   - \code
///     mapped_configuration.supercell().lattice().lat_column_mat()
///     \endcode
///   - \code mapped_child.lat_column_mat \endcode
///
///
///  ### Atomic coordinate and occupation mapping relations
///
/// The following relations hold between the mapped and unmapped atomic
/// coordinates:
///
/// \f[
///     \vec{r_1}(i) + \vec{d}(i) = V * Q * \vec{r_2}(p_i) + \vec{t}
/// \f]
///
/// - \f$\vec{r_1}(i)\f$: the Cartesian coordinate of the i-th site in the
/// parent superstructure. It is equal to:
///
///       xtal::SimpleStructure parent_superstructure = make_superstructure(
///           mapping.lattice_node.parent.transformation_matrix_to_super(),
///           config_mapper.struc_mapper().parent());
///       parent_superstructure.atom_info.coords.col(i);
///
/// - \f$\vec{d}(i)\f$: displacement associated with the molecule at the i-th
/// site in parent superstructure. If displacement ("disp") is a DoF of the
/// prim, then it is equal to:
///
///       mapped_configuration.configdof().local_dof("disp").standard_values().col(i)
///
///   Otherwise, displacement is a property and equal to:
///
///       mapped_properties.mol_info.properties["disp"].col(i)
///
/// - \f$\vec{r_2}(i)\f$: the Cartesian coordinate of the i-th site in the
/// child structure. It is equal to:
///
///       unmapped_child.atom_info.coords.col(i)
///
/// - \f$p_i\f$: A permutation vector, describes which atom in the unmapped
///   structure (\f$p_i\f$) is mapped to the i-th site of the mapped structure.
///   Values of \f$p_i\f$ greater than the number of atoms in the unmapped
///   structure indicate inferred vacancies. It is equal to:
///
///       mapping.atom_permutation
///
/// - \f$\vec{t}\f$: A translation vector. It is equal to:
///
///       mapping.atomic_node.translation
///
///
/// ### Global DoF and property mapping relations
///
/// Global properties, \f$v^{global}\f$, represented as vectors
/// (Eigen::VectorXd) in the \link xtal::DoFSet "standard basis"\endlink, are
/// related according to
///
/// \f[
///     v^{global, mapped} = M * v^{global, unmapped},
/// \f]
///
/// - \f$v^{global, unmapped}\f$: A global property of the child structure
///   being mapped. It is equal to
///
///       unmapped_child.properties[property_name]
///
/// - \f$v^{global, mapped}\f$: A global DoF or property of the result. It is
///   equal to
///
///       mapped_child.properties[property_name]
///
///   If it is a DoF of the prim, then it is also equal to
///
///       mapped_configuration.configdof().global_dof(property_name)
///           .standard_values()
///
///   Otherwise, it is a property and it is also equal to
///
///       mapped_properties.global[property_name]
///
/// - \f$M\f$: Matrix representation for the mapping transformation. It is
///   equal to
///
///       Eigen::MatrixXd M = AnisoValTraits(property_name).symop_to_matrix(
///           mapping.lattice_node.isometry,
///           mapping.lattice_node.stretch.inverse() *
///           mapping.atomic_node.translation,
///           mapping.atomic_node.time_reversal);
///
///
/// ### Local DoF and property mapping relations
///
/// Site properties, \f$v_i^{local}\f$, represented as vectors (columns of
/// Eigen::MatrixXd) in the \link xtal::DoFSet "standard basis"\endlink, are
/// transformed, and permuted according to
///
/// \f[
///     v_{i}^{local, mapped} = M * v_{p_{i}}^{local, unmapped},
/// \f]
///
/// - \f$v_{p_{i}}^{local, unmapped}\f$: An atomic property of the child
/// structure being mapped. \f$p_i\f$ is the atom permutation as described
/// previously. It is equal to
///
///       unmapped_child.atom_info.properties[property_name]
///           .col(mapping.atom_permutation[i])
///
/// - \f$v^{local, mapped}\f$: A molecular property of the resulting mapped
/// child structure. It is equal to
///
///       mapped_child.mol_info.properties[property_name].col(i)
///
/// - \f$M\f$: Matrix representation for the mapping transformation. It is
/// equal to
///
///       Eigen::MatrixXd M = AnisoValTraits(property_name).symop_to_matrix(
///           mapping.lattice_node.isometry,
///           mapping.lattice_node.stretch.inverse() *
///           mapping.atomic_node.translation,
///           mapping.atomic_node.time_reversal);
///
///

/// \class ConfigMapperSettings
/// \brief Settings for ConfigMapper
///
/// ConfigMapperSettings collects settings for the configuration mapping
/// method implemented by ConfigMapper.
///
/// There are 5 choices for the lattice mapping portion of the mapping method.
/// The choice of one of these methods determines which other parameters have
/// an effect. One and only one must be chosen. The methods are:
///
/// 1. fix_ideal = true: Force lattice mapping solutions of the form
///     \f$L_1 * T_1 = Q * L_2\f$, where
///
///   - \f$L_1\f$ = the prim lattice, as a column vector matrix
///   - \f$L_2\f$ = `child.lat_column_mat`, the unmapped child structure's
///     lattice vectors, as a column vector matrix
///   - \f$T_1\f$ = an integer transformation matrix, to be determined
///   - \f$Q\f$ = isometry matrix, constrained to be an elemenet of the prim
///   factor group, to be determined
///
/// 2. fix_lattice_mapping = true: Force lattice mapping solutions of the form
///     \f$S_1 = V * Q * L_2\f$, where
///
///   - \f$S_1\f$ = `this->configuration_lattice.value()`, lattice vectors, as
///     a column vector matrix, of an exact supercell of the prim lattice
///     vectors
///   - \f$L_2\f$ = `child.lat_column_mat`, the unmapped child structure's
///     lattice vectors, as a column vector matrix
///   - \f$V\f$, \f$Q\f$ = stretch and isometry matrices, to be determined
///
/// 3. fix_lattice = true: Force lattice mapping solutions of the form
///     \f$S_1 * N = V * Q * L_2\f$, where
///
///   - \f$S_1\f$ = `this->configuration_lattice.value()`, lattice vectors, as
///     a column vector matrix, of an exact supercell of the prim lattice
///     vectors
///   - \f$L_2\f$ = `child.lat_column_mat`, the unmapped child structure's
///     lattice vectors, as a column vector matrix
///   - \f$N\f$ = unimodular matrix, generates non-identical but equivalent
///     parent superlattices, to be determined
///   - \f$V\f$, \f$Q\f$ = stretch and isometry matrices, to be determined
///
/// 4. fix_lattice_volume_range = true: Force lattice mapping solutions of the
///     form \f$L_1 * T_1 * N = V * Q * L_2\f$, where
///
///         this->min_lattice_volume <=
///             (L_1 * T_1 * N).determinant() <= this->max_lattice_volume,
///
///   - \f$L_1\f$ = the prim lattice, as a column vector matrix
///   - \f$L_2\f$ = `child.lat_column_mat`, the unmapped child structure's
///     lattice vectors, as a column vector matrix
///   - \f$T_1\f$ = an integer transformation matrix, to be determined
///   - \f$L_2\f$ = `child.lat_column_mat`, the unmapped child structure's
///     lattice vectors, as a column vector matrix
///   - \f$N\f$ = unimodular matrix, generates non-identical but equivalent
///     parent superlattices, to be determined
///   - \f$V\f$, \f$Q\f$ = stretch and isometry matrices, to be determined
///
/// 5. auto_lattice_volume_range = true: The most general mapping approach, it
///   is equivalent `fix_lattice_volume_range` with a volume range that is
///   determined by constaints on the Va frac and change in volume between the
///   input structure and mapped structure.
///
/// If the mapping is succesful, one or more mapping solutions will be found
/// and the unmapped input structure will be mapped to the prim, creating a
/// `mapped_child` (xtal::Simpletructure), `mapped_configuration`
/// (Configuration), and `mapped_properties` (MappedProperties). The
/// combination of `mapped_configuration` and `mapped_properties` is equivalent
/// to the `mapped_child`.
///
/// As a post-processing step, the `mapped_configuration` is transformed to the
/// `final_configuration` which is (a) always in the canonical supercell and (b)
/// either:
/// 1. finalize_strict = true: (the 'strict' method) the configuration which
/// preserves the orientation of the original unmapped structure as much as
/// possible.
/// 2. otherwise, the finalized configuration is the canonical equivalent
/// configuration in the canonical supercell
///

namespace Local {
static int _permute_dist(xtal::MappingNode::MoleculeMap const &_perm,
                         xtal::MappingNode const &_node) {
  int result(0);
  for (int i = 0; i < _perm.size(); ++i) {
    result += std::abs(int(_node.atom_permutation[*(_perm[i].begin())]) - i);
  }
  return result;
}

//*******************************************************************************************

/// \brief Find symop (as PermuteIterator) that gives the most 'faithful'
/// equivalent mapping This means that
///   (1) the site permutation is as close to identity as possible (i.e.,
///   maximal character) (2) ties at (1) are broken by ensuring _node.isometry
///   is proper and close to zero rotation (i.e., maximal character*determinant)
///   (3) if (1) and (2) are ties, then we minimize _node.translation.norm()

template <typename IterType>
static IterType _strictest_equivalent(IterType begin, IterType end,
                                      xtal::MappingNode const &_node) {
  SymOp op(_node.isometry(), _node.translation(), false, _node.cost_tol());
  SymOp t_op;
  IterType best_it = begin;

  double tchar, tdist;
  int tdet, tpdist;
  double best_char = op.matrix().trace();
  double best_dist = op.tau().norm();
  int best_det = sgn(round(op.matrix().determinant()));
  int best_pdist = _permute_dist(_node.mol_map, _node);

  Coordinate tau(_node.lattice_node.parent.superlattice());
  while (begin != end) {
    t_op = begin->sym_op() * op;
    tdet = sgn(round(t_op.matrix().determinant()));
    bool skip_fg = false;
    if (tdet > best_det) {
      best_det = tdet;
      best_char = tdet * t_op.matrix().trace();
      best_pdist =
          _permute_dist(begin->combined_permute() * _node.mol_map, _node);
      tau.cart() = t_op.tau();
      tau.voronoi_within();
      best_dist = tau.const_cart().norm();
      best_it = begin;
    } else if (tdet == best_det) {
      tchar = tdet * t_op.matrix().trace();
      if (almost_equal(tchar, best_char)) {
        tpdist =
            _permute_dist(begin->combined_permute() * _node.mol_map, _node);
        if (tpdist > best_pdist) {
          best_det = tdet;
          best_char = tchar;
          best_pdist = tpdist;
          tau.cart() = t_op.tau();
          tau.voronoi_within();
          best_dist = tau.const_cart().norm();
          best_it = begin;
        } else if (tpdist == best_pdist) {
          tau.voronoi_within();
          tdist = tau.const_cart().norm();
          if (tdist < best_dist) {
            best_det = tdet;
            best_char = tchar;
            best_pdist = tpdist;
            best_dist = tdist;
            best_it = begin;
          }
        }
      } else if (tchar > best_char) {
        best_det = tdet;
        best_char = tchar;
        best_pdist =
            _permute_dist(begin->combined_permute() * _node.mol_map, _node);
        tau.cart() = t_op.tau();
        tau.voronoi_within();
        best_dist = tau.const_cart().norm();
        best_it = begin;
      } else {
        skip_fg = true;
      }
    } else {
      skip_fg = true;
    }

    if (skip_fg) {
      begin = begin_next_fg_op(begin, end);
    } else {
      ++begin;
    }
  }
  return best_it;
}
}  // namespace Local

namespace {

void check_equal(bool A, bool B, std::string message) {
  if (A != B) {
    throw std::runtime_error(message);
  }
}
void check_equal(Eigen::MatrixXd const &A, Eigen::MatrixXd const &B,
                 std::string message) {
  if (!almost_equal(A, B)) {
    throw std::runtime_error(message);
  }
}

Eigen::MatrixXd _scalar_as_matrix(double value) {
  Eigen::MatrixXd M(1, 1);
  M << value;
  return M;
}

/// \brief Reorders the permutation and compounds the spatial isometry (rotation
/// + translation) of _node with that of _it
xtal::MappingNode copy_apply(PermuteIterator const &_it,
                             xtal::MappingNode const &_node,
                             bool transform_cost_mat = true);

/// \brief Initializes configdof of Supercell '_scel' corresponding to an
/// idealized child structure (encoded by _child_struc) _child_struc is assumed
/// to have been idealized via structure-mapping or to be the result of
/// converting a configuration to a SimpleStructure. result.second gives list of
/// properties that were utilized in the course of building the configdof
std::pair<ConfigDoF, std::set<std::string> > to_configdof(
    SimpleStructure const &_child_struc, Supercell const &_scel);

xtal::StrucMapping::AllowedSpecies make_allowed_species(
    BasicStructure const &_prim, SimpleStructure::SpeciesMode _species_mode =
                                     SimpleStructure::SpeciesMode::ATOM) {
  xtal::StrucMapping::AllowedSpecies result(_prim.basis().size());
  Index i = 0;
  for (Site const &site : _prim.basis()) {
    for (Molecule const &mol : site.occupant_dof()) {
      if (_species_mode == SimpleStructure::SpeciesMode::MOL) {
        result[i].push_back(mol.name());
      } else if (_species_mode == SimpleStructure::SpeciesMode::ATOM) {
        if (mol.size() != 1) {
          throw std::runtime_error(
              "make_allowed_species may only be called on "
              "structures with single-atom species.");
        }
        result[i].push_back(mol.atom(0).name());
      }
    }
    ++i;
  }

  return result;
}

/// \brief Reorders the permutation and compounds the spatial isometry (rotation
/// + translation) of _node with that of _it
xtal::MappingNode copy_apply(PermuteIterator const &_it,
                             xtal::MappingNode const &_node,
                             bool transform_cost_mat) {
  xtal::MappingNode result(_node);
  SymOp op = _it.sym_op();

  // Apply symmetry to LatticeNode:
  // lattice_node.parent and lattice_node.child remain unchanged, but isometry
  // and stretch tensors are augmented
  //   parent superlattice Lp and child superlattice Lc are related via
  //      Lp = U * R * Lc  (where R is isometry and U is right stretch tensor of
  //      _node)
  //   'op' must be in point group of Lp, and invariance relation for Lp is
  //      op.matrix() * Lp * N = Lp  (where 'N' is integer fractional symop
  //      matrix)
  //   Substituting above expression and inserting Identity =
  //   op.matrix().transpose() * op.matrix() yields:
  //      Lp = (op.matrix() * U * op.matrix().transpose()) * (op.matrix() * R) *
  //      Lc * N
  //   where terms are grouped to reveal the transformation rules for matrices
  //   'U' and 'R' for the symmetry-related mapping The fractional matrix 'N' is
  //   not used in this routine, since _node records the child lattice in its
  //   undeformed state and all atomic coordinates are recorded in undeformed
  //   cartesian coordinates
  result.lattice_node.isometry = op.matrix() * _node.isometry();
  result.lattice_node.stretch =
      op.matrix() * _node.stretch() * op.matrix().transpose();

  // Apply symmetry to HungarianNode:
  // parent coordinates Cp (3xN) and child coordinates Cc (3xN) are related via
  // the mapping relation
  //    Cp = U * R * Cc * P.transpose() - D + T)]
  // where R is isometry, U is right stretch tensor, P is site permutation,
  // D is displacement field (3xN) and T is mapping translation (3x1 repeated
  // for N columns) 'op' must be in space group of Cp, and invariance relation
  // for Cp is
  //    Cp = op.matrix() * Cp * Ps.transpose() + op.tau()
  // where Ps = _it.combined_permute().matrix() and op.tau() is added col-wise.
  // The invariance is only valid to within a lattice translation of the basis
  // sites (which does not affect mapping score) Inserting the mapping relation
  // for Cp into the invariance relation yields
  //    Cp = [op.matrix() * U * op.matrix.transpose()] * [op.matrix * R] * Cc *
  //    [P.transpose() * Ps.transpose()] - [op.matrix() * D * Ps] + [op.matrix()
  //    * T + op.tau()]
  // where terms are grouped to reveal the transformation rules for mapping
  // permutation 'P', displacement field 'D', and mapping translation 'T' The
  // transormation rules for 'U' and 'R' are identical to those above and are
  // recorded in _node.lattice_node. These five translation rules specify all
  // considerations necessary to describe how a symop in the space group of the
  // parent crystal relates a mapping of the child crystal onto the parent
  // crystal to an equivalent mapping
  result.atomic_node.translation = op.matrix() * _node.translation() + op.tau();
  result.atom_displacement = op.matrix() * result.atom_displacement;
  for (Index i = 0; i < result.mol_map.size(); ++i) {
    result.mol_map[i] = _node.mol_map[_it.permute_ind(i)];
    result.mol_labels[i] = _node.mol_labels[_it.permute_ind(i)];
    // result.atom_displacement.col(i) = td.col(_it.permute_ind(i));
  }

  // Attempt to transfrom the constrained mapping problem and cost matrix
  // This is distinct from the relations described above, as the asignments are
  // not yet all known Instead of permuting the child indices, we will permute
  // the parent indices by the inverse permutation, which should have the same
  // effect in the end
  if (false) {  // Disabled due to changes related to molecule
    PermuteIterator inv_it = _it.inverse();
    for (Index i = 0; i < result.atomic_node.irow.size(); ++i)
      result.atomic_node.irow[i] =
          inv_it.permute_ind(_node.atomic_node.irow[i]);

    result.atomic_node.forced_on.clear();
    for (auto const &el : _node.atomic_node.forced_on)
      result.atomic_node.forced_on.emplace(inv_it.permute_ind(el.first),
                                           el.second);
    // No need to transform assignment vector or cost_mat, which are in terms of
    // the nominal indexing
  }
  return result;
}

xtal::SimpleStructure make_mapped_structure(
    xtal::MappingNode const &mapping_node,
    SimpleStructure const &unmapped_child,
    std::shared_ptr<Structure const> const &shared_prim) {
  throw std::runtime_error("TODO: make_mapped_structure");
}

/// \brief Copy child properties that are Prim DoF to create a Configuration
///
/// \param mapping_node Structure mapping information. The supercell of the
///     resulting configuration has lattice vectors equal to those of
///     `mapping_node.lattice_node.parent.superlattice()`.
/// \param mapped_child Child structure, already transformed so that it is
///     mapped to the parent superstructure. This means that
///     `mapped_child.mol_info.names[i]` specifies the molecule at site index i
///     in the resulting configuration, and for each prim DoF the
///     `DoFType::traits(dof_key).find_values` implementation is used to
///     construct DoF values from `mapped_child.properties` and
///     `mapped_child.mol_info.properties` without any further transformation.
/// \param shared_prim Shared prim structure that the child was mapped to.
///
Configuration make_mapped_configuration(
    xtal::MappingNode const &mapping_node, SimpleStructure const &mapped_child,
    std::shared_ptr<Structure const> const &shared_prim) {
  auto const &prim = *shared_prim;
  std::shared_ptr<Supercell const> mapped_shared_supercell =
      std::make_shared<Supercell const>(
          shared_prim, mapping_node.lattice_node.parent.superlattice());
  SimpleStructure::Info const &child_mol_info = mapped_child.mol_info;

  ConfigDoF configdof =
      make_configdof(*mapped_shared_supercell, prim.lattice().tol());

  Index i = 0;
  for (Index b = 0; b < prim.basis().size(); ++b) {
    for (Index l = 0; l < mapped_shared_supercell->volume(); ++l, ++i) {
      Index j = 0;
      for (; j < prim.basis()[b].occupant_dof().size(); ++j) {
        if (child_mol_info.names[i] == prim.basis()[b].occupant_dof()[j]) {
          configdof.occ(i) = j;
          break;
        }
      }
      if (j == prim.structure().basis()[b].occupant_dof().size())
        throw std::runtime_error(
            "Attempting to initialize ConfigDoF from SimpleStructure. Species "
            "'" +
            child_mol_info.names[i] + "' is not allowed on sublattice " +
            std::to_string(b));
    }
  }

  // note: the `find_values` methods throw if property is missing. is that the
  // right behavior?

  for (auto const &dof_key : xtal::global_dof_types(prim)) {
    auto val = DoFType::traits(dof_key).find_values(mapped_child.properties);
    configdof.global_dof(dof_key).from_standard_values(val.first);
  }

  for (auto const &dof_key : xtal::continuous_local_dof_types(prim)) {
    auto val = DoFType::traits(dof_key).find_values(child_mol_info.properties);
    configdof.local_dof(dof_key).from_standard_values(val.first);
  }

  return Configuration{mapped_shared_supercell, configdof};
}

}  // namespace

Index ConfigMapperResult::n_optimal(double tol /*=TOL*/) const {
  Index result = 0;
  auto it = maps.begin();
  while (it != maps.end() &&
         almost_equal(it->mapping.cost, maps.begin()->mapping.cost, tol)) {
    ++result;
    ++it;
  }
  return result;
}

ConfigMapper::ConfigMapper(std::shared_ptr<Structure const> const &_shared_prim,
                           ConfigMapperSettings const &_settings)
    : m_shared_prim(_shared_prim),
      m_struc_mapper(
          PrimStrucMapCalculator(
              *shared_prim(), adapter::Adapter<xtal::SymOpVector, SymGroup>()(
                                  shared_prim()->factor_group())),
          _settings.lattice_weight, _settings.max_volume_change,
          _settings.robust, _settings.soft_va_limit, _settings.cost_tol,
          _settings.min_va_frac, _settings.max_va_frac),
      m_settings(_settings) {
  if (settings().filter) {
    m_struc_mapper.set_filter(settings().filter);
  }

  for (Lattice const &lattice : settings().allowed_lattices) {
    m_struc_mapper.add_allowed_lattice(lattice);
  }
}

ConfigMapperResult ConfigMapper::import_structure(
    SimpleStructure const &child_struc, Configuration const *hint_ptr,
    std::vector<DoFKey> const &_hint_dofs) const {
  ConfigMapperResult result;
  double best_cost = xtal::StrucMapping::big_inf();
  Index k_best = settings().k_best;

  // If hint configuration is provided, score the mapping with parent=hint (w/
  // hint_dofs), child=child_struc
  double hint_cost;
  if (hint_ptr != nullptr) {
    xtal::StrucMapper tmapper(
        *struc_mapper().calculator().quasi_clone(
            make_simple_structure(*hint_ptr, _hint_dofs),
            make_point_group(
                hint_ptr->point_group(),
                hint_ptr->supercell().sym_info().supercell_lattice()),
            SimpleStructure::SpeciesMode::ATOM),
        struc_mapper().lattice_weight(), 0., struc_mapper().robust(),
        struc_mapper().soft_va_limit(), struc_mapper().cost_tol());

    auto config_maps = tmapper.map_deformed_struc_impose_lattice_node(
        child_struc,
        xtal::LatticeNode(hint_ptr->ideal_lattice(), hint_ptr->ideal_lattice(),
                          Lattice(child_struc.lat_column_mat),
                          Lattice(child_struc.lat_column_mat),
                          child_struc.atom_info.size()),
        k_best);

    if (!config_maps.empty()) {
      hint_cost = best_cost = config_maps.rbegin()->cost;
    }
  }

  // --- Mapped child to prim ---
  std::set<xtal::MappingNode> struc_maps;

  // Method 1: hint && 'ideal' option (exactly known lattice mapping)
  if (hint_ptr && settings().fix_lattice_mapping) {
    xtal::LatticeNode lattice_node(
        hint_ptr->prim().lattice(), hint_ptr->ideal_lattice(),
        Lattice(child_struc.lat_column_mat),
        Lattice(child_struc.lat_column_mat), child_struc.atom_info.size());
    struc_maps = struc_mapper().map_deformed_struc_impose_lattice_node(
        child_struc, lattice_node, k_best);

    // Method 2: hint && 'fix_lattice' option (exactly known lattice)
  } else if (hint_ptr && settings().fix_lattice) {
    struc_maps = struc_mapper().map_deformed_struc_impose_lattice(
        child_struc, hint_ptr->ideal_lattice(), k_best,
        best_cost + struc_mapper().cost_tol());
    if (struc_maps.empty())
      result.fail_msg = "Unable to map structure using same lattice as " +
                        hint_ptr->name() +
                        ". Try setting \"fix_lattice\" : false.";

    // Method 3: hint && 'fix_volume' option (exactly known volume)
  } else if (hint_ptr && settings().fix_lattice_volume_range) {
    Index vol = hint_ptr->supercell().volume();
    struc_maps = struc_mapper().map_deformed_struc_impose_lattice_vols(
        child_struc, vol, vol, k_best, best_cost + struc_mapper().cost_tol());
    if (struc_maps.empty())
      result.fail_msg =
          "Unable to map structure assuming volume = " + std::to_string(vol) +
          ". Try setting \"fix_volume\" : false.";

    // Method 4: no hint && 'ideal' option (ideal integer supercell of prim)
  } else if (settings().fix_ideal) {
    struc_maps = struc_mapper().map_ideal_struc(child_struc, k_best);
    if (struc_maps.empty())
      result.fail_msg =
          "Imported structure has lattice vectors that are not a perfect "
          "supercell of PRIM. Try setting \"ideal\" : false.";

    // Method 5: most general mapping
  } else {
    struc_maps = struc_mapper().map_deformed_struc(
        child_struc, k_best, best_cost + struc_mapper().cost_tol());
    if (struc_maps.empty())
      result.fail_msg =
          "Unable to map structure to prim. May be incompatible structure, or "
          "provided settings may be too restrictive.";
  }

  // For each mapping, construct mapped & final configurations & properties
  for (xtal::MappingNode const &mapping_node : struc_maps) {
    xtal::SimpleStructure const &unmapped_child = child_struc;

    // 1) construct mapped_child, mapped_configuration, mapped_properties ---
    // xtal::SimpleStructure mapped_child =
    //     make_mapped_structure(mapping_node, unmapped_child, shared_prim());
    xtal::SimpleStructure mapped_child =
        struc_mapper().calculator().resolve_setting(mapping_node,
                                                    unmapped_child);
    // TODO: ^ use standalone make_mapped_structure

    mapped_child.properties["lattice_deformation_cost"] =
        _scalar_as_matrix(mapping_node.lattice_node.cost);
    mapped_child.properties["atomic_deformation_cost"] =
        _scalar_as_matrix(mapping_node.atomic_node.cost);
    mapped_child.properties["total_cost"] =
        _scalar_as_matrix(mapping_node.cost);
    // TODO: move this ^ to `make_mapped_structure`

    Configuration mapped_configuration =
        make_mapped_configuration(mapping_node, mapped_child, shared_prim());
    MappedProperties mapped_properties =
        make_mapped_properties(mapped_child, *shared_prim());

    // 2) construct final_configuration, final_properties ---

    // 2a: put mapped_configuration in the canonical supercell
    xtal::Lattice canonical_superlattice = xtal::canonical::equivalent(
        mapped_configuration.ideal_lattice(), shared_prim()->point_group(),
        shared_prim()->lattice().tol());
    std::shared_ptr<Supercell> shared_canonical_supercell =
        std::make_shared<Supercell>(shared_prim(), canonical_superlattice);

    FillSupercell f{shared_canonical_supercell};
    Configuration configuration_in_canon_scel = f(mapped_configuration);
    SymOp symop_to_canon_scel{f.symop()};
    Permutation perm_to_canon_scel{f.permutation()};
    MappedProperties properties_in_canon_scel{sym::copy_apply(
        symop_to_canon_scel, perm_to_canon_scel, mapped_properties)};
    xtal::SimpleStructure structure_in_canon_scel = make_simple_structure(
        mapped_configuration.supercell(), mapped_configuration.configdof(),
        properties_in_canon_scel);

    // L_canon_scel = (g * L_mapped) * T
    Eigen::Matrix3d g = symop_to_canon_scel.matrix();
    Eigen::Matrix3d L_canon_scel = canonical_superlattice.lat_column_mat();
    Eigen::Matrix3d L_mapped =
        mapped_configuration.ideal_lattice().lat_column_mat();
    Eigen::Matrix3l transformation_matrix_to_canon_scel =
        lround((g * L_mapped).inverse() * L_canon_scel);

    // 2b: find the symop to either (i) the 'strict' configuration most similar
    // to the unmapped child, or (ii) the canonical configuration
    SupercellSymInfo canonical_supercell_sym_info =
        shared_canonical_supercell->sym_info();
    PermuteIterator perm_it = canonical_supercell_sym_info.permute_begin();
    if (settings().finalize_strict) {
      // Strictness transformation reduces permutation swaps, translation
      // magnitude, and isometry character
      perm_it = Local::_strictest_equivalent(
          canonical_supercell_sym_info.permute_begin(),
          canonical_supercell_sym_info.permute_end(), mapping_node);
    } else {
      perm_it = configuration_in_canon_scel.to_canonical();
    }
    SymOp symop_to_final = perm_it.sym_op();
    Permutation perm_to_final = perm_it.combined_permute();

    // 2c: construct the `final_configuration`, `final_properties` and
    // `final_structure`
    Configuration final_configuration = configuration_in_canon_scel;
    final_configuration.apply_sym(perm_it);

    MappedProperties final_properties =
        sym::copy_apply(perm_it, properties_in_canon_scel);

    xtal::SimpleStructure final_structure = make_simple_structure(
        final_configuration.supercell(), final_configuration.configdof(),
        final_properties);

    // 2e: check relationship between final_configuration and hint
    ConfigComparison hint_status =
        _make_hint_status(hint_ptr, final_configuration);

    // 3) save solution (a ConfigurationMapping)
    result.maps.emplace(
        unmapped_child, mapping_node, mapped_child, mapped_configuration,
        mapped_properties, symop_to_canon_scel, perm_to_canon_scel,
        transformation_matrix_to_canon_scel, structure_in_canon_scel,
        configuration_in_canon_scel, properties_in_canon_scel, symop_to_final,
        perm_to_final, final_structure, final_configuration, final_properties,
        hint_status, hint_cost);
  }

  return result;
}

ConfigComparison ConfigMapper::_make_hint_status(
    Configuration const *hint_ptr,
    Configuration const &final_configuration) const {
  if (hint_ptr == nullptr) {
    return ConfigComparison::None;
  }
  if (final_configuration.supercell() != hint_ptr->supercell()) {
    return ConfigComparison::NewScel;
  }

  ConfigIsEquivalent all_equiv(*hint_ptr);
  if (all_equiv(final_configuration)) {
    return ConfigComparison::Identical;
  }

  PermuteIterator perm_begin =
      final_configuration.supercell().sym_info().permute_begin();
  PermuteIterator perm_end =
      final_configuration.supercell().sym_info().permute_end();
  for (PermuteIterator it = perm_begin; it != perm_end; ++it) {
    if (all_equiv(it, final_configuration)) {
      return ConfigComparison::Equivalent;
    }
  }

  ConfigIsEquivalent occ_equiv(*hint_ptr, {"occ"});
  for (PermuteIterator it = perm_begin; it != perm_end; ++it) {
    if (occ_equiv(it, final_configuration)) {
      return ConfigComparison::Derivative;
    }
  }
  return ConfigComparison::NewOcc;
}

PrimStrucMapCalculator::PrimStrucMapCalculator(
    BasicStructure const &_prim, std::vector<xtal::SymOp> const &_symgroup,
    SimpleStructure::SpeciesMode _species_mode /*=StrucMapping::ATOM*/)
    : SimpleStrucMapCalculator(
          make_simple_structure(_prim),
          _symgroup.empty() ? xtal::make_factor_group(_prim) : _symgroup,
          _species_mode, make_allowed_species(_prim)),
      m_prim(_prim) {}

// / 3) For each solution, the "mapped configuration", "mapped properties", and
// / "mapped structure" are transformed to be set in the canonical equivalent
// / supercell. They are then transformed to an equivalent "final
// / configuration", "final properties", and "final structure" by a change of
// / supercell and/or rotation and permutation within the supercell, as
// / specified by method settings.
// /
// /
// / Mapping to the canonical supercell
// / ----------------------------------
// /
// / The configuration in the canonical supercell,
// `configuration_in_canon_scel`, / is constructed by finding the canonical
// equivalent supercell and the first / prim factor group operation,
// `symop_to_canon_scel` and transformation / matrix,
// `transformation_matrix_to_canon_scel`, that transforms the /
// `mapped_configuration` supercell to the canonical supercell lattice:
// /
// / \f[
// /        S_1^{canon} = R_g^{canon} * (L_1 * T_1 * N) * T_1^{canon},
// / \f]
// / where:
// /
// / - \f$S_1^{canon}\f$: The lattive vector column matrix for the canonical
// / equivalent supercell of the mapped configuration
// / (`configuration_in_canon_scel.ideal_lattice().lat_column_mat()`).
// / - \f$L_1 * T_1 * N\f$: The prim superlattice the input structure is mapped
// / to, as defined in StrucMapper
// / (`mapping.lattice_node.parent.superlattice().lat_column_mat()`).
// / - \f$R_g^{canon}\f$ = A prim factor group operation
// / (`symop_to_canon_scel.matrix()`),
// / - \f$T_1^{canon}\f$ = An integer transformation matrix
// / (`transformation_matrix_to_canon_scel`)
// /
// / Global properties and DoF are transformed according to:
// /
// /     matrix = AnisoValTraits(dofname).symop_to_matrix(
// /          g_canon.matrix, g_canon.tau, g_canon.time_reversal);
// /     v_global_in_canon_scel = matrix * v_global_mapped,
// /
// / By mapping site coordinates, the permutation, `permutation_to_canon_scel`
// / which maps sites to the canonical supercell can also be found. Then site
// / properties and DoF are transformed according to:
// /
// /     v_site_in_canon_scel[i] = matrix * v_site_unmapped[perm_canon[i]],
// /     matrix = AnisoValTraits(dofname).symop_to_matrix(g_canon.matrix,
// /                  g_canon.tau, g_canon.time_reversal),
// /     perm_canon = permutation_to_canon_scel, The permutation that maps sites
// /                  in `mapped_configuration` to sites in
// /                  `configuration_in_canon_scel`.
// /
// /
// / Selection of an equivalent configuration
// / ----------------------------------------
// /
// / As a post-processing step, the `mapped_configuration` is transformed to the
// / `final_configuration` which is either 1) (default method) the canonical
// / equivalent configuration in the canonical supercell, or 2) (the 'strict'
// / method) the configuration which preserves the orientation of the original
// / unmapped structure as much as possible.
// /
// / The relationship between the `configuration_in_canon_scel` to the
// / `final_configuration` is stored in the mapping results by `symop_to_final`
// / and `permutation_to_final` which transform properties and DoF as in the
// / previous step.
// /
// / Global properties and DoF are transformed according to:
// /     v_global_final = matrix * v_global_in_canon_scel,
// /     matrix = AnisoValTraits(dofname).symop_to_matrix(g_final.matrix,
// /                  g_final.tau, g_final.time_reversal)
// /     g_final = symop_to_final, a prim factor group operation
// /
// / By mapping site coordinates, the permutation, `permutation_to_canon_scel`
// / which maps sites to the canonical supercell can also be found. Then site
// / properties and DoF are transformed according to:
// /     v_site_final[i] = matrix * v_site_in_canon_scel[perm_final[i]],
// /     matrix = AnisoValTraits(dofname).symop_to_matrix(g_final.matrix,
// /                  g_final.tau, g_final.time_reversal),
// /     perm_final = permutation_to_final, The permutation that maps sites
// /                  in `configuration_in_canon_scel` to sites in
// /                  `final_configuration`.
// /

}  // namespace CASM
