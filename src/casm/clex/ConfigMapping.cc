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

/// \class ConfigMapper
/// \brief Implements a method for mapping structures to configurations
/// \ingroup Configuration
///
/// The ConfigMapper class implements a method mapping an arbitrary crystal
/// structure to a configuration of a crystal template as described by a
/// xtal::BasicStructure.
///
///
/// Configuration mapping method summary
/// ------------------------------------
///
/// Reminders:
/// - parent: the reference structure
/// - child: the structure being mapped to parent
///
/// Steps:
/// 1) Generate and score potential mappings of the "child structure" to the
/// "parent structure" (reference). Lattice mapping is done first, then atomic
/// assignment for each potential lattice mapping solution. A variety of
/// settings control which lattice mappings are considered, how many solutions
/// are stored, and how they are scored. See `StrucMapper` for details.
///
/// 2) The best scoring mappings are used to create "mapped structures",
/// which must be de-rotated according to the mapping, expressed in terms of
/// molecules, ordered according to the atomic assignment and to match the
/// ordering of sites in the appropriate supercell, and have displacements,
/// strains, and other DoF and properties appropriately de-rotated, permuted,
/// scaled (if extensive), and expressed in terms of the reference structure,
/// so that equivalent Configurations ("mapped configurations") may be
/// constructed.
///
/// 3) The "mapped configuration" may be transformed to an equivalent "final
/// configuration" by a change of supercell and/or rotation and permutation
/// within the supercell, as specified by method settings.
///
/// --- Making the mapped child structure: ---
///
/// After a total mapping solution is obtained, the `mapped_child` is
/// constructed. The `mapped_child` is equivalent to the original unmapped child
/// structure, with the following transformations:
/// - Atoms and inferred vacancies are resolved to molecules
/// - Molecules are re-ordered to match the order of sites in the Supercell with
/// lattice vectors equal to those of the mapped structure (TODO: check this)
/// - Molecule displacements are set to the average of the coordinates of atoms
/// in the molecule, and then molecules are translated such that average
/// displacement of molecules relative to the ideal sites is 0
/// - Isometry and permutation are applied to properties of the unmapped child
/// structure
///
/// For global properties this means:
///     v_global_mapped = matrix * v_global_unmapped,
///     matrix = AnisoValTraits(property_type).symop_to_matrix(
///                  Q^{N}, (V^{N}).inverse() * trans, time_reversal)
///
/// For site properties this means:
///     v_site_mapped[i] = matrix * v_global_unmapped[perm[i]]
///     matrix = AnisoValTraits(property_type).symop_to_matrix(
///                  Q^{N}, (V^{N}).inverse() * trans, time_reversal)
///
///
/// --- Constructing the mapped configuration and mapped properties: ---
///
/// Given the mapped, molecule-ized, child structure, `mapped_child`, the
/// `mapped_configuration` are `mapped_properties` are constructed.
///
///
/// --- Selection of an equivalent configuration: ---
///
/// As a post-processing step, the `mapped_configuration` is transformed to the
/// `final_configuration` which is either 1) (default method) the canonical
/// equivalent configuration in the canonical supercell, or 2) (the 'strict'
/// method) the configuration which preserves the orientation of the original
/// unmapped structure as much as possible.
///
/// The relationship between the `mapped_configuration` and the
/// `final_configuration` is stored in the mapping results by `symop_to_final`
/// and `transformation_matrix_to_final` via:
/// - final_configuration.ideal_lattice() = (matrix_to_final *
/// mapped_configuration.ideal_lattice()) * transformation_matrix_to_final
/// - final_configuration.global_dofs(dofname) = matrix *
/// mapped_configuration.global_dofs(dofname)
/// - final_configuration.local_dofs(dofname)[i] = matrix *
/// mapped_configuration.local_dofs(dofname)[permutation_to_final[i]],
///
/// where matrix =
/// AnisoValTraits(dofname).symop_to_matrix(symop_to_final.matrix,
/// symop_to_final.tau, symop_to_final.time_reversal)
///
/// A similar transformation is used to construct `final_properties` from
/// `mapped_properties`. The relationship between mapped and final properties is
/// the same as the relationship between mapped and final DoF.
///

namespace CASM {

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
    std::cout << "A: \n" << A << std::endl;
    std::cout << "B: \n" << B << std::endl << std::endl;
    throw std::runtime_error(message);
  }
}
void check_equal(Eigen::MatrixXd const &A, Eigen::MatrixXd const &B,
                 std::string message) {
  if (!almost_equal(A, B)) {
    std::cout << "A: \n" << A << std::endl;
    std::cout << "B: \n" << B << std::endl << std::endl;
    throw std::runtime_error(message);
  }
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
/// \param mapped_child Child structure, already transformed so that it is
///     mapped to the parent superstructure. This means that
///     `mapped_child.mol_info.names[i]` specifies the molecule at site index i
///     in the resulting configuration, and mapped_child.properties and
///     mapped_child.mol_info.properties are DoF values in the standard basis
///     which can be directly written to global and local Configuration DoF
///     without any further rotation, etc.
/// \param shared_prim Shared prim structure that the child was mapped to.
///
Configuration make_mapped_configuration(
    SimpleStructure const &mapped_child,
    std::shared_ptr<Structure const> const &shared_prim) {
  auto const &prim = *shared_prim;
  Lattice mapped_lattice{mapped_child.lat_column_mat, prim.lattice().tol()};
  std::shared_ptr<Supercell const> mapped_shared_supercell =
      std::make_shared<Supercell const>(shared_prim, mapped_lattice);
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
                           ConfigMapping::Settings const &_settings,
                           double _tol)
    : m_shared_prim(_shared_prim),
      m_struc_mapper(
          PrimStrucMapCalculator(
              *shared_prim(), adapter::Adapter<xtal::SymOpVector, SymGroup>()(
                                  shared_prim()->factor_group())),
          _settings.lattice_weight, _settings.max_vol_change,
          _settings.options(),
          _tol > 0. ? _tol : shared_prim()->lattice().tol(),
          _settings.min_va_frac, _settings.max_va_frac),
      m_settings(_settings) {
  if (settings().filter) {
    m_struc_mapper.set_filter(settings().filter);
  }

  for (Lattice const &lattice : settings().forced_lattices) {
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
        struc_mapper().lattice_weight(), 0., struc_mapper().options(),
        struc_mapper().cost_tol());

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
  if (hint_ptr && settings().ideal) {
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
  } else if (hint_ptr && settings().fix_volume) {
    Index vol = hint_ptr->supercell().volume();
    struc_maps = struc_mapper().map_deformed_struc_impose_lattice_vols(
        child_struc, vol, vol, k_best, best_cost + struc_mapper().cost_tol());
    if (struc_maps.empty())
      result.fail_msg =
          "Unable to map structure assuming volume = " + std::to_string(vol) +
          ". Try setting \"fix_volume\" : false.";

    // Method 4: no hint && 'ideal' option (ideal integer supercell of prim)
  } else if (settings().ideal) {
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

    // TODO: move this to `make_mapped_structure`
    // mapped_structure.properties["lattice_deformation_cost"] =
    //     mapping_node.lattice_node.cost;
    // mapped_structure.properties["atomic_deformation_cost"] =
    //     mapping_node.atomic_node.cost;
    // mapped_structure.properties["total_cost"] = mapping_node.cost;

    xtal::SimpleStructure mapped_child =
        make_mapped_structure(mapping_node, unmapped_child, shared_prim());
    Configuration mapped_configuration =
        make_mapped_configuration(mapped_child, shared_prim());
    MappedProperties mapped_properties =
        make_mapped_properties(mapped_child, *shared_prim());

    // 2) construct final_configuration, final_properties ---

    // 2a: put mapped_configuration in the canonical supercell
    xtal::Lattice canonical_superlattice = xtal::canonical::equivalent(
        mapped_configuration.ideal_lattice(), shared_prim()->point_group(),
        shared_prim()->lattice().tol());
    std::shared_ptr<Supercell> shared_canonical_supercell =
        std::make_shared<Supercell>(shared_prim(), canonical_superlattice);
    Configuration mapped_configuration_in_canonical_supercell =
        fill_supercell(mapped_configuration, shared_canonical_supercell);
    MappedProperties mapped_properties_in_canonical_supercell =
        mapped_properties;
    // TODO: ^ transform properties also, and record symop to canonical

    // 2b: find the symop to either (i) the 'strict' configuration most similar
    // to the unmapped child, or (ii) the canonical configuration
    SupercellSymInfo canonical_supercell_sym_info =
        shared_canonical_supercell->sym_info();
    PermuteIterator perm_it = canonical_supercell_sym_info.permute_begin();
    if (settings().strict) {
      // Strictness transformation reduces permutation swaps, translation
      // magnitude, and isometry character
      perm_it = Local::_strictest_equivalent(
          canonical_supercell_sym_info.permute_begin(),
          canonical_supercell_sym_info.permute_end(), mapping_node);
    } else {
      perm_it = mapped_configuration_in_canonical_supercell.to_canonical();
    }

    // 2c: record the symop_to_final and transformation_matrix_to_final
    xtal::SymOp symop_to_final =
        adapter::Adapter<xtal::SymOp, SymOp>()(perm_it.sym_op());

    // L_final = matrix * L_mapped * T
    Eigen::Matrix3d M = symop_to_final.matrix;
    Eigen::Matrix3d L_final = canonical_superlattice.lat_column_mat();
    Eigen::Matrix3d L_mapped =
        mapped_configuration.ideal_lattice().lat_column_mat();
    Eigen::Matrix3l transformation_matrix_to_final =
        lround((M * L_mapped).inverse() * L_final);

    // 2d: construct the `final_configuration`
    Configuration final_configuration =
        mapped_configuration_in_canonical_supercell;
    final_configuration.apply_sym(perm_it);

    // 2e: construct the `final_properties`
    MappedProperties final_properties =
        copy_apply(perm_it, mapped_properties_in_canonical_supercell);

    // 2e: check relationship between final_configuration and hint
    HintStatus hint_status = _make_hint_status(hint_ptr, final_configuration);

    // 3) save solution (a ConfigMapperResult::ConfigurationMapping)
    result.maps.emplace(unmapped_child, mapping_node, mapped_child,
                        mapped_configuration, mapped_properties, symop_to_final,
                        transformation_matrix_to_final, final_configuration,
                        final_properties, hint_status, hint_cost);
  }

  return result;
}

ConfigMapperResult::HintStatus ConfigMapper::_make_hint_status(
    Configuration const *hint_ptr,
    Configuration const &final_configuration) const {
  if (hint_ptr == nullptr) {
    return HintStatus::None;
  }
  if (final_configuration.supercell() != hint_ptr->supercell()) {
    return HintStatus::NewScel;
  }

  ConfigIsEquivalent all_equiv(*hint_ptr);
  if (all_equiv(final_configuration)) {
    return HintStatus::Identical;
  }

  PermuteIterator perm_begin =
      final_configuration.supercell().sym_info().permute_begin();
  PermuteIterator perm_end =
      final_configuration.supercell().sym_info().permute_end();
  for (PermuteIterator it = perm_begin; it != perm_end; ++it) {
    if (all_equiv(it, final_configuration)) {
      return HintStatus::Equivalent;
    }
  }

  ConfigIsEquivalent occ_equiv(*hint_ptr, {"occ"});
  for (PermuteIterator it = perm_begin; it != perm_end; ++it) {
    if (occ_equiv(it, final_configuration)) {
      return HintStatus::Derivative;
    }
  }
  return HintStatus::NewOcc;
}

PrimStrucMapCalculator::PrimStrucMapCalculator(
    BasicStructure const &_prim, std::vector<xtal::SymOp> const &_symgroup,
    SimpleStructure::SpeciesMode _species_mode /*=StrucMapping::ATOM*/)
    : SimpleStrucMapCalculator(
          make_simple_structure(_prim),
          _symgroup.empty() ? xtal::make_factor_group(_prim) : _symgroup,
          _species_mode, make_allowed_species(_prim)),
      m_prim(_prim) {}

}  // namespace CASM
