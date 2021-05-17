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
#include "casm/clex/ParamComposition.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/SimpleStructureTools.hh"
#include "casm/clex/Supercell.hh"
#include "casm/completer/Handlers.hh"
#include "casm/crystallography/Adapter.hh"
#include "casm/crystallography/BasicStructureTools.hh"
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

namespace Local {
static int _permute_dist(MappingNode::MoleculeMap const &_perm,
                         MappingNode const &_node) {
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
                                      MappingNode const &_node) {
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

//*******************************************************************************************
Index ConfigMapperResult::n_optimal(double tol /*=TOL*/) const {
  Index result = 0;
  auto it = maps.begin();
  while (it != maps.end() &&
         almost_equal((it->first).cost, (maps.begin()->first).cost, tol)) {
    ++result;
    ++it;
  }
  return result;
}

//*******************************************************************************************
namespace ConfigMapping {

xtal::StrucMapping::AllowedSpecies _allowed_species(
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
              "ConfigMapping::_allowed_species may only be called on "
              "structures with single-atom species.");
        }
        result[i].push_back(mol.atom(0).name());
      }
    }
    ++i;
  }

  return result;
}
}  // namespace ConfigMapping

//*******************************************************************************************

/// \brief Reorders the permutation and compounds the spatial isometry (rotation
/// + translation) of _node with that of _it
MappingNode copy_apply(PermuteIterator const &_it, MappingNode const &_node,
                       bool transform_cost_mat) {
  MappingNode result(_node);
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

//*******************************************************************************************

/// \brief Initializes configdof of Supercell '_scel' corresponding to an
/// idealized child structure (encoded by _child_struc) _child_struc is assumed
/// to have been idealized via structure-mapping or to be the result of
/// converting a configuration to a SimpleStructure. result.second gives list of
/// properties that were utilized in the course of building the configdof
std::pair<ConfigDoF, std::set<std::string> > to_configdof(
    SimpleStructure const &_child_struc, Supercell const &_scel) {
  SimpleStructure::Info const &c_info(_child_struc.mol_info);
  std::pair<ConfigDoF, std::set<std::string> > result(
      make_configdof(_scel, TOL), {});
  PrimClex::PrimType const &prim(_scel.prim());
  Index i = 0;
  for (Index b = 0; b < prim.basis().size(); ++b) {
    for (Index l = 0; l < _scel.volume(); ++l, ++i) {
      Index j = 0;
      for (; j < prim.basis()[b].occupant_dof().size(); ++j) {
        if (c_info.names[i] == prim.basis()[b].occupant_dof()[j]) {
          result.first.occ(i) = j;
          break;
        }
      }
      if (j == prim.structure().basis()[b].occupant_dof().size())
        throw std::runtime_error(
            "Attempting to initialize ConfigDoF from SimpleStructure. Species "
            "'" +
            c_info.names[i] + "' is not allowed on sublattice " +
            std::to_string(b));
    }
  }

  for (auto const &dof : result.first.global_dofs()) {
    auto val = DoFType::traits(dof.first).find_values(_child_struc.properties);
    result.first.global_dof(dof.first).from_standard_values(val.first);
    result.second.insert(val.second.begin(), val.second.end());
  }

  for (auto const &dof : result.first.local_dofs()) {
    auto val = DoFType::traits(dof.first).find_values(c_info.properties);
    result.first.local_dof(dof.first).from_standard_values(val.first);
    result.second.insert(val.second.begin(), val.second.end());
  }

  return result;
}
//*******************************************************************************************

PrimStrucMapCalculator::PrimStrucMapCalculator(
    BasicStructure const &_prim, std::vector<xtal::SymOp> const &_symgroup,
    SimpleStructure::SpeciesMode _species_mode /*=StrucMapping::ATOM*/)
    : SimpleStrucMapCalculator(
          make_simple_structure(_prim),
          _symgroup.empty() ? xtal::make_factor_group(_prim) : _symgroup,
          _species_mode, ConfigMapping::_allowed_species(_prim)),
      m_prim(_prim) {}

//*******************************************************************************************

// void ConfigMapper::add_allowed_lattices(
//     std::vector<std::string> const &_lattice_names) {
//   for (std::string const &_name : _lattice_names) {
//     auto it = primclex().template db<Supercell>().find(_name);
//     if (it == primclex().template db<Supercell>().end())
//       throw std::runtime_error(
//           "Could not add mapping lattice constraint " + _name +
//           " because no supercell having that name exists in database.\n");
//     m_struc_mapper.add_allowed_lattice(it->lattice());
//   }
// }

//*******************************************************************************************

void ConfigMapper::clear_allowed_lattices() {
  m_struc_mapper.clear_allowed_lattices();
}

//*******************************************************************************************

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
  // if (!settings().filter.empty()) {
  //   /// If a filter string is specified, construct the Supercell query filter
  //   /// for restricting potential supercells for mapping
  //   DataFormatter<Supercell> formatter =
  //       _pclex.settings().query_handler<Supercell>().dict().parse(
  //           settings().filter);
  //   auto filter = [formatter, &_pclex](Lattice const &parent,
  //                                      Lattice const &child) -> bool {
  //     ValueDataStream<bool> check_stream;
  //     check_stream << formatter(Supercell(&_pclex, parent));
  //     return check_stream.value();
  //   };
  //
  //   m_struc_mapper.set_filter(filter);
  // }
  if (settings().filter) {
    m_struc_mapper.set_filter(settings().filter);
  }

  for (Lattice const &lattice : settings().forced_lattices) {
    // auto it = _pclex.db<Supercell>().find(scel);
    // if (it == _pclex.db<Supercell>().end())
    //   throw std::runtime_error("Cannot restrict mapping to lattice " + scel +
    //                            ". Superlattice does not exist in project.");
    m_struc_mapper.add_allowed_lattice(lattice);
  }
}

//*******************************************************************************************
ConfigMapperResult ConfigMapper::import_structure(
    SimpleStructure const &_child_struc, Configuration const *hint_ptr,
    std::vector<DoFKey> const &_hint_dofs) const {
  return import_structure(_child_struc, settings().k_best, hint_ptr,
                          _hint_dofs);
}

//*******************************************************************************************

ConfigMapperResult ConfigMapper::import_structure(
    SimpleStructure const &child_struc, Index k, Configuration const *hint_ptr,
    std::vector<DoFKey> const &_hint_dofs) const {
  ConfigMapperResult result;
  double best_cost = xtal::StrucMapping::big_inf();

  // bool is_new_config(true);
  double hint_cost;
  if (hint_ptr != nullptr) {
    StrucMapper tmapper(
        *struc_mapper().calculator().quasi_clone(
            make_simple_structure(*hint_ptr, _hint_dofs),
            make_point_group(
                hint_ptr->point_group(),
                hint_ptr->supercell().sym_info().supercell_lattice()),
            SimpleStructure::SpeciesMode::ATOM),
        struc_mapper().lattice_weight(), 0., struc_mapper().options(),
        struc_mapper().cost_tol());

    /*
    auto config_maps = tmapper.map_deformed_struc_impose_lattice(child_struc,
                                                                 hint_ptr->ideal_lattice(),
                                                                 1);
    */
    auto config_maps = tmapper.map_deformed_struc_impose_lattice_node(
        child_struc,
        xtal::LatticeNode(hint_ptr->ideal_lattice(), hint_ptr->ideal_lattice(),
                          Lattice(child_struc.lat_column_mat),
                          Lattice(child_struc.lat_column_mat),
                          child_struc.atom_info.size()),
        k);

    // Refactor into external routine A. This is too annoying with the current
    // way that supercells are managed
    if (!config_maps.empty()) {
      hint_cost = best_cost = config_maps.rbegin()->cost;
    }
    //\End routine A
    // TODO: Check equivalence with hint_ptr
  }

  std::set<MappingNode> struc_maps;

  if (hint_ptr && settings().ideal) {
    xtal::LatticeNode lattice_node(
        hint_ptr->prim().lattice(), hint_ptr->ideal_lattice(),
        Lattice(child_struc.lat_column_mat),
        Lattice(child_struc.lat_column_mat), child_struc.atom_info.size());
    struc_maps = struc_mapper().map_deformed_struc_impose_lattice_node(
        child_struc, lattice_node, k);
  } else if (hint_ptr && settings().fix_lattice) {
    struc_maps = struc_mapper().map_deformed_struc_impose_lattice(
        child_struc, hint_ptr->ideal_lattice(), k,
        best_cost + struc_mapper().cost_tol());
    if (struc_maps.empty())
      result.fail_msg = "Unable to map structure using same lattice as " +
                        hint_ptr->name() +
                        ". Try setting \"fix_lattice\" : false.";

  } else if (hint_ptr && settings().fix_volume) {
    Index vol = hint_ptr->supercell().volume();
    struc_maps = struc_mapper().map_deformed_struc_impose_lattice_vols(
        child_struc, vol, vol, k, best_cost + struc_mapper().cost_tol());
    if (struc_maps.empty())
      result.fail_msg =
          "Unable to map structure assuming volume = " + std::to_string(vol) +
          ". Try setting \"fix_volume\" : false.";

  } else if (settings().ideal) {
    struc_maps = struc_mapper().map_ideal_struc(child_struc, k);
    if (struc_maps.empty())
      result.fail_msg =
          "Imported structure has lattice vectors that are not a perfect "
          "supercell of PRIM. Try setting \"ideal\" : false.";
  } else {
    struc_maps = struc_mapper().map_deformed_struc(
        child_struc, k, best_cost + struc_mapper().cost_tol());
    if (struc_maps.empty())
      result.fail_msg =
          "Unable to map structure to prim. May be incompatible structure, or "
          "provided settings may be too restrictive.";
  }

  // Refactor into external routine A. This is too annoying with the current way
  // that supercells are managed
  for (auto const &map : struc_maps) {
    std::shared_ptr<Supercell> shared_scel = std::make_shared<Supercell>(
        shared_prim(), map.lattice_node.parent.superlattice());
    SimpleStructure resolved_struc =
        struc_mapper().calculator().resolve_setting(map, child_struc);
    std::pair<ConfigDoF, std::set<std::string> > tdof =
        to_configdof(resolved_struc, *shared_scel);
    Configuration tconfig(shared_scel, tdof.first);
    PermuteIterator perm_it = shared_scel->sym_info().permute_begin();
    if (settings().strict) {
      // Strictness transformation reduces permutation swaps, translation
      // magnitude, and isometry character
      perm_it = Local::_strictest_equivalent(
          shared_scel->sym_info().permute_begin(),
          shared_scel->sym_info().permute_end(), map);
    } else {
      perm_it = tconfig.to_canonical();
    }
    tconfig.apply_sym(perm_it);
    MappingNode resolved_node = copy_apply(perm_it, map);
    resolved_struc =
        struc_mapper().calculator().resolve_setting(resolved_node, child_struc);
    result.maps.emplace(
        resolved_node, ConfigMapperResult::Individual(std::move(tconfig),
                                                      std::move(resolved_struc),
                                                      std::move(tdof.second)));
  }
  //\End routine A

  if (hint_ptr != nullptr) {
    ConfigIsEquivalent all_equiv(*hint_ptr);

    ConfigIsEquivalent occ_equiv(*hint_ptr, {"occ"});

    for (auto &map : result.maps) {
      map.second.hint_cost = hint_cost;

      if (map.second.config.supercell() != hint_ptr->supercell()) {
        map.second.hint_status = HintStatus::NewScel;
        continue;
      }
      if (all_equiv(map.second.config)) {
        map.second.hint_status = HintStatus::Identical;
        continue;
      }
      PermuteIterator perm_begin =
          map.second.config.supercell().sym_info().permute_begin();
      PermuteIterator perm_end =
          map.second.config.supercell().sym_info().permute_end();
      for (PermuteIterator it = perm_begin; it != perm_end; ++it) {
        if (all_equiv(it, map.second.config)) {
          map.second.hint_status = HintStatus::Equivalent;
          break;
        }
      }
      if (map.second.hint_status == HintStatus::Equivalent) continue;
      map.second.hint_status = HintStatus::NewOcc;
      for (PermuteIterator it = perm_begin; it != perm_end; ++it) {
        if (occ_equiv(it, map.second.config)) {
          map.second.hint_status = HintStatus::Derivative;
          break;
        }
      }
    }
  }

  return result;
}
}  // namespace CASM
