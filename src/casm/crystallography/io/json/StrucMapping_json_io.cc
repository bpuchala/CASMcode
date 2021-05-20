#include "casm/crystallography/io/json/StrucMapping_json_io.hh"

#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/crystallography/StrucMapping.hh"

namespace CASM {

jsonParser &to_json(xtal::MappingNode const &mapping_node, jsonParser &json) {
  // json["descriptions"]["lattice_mapping"] = "Defined such that
  // `parent_supercell_lattice_column_matrix == stretch * isometry *
  // unmapped_child_supercell_lattice_column_matrix`, where
  // `parent_supercell_lattice_column_matrix` is an ideal supercell of the
  // parent structure being mapped to (usually a prim), and
  // `unmapped_child_supercell_lattice_column_matrix` is a supercell of the
  // structure being mapped.";
  //
  // json["descriptions"]["forced_assignments"] = "Listed as an array of pairs
  // [parent_site_index, child_site_index], which were enforced such that the
  // assignment maps child_site_index onto parent_site_index.";
  //
  // json["description"]["translation"] =
  // "translation=parent_coord.col(i)-mapped_child_coord.col(permutation[i])+displacement.col(i),
  // for each site index, i. `mapped_child_coord` are the coordinates in the
  // unrotated, undeformed (i.e. mapped) child structure."
  //
  // json["description"]["atom_displacement"] = "Nx3 matrix of displacements for
  // all sites in parent supercell. Sites with vacancies are included, but set
  // to zero displacement."
  //
  // json["description"]["atom_permutation"] = "The value, j =
  // atom_permutation[i], indicates that j-th atom of the child superstructure
  // maps onto the i-th site of parent superstructure. If the parent has N sites
  // and the child has M<N atoms, vacancies are designated by values j>=M."
  //
  // json["description"]["mol_map"] = "mol_map[j] is a set of indices of atoms
  // in the parent structure that comprise the molecule at j-th molecular site
  // in the parent supercell."
  //
  // json["description"]["mol_labels"] = "Labels assigned to molecules on each
  // site, as [name, occupant_index]."

  // lattice mapping data:
  xtal::LatticeNode const &lattice_node = mapping_node.lattice_node;
  json["stretch"] = lattice_node.stretch;
  json["isometry"] = lattice_node.isometry;
  json["parent_supercell_lattice_column_matrix"] =
      lattice_node.parent.superlattice().lat_column_mat();
  json["parent_transformation_matrix_to_super"] =
      lattice_node.parent.transformation_matrix_to_super();
  json["lattice_deformation_cost"] = lattice_node.cost;
  json["lattice_deformation_cost_method"] = lattice_node.cost_method;

  // occupation mapping (assignment) data:
  xtal::AssignmentNode const &assignment_node = mapping_node.atomic_node;
  json["translation"] = assignment_node.translation;
  json["time_reversal"] = assignment_node.time_reversal;
  json["forced_assignments"] = assignment_node.forced_on;
  json["atomic_deformation_cost"] = assignment_node.cost;
  json["atomic_deformation_cost_method"] = assignment_node.cost_method;

  json["atom_displacement"] = mapping_node.atom_displacement.transpose();
  json["atom_permutation"] = mapping_node.atom_permutation;
  json["mol_map"] = mapping_node.mol_map;
  json["mol_labels"] = mapping_node.mol_labels;

  // total mapping data:
  json["lattice_weight"] = mapping_node.lattice_weight;
  json["atomic_weight"] = mapping_node.atomic_weight;
  json["total_cost"] = mapping_node.cost;
  json["total_cost_method"] = mapping_node.cost_method;

  return json;
}

}  // namespace CASM
