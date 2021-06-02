#include "autotools.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clex/SimpleStructureTools.hh"
#include "casm/crystallography/Adapter.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/SimpleStrucMapCalculator.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/StrucMapping.hh"
#include "casm/crystallography/io/json/StrucMapping_json_io.hh"
#include "gtest/gtest.h"

// BCC mapping tests

using namespace CASM;

// \brief Confirm that the MappingNode maps the unmapped child to the mapped
// child constructed by `mapper.resolve_setting` as expected
void assert_mapping_relations(std::set<xtal::MappingNode> const &mappings,
                              xtal::StrucMapper const &mapper,
                              xtal::SimpleStructure const &unmapped_child);

xtal::BasicStructure make_bcc_basicstructure() {
  // binary {A, B} conventional BCC base structure

  using namespace xtal;

  Molecule A = Molecule::make_atom("A");
  Molecule B = Molecule::make_atom("B");

  Lattice lat{Eigen::Vector3d{4.0, 0.0, 0.0}, Eigen::Vector3d{0.0, 4.0, 0.0},
              Eigen::Vector3d{0.0, 0.0, 4.0}};

  BasicStructure struc{lat};
  struc.set_basis({Site{Coordinate{0.0, 0.0, 0.0, lat, CART}, {A, B}},
                   Site{Coordinate{2.0, 2.0, 2.0, lat, CART}, {A, B}}});
  return struc;
}

xtal::SimpleStructure make_bcc_simplestructure_1x1x1() {
  xtal::SimpleStructure simple;

  // clang-format off
  simple.lat_column_mat <<
    4., 0., 0.,
    0., 4., 0.,
    0., 0., 4.1;

  simple.atom_info.coords = Eigen::MatrixXd(3, 2);
  simple.atom_info.coords <<
    0., 2.1,
    0., 2.,
    0., 2.;

  simple.atom_info.names = {"A", "A"};
  // clang-format on

  return simple;
}

xtal::SimpleStructure make_bcc_simplestructure_2x1x1() {
  xtal::SimpleStructure simple;

  // clang-format off
  simple.lat_column_mat <<
    8., 0., 0.,
    0., 4., 0.,
    0., 0., 4.;

  simple.atom_info.coords = Eigen::MatrixXd(3, 4);
  simple.atom_info.coords <<
    0., 2., 4., 6.,
    0., 2., 0., 2.,
    0., 2., 0., 2.;

  simple.atom_info.names = {"A", "A", "A", "A"};
  // clang-format on

  return simple;
}

xtal::SymOpVector identity_group() { return {CASM::xtal::SymOp::identity()}; }

// This fixture lays out all the StrucMapper parameters

class StrucMapperTest : public testing::Test {
 protected:
  // parent: (reference structure)
  xtal::BasicStructure basicstructure;
  xtal::SimpleStructure parent;
  xtal::SymOpVector parent_factor_group;

  // child: (what is being mapped to parent)
  xtal::SimpleStructure child;

  // --- StrucMapper constructor parameters: ---

  // Specifies parent structure, factor group, allowed species on each
  // structure site (default: SimpleStrucMapCalculator(parent));
  xtal::SimpleStrucMapCalculator calculator;

  // Total cost weighting, lattice deformation component (default: 0.5)
  // total_cost = lattice_weight * lattice_deformation_cost +
  //              (1-lattice_weight) * atomic_deformation_cost
  double lattice_weight;

  // Potential parent superlattice volume constraint (default: 0.5)
  double max_volume_change;

  // Robust, soft_va_limit options (default: none)
  //   (bitwise-OR of StrucMapper::Options)
  xtal::StrucMapper::Options struc_mapper_options;

  // Mapping cost comparison tolerance (default: TOL)
  double cost_tol;

  // Potential parent superlattice volume constraint (default: 0.)
  double min_va_frac;

  // Potential parent superlattice volume constraint (default: 1.)
  double max_va_frac;

  // --- Standard mapping method parameters: ---

  // default: 1
  Index k_best = 1;

  // default: xtal::StrucMapping::big_inf()
  double max_cost = xtal::StrucMapping::big_inf();

  // default: -TOL
  double min_cost = -TOL;

  // default: false
  bool keep_invalid = false;

  // default: {xtal::SymOp::identity()}
  xtal::SymOpVector child_factor_group = identity_group();

  StrucMapperTest()
      :

        // parent
        basicstructure(make_bcc_basicstructure()),
        parent(make_simple_structure(basicstructure)),
        parent_factor_group(xtal::make_factor_group(basicstructure)),

        // child
        child(make_bcc_simplestructure_2x1x1()),

        // StrucMapper constructor args:
        calculator(parent, parent_factor_group,
                   CASM::xtal::SimpleStructure::SpeciesMode::ATOM,
                   allowed_molecule_names(basicstructure)),
        lattice_weight(0.5),
        max_volume_change(0.5),
        struc_mapper_options(xtal::StrucMapper::Options::none),
        cost_tol(TOL),
        min_va_frac(0.),
        max_va_frac(1.),

        // Standard mapping method parameters:
        k_best(1),
        max_cost(xtal::StrucMapping::big_inf()),
        min_cost(-TOL),
        keep_invalid(false),
        child_factor_group(identity_group()) {
    // Additional mapping options: (for reference or testing)

    // LatticeMap unimodular matrix N element range (default: 1)
    // mapper.set_lattice_transformation_range(1);

    // Set whether to use symmetry-breaking strain cost (default: false)
    // mapper.set_symmetrize_lattice_cost(false);

    // Set whether to use symmetry-breaking strain cost (default: false)
    // auto parent_permute_group =
    //     xtal::make_permutation_representation(basic_struc, parent_fg);
    // mapper.set_symmetrize_atomic_cost(true, parent_factor_group,
    // parent_permute_group);

    // Specify particular parent superlattices to consider (default none)
    // // clang-format off
    // Eigen::Matrix3d L;
    // L <<
    //   8., 0., 0.,
    //   0., 4., 0.,
    //   0., 0., 4.;
    // // clang-format on
    // mapper.add_allowed_lattice(Lattice(L));

    // Set a lattice filter function (default none)
    // auto filter_f = [](Lattice const &parent_prim_lattice,
    //                    Lattice const &proposed_parent_superlattice) {
    //   return volume(proposed_parent_superlattice)/
    //   volume(parent_prim_superlattice) < 2. + TOL);
    // }
    // mapper.set_filter(filter_f);
  }

  // void print_mappings_json(std::set<xtal::MappingNode> const &mappings) {
  //   jsonParser json;
  //   to_json(mappings, json["mappings"]);
  //   std::cout << json << std::endl;
  // }
};

/// The `StrucMapper::map_X_struc[_Y]` methods run the structure mapping
/// algorithm in general or under different simplifying assumptions:
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
/// The following examples all map a 4-atom bcc structure to a 2-atom
/// conventional bcc parent structure, resulting in 2 equivalent mappings each,
/// and then call `assert_mapping_relations` to check / demonstrate how
/// the "unmapped child", "mapped child", and "parent" structures are related
/// via the xtal::MappingNode mapping results.

TEST_F(StrucMapperTest, MapDeformedStruc0) {
  xtal::StrucMapper mapper(calculator, lattice_weight, max_volume_change,
                           struc_mapper_options, cost_tol, min_va_frac,
                           max_va_frac);

  std::set<xtal::MappingNode> mappings = mapper.map_deformed_struc(
      child, k_best, max_cost, min_cost, keep_invalid, child_factor_group);

  EXPECT_EQ(mappings.size(), 2);
  assert_mapping_relations(mappings, mapper, child);
}

TEST_F(StrucMapperTest, MapDeformedStrucImposeLatticeVols0) {
  xtal::StrucMapper mapper(calculator, lattice_weight, max_volume_change,
                           struc_mapper_options, cost_tol, min_va_frac,
                           max_va_frac);

  Index min_vol = 1;
  Index max_vol = 4;
  std::set<xtal::MappingNode> mappings =
      mapper.map_deformed_struc_impose_lattice_vols(
          child, min_vol, max_vol, k_best, max_cost, min_cost, keep_invalid,
          child_factor_group);

  EXPECT_EQ(mappings.size(), 2);
  assert_mapping_relations(mappings, mapper, child);
}

TEST_F(StrucMapperTest, MapDeformedStrucImposeLattice0) {
  xtal::StrucMapper mapper(calculator, lattice_weight, max_volume_change,
                           struc_mapper_options, cost_tol, min_va_frac,
                           max_va_frac);

  // clang-format off
  Eigen::Matrix3d L;
  L <<
    8., 0., 0.,
    0., 4., 0.,
    0., 0., 4.;
  // clang-format on
  std::set<xtal::MappingNode> mappings =
      mapper.map_deformed_struc_impose_lattice(child, xtal::Lattice(L), k_best,
                                               max_cost, min_cost, keep_invalid,
                                               child_factor_group);

  EXPECT_EQ(mappings.size(), 2);
  assert_mapping_relations(mappings, mapper, child);
}

TEST_F(StrucMapperTest, MapDeformedStrucImposeLatticeNode0) {
  xtal::StrucMapper mapper(calculator, lattice_weight, max_volume_change,
                           struc_mapper_options, cost_tol, min_va_frac,
                           max_va_frac);

  // clang-format off
  Eigen::Matrix3d L;
  L <<
    8., 0., 0.,
    0., 4., 0.,
    0., 0., 4.;
  // clang-format on
  xtal::Lattice parent_prim_lattice{parent.lat_column_mat};
  xtal::Lattice parent_superlattice{L};
  xtal::Lattice unmapped_child_prim_lattice{child.lat_column_mat};
  xtal::Lattice unmapped_child_superlattice{child.lat_column_mat};
  Index child_N_atom = 0;  // no longer used
  xtal::LatticeNode lattice_node(parent_prim_lattice, parent_superlattice,
                                 unmapped_child_prim_lattice,
                                 unmapped_child_superlattice, child_N_atom);
  std::set<xtal::MappingNode> mappings =
      mapper.map_deformed_struc_impose_lattice_node(
          child, lattice_node, k_best, max_cost, min_cost, keep_invalid);

  EXPECT_EQ(mappings.size(), 2);
  assert_mapping_relations(mappings, mapper, child);
}

TEST_F(StrucMapperTest, MapIdealStruc0) {
  xtal::StrucMapper mapper(calculator, lattice_weight, max_volume_change,
                           struc_mapper_options, cost_tol, min_va_frac,
                           max_va_frac);

  std::set<xtal::MappingNode> mappings =
      mapper.map_ideal_struc(child, k_best, max_cost, min_cost, keep_invalid);

  EXPECT_EQ(mappings.size(), 2);
  assert_mapping_relations(mappings, mapper, child);
}

// \brief Confirm that the MappingNode maps the unmapped child to the mapped
// child constructed by `mapper.resolve_setting` as expected
void assert_mapping_relations(std::set<xtal::MappingNode> const &mappings,
                              xtal::StrucMapper const &mapper,
                              xtal::SimpleStructure const &unmapped_child) {
  xtal::SimpleStructure const &parent = mapper.calculator().parent();

  for (xtal::MappingNode const &mapping : mappings) {
    xtal::SimpleStructure mapped_child =
        mapper.calculator().resolve_setting(mapping, unmapped_child);
    xtal::Lattice mapped_child_lattice(mapped_child.lat_column_mat);

    // lattice mapping relations:
    //
    //     L1 * T1 * N = V^{N} * Q^{N} * L2 * T2
    //
    // lattice_node.parent.prim_lattice() = L1
    // lattice_node.parent.superlattice() = L1 * T1 * N
    // lattice_node.child.prim_lattice() = V^{N} * Q^{N} * L2
    // lattice_node.child.superlattice() = V^{N} * Q^{N} * L2 * T2
    // S1 = L1 * T1
    // S2 = L2 * T2
    // mapped_child.lat_column_mat = Q^{N} * L2 * T2

    xtal::LatticeNode const &lattice_node = mapping.lattice_node;
    Eigen::Matrix3d S1xN = lattice_node.parent.superlattice().lat_column_mat();
    Eigen::Matrix3l T1xN = lattice_node.parent.transformation_matrix_to_super();
    Eigen::Matrix3d L2 = unmapped_child.lat_column_mat;
    Eigen::Matrix3l T2 = lattice_node.child.transformation_matrix_to_super();
    Eigen::Matrix3d S2 = L2 * T2.cast<double>();
    Eigen::Matrix3d V = lattice_node.stretch;
    Eigen::Matrix3d Q = lattice_node.isometry;
    Eigen::Matrix3d U_reverse = V.inverse();

    ASSERT_TRUE(almost_equal(S1xN, V * Q * S2));
    ASSERT_TRUE(almost_equal(mapped_child.lat_column_mat, Q * S2));

    xtal::SimpleStructure parent_superstructure =
        make_superstructure(T1xN.cast<int>(), parent);
    parent_superstructure.within();

    // atomic assignment relations:
    //
    //     mapped_child.mol_info.names[i] =
    //         unmapped_child.atom_info.names[perm[i]],
    //
    //     mapped_child.atom_info.coords.col(i) =
    //         parent_superstructure.coords(i) +
    //         mapping.atom_displacement.col(i) =
    //         mapping.lattice_node.stretch * mapping.lattice_node.isometry *
    //         unmapped_child.atom_info.coords.col(perm[i]) +
    //             mapping.atomic_node.translation,
    //
    //     where perm = mapping.atom_permutation.

    auto const &perm = mapping.atom_permutation;
    auto const &trans = mapping.atomic_node.translation;

    ASSERT_EQ(mapped_child.mol_info.size(),
              unmapped_child.atom_info.size() * T2.determinant());
    // TODO: update this ^ when mapping molecules

    for (Index i = 0; i < mapped_child.mol_info.size(); ++i) {
      // assert permutation maps atom names -> mol names
      ASSERT_EQ(mapped_child.mol_info.names[i],
                unmapped_child.atom_info.names[perm[i]]);
      // TODO: update this ^ when mapping molecules

      // compare coordinates using parent superlattice:
      {
        // r1[i] + disp[i] = V * Q * r2[perm[i]] + trans

        // mapped child coordinates, stretched to match parent superlattice
        // r = V * r[i]
        xtal::Coordinate mapped_child_coord(
            V * mapped_child.mol_info.coords.col(i),
            lattice_node.parent.superlattice(), CART);

        // displaced, but not deformed parent coordinate, in parent superlattice
        // r = r1[i] + disp[i]
        xtal::Coordinate displaced_parent_coord(
            parent_superstructure.atom_info.coords.col(i) +
                mapping.atom_displacement.col(i),
            lattice_node.parent.superlattice(), CART);

        // unmapped child coordinates, transformed to match parent superlattice:
        // r = V * Q * r2[perm[i]] + trans
        xtal::Coordinate transformed_child_coord(
            V * Q * unmapped_child.atom_info.coords.col(perm[i]) + trans,
            lattice_node.parent.superlattice(), CART);

        ASSERT_TRUE(transformed_child_coord.robust_min_dist(
                        displaced_parent_coord) < TOL);
        ASSERT_TRUE(mapped_child_coord.robust_min_dist(displaced_parent_coord) <
                    TOL);
        ASSERT_TRUE(
            mapped_child_coord.robust_min_dist(transformed_child_coord) < TOL);
      }

      // compare coordinates using mapped_child_lattice
      {
        // U_reverse * (r1[i] + disp[i]) = Q * r2[perm[i]] + U_reverse * trans

        // mapped child coordinates, stretched to match parent superlattice
        // r = r[i]
        xtal::Coordinate mapped_child_coord(mapped_child.mol_info.coords.col(i),
                                            mapped_child_lattice, CART);

        // displaced, but not deformed parent coordinate, in parent superlattice
        // r = U_reverse * (r1[i] + disp[i])
        xtal::Coordinate displaced_parent_coord(
            U_reverse * (parent_superstructure.atom_info.coords.col(i) +
                         mapping.atom_displacement.col(i)),
            mapped_child_lattice, CART);

        // unmapped child coordinates, transformed to match parent superlattice:
        // r = Q * r2[perm[i]] + U_reverse * trans
        xtal::Coordinate transformed_child_coord(
            Q * unmapped_child.atom_info.coords.col(perm[i]) +
                U_reverse * trans,
            mapped_child_lattice, CART);

        ASSERT_TRUE(transformed_child_coord.robust_min_dist(
                        displaced_parent_coord) < TOL);
        ASSERT_TRUE(mapped_child_coord.robust_min_dist(displaced_parent_coord) <
                    TOL);
        ASSERT_TRUE(
            mapped_child_coord.robust_min_dist(transformed_child_coord) < TOL);
      }
    }
  }
}
