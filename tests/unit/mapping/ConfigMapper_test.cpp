#include "autotools.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/container/stream_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clex/ConfigMapping.hh"
#include "casm/clex/SimpleStructureTools.hh"
#include "casm/clex/Supercell.hh"
#include "casm/clex/io/json/ConfigMapping_json_io.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/CanonicalForm.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/io/SimpleStructureIO.hh"
#include "casm/crystallography/io/UnitCellCoordIO.hh"
#include "casm/misc/algorithm.hh"
#include "casm/strain/StrainConverter.hh"
#include "gtest/gtest.h"

// BCC mapping tests

using namespace CASM;

// \brief Confirm that the MappingNode maps the unmapped child to the mapped
// child
//
// From StrucMapper_test.cpp
void assert_mapping_relations(xtal::MappingNode const &mapping,
                              xtal::SimpleStructure const &parent,
                              xtal::SimpleStructure const &mapped_child,
                              xtal::SimpleStructure const &unmapped_child);

// \brief Confirm that the ConfigMapperResult results map as expected
void assert_mapping_relations(
    ConfigMapperResult const &config_mapper_result,
    ConfigMapperSettings const &config_mapper_settings,
    ConfigMapper const &mapper, xtal::SimpleStructure const &unmapped_child);

// \brief Confirm that a ConfigurationMapping maps as expected
void assert_mapping_relations(ConfigurationMapping const &map,
                              xtal::SimpleStructure const &parent,
                              xtal::SimpleStructure const &unmapped_child);

// Check that lattice and strain are consistent among structure, configuration,
// and properties
void assert_strain_consistency(xtal::SimpleStructure const &structure,
                               Configuration const &configuration,
                               MappedProperties const &properties,
                               Eigen::Matrix3d &deformation_gradient);

// Check that coordinates and displacements are consistent among structure,
// configuration, and properties
void assert_disp_consistency(xtal::SimpleStructure const &structure,
                             Configuration const &configuration,
                             MappedProperties const &properties,
                             Eigen::Matrix3d const &deformation_gradient);

// Check that structure, configuration, and properties are consistent
void assert_consistency(xtal::SimpleStructure const &structure,
                        Configuration const &configuration,
                        MappedProperties const &properties);

// Check structure transformation
void assert_structure_transformation(
    SymOp const &op, Permutation const &mol_permutation,
    Eigen::Matrix3l const &transformation_matrix_to_super,
    xtal::SimpleStructure const &final_child,
    xtal::SimpleStructure const &init_child);

// Check configuration transformation
void assert_configuration_transformation(
    SymOp const &op, Permutation const &permutation,
    Eigen::Matrix3l const &transformation_matrix_to_super,
    Configuration const &final_configuration,
    Configuration const &init_configuration);

class ConfigMapperBCCTest : public testing::Test {
 protected:
  // prim: (reference structure)
  std::shared_ptr<Structure const> shared_prim;

  // child: (what is being mapped to prim)
  xtal::SimpleStructure child;

  ConfigMapperSettings config_mapper_settings;

  ConfigMapperBCCTest()
      : shared_prim(std::make_shared<Structure const>(make_bcc_prim())),
        child(make_bcc_simplestructure_2x1x1()) {}

  // binary {A, B} BCC prim structure
  static xtal::BasicStructure make_bcc_prim();

  // conventional 1x1x1 bcc cell
  static xtal::SimpleStructure make_bcc_simplestructure_1x1x1();

  // conventional 2x1x1 bcc cell
  static xtal::SimpleStructure make_bcc_simplestructure_2x1x1();
};

TEST_F(ConfigMapperBCCTest, Test0) {
  ConfigMapper mapper(shared_prim, config_mapper_settings);
  ConfigMapperResult result = mapper.import_structure(child);

  EXPECT_EQ(result.maps.size(), 1);

  jsonParser json;
  json.put_obj();
  to_json(child, json["child"], {}, CART);
  std::cout << json << std::endl;

  json.put_obj();
  to_json(result, json, CART);
  std::cout << json << std::endl;

  assert_mapping_relations(result, config_mapper_settings, mapper, child);
}

// \brief Confirm that the ConfigMapperResult results map as expected
void assert_mapping_relations(
    ConfigMapperResult const &config_mapper_result,
    ConfigMapperSettings const &config_mapper_settings,
    ConfigMapper const &config_mapper,
    xtal::SimpleStructure const &unmapped_child) {
  auto const &parent = config_mapper.struc_mapper().parent();
  for (ConfigurationMapping const &map : config_mapper_result.maps) {
    assert_mapping_relations(map, parent, unmapped_child);

    if (config_mapper_settings.finalize_strict == false) {
      ASSERT_TRUE(map.final_configuration.is_canonical());
    }
  }
}

// \brief Confirm that a ConfigurationMapping maps as expected
void assert_mapping_relations(ConfigurationMapping const &map,
                              xtal::SimpleStructure const &parent,
                              xtal::SimpleStructure const &unmapped_child) {
  // unmapped_child -> mapped_child
  assert_mapping_relations(map.mapping, parent, map.mapped_child,
                           unmapped_child);

  // -> mapped_configuration, mapped_properties
  // std::cout << "\nCHECK MAPPED:\n" << std::endl;
  assert_consistency(map.mapped_child, map.mapped_configuration,
                     map.mapped_properties);

  // // apply symop_to_canon_scel, permutation_to_canon_scel,
  // // transformation_matrix_to_canon_scel
  // // -> configuration_in_canon_scel, properties_in_canon_scel,
  // // structure_in_canon_scel
  // std::cout << "\nCHECK IN_CANON_SCEL:\n" << std::endl;
  ASSERT_TRUE(map.configuration_in_canon_scel.supercell().is_canonical());
  assert_structure_transformation(
      map.symop_to_canon_scel, map.permutation_to_canon_scel,
      map.transformation_matrix_to_canon_scel, map.structure_in_canon_scel,
      map.mapped_child);
  assert_consistency(map.structure_in_canon_scel,
                     map.configuration_in_canon_scel,
                     map.properties_in_canon_scel);

  // // apply symop_to_final, permutation_to_final,
  // // -> final_configuration, final_properties, final_structure
  // std::cout << "\nCHECK FINAL:\n" << std::endl;
  assert_configuration_transformation(
      map.symop_to_final, map.permutation_to_final,
      map.transformation_matrix_to_final, map.final_configuration,
      map.configuration_in_canon_scel);
  // assert_transformation(symop_to_final, permutation_to_final,
  // final_configuration, configuration_in_canon_scel);
  assert_consistency(map.final_structure, map.final_configuration,
                     map.final_properties);
}

// Check that lattice and strain are consistent among structure, configuration,
// and properties
//
// \param deformation gradient Gets set from DoF or properties
//
void assert_strain_consistency(xtal::SimpleStructure const &structure,
                               Configuration const &configuration,
                               MappedProperties const &properties,
                               Eigen::Matrix3d &deformation_gradient) {
  auto const &prim = configuration.supercell().prim();

  // lattice:
  std::string strain_key;  // "Ustrain", "Bstrain", etc.
  Eigen::VectorXd unrolled_metric;

  // strain is in one of configuration or properties:
  ASSERT_NE(xtal::has_strain_dof(prim), has_strain_property(properties));

  if (xtal::has_strain_dof(prim)) {
    strain_key = xtal::get_strain_dof_key(prim);
    unrolled_metric =
        configuration.configdof().global_dof(strain_key).standard_values();
  } else {
    strain_key = get_strain_property_key(properties);
    unrolled_metric = properties.global.at(strain_key);
  }

  StrainConverter c(xtal::get_strain_metric(strain_key));
  Eigen::Matrix3d F = c.unrolled_strain_metric_to_F(unrolled_metric);
  Eigen::Matrix3d lattice_vectors_with_strain =
      F * configuration.ideal_lattice().lat_column_mat();

  ASSERT_TRUE(
      almost_equal(structure.lat_column_mat, lattice_vectors_with_strain));

  deformation_gradient = F;
}

// Check that coordinates and displacements are consistent among structure,
// configuration, and properties
//
// - checks structure.mol_info, not atom_info
void assert_disp_consistency(xtal::SimpleStructure const &structure,
                             Configuration const &configuration,
                             MappedProperties const &properties,
                             Eigen::Matrix3d const &deformation_gradient) {
  auto const &prim = configuration.supercell().prim();

  // F * (r_ideal + disp) = r_structure

  // coordinates:

  // ideal coords:
  Index n_sites = configuration.size();
  auto f = configuration.supercell().sym_info().unitcellcoord_index_converter();
  Eigen::MatrixXd r_ideal(3, n_sites);
  for (Index l = 0; l < n_sites; l++) {
    // jsonParser json;
    // json = f(l);
    // std::cout << "l: " << l << "  f(l): " << json << "  coord: " <<
    // f(l).coordinate(prim).const_cart().transpose() << std::endl;
    r_ideal.col(l) = f(l).coordinate(prim).const_cart();
  }

  // structure coordinates:
  ASSERT_EQ(structure.mol_info.size(), configuration.size());
  Eigen::MatrixXd r_structure = structure.mol_info.coords;

  // displacements:

  // disp is in one of configuration or properties:
  auto local_dof_types = xtal::all_local_dof_types(prim);
  bool has_disp_dof = contains(local_dof_types, "disp");
  bool has_disp_property = properties.site.count("disp");
  ASSERT_NE(has_disp_dof, has_disp_property);

  Eigen::MatrixXd disp;
  if (has_disp_dof) {
    disp = configuration.configdof().local_dof("disp").standard_values();
  } else {
    disp = properties.site.at("disp");
  }

  // check coord equivalence, accounting for periodic boundaries
  Eigen::MatrixXd r_lhs = deformation_gradient * (r_ideal + disp);
  xtal::Lattice lattice{structure.lat_column_mat};
  // std::cout << "lat_column_mat: \n" << structure.lat_column_mat << std::endl;
  for (Index l = 0; l < n_sites; l++) {
    xtal::Coordinate lhs_coord(r_lhs.col(l), lattice, CART);
    xtal::Coordinate rhs_coord(r_structure.col(l), lattice, CART);
    // std::cout << "l: " << l << "  lhs: " <<
    // lhs_coord.const_cart().transpose() << "  rhs: " <<
    // rhs_coord.const_cart().transpose() << std::endl;
    ASSERT_TRUE(lhs_coord.robust_min_dist(rhs_coord) < TOL);
  }
}

// Check that structure, configuration, and properties are consistent
void assert_consistency(xtal::SimpleStructure const &structure,
                        Configuration const &configuration,
                        MappedProperties const &properties) {
  Eigen::Matrix3d F;
  assert_strain_consistency(structure, configuration, properties, F);
  assert_disp_consistency(structure, configuration, properties, F);
}

// Check structure transformation
//
// Checks:
// - lattice transformation
// - coordinate transformation
// - global properties transformation
// - site properties (mol_info) transformation
//
// Does not check:
// - atom properties (atom_info) transformation
// - anisotropic occ indices (yet) // TODO: anisotropic occ
void assert_structure_transformation(
    SymOp const &op, Permutation const &mol_permutation,
    Eigen::Matrix3l const &transformation_matrix_to_super,
    xtal::SimpleStructure const &final_child,
    xtal::SimpleStructure const &init_child) {
  // check lattice transformation
  // L_final = op.matrix() * L_init * T
  // std::cout << "\nCHECK LATTICE:\n" << std::endl;
  // std::cout << "L_final:\n" << final_child.lat_column_mat << std::endl;
  // std::cout << "L_init:\n" << init_child.lat_column_mat << std::endl;
  // std::cout << "op.matrix():\n" << op.matrix() << std::endl;
  // std::cout << "T:\n" << transformation_matrix_to_super.cast<double>() <<
  // std::endl; std::cout << "L_transformed:\n" << op.matrix() *
  // init_child.lat_column_mat * transformation_matrix_to_super.cast<double>()
  // << std::endl;
  ASSERT_TRUE(almost_equal(final_child.lat_column_mat,
                           op.matrix() * init_child.lat_column_mat *
                               transformation_matrix_to_super.cast<double>()));

  // check coordinate transformation (accounting for periodic boundaries)
  // r_final[l] = op.matrix() * r_init[mol_permutation[l]]
  // std::cout << "\nCHECK COORDS:\n" << std::endl;
  // std::cout << "perm: " << mol_permutation << std::endl;
  xtal::Lattice final_lattice{final_child.lat_column_mat};
  for (Index l = 0; l < final_child.mol_info.size(); l++) {
    xtal::Coordinate r_final{final_child.mol_info.coords.col(l), final_lattice,
                             CART};
    xtal::Coordinate r_transformed{
        op.matrix() * init_child.mol_info.coords.col(mol_permutation[l]),
        final_lattice, CART};
    // std::cout << "r_final:" << r_final.const_cart().transpose() << std::endl;
    // std::cout << "r_transformed:" << r_transformed.const_cart().transpose()
    // << std::endl;
    ASSERT_TRUE(r_final.robust_min_dist(r_transformed) < TOL);
  }

  // check transformation of global properties
  // std::cout << "\nCHECK GLOBAL PROPERTIES:\n" << std::endl;
  for (auto const &pair : init_child.properties) {
    // std::cout << "checking: " << pair.first << std::endl;
    ASSERT_TRUE(final_child.properties.count(pair.first));

    Eigen::MatrixXd final_value = final_child.properties.at(pair.first);
    Eigen::MatrixXd init_value = init_child.properties.at(pair.first);
    Eigen::MatrixXd matrix =
        AnisoValTraits(pair.first)
            .symop_to_matrix(op.matrix(), op.tau(), op.time_reversal());
    // std::cout << "final_value:" << final_value.transpose() << std::endl;
    // std::cout << "init_value:" << init_value.transpose() << std::endl;
    // std::cout << "matrix:" << init_value.transpose() << std::endl;
    // std::cout << "transformed_value:" << (matrix * init_value).transpose() <<
    // std::endl;
    ASSERT_TRUE(almost_equal(final_value, matrix * init_value));
  }

  // check transformation and permutation of local properties
  // std::cout << "\nCHECK SITE PROPERTIES:\n" << std::endl;
  ASSERT_EQ(final_child.mol_info.size(), init_child.mol_info.size());
  ASSERT_EQ(final_child.mol_info.size(), mol_permutation.size());
  for (auto const &pair : init_child.mol_info.properties) {
    // std::cout << "checking: " << pair.first << std::endl;
    ASSERT_TRUE(final_child.mol_info.properties.count(pair.first));

    for (Index l = 0; l < final_child.mol_info.size(); l++) {
      Eigen::MatrixXd final_value =
          final_child.mol_info.properties.at(pair.first).col(l);
      Eigen::MatrixXd init_value =
          init_child.mol_info.properties.at(pair.first).col(mol_permutation[l]);
      Eigen::MatrixXd matrix =
          AnisoValTraits(pair.first)
              .symop_to_matrix(op.matrix(), op.tau(), op.time_reversal());
      // std::cout << "\nl: " << l << std::endl;
      // std::cout << "final_value:" << final_value.transpose() << std::endl;
      // std::cout << "init_value:" << init_value.transpose() << std::endl;
      // std::cout << "matrix:" << init_value.transpose() << std::endl;
      // std::cout << "transformed_value:" << (matrix * init_value).transpose()
      // << std::endl;
      ASSERT_TRUE(almost_equal(final_value, matrix * init_value));
    }
  }
}

// Check configuration transformation
//
// Checks:
// - lattice
// - occ (not valid for anisotropic yet) TODO: anisotropic occ
// - global dof
// - site dof
//
void assert_configuration_transformation(
    SymOp const &op, Permutation const &permutation,
    Eigen::Matrix3l const &transformation_matrix_to_super,
    Configuration const &final_configuration,
    Configuration const &init_configuration) {
  // check lattice transformation
  // L_final = op.matrix() * L_init * T
  // std::cout << "\nCHECK LATTICE:\n" << std::endl;
  // std::cout << "L_final:\n" <<
  // final_configuration.ideal_lattice().lat_column_mat() << std::endl;
  // std::cout << "L_init:\n" <<
  // init_configuration.ideal_lattice().lat_column_mat() << std::endl; std::cout
  // << "op.matrix():\n" << op.matrix() << std::endl; std::cout << "T:\n" <<
  // transformation_matrix_to_super.cast<double>() << std::endl; std::cout <<
  // "L_transformed:\n" << op.matrix() *
  // init_configuration.ideal_lattice().lat_column_mat() *
  // transformation_matrix_to_super.cast<double>() << std::endl;
  ASSERT_TRUE(almost_equal(
      final_configuration.ideal_lattice().lat_column_mat(),
      op.matrix() * init_configuration.ideal_lattice().lat_column_mat() *
          transformation_matrix_to_super.cast<double>()));

  // check occ // TODO: check anisotropic occ
  for (Index l = 0; l < permutation.size(); l++) {
    ASSERT_EQ(final_configuration.occ(l),
              init_configuration.occ(permutation[l]));
  }

  // check transformation of global dof
  // std::cout << "\nCHECK GLOBAL DOF:\n" << std::endl;
  for (auto const &pair : init_configuration.configdof().global_dofs()) {
    // std::cout << "checking: " << pair.first << std::endl;
    ASSERT_TRUE(
        final_configuration.configdof().global_dofs().count(pair.first));

    Eigen::MatrixXd final_value = final_configuration.configdof()
                                      .global_dof(pair.first)
                                      .standard_values();
    Eigen::MatrixXd init_value =
        init_configuration.configdof().global_dof(pair.first).standard_values();
    Eigen::MatrixXd matrix =
        AnisoValTraits(pair.first)
            .symop_to_matrix(op.matrix(), op.tau(), op.time_reversal());
    // std::cout << "final_value:" << final_value.transpose() << std::endl;
    // std::cout << "init_value:" << init_value.transpose() << std::endl;
    // std::cout << "matrix:" << init_value.transpose() << std::endl;
    // std::cout << "transformed_value:" << (matrix * init_value).transpose() <<
    // std::endl;
    ASSERT_TRUE(almost_equal(final_value, matrix * init_value));
  }

  // check transformation and permutation of local dof
  // std::cout << "\nCHECK SITE DOF:\n" << std::endl;
  auto const &final_local_dofs = final_configuration.configdof().local_dofs();
  auto const &init_local_dofs = init_configuration.configdof().local_dofs();
  ASSERT_EQ(final_local_dofs.size(), init_local_dofs.size());
  ASSERT_EQ(final_configuration.size(), permutation.size());
  for (auto const &pair : init_local_dofs) {
    // std::cout << "checking: " << pair.first << std::endl;
    ASSERT_TRUE(final_local_dofs.count(pair.first));

    for (Index l = 0; l < permutation.size(); l++) {
      Eigen::MatrixXd final_value =
          final_local_dofs.at(pair.first).standard_values().col(l);
      Eigen::MatrixXd init_value =
          init_local_dofs.at(pair.first).standard_values().col(permutation[l]);
      Eigen::MatrixXd matrix =
          AnisoValTraits(pair.first)
              .symop_to_matrix(op.matrix(), op.tau(), op.time_reversal());
      // std::cout << "\nl: " << l << std::endl;
      // std::cout << "final_value:" << final_value.transpose() << std::endl;
      // std::cout << "init_value:" << init_value.transpose() << std::endl;
      // std::cout << "matrix:" << init_value.transpose() << std::endl;
      // std::cout << "transformed_value:" << (matrix * init_value).transpose()
      // << std::endl;
      ASSERT_TRUE(almost_equal(final_value, matrix * init_value));
    }
  }
}

// binary {A, B} BCC prim structure
xtal::BasicStructure ConfigMapperBCCTest::make_bcc_prim() {
  using namespace xtal;

  Molecule A = Molecule::make_atom("A");
  Molecule B = Molecule::make_atom("B");

  Lattice lat{Eigen::Vector3d{4.0, 0.0, 0.0}, Eigen::Vector3d{0.0, 4.0, 0.0},
              Eigen::Vector3d{0.0, 0.0, 4.0}};

  SiteDoFSet disp{AnisoValTraits::disp()};
  BasicStructure struc{lat};
  struc.set_basis({Site{Coordinate{0.0, 0.0, 0.0, lat, CART}, {A, B}, {disp}},
                   Site{Coordinate{2.0, 2.0, 2.0, lat, CART}, {A, B}, {disp}}});

  // Add global DoF
  // GLstrain: Green-Lagrange strain
  struc.set_global_dofs({AnisoValTraits::strain("GL")});

  jsonParser json;
  json = jsonParser::object();
  std::cout << "struc: " << to_json(xtal::global_dof_types(struc), json)
            << std::endl;

  BasicStructure prim = make_primitive(struc);
  prim.set_lattice(xtal::canonical::equivalent(prim.lattice()), CART);

  json = jsonParser::object();
  std::cout << "prim: " << to_json(xtal::global_dof_types(prim), json)
            << std::endl;
  return prim;
}

// conventional 1x1x1 bcc cell
xtal::SimpleStructure ConfigMapperBCCTest::make_bcc_simplestructure_1x1x1() {
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

// conventional 2x1x1 bcc cell
xtal::SimpleStructure ConfigMapperBCCTest::make_bcc_simplestructure_2x1x1() {
  xtal::SimpleStructure simple;

  // clang-format off
  simple.lat_column_mat <<
    8., 0., 0.,
    0., 4., 0.,
    0., 0., 4.1;

  simple.atom_info.coords = Eigen::MatrixXd(3, 4);
  simple.atom_info.coords <<
    0., 2., 4.01, 6.01,
    0., 2., 0., 2.,
    0., 2., 0., 2.;

  simple.atom_info.names = {"A", "A", "A", "B"};
  // clang-format on

  return simple;
}
