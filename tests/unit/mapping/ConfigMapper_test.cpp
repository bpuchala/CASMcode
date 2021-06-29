#include "autotools.hh"
#include "casm/casm_io/container/json_io.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/clex/ConfigMapping.hh"
#include "casm/clex/SimpleStructureTools.hh"
#include "casm/clex/io/json/ConfigMapping_json_io.hh"
#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/BasicStructureTools.hh"
#include "casm/crystallography/CanonicalForm.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/SimpleStructureTools.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/crystallography/io/SimpleStructureIO.hh"
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
void assert_mapping_relations(ConfigMapperResult const &config_mapper_result,
                              xtal::StrucMapper const &mapper,
                              xtal::SimpleStructure const &unmapped_child);

// \brief Confirm that a ConfigurationMapping maps as expected
void assert_mapping_relations(ConfigurationMapping const &map,
                              xtal::SimpleStructure const &parent,
                              xtal::SimpleStructure const &unmapped_child);

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

  EXPECT_EQ(result.maps.size(), 1);
}

// \brief Confirm that the ConfigMapperResult results map as expected
void assert_mapping_relations(ConfigMapperResult const &config_mapper_result,
                              ConfigMapper const &config_mapper,
                              xtal::SimpleStructure const &unmapped_child) {
  auto const &parent = config_mapper.struc_mapper().parent();
  for (ConfigurationMapping const &map : config_mapper_result.maps) {
    assert_mapping_relations(map, parent, unmapped_child);
  }
}

// \brief Confirm that a ConfigurationMapping maps as expected
void assert_mapping_relations(ConfigurationMapping const &map,
                              xtal::SimpleStructure const &parent,
                              xtal::SimpleStructure const &unmapped_child) {
  // unmapped_child -> mapped_child
  assert_mapping_relations(map.mapping, parent, map.mapped_child,
                           unmapped_child);

  // // -> mapped_configuration, mapped_properties
  // assert_consistency(mapped_child, mapped_configuration, mapped_properties);

  // // apply symop_to_canon_scel, permutation_to_canon_scel,
  // // transformation_matrix_to_canon_scel
  // // -> configuration_in_canon_scel, properties_in_canon_scel,
  // // structure_in_canon_scel
  // assert_transformation(symop_to_canon_scel, permutation_to_canon_scel,
  //                       transformation_matrix_to_canon_scel,
  //                       structure_in_canon_scel, mapped_child);
  // assert_consistency(structure_in_canon_scel, configuration_in_canon_scel,
  //                    properties_in_canon_scel);
  //
  // // apply symop_to_final, permutation_to_final,
  // // -> final_configuration, final_properties, final_structure
  // assert_transformation(configuration_in_canon_scel, symop_to_final,
  //                       permutation_to_final, final_configuration);
  // assert_consistency(final_structure, final_configuration, final_properties);
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
