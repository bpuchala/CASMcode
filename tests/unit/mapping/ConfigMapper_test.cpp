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
#include "gtest/gtest.h"

// BCC mapping tests

using namespace CASM;

// // \brief Confirm that the ConfigMapperResult results map as expected
// void assert_mapping_relations(ConfigMapperResult const &config_mapper_result,
//                               xtal::StrucMapper const &mapper,
//                               xtal::SimpleStructure const &unmapped_child);

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

  jsonParser json;
  to_json(result, json, CART);
  std::cout << json << std::endl;

  EXPECT_EQ(result.maps.size(), 1);
}

// binary {A, B} BCC prim structure
xtal::BasicStructure ConfigMapperBCCTest::make_bcc_prim() {
  using namespace xtal;

  Molecule A = Molecule::make_atom("A");
  Molecule B = Molecule::make_atom("B");

  Lattice lat{Eigen::Vector3d{4.0, 0.0, 0.0}, Eigen::Vector3d{0.0, 4.0, 0.0},
              Eigen::Vector3d{0.0, 0.0, 4.0}};

  BasicStructure struc{lat};
  struc.set_basis({Site{Coordinate{0.0, 0.0, 0.0, lat, CART}, {A, B}},
                   Site{Coordinate{2.0, 2.0, 2.0, lat, CART}, {A, B}}});

  BasicStructure prim = make_primitive(struc);
  prim.set_lattice(xtal::canonical::equivalent(prim.lattice()), CART);
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
    0., 0., 4.;

  simple.atom_info.coords = Eigen::MatrixXd(3, 4);
  simple.atom_info.coords <<
    7.98, 1.99, 4.01, 6.02,
    0., 2., 0., 2.,
    0., 2., 0., 2.;

  simple.atom_info.names = {"A", "A", "A", "A"};
  // clang-format on

  return simple;
}
