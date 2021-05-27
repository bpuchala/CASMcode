#include "casm/crystallography/LatticeMap.hh"

#include "autotools.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "gtest/gtest.h"

using namespace CASM;

TEST(LatticeMapTest, Test1) {
  // Eigen::Matrix3d V;
  // V << 0.998559791924, -0.000000000970, 0.000000000000,  //
  //     -0.000000000970, 0.998559788571, 0.000000000000,   //
  //     0.000000000000, 0.000000000000, 0.997675852223;    //
  //
  // Eigen::Matrix3d Q;
  // Q << -0.500000000000, -0.866025403784, 0.000000000000,  //
  //     -0.866025403784, 0.500000000000, 0.000000000000,    //
  //     0.000000000000, 0.000000000000, -1.000000000000;    //

  // parent
  Eigen::Matrix3d L1;
  // L1 << 3.233986860000, -1.616993430000, 0.000000000000,  //
  //     0.000000000000, 2.800714770000, 0.000000000000,     //
  //     0.000000000000, 0.000000000000, 5.168678340000;     //
  L1 << 3.233986860000, -1.616993430000, 0.000000000000,  //
      0.000000000000, 2.800714770000, 0.000000000000,     //
      0.000000000000, 0.000000000000, 10.337356680000;    //

  xtal::Lattice L1_lattice(L1);
  auto L1_point_group = xtal::make_point_group(L1_lattice);

  // child
  Eigen::Matrix3d L2;
  // L2 << 3.238651197049, 0.000000000000, 0.000000000000,  //
  //     -1.619325598525, 2.804754204349, 0.000000000000,   //
  //     0.000000000000, 0.000000000000, 5.180719096772;    //
  L2 << 3.269930775653, 0.000000000000, 0.000000000000,  //
      -1.634965387827, 2.831843113861, 0.000000000000,   //
      0.000000000000, 0.000000000000, 10.464806115486;   //

  xtal::Lattice L2_lattice(L2);
  auto L2_point_group = xtal::make_point_group(L2_lattice);

  std::cout << "L1: \n" << L1 << std::endl;
  std::cout << "L2: \n" << L2 << std::endl;

  // this is no longer used by LatticeMap
  Index n_atoms = 0;

  // range of elements of N matrix (controls number of potential mappings to be
  /// considered... larger is more)
  int unimodular_element_range = 1;

  // default
  Eigen::MatrixXd strain_gram_mat = Eigen::MatrixXd::Identity(9, 9);

  double max_lattice_cost = 1e20;

  bool use_symmetry_breaking_strain_cost = false;

  xtal::LatticeMap lattice_map{L1_lattice,
                               L2_lattice,
                               n_atoms,
                               unimodular_element_range,
                               L1_point_group,
                               L2_point_group,
                               strain_gram_mat,
                               max_lattice_cost,
                               use_symmetry_breaking_strain_cost};

  do {
    std::cout << "---" << std::endl << std::endl;
    std::cout << "strain_cost: " << lattice_map.strain_cost() << std::endl;
    std::cout << "F: \n" << lattice_map.deformation_gradient() << std::endl;
    std::cout << "N: \n" << lattice_map.matrixN() << std::endl;
    std::cout << "reduced_parent: \n"
              << lattice_map.reduced_parent_matrix() << std::endl;
    std::cout << "reduced_child: \n"
              << lattice_map.reduced_child_matrix() << std::endl
              << std::endl;

    // child = deformation_gradient*parent*N
    std::cout << lattice_map.deformation_gradient() * L1 * lattice_map.matrixN()
              << std::endl
              << std::endl;

    lattice_map.next_mapping_better_than(max_lattice_cost);
  } while (lattice_map.strain_cost() < (max_lattice_cost + TOL));

  // EXPECT_TRUE(almost_equal(L1, V * Q * L2));
}
