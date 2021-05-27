#include "casm/crystallography/LatticeMap.hh"

#include "autotools.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "gtest/gtest.h"

using namespace CASM;

TEST(LatticeMapTest, Test1) {
  // parent
  Eigen::Matrix3d L1;
  L1 << 3.233986860000, -1.616993430000, 0.000000000000,  //
      0.000000000000, 2.800714770000, 0.000000000000,     //
      0.000000000000, 0.000000000000, 10.337356680000;    //

  xtal::Lattice L1_lattice(L1);
  auto L1_point_group = xtal::make_point_group(L1_lattice);

  // child
  Eigen::Matrix3d L2;
  L2 << 3.269930775653, 0.000000000000, 0.000000000000,  //
      -1.634965387827, 2.831843113861, 0.000000000000,   //
      0.000000000000, 0.000000000000, 10.464806115486;   //

  xtal::Lattice L2_lattice(L2);
  auto L2_point_group = xtal::make_point_group(L2_lattice);

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

  EXPECT_TRUE(
      almost_equal(lattice_map.parent_matrix(), L1_lattice.lat_column_mat()));
  EXPECT_TRUE(
      almost_equal(lattice_map.child_matrix(), L2_lattice.lat_column_mat()));

  Index count = 0;
  do {
    count++;
    EXPECT_TRUE(almost_equal(
        lattice_map.child_matrix(),
        lattice_map.deformation_gradient() * L1 * lattice_map.matrixN()));

    lattice_map.next_mapping_better_than(lattice_map.strain_cost());
  } while (lattice_map.strain_cost() < (lattice_map.strain_cost() + TOL));

  // std::cout << "count: " << count << std::endl;
  // EXPECT_EQ(count, 12);

  Eigen::Matrix3d V_expected;
  V_expected << 1.144098197277, -0.032761728177, 0.000000000000,  //
      -0.032761728177, 0.855879024393, 0.000000000000,            //
      0.000000000000, 0.000000000000, 0.987821137432;             //

  Eigen::Matrix3d Q_expected;
  Q_expected << 0.0382504439139, 0.999268183993, 0.000000000000,  //
      -0.999268183993, 0.038250443914, 0.000000000000,            //
      0.000000000000, 0.000000000000, 1.000000000000;             //

  // using LatticeNode defintions for stretch=V, isometry=Q
  Eigen::Matrix3d F_reverse = lattice_map.deformation_gradient();
  Eigen::Matrix3d N = lattice_map.matrixN();
  Eigen::Matrix3d stretch = polar_decomposition(F_reverse).inverse();
  Eigen::Matrix3d isometry = (F_reverse * stretch).transpose();

  EXPECT_TRUE(almost_equal(L1 * N, stretch * isometry * L2));

  EXPECT_TRUE(almost_equal(stretch, V_expected))
      << "Unexpected stretch. Found stretch:\n"
      << stretch;
  EXPECT_TRUE(almost_equal(isometry, Q_expected))
      << "Unexpected isometry. Found isometry:\n"
      << isometry;
}
