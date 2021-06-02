#include "casm/crystallography/LatticeMap.hh"

#include "autotools.hh"
#include "casm/crystallography/Strain.hh"
#include "casm/crystallography/SymTools.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "gtest/gtest.h"

using namespace CASM;

class StrainCostTest : public testing::Test {
 protected:
  xtal::SymOpVector cubic_point_group;
  xtal::SymOpVector tetragonal_point_group;
  xtal::SymOpVector hexagonal_point_group;

  StrainCostTest();

  static Eigen::Matrix3d identity() { return Eigen::Matrix3d::Identity(); }

  static Eigen::Matrix3d tetragonal_stretch() {
    Eigen::Matrix3d V;
    V << 1.01, 0.0, 0.0,  //
        0.0, 1.0, 0.0,    //
        0.0, 0.0, 1.0;    //
    return V;
  }

  static Eigen::Matrix3d shear_stretch() {
    Eigen::Matrix3d V;
    V << 1.0, 0.01, 0.0,  //
        0.01, 1.0, 0.0,   //
        0.0, 0.0, 1.0;    //
    return V;
  }

  static xtal::Lattice hexagonal_lattice() {
    Eigen::Matrix3d L;
    double a = 1.;
    double c = 2. * sqrt(2.) / sqrt(3.);
    L << a, a / 2., 0.0,              //
        0.0, a * sqrt(3.) / 2., 0.0,  //
        0.0, 0.0, c;                  //
    return xtal::Lattice(L);
  }
};

StrainCostTest::StrainCostTest()
    : cubic_point_group(xtal::make_point_group(xtal::Lattice(identity()))),
      tetragonal_point_group(
          xtal::make_point_group(xtal::Lattice(tetragonal_stretch()))),
      hexagonal_point_group(xtal::make_point_group(hexagonal_lattice())) {}

TEST_F(StrainCostTest, IdentityTest) {
  // no strain
  Eigen::Matrix3d F = identity();
  double cost;

  cost = xtal::isotropic_strain_cost(F);
  EXPECT_TRUE(almost_equal(cost, 0.0)) << "isotropic cost: " << cost;

  cost = xtal::symmetry_breaking_strain_cost(F, cubic_point_group);
  EXPECT_TRUE(almost_equal(cost, 0.0))
      << "cubic symmetry-breaking cost: " << cost;

  cost = xtal::symmetry_breaking_strain_cost(F, tetragonal_point_group);
  EXPECT_TRUE(almost_equal(cost, 0.0))
      << "tetragonal symmetry-breaking cost: " << cost;

  cost = xtal::symmetry_breaking_strain_cost(F, hexagonal_point_group);
  EXPECT_TRUE(almost_equal(cost, 0.0))
      << "hexagonal symmetry-breaking cost: " << cost;
}

TEST_F(StrainCostTest, IsotropicTest) {
  // volume scaled
  Eigen::Matrix3d F = 1.1 * identity();
  double cost;

  cost = xtal::isotropic_strain_cost(F);
  EXPECT_TRUE(almost_equal(cost, 0.0)) << "isotropic cost: " << cost;

  cost = xtal::symmetry_breaking_strain_cost(F, cubic_point_group);
  EXPECT_TRUE(almost_equal(cost, 0.0))
      << "cubic symmetry-breaking cost: " << cost;

  cost = xtal::symmetry_breaking_strain_cost(F, tetragonal_point_group);
  EXPECT_TRUE(almost_equal(cost, 0.0))
      << "tetragonal symmetry-breaking cost: " << cost;

  cost = xtal::symmetry_breaking_strain_cost(F, hexagonal_point_group);
  EXPECT_TRUE(almost_equal(cost, 0.0))
      << "hexagonal symmetry-breaking cost: " << cost;
}

TEST_F(StrainCostTest, TetragonalTest) {
  // tetragonal stretch
  Eigen::Matrix3d F = tetragonal_stretch();
  double cost;

  cost = xtal::isotropic_strain_cost(F);
  EXPECT_FALSE(almost_equal(cost, 0.0)) << "isotropic cost: " << cost;

  cost = xtal::symmetry_breaking_strain_cost(F, cubic_point_group);
  EXPECT_FALSE(almost_equal(cost, 0.0))
      << "cubic symmetry-breaking cost: " << cost;

  cost = xtal::symmetry_breaking_strain_cost(F, tetragonal_point_group);
  EXPECT_TRUE(almost_equal(cost, 0.0))
      << "tetragonal symmetry-breaking cost: " << cost;

  cost = xtal::symmetry_breaking_strain_cost(F, hexagonal_point_group);
  EXPECT_FALSE(almost_equal(cost, 0.0))
      << "hexagonal symmetry-breaking cost: " << cost;
}

TEST_F(StrainCostTest, TetragonalTest2) {
  // tetragonal stretch
  Eigen::Matrix3d Vx;
  Vx << 1.01, 0.0, 0.0,  //
      0.0, 1.0, 0.0,     //
      0.0, 0.0, 1.0;     //

  Eigen::Matrix3d Vy;
  Vy << 1.0, 0.0, 0.0,  //
      0.0, 1.01, 0.0,   //
      0.0, 0.0, 1.0;    //

  Eigen::Matrix3d Vz;
  Vz << 1.0, 0.0, 0.0,  //
      0.0, 1.0, 0.0,    //
      0.0, 0.0, 1.01;   //

  double cost_x, cost_y, cost_z;

  cost_x = xtal::isotropic_strain_cost(Vx);
  cost_y = xtal::isotropic_strain_cost(Vy);
  cost_z = xtal::isotropic_strain_cost(Vz);
  EXPECT_TRUE(almost_equal(cost_x, cost_y))
      << "isotropic cost_x: " << cost_x << " cost_y: " << cost_y;
  EXPECT_TRUE(almost_equal(cost_y, cost_z))
      << "isotropic cost_y: " << cost_y << " cost_z: " << cost_z;

  cost_x = xtal::symmetry_breaking_strain_cost(Vx, cubic_point_group);
  cost_y = xtal::symmetry_breaking_strain_cost(Vy, cubic_point_group);
  cost_z = xtal::symmetry_breaking_strain_cost(Vz, cubic_point_group);
  EXPECT_TRUE(almost_equal(cost_x, cost_y))
      << "cubic symmetry-breaking cost_x: " << cost_x << " cost_y: " << cost_y;
  EXPECT_TRUE(almost_equal(cost_y, cost_z))
      << "cubic symmetry-breaking cost_y: " << cost_y << " cost_z: " << cost_z;

  cost_x = xtal::symmetry_breaking_strain_cost(Vx, tetragonal_point_group);
  cost_y = xtal::symmetry_breaking_strain_cost(Vy, tetragonal_point_group);
  cost_z = xtal::symmetry_breaking_strain_cost(Vz, tetragonal_point_group);
  EXPECT_FALSE(almost_equal(cost_x, cost_y))
      << "tetragonal symmetry-breaking cost_x: " << cost_x
      << " cost_y: " << cost_y;
  EXPECT_TRUE(almost_equal(cost_y, cost_z))
      << "tetragonal symmetry-breaking cost_y: " << cost_y
      << " cost_z: " << cost_z;
}

TEST_F(StrainCostTest, ShearTest) {
  // shear stretch
  Eigen::Matrix3d F = shear_stretch();
  double cost;

  cost = xtal::isotropic_strain_cost(F);
  EXPECT_FALSE(almost_equal(cost, 0.0)) << "isotropic cost: " << cost;

  cost = xtal::symmetry_breaking_strain_cost(F, cubic_point_group);
  EXPECT_FALSE(almost_equal(cost, 0.0))
      << "cubic symmetry-breaking cost: " << cost;

  cost = xtal::symmetry_breaking_strain_cost(F, tetragonal_point_group);
  EXPECT_FALSE(almost_equal(cost, 0.0))
      << "tetragonal symmetry-breaking cost: " << cost;

  cost = xtal::symmetry_breaking_strain_cost(F, hexagonal_point_group);
  EXPECT_FALSE(almost_equal(cost, 0.0))
      << "hexagonal symmetry-breaking cost: " << cost;
}

TEST_F(StrainCostTest, ShearTest2) {
  // shear stretch
  Eigen::Matrix3d Vxy;
  Vxy << 1.0, 0.01, 0.0,  //
      0.01, 1.0, 0.0,     //
      0.0, 0.0, 1.0;      //

  Eigen::Matrix3d Vxz;
  Vxz << 1.0, 0.0, 0.01,  //
      0.0, 1.0, 0.0,      //
      0.01, 0.0, 1.0;     //

  Eigen::Matrix3d Vyz;
  Vyz << 1.0, 0.0, 0.0,  //
      0.0, 1.0, 0.01,    //
      0.0, 0.01, 1.0;    //

  double cost_xy, cost_xz, cost_yz;

  cost_xy = xtal::isotropic_strain_cost(Vxy);
  cost_xz = xtal::isotropic_strain_cost(Vxz);
  cost_yz = xtal::isotropic_strain_cost(Vyz);
  EXPECT_TRUE(almost_equal(cost_xy, cost_xz))
      << "isotropic cost_xy: " << cost_xy << " cost_xz: " << cost_xz;
  EXPECT_TRUE(almost_equal(cost_xy, cost_yz))
      << "isotropic cost_xy: " << cost_xy << " cost_yz: " << cost_yz;

  cost_xy = xtal::symmetry_breaking_strain_cost(Vxy, cubic_point_group);
  cost_xz = xtal::symmetry_breaking_strain_cost(Vxz, cubic_point_group);
  cost_yz = xtal::symmetry_breaking_strain_cost(Vyz, cubic_point_group);
  EXPECT_TRUE(almost_equal(cost_xy, cost_xz))
      << "cubic symmetry-breaking cost_xy: " << cost_xy
      << " cost_xz: " << cost_xz;
  EXPECT_TRUE(almost_equal(cost_xy, cost_yz))
      << "cubic symmetry-breaking cost_xy: " << cost_xy
      << " cost_yz: " << cost_yz;

  cost_xy = xtal::symmetry_breaking_strain_cost(Vxy, tetragonal_point_group);
  cost_xz = xtal::symmetry_breaking_strain_cost(Vxz, tetragonal_point_group);
  cost_yz = xtal::symmetry_breaking_strain_cost(Vyz, tetragonal_point_group);
  EXPECT_TRUE(almost_equal(cost_xy, cost_xz))
      << "tetragonal symmetry-breaking cost_xy: " << cost_xy
      << " cost_xz: " << cost_xz;
  EXPECT_TRUE(almost_equal(cost_xy, cost_yz))
      << "tetragonal symmetry-breaking cost_xy: " << cost_xy
      << " cost_yz: " << cost_yz;
}

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

  // range of elements of N matrix (controls number of potential mappings to be
  /// considered... larger is more)
  int unimodular_element_range = 1;

  double max_lattice_cost = 1e20;

  bool use_symmetry_breaking_strain_cost = false;

  xtal::LatticeMap lattice_map{L1_lattice,
                               L2_lattice,
                               unimodular_element_range,
                               L1_point_group,
                               L2_point_group,
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
  Eigen::Matrix3d stretch = strain::right_stretch_tensor(F_reverse).inverse();
  Eigen::Matrix3d isometry = (F_reverse * stretch).transpose();

  EXPECT_TRUE(almost_equal(L1 * N, stretch * isometry * L2));

  EXPECT_TRUE(almost_equal(stretch, V_expected))
      << "Unexpected stretch. Found stretch:\n"
      << stretch;
  EXPECT_TRUE(almost_equal(isometry, Q_expected))
      << "Unexpected isometry. Found isometry:\n"
      << isometry;
}
