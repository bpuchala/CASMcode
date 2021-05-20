#include "casm/crystallography/LatticeMap.hh"

#include "autotools.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "gtest/gtest.h"

TEST(LatticeMapTest, Test1) {
  Eigen::Matrix3d V;
  V << 0.998559791924, -0.000000000970, 0.000000000000,  //
      -0.000000000970, 0.998559788571, 0.000000000000,   //
      0.000000000000, 0.000000000000, 0.997675852223;    //

  Eigen::Matrix3d Q;
  Q << -0.500000000000, -0.866025403784, 0.000000000000,  //
      -0.866025403784, 0.500000000000, 0.000000000000,    //
      0.000000000000, 0.000000000000, -1.000000000000;    //

  // parent lat_column_mat
  Eigen::Matrix3d L1;
  L1 << 3.233986860000, -1.616993430000, 0.000000000000,  //
      0.000000000000, 2.800714770000, 0.000000000000,     //
      0.000000000000, 0.000000000000, 5.168678340000;     //

  // child lat_column_mat
  Eigen::Matrix3d L2;
  L2 << 3.238651197049, 0.000000000000, 0.000000000000,  //
      -1.619325598525, 2.804754204349, 0.000000000000,   //
      0.000000000000, 0.000000000000, 5.180719096772;    //

  std::cout << "L1: \n" << std::endl;

  EXPECT_TRUE(almost_equal(L1, V * Q * L2));
}
