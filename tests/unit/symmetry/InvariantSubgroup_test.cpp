#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

/// What is being tested:
#include "casm/symmetry/InvariantSubgroup.hh"
#include "casm/symmetry/InvariantSubgroup_impl.hh"

/// What is being used to test it:
#include "Common.hh"
#include "casm/app/AppIO_impl.hh"
#include "casm/casm_io/jsonFile.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/symmetry/Orbit_impl.hh"

using namespace CASM;

BOOST_AUTO_TEST_SUITE(InvariantSubgroupTest)

BOOST_AUTO_TEST_CASE(Test0) {
  test::ZrOProj proj;
  proj.check_init();

  Logging logging = Logging::null();
  PrimClex primclex(proj.dir, logging);
  const Structure &prim = primclex.prim();
  const SymGroup &prim_fg = primclex.prim().factor_group();
  const Lattice &lat = prim.lattice();
  Supercell prim_scel(&primclex, Eigen::Matrix3i::Identity());

  Log &log = primclex.log();

  BOOST_CHECK_EQUAL(true, true);

  // Make PrimPeriodicIntegralClusterOrbit
  jsonFile bspecs {"tests/unit/kinetics/ZrO_bspecs_0.json"};

  std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
  make_prim_periodic_orbits(
    primclex.prim(),
    bspecs,
    alloy_sites_filter,
    primclex.crystallography_tol(),
    std::back_inserter(orbits),
    primclex.log());
  BOOST_CHECK_EQUAL(true, true);

  PrototypePrinter<IntegralCluster> printer;
  printer.opt.coord_type = INTEGRAL;
  // print_clust(orbits.begin(), orbits.end(), primclex.log(), printer);

  // Make cluster groups & check size, based on prim.factor_group symmetry
  {
    for(const auto &orbit : orbits) {

      // Test make_invariant_subgroup using orbit generators
      SymGroup cluster_group_a = make_invariant_subgroup(orbit.prototype(), prim_fg, orbit.sym_compare());

      // Test make_invariant_subgroup using the orbit equivalence map & orbit.prototype()
      SymGroup cluster_group_b = make_invariant_subgroup(orbit);
      BOOST_CHECK_EQUAL(cluster_group_a.size(), cluster_group_b.size());

      // Test make_invariant_subgroup using Supercell
      std::vector<PermuteIterator> cluster_group_c = make_invariant_subgroup(orbit.prototype(), prim_scel);
      BOOST_CHECK_EQUAL(cluster_group_a.size(), cluster_group_c.size());

    }
    //test::print_computed_result(std::cout, "cluster_group_size", cluster_group_size);
    BOOST_CHECK_EQUAL(true, true);
  }

  // Individual cluster test cases
  {
    // Make vol 2 supercell & and background configuration
    Eigen::Vector3d a, b, c;
    std::tie(a, b, c) = primclex.prim().lattice().vectors();
    Supercell scel_vol2 {&primclex, Lattice(2 * a, 1 * b, 1 * c)};
    Configuration config(scel_vol2);
    config.init_occupation();

    // Scel sym_compare
    ScelPeriodicSymCompare<IntegralCluster> scel_sym_compare(scel_vol2);

    // Get the config factor group (should just be all Supercell operations)
    std::vector<PermuteIterator> _config_fg = config.factor_group();
    SymGroup config_fg = make_sym_group(_config_fg.begin(), _config_fg.end());

    BOOST_CHECK_EQUAL(scel_vol2.factor_group().size(), 8);
    BOOST_CHECK_EQUAL(config_fg.size(), 16);

    // a couple test cases:

    // null cluster
    {
      IntegralCluster clust {prim};
      SymGroup cluster_group = make_invariant_subgroup(clust, config_fg, scel_sym_compare);
      BOOST_CHECK_EQUAL(cluster_group.size(), 16);
    }

    // point cluster
    {
      IntegralCluster clust {prim};
      clust.elements().push_back(UnitCellCoord {prim, 2, 0, 0, 0});
      SymGroup cluster_group = make_invariant_subgroup(clust, config_fg, scel_sym_compare);
      BOOST_CHECK_EQUAL(cluster_group.size(), 4);
    }

    // pair cluster
    {
      IntegralCluster clust {prim};
      clust.elements().push_back(UnitCellCoord {prim, 3, 1, 0, 0});
      clust.elements().push_back(UnitCellCoord {prim, 2, 1, 1, 1});
      SymGroup cluster_group = make_invariant_subgroup(clust, config_fg, scel_sym_compare);
      BOOST_CHECK_EQUAL(cluster_group.size(), 2);
    }

    // pair cluster - equivalent to previous cluster by prim symmetry, different by scel symmetry
    {
      IntegralCluster clust {prim};
      clust.elements().push_back(UnitCellCoord {prim, 3, 1, 0, 0});
      clust.elements().push_back(UnitCellCoord {prim, 2, 2, 1, 1});
      SymGroup cluster_group = make_invariant_subgroup(clust, config_fg, scel_sym_compare);
      BOOST_CHECK_EQUAL(cluster_group.size(), 1);
    }
  }

  // Using a reduced symmetry vol 2 Supercell, check three different methods of generating the invariant subgroup
  {
    // Make vol 2 supercell & and background configuration
    Eigen::Vector3d a, b, c;
    std::tie(a, b, c) = primclex.prim().lattice().vectors();
    Supercell scel_vol2 {&primclex, Lattice(2 * a, 1 * b, 1 * c)};
    Configuration config(scel_vol2);
    config.init_occupation();

    // Scel sym_compare
    ScelPeriodicSymCompare<IntegralCluster> scel_sym_compare(scel_vol2);

    // Get the config factor group (should just be all Supercell operations)
    std::vector<PermuteIterator> _config_fg = config.factor_group();
    SymGroup config_fg = make_sym_group(_config_fg.begin(), _config_fg.end());

    BOOST_CHECK_EQUAL(scel_vol2.factor_group().size() * 2, _config_fg.size());
    BOOST_CHECK_EQUAL(config_fg.size(), _config_fg.size());

    for(const auto &orbit : orbits) {

      OrbitGenerators<ScelPeriodicIntegralClusterOrbit> generators {config_fg, scel_sym_compare};
      for(const auto &eq : orbit) {
        for(auto it = scel_vol2.translate_begin(); it != scel_vol2.translate_end(); ++it) {
          generators.insert(copy_apply(it, eq));
        }
      }

      Index el_sum = 0;
      for(const auto &el : generators.elements) {

        ScelPeriodicIntegralClusterOrbit suborbit {el, config_fg, scel_sym_compare};
        el_sum += suborbit.size();

        // Test make_invariant_subgroup using orbit generators
        SymGroup cluster_group_a = make_invariant_subgroup(el, config_fg, scel_sym_compare);

        // group size must equal number of elements * invariant group size
        BOOST_CHECK_EQUAL(config_fg.size(), suborbit.size()*cluster_group_a.size());

        // Test make_invariant_subgroup using the orbit equivalence map & orbit.prototype()
        SymGroup cluster_group_b = make_invariant_subgroup(suborbit);
        BOOST_CHECK_EQUAL(cluster_group_a.size(), cluster_group_b.size());

        // Test make_invariant_subgroup using Supercell
        std::vector<PermuteIterator> cluster_group_c = make_invariant_subgroup(el, scel_vol2);
        BOOST_CHECK_EQUAL(cluster_group_a.size(), cluster_group_c.size());
      }

      if(!orbit.prototype().size()) {
        BOOST_CHECK_EQUAL(orbit.size(), 1);
        BOOST_CHECK_EQUAL(el_sum, 1);
      }
      else {
        BOOST_CHECK_EQUAL(2 * orbit.size(), el_sum);
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
