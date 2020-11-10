#include "gtest/gtest.h"

#include "casm/casm_io/Log.hh"
#include "casm/clusterography/ClusterOrbits_impl.hh"
#include "casm/clusterography/ClusterSpecs_impl.hh"

#include "crystallography/TestStructures.hh" // for test::ZrO_prim

namespace example {

  std::vector<CASM::OrbitBranchSpecs<CASM::PrimPeriodicIntegralClusterOrbit> >
  make_orbit_branch_specs(
    std::shared_ptr<CASM::Structure const> const &shared_prim);

  std::vector<CASM::IntegralClusterOrbitGenerator> make_custom_generators(
    std::shared_ptr<CASM::Structure const> const &shared_prim);

  void check_orbits(
    std::shared_ptr<CASM::Structure const> const &shared_prim,
    std::vector<CASM::PrimPeriodicIntegralClusterOrbit> const &orbits);
}

// While CASM provides a specific function for generating the asymmetric unit, it also provides
//   methods to generate clusters in more general cases using:
//   - the OrbitBranchSpecs class
//   - the IntegralClusterOrbitGenerator class
//   - the `make_orbits` function
//   - classes derived from the ClusterSpecs class
//
// The standard CASM cluster generation methods work recursively, meaning 2-point clusters are
//   generated by attempting to add sites to 1-point clusters, 3-point clusters are generated by
//   attempting to add sites to 2-point clusters, etc.
//
// The "orbit tree" of all cluster orbits is made up of "orbit branches" where each "branch" of the
//   "orbit tree" is made up of all clusters with the same number of sites. For example, the n=2
//   "orbit branch" consists of all orbits of 2-point clusters, the n=3 "orbit branch" consists
//   of all orbits of 3-point clusters, etc.
//
// The OrbitBranchSpecs class specifies how one orbit branch should be constructed, given the
//   previous branch as input. OrbitBranchSpecs has:
//   - A pointer to the prim Structure
//   - A pointer to the generating group
//   - A "SymCompare" functor
//   - A vector of "candidate sites": UnitCellCoord indicating sites that should be added to
//     to clusters from the previous branch in an attempt to find new unique clusters
//   - A filter function (std::function<bool (ClusterType)>) which returns false for clusters that
//     should not be used to construct an Orbit. This can be used to specify truncation criteria
//     for orbits in a branch. The most common criteria is a site-to-site distance cutoff.
//
// Sometimes it is known or suspected that a particular cluster orbit is necessary to include in the
//   cluster expansion basis, perhaps due to a known ground state ordering. The
//   IntegralClusterOrbitGenerator class provides a way to specify that particular cluster orbits
//   (and orbits of subclusters) be generated. It consists of:
//   - an IntegralCluster prototype
//   - a bool "include_subclusters" option, which if true tells the orbit generating method to also
//     generate all orbits that are subclusters of the prototype cluster. For example, if prototype
//     is a triplet cluster, each pair of sites in the triplet would be used to generate pair
//     cluster orbits, and each single site in the triplet would be used to generate single point
//     cluster orbits.
//
// The classes derived from the ClusterSpecs class specify how to construct all the orbit branches.
//   They provide:
//   - A method to generate all orbits, as determined by implementation-specific parameters
//   - A method to generate orbits given a vector of generating elements
//
// Depending on the particular method being implemented, they will typically do this by specifying
//   how many orbit branches to construct and how to generate an OrbitBranchSpecs instance for each
//   branch by providing:
//   - A pointer to the prim Structure
//   - A pointer to the generating group
//   - A "SymCompare" functor
//   - A method to generate "candidate sites" for each orbit branch
//   - A method to generate a cluster filter function for each orbit branch
//   - A method to specify custom clusters, as with a vector of IntegralClusterOrbitGenerator.
//
// CASM provides some functions that can help build the "Specs" classes:
// - SiteFilterFunction (typedef std::function<bool (xtal::Site)>):
//   - all_sites_filter: Generate clusters using all Sites
//   - alloy_sites_filter: Generate clusters using Site with site_occupant.size() > 1
//   - dof_sites_filter(const std::vector<DoFKey> &dofs = {}): Generate clusters using Site with
//     specified DoF
// - ClusterFilterFunction (typedef std::function<bool (IntegralCluster)>):
//   - all_clusters_filter(): Accept all clusters
//   - max_length_cluster_filter(double max_length): Accept clusters with max pair distance less
//     than max_length
//   - within_scel_max_length_cluster_filter(double max_length, Eigen::Matrix3l const &
//     superlattice_matrix): Accept clusters with max pair distance (using closest images) less than
//     max_length
// - CandidateSitesFunction:
//   (typedef std::function<std::vector<xtal::UnitCellCoord> (Structure const &, SiteFilterFunction)>)
//   - empty_neighborhood(): No sites (for null orbit, or global dof only)
//   - origin_neighborhood(): Only sites in the origin unit cell
//   - scel_neighborhood(Eigen::Matrix3l const &superlattice_matrix): All Sites in the supercell
//     defined by the superlattice_matrix
//   - max_length_neighborhood(double max_length): Sites within max_length distance to any site in
//     the origin unit cell.
//   - cutoff_radius_neighborhood(IntegralCluster const &phenomenal, double cutoff_radius): Sites
//     within cutoff_radius distance to any site in the phenomenal cluster
//   - within_scel_cutoff_radius_neighborhood(IntegralCluster const &phenomenal,
//     double cutoff_radius, Eigen::Matrix3l const &superlattice_matrix): Sites within cutoff_radius
//     distance (using closest images) to any site in the phenomenal cluster

// Note:
// To print orbit results generated by this example, use this:
// #include "casm/app/AppIO_impl.hh"
//
// void print_orbits(std::vector<PrimPeriodicIntegralClusterOrbit> const &orbits) {
//   OrbitPrinterOptions printer_options;
//   printer_options.coord_type = CASM::INTEGRAL;
//   ProtoSitesPrinter printer {printer_options};
//   print_clust(orbits.begin(), orbits.end(), CASM::log(), printer);
// }

TEST(ExampleClusterographyClusterSpecs, OrbitBranchSpecs) {

  // An example using OrbitBranchSpecs, IntegralClusterOrbitGenerator for custom cluster orbits,
  //   and the `make_orbits` function directly to generate the n=0, n=1, and n=2 branches, along
  //   with custom n=3 orbits.

  // Construct a ZrO prim
  auto shared_prim = std::make_shared<CASM::Structure const>(test::ZrO_prim());

  // PrimPeriodicIntegralClusterOrbit is a typedef for:
  // - Orbit<PrimPeriodicSymCompare<IntegralCluster>>
  typedef CASM::PrimPeriodicIntegralClusterOrbit orbit_type;

  // Local function to construct OrbitBranchSpecs for:
  // - null orbit branch (cluster expansion constant term)
  // - point cluster orbit branch
  // - pair clusters up to max length 5.17=c-axis length + a little
  std::vector<CASM::OrbitBranchSpecs<orbit_type>> orbit_branch_specs;
  orbit_branch_specs = example::make_orbit_branch_specs(shared_prim);

  // --- Make orbits using OrbitBranchSpecs and custom generators

  // Store generated orbits in this container
  std::vector<orbit_type> orbits;

  // Print messages about orbit generation progress to `status` stream, in this case do not print.
  std::ostream &status = CASM::null_log();

  // Construct orbits
  make_orbits(
    orbit_branch_specs.begin(), // range of OrbitBranchSpecs
    orbit_branch_specs.end(),
    example::make_custom_generators(shared_prim), // custom triplet clusters (and subclusters)
    std::back_inserter(orbits), // orbit vector inserter
    status
  );

  example::check_orbits(shared_prim, orbits);
}

TEST(ExampleClusterographyClusterSpecs, PeriodicMaxLengthClusterSpecs) {

  // An example using PeriodicMaxLengthClusterSpecs to generate the same cluster orbits as the
  //   previous example.

  // Construct a ZrO prim (same as above)
  auto shared_prim = std::make_shared<CASM::Structure const>(test::ZrO_prim());

  // Use the complete factor group as the generating group
  CASM::SymGroup const &generating_group = shared_prim->factor_group();

  // Cluster "max_length" cutoff values: max_length[b] is cutoff value for orbit branch b
  // - null and point cluster values are ignored
  // - pair clusters up to max length 5.17=c-axis length + a little
  std::vector<double> max_length = {0, 0, 5.17};

  // PeriodicMaxLengthClusterSpecs is constructed with:
  CASM::PeriodicMaxLengthClusterSpecs cluster_specs {
    shared_prim,                        // the prim Structure
    generating_group,                   // the orbit generating group
    CASM::alloy_sites_filter,           // filter to include sites with >1 occupant in cluster orbits
    max_length,                         // cluster "max_length" cutoff values
    example::make_custom_generators(shared_prim) // custom triplet clusters (and subclusters)
  };

  // Print messages about orbit generation progress to `status` stream, in this case do not print.
  std::ostream &status = CASM::null_log();

  // ClusterSpecs::PeriodicOrbitVec is typedef of
  //   std::vector<Orbit<PrimPeriodicSymCompare<IntegralCluster>>>
  CASM::ClusterSpecs::PeriodicOrbitVec orbits = cluster_specs.make_periodic_orbits(status);

  example::check_orbits(shared_prim, orbits);
}


namespace example {

  std::vector<CASM::OrbitBranchSpecs<CASM::PrimPeriodicIntegralClusterOrbit> >
  make_orbit_branch_specs(std::shared_ptr<CASM::Structure const> const &shared_prim) {

    typedef CASM::PrimPeriodicSymCompare<CASM::IntegralCluster> sym_compare_type;
    typedef CASM::Orbit<sym_compare_type> orbit_type;

    // Use the complete factor group as the generating group
    CASM::SymGroup const &generating_group = shared_prim->factor_group();

    // Use a PrimPeriodicSymCompare<IntegralCluster> "SymCompare" functor:
    sym_compare_type sym_compare {shared_prim, shared_prim->lattice().tol()};

    // A SiteFilterFunction returns true if a Site should be included and false if it should be
    // excluded. It is a typedef for:
    // - std::function<bool (xtal::Site)>
    // In this case, alloy_sites_filter includes only sites with > 1 allowed occupant
    CASM::SiteFilterFunction site_filter = CASM::alloy_sites_filter;

    // A vector that will be filled with candidate sites for constructing clusters
    std::vector<CASM::xtal::UnitCellCoord> candidate_sites;

    // A CandidateSitesFunction generates a vector of UnitCellCoord from a Structure and
    // SiteFilterFuntion. It is typedef for:
    // - std::function<std::vector<xtal::UnitCellCoord> (Structure const &, SiteFilterFunction)>
    CASM::CandidateSitesFunction candidate_sites_f;

    // A ClusterFilterFunction returns true if an IntegralCluster should be included and false if it
    // should be excluded. It is a typedef for:
    // - std::function<std::vector<xtal::UnitCellCoord> (Structure const &, SiteFilterFunction)>
    CASM::ClusterFilterFunction cluster_filter;

    // A vector of OrbitBranchSpecs for our orbit_type. Include one for each orbit branch that should
    //   be generated, including the null orbit branch.
    std::vector<CASM::OrbitBranchSpecs<orbit_type> > orbit_branch_specs;

    // --- Construct n=0 OrbitBranchSpecs ---
    // For the null cluster orbit, we include no candidate sites. The
    //   cluster_filter value should always return true so we use `all_clusters_filter`.

    // Using `empty_neighborhood` is not strictly necessary, the empty `candidate_sites` vector would
    //   suffice. It is used here to illustrate its use.
    candidate_sites_f = CASM::empty_neighborhood();
    candidate_sites = candidate_sites_f(*shared_prim, site_filter);
    cluster_filter = CASM::all_clusters_filter();
    orbit_branch_specs.emplace_back(*shared_prim, candidate_sites.begin(), candidate_sites.end(),
                                    generating_group, cluster_filter, sym_compare);

    // --- Construct n=1 OrbitBranchSpecs ---
    // For the single point cluster orbit branch, include all origin unit cell sites in candidate
    //   sites. The cluster_filter value should always return true so we use `all_clusters_filter`.
    candidate_sites_f = CASM::origin_neighborhood();
    candidate_sites = candidate_sites_f(*shared_prim, site_filter);
    cluster_filter = CASM::all_clusters_filter();
    orbit_branch_specs.emplace_back(*shared_prim, candidate_sites.begin(), candidate_sites.end(),
                                    generating_group, cluster_filter, sym_compare);

    // --- Construct n=2 OrbitBranchSpecs, using pairs with max_length<5.17 ---
    // For the pair cluster orbit branch, include all origin unit cell sites in candidate
    //   sites and sites in neighboring unit cells as necessary to be sure pair clusters up to length
    //   `max_length` are included. The cluster_filter will return true only for clusters whose pair
    //   distance is less `max_length`.
    double max_length = 5.17;
    candidate_sites_f = CASM::max_length_neighborhood(max_length);
    candidate_sites = candidate_sites_f(*shared_prim, site_filter);
    cluster_filter = CASM::max_length_cluster_filter(max_length);
    orbit_branch_specs.emplace_back(*shared_prim, candidate_sites.begin(), candidate_sites.end(),
                                    generating_group, cluster_filter, sym_compare);

    return orbit_branch_specs;
  }

  std::vector<CASM::IntegralClusterOrbitGenerator> make_custom_generators(
    std::shared_ptr<CASM::Structure const> const &shared_prim) {

    // --- Custom cluster orbit generators for triplet clusters ---
    // For the custom triplet cluster orbits we directly specify the cluster prototype sites and set
    //   the `include_subclusters` flag to `true`. In this example all the subclusters are included
    //   already for the first two custom orbits.

    CASM::IntegralCluster triplet_cluster {*shared_prim};
    bool include_subclusters = true;
    std::vector<CASM::IntegralClusterOrbitGenerator> custom_generators;

    // triplet cluster made of [a, 0, 0], [a, a, 0] (vectors from reference site {2,0,0,0})
    triplet_cluster.elements().clear();
    triplet_cluster.elements().emplace_back(2, 0, 0, 0);
    triplet_cluster.elements().emplace_back(2, 1, 0, 0);
    triplet_cluster.elements().emplace_back(2, 1, 1, 0);
    custom_generators.emplace_back(triplet_cluster, include_subclusters);

    // triplet cluster made of [a, 0, 0], [0, 0, c/2]
    triplet_cluster.elements().clear();
    triplet_cluster.elements().emplace_back(2, 0, 0, 0);
    triplet_cluster.elements().emplace_back(3, 0, 0, 0);
    triplet_cluster.elements().emplace_back(2, 1, 0, 0);
    custom_generators.emplace_back(triplet_cluster, include_subclusters);

    // triplet cluster made of [0, 0, c], [0, 0, 2*c]
    triplet_cluster.elements().clear();
    triplet_cluster.elements().emplace_back(2, 0, 0, 0);
    triplet_cluster.elements().emplace_back(2, 0, 0, 1);
    triplet_cluster.elements().emplace_back(2, 0, 0, 2);
    custom_generators.emplace_back(triplet_cluster, include_subclusters);

    return custom_generators;
  }

  // Check the generated PrimPeriodicIntegralClusterOrbit for this example:
  // - ZrO prim (a=3.23398686, c=5.16867834)
  // - alloy_sites_filter
  // - n=2 orbit branch: max_length < 5.17 (= c + eps)
  // - three custom n=3 orbits plus subclusters (from make_custom_generators)
  void check_orbits(
    std::shared_ptr<CASM::Structure const> const &shared_prim,
    std::vector<CASM::PrimPeriodicIntegralClusterOrbit> const &orbits) {

    // Lattice vectors, for checking orbit results
    Eigen::Vector3d a, b, c;
    std::tie(a, b, c) = shared_prim->lattice().vectors();

    // Expected orbits:
    EXPECT_EQ(orbits.size(), 10);

    EXPECT_EQ(orbits[0].size(), 1); // null orbit
    EXPECT_EQ(orbits[0].prototype().size(), 0);

    EXPECT_EQ(orbits[1].size(), 2); // [Va,O] site orbit
    EXPECT_EQ(orbits[1].prototype().size(), 1);

    EXPECT_EQ(orbits[2].size(), 2); // 2NN pair [0, 0, c/2]
    EXPECT_EQ(orbits[2].prototype().size(), 2);
    EXPECT_TRUE(CASM::almost_equal(orbits[2].invariants().displacement().back(), (c / 2.).norm()));

    EXPECT_EQ(orbits[3].size(), 6); // 1NN pair [a, 0, 0]
    EXPECT_EQ(orbits[3].prototype().size(), 2);
    EXPECT_TRUE(CASM::almost_equal(orbits[3].invariants().displacement().back(), a.norm()));

    EXPECT_EQ(orbits[4].size(), 12); // 3NN pair [a, a, c/2]
    EXPECT_EQ(orbits[4].prototype().size(), 2);
    EXPECT_TRUE(CASM::almost_equal(orbits[4].invariants().displacement().back(), (a + b + c / 2.).norm()));

    EXPECT_EQ(orbits[5].size(), 2); // 4NN pair [0, 0, c]
    EXPECT_EQ(orbits[5].prototype().size(), 2);
    EXPECT_TRUE(CASM::almost_equal(orbits[5].invariants().displacement().back(), c.norm()));

    EXPECT_EQ(orbits[6].size(), 2); // 4NN pair [0, 0, 2*c] (include as subcluster of orbits[9])
    EXPECT_EQ(orbits[6].prototype().size(), 2);
    EXPECT_TRUE(CASM::almost_equal(orbits[6].invariants().displacement().back(), (2.*c).norm()));

    EXPECT_EQ(orbits[7].size(), 4); // triplet cluster [a, 0, 0], [a, a, 0]
    EXPECT_EQ(orbits[7].prototype().size(), 3);
    EXPECT_TRUE(CASM::almost_equal(orbits[7].invariants().displacement()[0], a.norm()));
    EXPECT_TRUE(CASM::almost_equal(orbits[7].invariants().displacement()[1], a.norm()));
    EXPECT_TRUE(CASM::almost_equal(orbits[7].invariants().displacement()[2], a.norm()));

    EXPECT_EQ(orbits[8].size(), 24); // triplet cluster [a, 0, 0], [0, 0, c/2]
    EXPECT_EQ(orbits[8].prototype().size(), 3);
    EXPECT_TRUE(CASM::almost_equal(orbits[8].invariants().displacement()[0], (c / 2.).norm()));
    EXPECT_TRUE(CASM::almost_equal(orbits[8].invariants().displacement()[1], a.norm()));
    EXPECT_TRUE(CASM::almost_equal(orbits[8].invariants().displacement()[2], (a + c / 2.).norm()));

    EXPECT_EQ(orbits[9].size(), 2); // triplet cluster [0, 0, c], [0, 0, 2*c]
    EXPECT_EQ(orbits[9].prototype().size(), 3);
    EXPECT_TRUE(CASM::almost_equal(orbits[9].invariants().displacement()[0], c.norm()));
    EXPECT_TRUE(CASM::almost_equal(orbits[9].invariants().displacement()[1], c.norm()));
    EXPECT_TRUE(CASM::almost_equal(orbits[9].invariants().displacement()[2], (2.*c).norm()));
  }

}
