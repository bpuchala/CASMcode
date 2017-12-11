#include "casm/crystallography/Site.hh"
#include "casm/clusterography/ClusterOrbits.hh"
#include "casm/clex/PrimClex.hh"
#include "casm/clex/HasPrimClex_impl.hh"
#include "casm/app/AppIO.hh"
#include "casm/database/DiffTransOrbitDatabase.hh"

#include "casm/kinetics/DiffusionTransformationEnum_impl.hh"
#include "casm/clusterography/ClusterSymCompare_impl.hh"
#include "casm/app/AppIO_impl.hh"
#include "casm/casm_io/InputParser_impl.hh"
#include "casm/clusterography/ClusterSpecsParser_impl.hh"

extern "C" {
  CASM::EnumInterfaceBase *make_DiffusionTransformationEnum_interface() {
    return new CASM::EnumInterface<CASM::Kinetics::DiffusionTransformationEnum>();
  }
}

namespace CASM {

  namespace Kinetics {

    namespace {

      Log &operator<<(Log &out, const IntegralCluster &clust) {
        SitesPrinter printer;
        printer.print(clust, out);
        out << std::endl;
        return out;
      }
    }

    /// \brief Construct with an IntegralCluster
    DiffusionTransformationEnum::DiffusionTransformationEnum(const IntegralCluster &clust) :
      m_cluster(clust),
      m_current(new DiffusionTransformation(clust.prim())) {

      // initialize to/from counter
      _init_occ_counter();

      // initialize the specie trajectory for the first valid set of to/from occupation values
      m_from_loc = _init_from_loc(m_occ_counter());
      m_to_loc = _init_to_loc(m_occ_counter());

      // set initial DiffTrans
      _set_current();

      if(!m_current->is_valid()) {
        increment();
      }

      if(!m_occ_counter.valid()) {
        _invalidate();
      }
      else {
        this->_initialize(&(*m_current));
        _set_step(0);
      }
    }

    const std::string DiffusionTransformationEnum::enumerator_name = "DiffusionTransformationEnum";
    const std::string DiffusionTransformationEnum::interface_help =
      "DiffusionTransformationEnum: \n\n"

      "  cspecs: JSON object \n"
      "    Indicate clusters to enumerate all occupational diffusion transformations. The \n"
      "    JSON item \"cspecs\" should be a cspecs style initialization of cluster number and sizes.\n"
      "    See below.          \n\n"

      "  require: JSON array of strings (optional,default=[]) \n "
      "    Indicate required species (atom or molecules names) to enforce that a given species \n"
      "    must be a part of the diffusion transformation. The JSON array \"require\" should be \n"
      "    an array of species names. i.e. \"require\": [\"Va\",\"O\"] \n\n"

      "  exclude: JSON array of strings (optional,default=[]) \n "
      "    Indicate excluded species (atom or molecules names) to enforce that a given species \n"
      "    must not be a part of the diffusion transformation. The JSON array \"exclude\" should \n"
      "    be an array of species names. i.e. \"exclude\": [\"Al\",\"Ti\"] \n\n"

      "  dry_run: bool (optional, default=false)\n"
      "    Perform dry run.\n\n"

      "  coordinate_mode: string (optional, default=FRAC)\n"
      "    Coordinate mode (FRAC, CART, INTEGRAL) for printing orbits.\n\n"

      "  orbit_print_mode: string (optional, default=\"PROTO\")\n"
      "    Mode (FULL, PROTO) to select printing full orbits or just orbit prototypes.\n\n"

      "  Example:\n"
      "  {\n"
      "   \"require\":[\"Va\"],\n"
      "   \"exclude\":[],\n"
      "   \"cspecs\":{\n"
      "      \"orbit_branch_specs\" : { \n"
      "       \"2\" : {\"max_length\" : 5.01},\n"
      "       \"3\" : {\"max_length\" : 5.01}\n"
      "      }\n"
      "    }\n"
      "  }\n\n"
      ;

    template<typename T, typename...Args>
    T get_else(const jsonParser &json, const std::string &key, const T &default_value, Args &&... args) {
      auto it = json.find(key);
      if(it != json.end()) {
        return it->get<T>(std::forward<Args>(args)...);
      }
      else {
        return default_value;
      }
    }


    class DiffTransEnumParser : InputParser, HasPrimClex<CRTPBase<DiffTransEnumParser>> {

    public:

      DiffTransEnumParser(
        const PrimClex &_primclex,
        jsonParser &_input,
        fs::path _path,
        bool _required) :
        InputParser(_input, _path, _required),
        m_primclex(_primclex) {

        auto _relpath = Relpath(_path);

        this->kwargs["cspecs"] =
          std::make_shared<PrimPeriodicClustersByMaxLength>(
            input, _relpath("cspecs"), true);

        // check that "require" and "exclude" species names are valid
        {

        }
      }

      //      void require_valid_species(std::string option) {
      //        auto ptr = this->optional<std::set<std::string>>("require");
      //        if(ptr) {
      //          auto struc_specie = prim().struc_specie();
      //          auto struc_mol = prim().struc_molecule();
      //
      //          auto res = std::all_of(ptr->begin(), ptr->end(), [&](std::string val) {});
      //          require_valid_species("require");
      //        }
      //      }

      std::set<std::string> required_species() const {
        return self.get_else<std::set<std::string>>("require", std::set<std::string>());
      }

      std::set<std::string> excluded_species() const {
        return self.get_else<std::set<std::string>>("exclude", std::set<std::string>());
      }

      bool dry_run() const {
        return self.get_else<bool>("dry_run", false);
      }

      COORD_TYPE coordinate_mode() const {
        return self.get_else<COORD_TYPE>(traits<COORD_TYPE>::name, COORD_TYPE::FRAC);
      }

      ORBIT_PRINT_MODE orbit_print_mode() const {
        return self.get_else<ORBIT_PRINT_MODE>(traits<ORBIT_PRINT_MODE>::name, ORBIT_PRINT_MODE::PROTO);
      }

      const PrimPeriodicClustersByMaxLength &cspecs() const {
        return *m_cspecs_parser;
      }

      const PrimClex &primclex() const {
        return m_primclex;
      }

    private:
      const PrimClex &m_primclex;
      std::shared_ptr<PrimPeriodicClustersByMaxLength> m_cspecs_parser;

    };

    /// Implements increment
    void DiffusionTransformationEnum::increment() {

      // get the next valid DiffTrans
      // to do so, get the next valid specie trajectory
      do {

        // by taking a permutation of possible 'to' specie position
        bool valid_perm = std::next_permutation(m_to_loc.begin(), m_to_loc.end());

        // if no more possible specie trajectory,
        if(!valid_perm) {

          // get next valid from/to occupation values
          do {
            m_occ_counter++;
            _update_current_occ_transform();
          }
          while(m_occ_counter.valid() && !m_current->is_valid_occ_transform());

          // if no more possible from/to occupation values, return
          if(!m_occ_counter.valid()) {
            _invalidate();
            return;
          }
          else {
            m_from_loc = _init_from_loc(m_occ_counter());
            m_to_loc = _init_to_loc(m_occ_counter());
            _update_current_occ_transform();
            _set_current_loc();
          }
        }
        _update_current_to_loc();
      }
      while(!m_current->is_valid_specie_traj());

      _increment_step();
    }

    /// Implements run
    template<typename DatabaseType>
    int DiffusionTransformationEnum::run(
      const PrimClex &primclex,
      const jsonParser &_kwargs,
      const Completer::EnumOption &enum_opt,
      DatabaseType &db) {

      bool dry_run = CASM::dry_run(_kwargs, enum_opt);
      std::string dry_run_msg = CASM::dry_run_msg(dry_run);
      std::string lead;

      jsonParser kwargs;
      if(!_kwargs.contains("cspecs")) {
        primclex.err_log() << "DiffusionTransformationEnum currently has no default and requires a correct JSON with a bspecs tag within it" << std::endl;
        throw std::runtime_error("Error in DiffusionTransformationEnum: cspecs not found");
      }

      std::vector<std::string> require;
      std::vector<std::string> exclude;
      if(_kwargs.contains("require")) {
        for(auto it = _kwargs["require"].begin(); it != _kwargs["require"].end(); ++it) {
          require.push_back(from_json<std::string>(*it));
        }
      }
      if(_kwargs.contains("exclude")) {
        for(auto it = _kwargs["exclude"].begin(); it != _kwargs["exclude"].end(); ++it) {
          exclude.push_back(from_json<std::string>(*it));
        }
      }
      std::vector<std::string> filter_expr = make_enumerator_filter_expr(_kwargs, enum_opt);

      COORD_TYPE coord_mode;
      _kwargs.get_else(coord_mode, traits<COORD_TYPE>::name, COORD_TYPE::FRAC);
      ORBIT_PRINT_MODE orbit_print_mode;
      _kwargs.get_else(orbit_print_mode, traits<ORBIT_PRINT_MODE>::name, ORBIT_PRINT_MODE::PROTO);
      Printer<Kinetics::DiffusionTransformation> dt_printer(6, '\n', coord_mode);

      Index Ninit = db.size();
      Log &log = primclex.log();
      lead = dry_run_msg;
      log << lead << "# diffusion transformations in this project: " << Ninit << "\n" << std::endl;

      log.begin(enumerator_name);
      log.increase_indent();
      log << std::endl;

      log.begin<Log::verbose>("Calculate cluster orbits");
      std::vector<PrimPeriodicIntegralClusterOrbit> orbits;
      auto end = make_prim_periodic_orbits(
                   primclex.prim(),
                   _kwargs["cspecs"],
                   alloy_sites_filter,
                   primclex.crystallography_tol(),
                   std::back_inserter(orbits),
                   primclex.log());
      print_clust(orbits.begin(), orbits.end(), log, orbit_print_mode, coord_mode);
      log << std::endl;

      log.begin<Log::verbose>("Calculate diff_trans orbits");
      std::vector< PrimPeriodicDiffTransOrbit > diff_trans_orbits;
      auto end2 = make_prim_periodic_diff_trans_orbits(
                    orbits.begin(),
                    orbits.end(),
                    primclex.crystallography_tol(),
                    std::back_inserter(diff_trans_orbits),
                    &primclex);
      print_clust(diff_trans_orbits.begin(), diff_trans_orbits.end(), log, orbit_print_mode, coord_mode);
      log << std::endl;

      log.begin<Log::verbose>("Check diff_trans orbits");
      log << lead << "COORD_MODE = " << coord_mode << std::endl << std::endl;
      lead = log.indent_str() + dry_run_msg;
      for(auto &diff_trans_orbit : diff_trans_orbits) {
        log << lead << "Checking: \n";
        log.increase_indent();
        dt_printer.print(diff_trans_orbit.prototype(), log);
        auto specie_count = diff_trans_orbit.prototype().specie_count();

        if(!includes_all(specie_count, require.begin(), require.end())) {
          log << lead << "- Missing required species, do not insert" << std::endl;
          log.decrease_indent();
          log << std::endl;
          continue;
        }
        if(!excludes_all(specie_count, exclude.begin(), exclude.end())) {
          log << lead << "- Includes excluded species, do not insert" << std::endl;
          log.decrease_indent();
          log << std::endl;
          continue;
        }

        //insert current into database
        auto res = db.insert(diff_trans_orbit);
        if(res.second) {
          log << lead << "- Inserted as: " << res.first->name() << std::endl;
        }
        else {
          log << lead << "- Already exists: " << res.first->name() << std::endl;
        }
        log << std::endl;
        log.decrease_indent();

      }

      log << lead << "  DONE." << std::endl << std::endl;
      log.decrease_indent();
      lead = log.indent_str() + dry_run_msg;

      Index Nfinal = db.size();

      log << lead << "# new diffusion transformations: " << Nfinal - Ninit << "\n";
      log << lead << "# diffusion transformations in this project: " << Nfinal << "\n" << std::endl;

      return 0;
    }

    /// Implements run
    int DiffusionTransformationEnum::run(
      const PrimClex &primclex,
      const jsonParser &_kwargs,
      const Completer::EnumOption &enum_opt) {

      auto &db = primclex.db<PrimPeriodicDiffTransOrbit>();
      bool dry_run = CASM::dry_run(_kwargs, enum_opt);

      DiffusionTransformationEnum::run(primclex, _kwargs, enum_opt, db);

      if(!dry_run) {
        primclex.log() << "Writing diffusion transformation database..." << std::endl;
        db.commit();
        primclex.log() << "  DONE" << std::endl;
      }

      return 0;
    }

    template int DiffusionTransformationEnum::run<std::set<PrimPeriodicDiffTransOrbit>>(
      const PrimClex &,
      const jsonParser &,
      const Completer::EnumOption &,
      std::set<PrimPeriodicDiffTransOrbit> &);

    // -- Unique -------------------

    const Structure &DiffusionTransformationEnum::prim() const {
      return m_cluster.prim();
    }

    const IntegralCluster &DiffusionTransformationEnum::cluster() const {
      return m_cluster;
    }

    /// \brief The occ_counter contains the from/to occupation values for each site
    ///
    /// - Layout is: [from values | to values ]
    ///
    void DiffusionTransformationEnum::_init_occ_counter() {
      Index N = cluster().size();
      std::vector<Index> init_occ(N * 2, 0);
      std::vector<Index> final_occ(N * 2);
      for(int i = 0; i < N; i++) {
        final_occ[i] = cluster()[i].site().site_occupant().size() - 1;
        final_occ[i + N] = final_occ[i];
      }
      std::vector<Index> incr(N * 2, 1);

      m_occ_counter = Counter<std::vector<Index> >(init_occ, final_occ, incr);
    }

    /// \brief Returns container of 'from' specie locations
    std::vector<SpecieLocation> DiffusionTransformationEnum::_init_from_loc(const std::vector<Index> &occ_values) {
      return _init_loc(occ_values, 0);
    }

    /// \brief Returns container of 'to' specie locations
    std::vector<SpecieLocation> DiffusionTransformationEnum::_init_to_loc(const std::vector<Index> &occ_values) {
      return _init_loc(occ_values, cluster().size());
    }

    /// \brief Returns container of 'from' or 'to' specie locations
    ///
    /// - offset == 0 for 'from', N for 'to' specie locations
    ///
    std::vector<SpecieLocation> DiffusionTransformationEnum::_init_loc(const std::vector<Index> &occ_values, Index offset) {

      Index N = cluster().size();
      std::vector<SpecieLocation> loc;
      // for each 'from' occupant
      for(Index i = 0; i < N; ++i) {
        Index occ = occ_values[i + offset];
        UnitCellCoord uccoord = cluster()[i];
        Index mol_size = uccoord.site().site_occupant()[occ].size();
        // for each specie
        for(Index j = 0; j < mol_size; ++j) {
          loc.emplace_back(uccoord, occ, j);
        }
      }
      return loc;
    }

    /// \brief Uses m_cluster, m_occ_counter, m_from_loc, and m_to_loc to set m_current
    void DiffusionTransformationEnum::_set_current() {
      m_current->occ_transform().clear();
      m_current->occ_transform().reserve(cluster().size());
      for(const auto &uccoord : cluster()) {
        m_current->occ_transform().emplace_back(uccoord, 0, 0);
      }
      _update_current_occ_transform();
      _set_current_loc();
      _update_current_to_loc();
    }

    void DiffusionTransformationEnum::_update_current_occ_transform() {
      auto N = cluster().size();
      Index i = 0;
      for(auto &t : m_current->occ_transform()) {
        t.from_value = m_occ_counter()[i];
        t.to_value = m_occ_counter()[i + N];
        ++i;
      }
    }

    void DiffusionTransformationEnum::_set_current_loc() {
      m_current->specie_traj().clear();
      m_current->specie_traj().reserve(m_from_loc.size());
      for(const auto &t : m_from_loc) {
        m_current->specie_traj().emplace_back(t, t);
      }
    }

    void DiffusionTransformationEnum::_update_current_to_loc() {
      auto it = m_to_loc.begin();
      for(auto &t : m_current->specie_traj()) {
        t.to = *it++;
      }
    }

#define  PRIM_PERIODIC_DIFF_TRANS_ORBITS_INST(INSERTER, CLUSTER_IT) \
    \
    template INSERTER make_prim_periodic_diff_trans_orbits<INSERTER, CLUSTER_IT>( \
      CLUSTER_IT begin, \
      CLUSTER_IT end, \
      double xtal_tol, \
      INSERTER result, \
      const PrimClex*); \
    \

#define _VECTOR_IT(ORBIT) std::vector<ORBIT>::iterator
#define _VECTOR_INSERTER(ORBIT) std::back_insert_iterator<std::vector<ORBIT> >

#define  PRIM_PERIODIC_DIFF_TRANS_ORBITS_VECTOR_INST(DIFF_TRANS_ORBIT, CLUSTER_ORBIT) \
      PRIM_PERIODIC_DIFF_TRANS_ORBITS_INST( \
        _VECTOR_INSERTER(DIFF_TRANS_ORBIT), \
        _VECTOR_IT(CLUSTER_ORBIT))

    PRIM_PERIODIC_DIFF_TRANS_ORBITS_VECTOR_INST(
      PrimPeriodicDiffTransOrbit,
      PrimPeriodicIntegralClusterOrbit)

  }
}
