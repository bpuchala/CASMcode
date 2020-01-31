#include <cstring>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include "casm/global/definitions.hh"
#include "casm/misc/CASM_Eigen_math.hh"
#include "casm/casm_io/json/jsonParser.hh"
#include "casm/app/DirectoryStructure.hh"
#include "casm/app/HamiltonianModules.hh"
#include "casm/app/ProjectBuilder.hh"
#include "casm/app/casm_functions.hh"
#include "casm/app/AppIO.hh"
#include "casm/completer/Handlers.hh"
#include "casm/crystallography/CanonicalForm.hh"
#include "casm/crystallography/Structure.hh"

namespace CASM {

  namespace Completer {
    InitOption::InitOption(): OptionHandlerBase("init") {}

    void InitOption::initialize() {
      add_help_suboption();
      add_file_path_suboption();
      add_configlist_suboption();
      add_confignames_suboption();
      add_dofs_suboption();
      m_desc.add_options()
      ("sub,s", "Initialize a project in sub-directory of an existing project. After initialization, commands executed below the sub-directory will act on the sub-directory project; commmands executed above the sub-directory will act on the original project.")
      ("generate-prim", "Create prim.json file for specified configuration(s). Each prim.json will be written to the training_data directory of the corresponding configuration.")
      ("force,f", "Force using a non-reduced, non-primitive, or left-handed PRIM.");
      return;
    }
  }

  // ///////////////////////////////////////
  // 'init' function for casm
  //    (add an 'if-else' statement in casm.cpp to call this)

  int init_command(const CommandArgs &args) {

    std::string name;
    po::variables_map vm;
    HamiltonianModules modules;

    /// Set command line options using boost program_options
    Completer::InitOption init_opt;

    try {
      po::store(po::parse_command_line(args.argc(), args.argv(), init_opt.desc()), vm); // can throw

      /** --help option
      */
      if(vm.count("help")) {
        args.log() << "\n";
        args.log() << init_opt.desc() << std::endl;

        return 0;
      }

      if(vm.count("desc")) {
        args.log() << "\n";
        args.log() << init_opt.desc() << std::endl;

        args.log() << "DESCRIPTION                                                \n" <<
                   "    Initialize a new CASM project in the current directory.\n" <<
                   "    - Expects a prim.json file in the current directory    \n" <<
                   "    - If not found, looks for a PRIM file in the current   \n" <<
                   "      directory and creates prim.json.                     \n\n";

        return 0;
      }

      po::notify(vm); // throws on error, so do after help in case
      // there are any problems
    }
    catch(po::error &e) {
      args.err_log() << "ERROR: " << e.what() << std::endl << std::endl;
      args.err_log() << init_opt.desc() << std::endl;
      return ERR_INVALID_ARG;
    }
    catch(std::exception &e) {
      args.err_log() << "Unhandled Exception reached the top of main: "
                     << e.what() << ", application will now exit" << std::endl;
      return ERR_UNKNOWN;

    }


    fs::path root = fs::current_path();
    if(!init_opt.file_path().empty()) {
      root = init_opt.file_path();
    }
    else if(!args.root.empty()) {
      root = args.root;
    }

    fs::path existing = find_casmroot(root);
    if(vm.count("sub")) {
      if(existing == root) {
        args.log() << "Directory '" << root << "' is already the head directory of a casm project." << std::endl;
        return ERR_OTHER_PROJ;
      }
      if(!root.empty()) {
        args.log() << "No existing project found. Cannot create sub-directory project at '" << root << "'." << std::endl;
        return ERR_OTHER_PROJ;

      }
    }
    else if(!existing.empty()) {
      args.log() << "Already in a casm project." << std::endl;
      return ERR_OTHER_PROJ;
    }

    args.log() << "\n***************************\n" << std::endl;

    DirectoryStructure dir(root);
    Structure prim;

    // if prim.json does not exist, try to read PRIM and create prim.json
    if(!fs::is_regular_file(dir.prim())) {

      if(!fs::is_regular_file(dir.PRIM())) {
        args.log() << "Error in 'casm init': Neither 'prim.json' nor 'PRIM' found.\n\n";

        args.log() << "Run 'casm format --prim' for the format of the 'prim.json' file.\n\n";

        args.log() << "For step by step help use: 'casm status -n'\n\n";

        return ERR_MISSING_INPUT_FILE;
      }

      try {
        fs::ifstream poscar_prim;
        poscar_prim.open(dir.PRIM());
        prim.read(poscar_prim);
        poscar_prim.close();
      }
      catch(std::runtime_error &e) {

        args.err_log() << "ERROR: No prim.json exists. PRIM exists, but it could not be read.\n";
        args.err_log() << e.what() << std::endl;
        return ERR_INVALID_INPUT_FILE;
      }


      std::string poscar_prim_title = prim.title();
      args.log() << "Converting 'PRIM' to 'prim.json'.\n\n" << std::endl;

      args.log() << "Please enter a short title for this project.\n";
      args.log() << "  Use something suitable as a prefix for files specific to this project, such as 'ZrO' or 'TiAl'.\n\n";

      std::string ttitle;
      args.log() << "Title: ";
      std::cin >> ttitle;
      prim.set_title(ttitle);
      args.log() << "\n\n";

      jsonParser json;
      write_prim(prim, json, FRAC);
      json["description"] = poscar_prim_title;
      fs::ofstream primfile(dir.prim());
      json.print(primfile);
      primfile.close();

    }

    jsonParser prim_json;

    try {

      prim_json = jsonParser(dir.prim());

      prim = Structure(read_prim(prim_json, TOL, &modules));
    }
    catch(std::runtime_error &e) {
      args.err_log() << e.what() << std::endl;

      return ERR_INVALID_INPUT_FILE;
    }

    /// Check if PRIM is primitive
    BasicStructure<Site> true_prim;
    true_prim.set_title(prim.title());
    if(!prim.is_primitive(true_prim)) {
      if(!vm.count("force")) {
        args.err_log() << "ERROR: The structure in the prim.json file is not primitive. Writing the most       \n"
                       << "       primitive structure to file 'prim.true.json'.\n\n";

        Structure tmp(true_prim);
        Lattice lat_niggli = xtal::canonical::equivalent(true_prim.lattice(), tmp.point_group(), TOL);
        tmp.set_lattice(lat_niggli, CART);

        fs::ofstream primfile(root / "prim.true.json");
        jsonParser json;
        write_prim(tmp, json, FRAC);
        json["description"] = prim_json["description"];
        json.print(primfile);
        primfile.close();

        args.err_log() << "If you want to use the current prim.json anyway, re-run with the --force option. Some\n"
                       << "CASM features cannot be used with a non-primitive starting structure.\n";
        return ERR_INVALID_INPUT_FILE;
      }
      else {
        args.err_log() << "WARNING: The structure in the prim.json file is not primitive. Continuing anyway    \n"
                       << "         because the --force option is on.\n\n";
      }
    }

    /// Check that the PRIM is in reduced form:
    Lattice niggli_lat = xtal::canonical::equivalent(prim.lattice(), prim.point_group(), TOL);

    bool is_standard_niggli = almost_equal(niggli_lat.lat_column_mat(), prim.lattice().lat_column_mat());

    if(!is_standard_niggli) {
      if(!vm.count("force")) {
        if(!is_standard_niggli) {
          args.err_log() << "ERROR: The structure in the prim.json file is not the niggli cell in the CASM standard\n"
                         << "       orientation. Writing the suggested structure to 'prim.niggli.json'.\n\n"
                         << "       If you want to use the current prim.json anyway, re-run with the --force option.\n";
        }

        Structure tmp(true_prim);
        Lattice lat_niggli = xtal::canonical::equivalent(true_prim.lattice(), tmp.point_group(), TOL);
        tmp.set_lattice(lat_niggli, CART);

        fs::ofstream primfile(root / "prim.niggli.json");
        jsonParser json;
        write_prim(tmp, json, FRAC);
        json["description"] = prim_json["description"];
        json.print(primfile);
        primfile.close();

        return ERR_INVALID_INPUT_FILE;
      }
      else {
        args.err_log() << "WARNING: The structure in the prim.json file is not the standard orientation Niggli\n"
                       << "         cell. Continuing anyway because the --force option is on.\n\n";
        //return 1;
      }

    }

    //Check if the lattice is right handed, and if not print PRIM.right_handed.json
    if(!prim.lattice().is_right_handed()) {
      if(!vm.count("force")) {
        args.err_log() << "ERROR: The structure in prim.json is not right-handed. Some electronic-"
                       << "structure codes will not accept this input. If you would like to "
                       << "keep this PRIM, re-run with the --force option. Writing the "
                       << "right-handed structure to PRIM.right_handed.json" << std::endl;

        prim.set_lattice(Lattice(prim.lattice()).make_right_handed(), CART);
        prim.within();

        fs::ofstream primfile(root / "prim.right_handed.json");
        jsonParser json;
        write_prim(prim, json, FRAC);
        json["description"] = prim_json["description"];
        json.print(primfile);
        primfile.close();
        return ERR_INVALID_INPUT_FILE;
      }
      else {
        args.err_log() << "WARNING: The structure in the prim.json file is not right-handed. Continuing anyway    \n"
                       << "         because the --force option is on.\n\n";
      }
    }

    try {
      args.log() << "Initializing CASM project '" << prim.title() << "'" << std::endl;
      ProjectBuilder builder(root, prim.title(), "formation_energy");
      builder.build();
    }
    catch(std::runtime_error &e) {
      args.err_log() << "ERROR: Could not build CASM project.\n";
      args.err_log() << e.what() << std::endl;
      return ERR_INVALID_INPUT_FILE;
    }

    args.log() << "  DONE" << std::endl;
    args.log() << std::endl;

    return 0;
  };

}
