import os, math, sys, json, re, warnings, shutil
import pbs
import vasp
import casm
import casm.project
from casm.project import Project, Selection
import vaspwrapper
import configproperties

class ConfigProperties(ConfigPropertiesBase):
    """ Derived class with additional configuration properties and errors

    Additional Attributes
    --------------------

      n_images: int
        number of images in the calculation

    """

    def __init__(self, config_data):
        """ 
             Add additional properties and error messages
        """
        ConfigPropertiesBase.__init__(self, config_data)

        ## Error message if "n_images" not present in settings
        if not "n_images" in self.settings:
            raise vaspwrapper.VaspWrapperError("Could not find \"n_images\" in \"calc.json\" in an appropriate \"settings\" directory")
            sys.stdout.flush()

        self.n_images = self.settings["n_images"]
        # append the n_images to calcdir
        self.calcdir = os.path.join(self.calcdir, "N_images_{}".format(self.n_images))

class Neb(object):
    """The Relax class contains functions for setting up, executing, and parsing a VASP relaxation.

        The relaxation creates the following directory structure:
        config/
          calctype.name/
              run.0/
              ....

        'run.i' directories are only created when ready.

        This automatically looks for VASP settings files using:
          casm.project.DirectoryStructure.settings_path_crawl

    """
    def __init__(self, selection, auto=True, sort=True):
        """
        Construct a VASP neb job object.

        Arguments
        ----------

            selection: casm.project.Selection object, default= yet to be implemented #todo
              Selection of all DiffTransConfigurations to submit a NEB calculation. 
              default should be MASTER selection and yet to be implemented

            auto: boolean, optional, default=True,
              Use True to use the pbs module's JobDB to manage pbs jobs

            sort: boolean, optional, default=True,
              Use True to sort atoms in POSCAR by type

        """
        print "Construct a casm.vaspwrapper.Relax instance:"

        self.selection = selection
        self.auto = auto
        self.sort = sort

    def _pre_setup(self):
        for config_data in self.selection.data:
            config_obj = ConfigProperties(config_data)
            try:
                os.makedirs(config_obj.calcdir)
                for i in range(config_obj.n_images + 2):
                    os.makedirs(os.path.join(config_obj.calcdir, str(i).zfill(2)))
            except:
                pass

    def setup(self):
        """ Setup initial relaxation run

            Uses the following files from the most local .../settings/calctype.name directory:
                INCAR: VASP input settings
                KPOINTS: VASP kpoints settings
                POSCAR: reference for KPOINTS if KPOINTS mode is not A/AUTO/Automatic
                SPECIES: info for each species such as which POTCAR files to use, MAGMOM, GGA+U, etc.

            Uses the following files from the .../config directory:
                POS: structure of the configuration to be relaxed

        """
        self._pre_setup()
        ## write the selection Interpolation command
        proj = self.selection.proj
        dict = {}
        for config_data in self.selection.data:
            config_obj = ConfigProperties(config_data)
            conf_dict = {"n_images" : self.n_images, "calctype" : self.calctype}
            dict[config_obj.configname] = conf_dict
        args = "enum --method DiffTransConfigInterpolation -j neb_tmp.json"
        args += json.dumps(dict)
        output = proj.command(args)
        #remove the tmp file

        for config_data in self.selection.data:
            config_obj = ConfigProperties(config_data)
            # Find required input files in CASM project directory tree
            vaspfiles = casm.vaspwrapper.vasp_input_file_names(config_obj.casm_directories, config_obj.configname,
                                                               config_obj.clex, config_obj.calc_subdir)
            incarfile, prim_kpointsfile, prim_poscarfile, temp_poscarfile, speciesfile = vaspfiles
            # Find optional input files
            extra_input_files = []
            for s in config_obj.settings["extra_input_files"]:
                extra_input_files.append(config_obj.casm_directories.settings_path_crawl(s, config_obj.configname,
                                                                                         config_obj.clex, config_obj.calc_subdir))
                if extra_input_files[-1] is None:
                    raise vasp.VaspError("Neb.setup failed. Extra input file " + s + " not found in CASM project.")
            if config_obj.settings["initial"]:
                extra_input_files += [config_obj.casm_directories.settings_path_crawl(config_obj.settings["initial"],
                                                                                      config_obj.configname, config_obj.clex,
                                                                                      config_obj.calc_subdir)]
                if extra_input_files[-1] is None:
                    raise vasp.VaspError("Neb.setup failed. No initial INCAR file " + config_obj.settings["initial"] + " found in CASM project.")
            if config_obj.settings["final"]:
                extra_input_files += [config_obj.casm_directories.settings_path_crawl(config_obj.settings["final"],
                                                                                      config_obj.configname, config_obj.clex,
                                                                                      config_obj.calc_subdir)]
                if extra_input_files[-1] is None:
                    raise vasp.VaspError("Neb.setup failed. No final INCAR file " + config_obj.settings["final"] + " found in CASM project.")
            sys.stdout.flush()

            #make vasp input files
            sample_super_poscarfile = os.path.join(config_obj.calcdir, "00", "POSCAR")
            vasp.io.write_vasp_input(config_obj.calcdir, incarfile, prim_kpointsfile, prim_poscarfile,
                                     sample_super_poscarfile, speciesfile, self.sort, extra_input_files,
                                     config_obj.settings["strict_kpoints"])
            ## settings the images tag in incar file
            tmp_dict = {"images": config_obj.settings["n_images"]}
            vasp.io.set_incar_tag(tmp_dict, config_obj.calcdir)



######################

            
class submit(object):

    def __init__(self, selection = None):
        self.selection = selection

    
    """Submit a PBS job for this VASP relaxation"""

        print "Submitting..."
        print "Configuration:", self.configname
        # first, check if the job has already been submitted and is not completed
        db = pbs.JobDB()
        print "Calculation directory:", self.calcdir
        id = db.select_regex_id("rundir", self.calcdir)
        print "JobID:", id
        sys.stdout.flush()
        if id != []:
            for j in id:
                job = db.select_job(j)
                # taskstatus = ["Incomplete","Complete","Continued","Check","Error:.*","Aborted"]
                # jobstatus = ["C","Q","R","E","W","H","M"]
                if job["jobstatus"] != "C":
                    print "JobID:", job["jobid"], "  Jobstatus:", job["jobstatus"], "  Not submitting."
                    sys.stdout.flush()
                    return
                #elif job["taskstatus"] in ["Complete", "Check"] or re.match( "Error:.*", job["taskstatus"]):
                #    print "JobID:", job["jobid"], "  Taskstatus:", job["taskstatus"], "  Not submitting."
                #    sys.stdout.flush()
                #    return


        # second, only submit a job if relaxation status is "incomplete"

        # construct the Relax object
        calculation = vasp.Neb(self.calcdir, self.run_settings())

        # check the current status
        (status, task) = calculation.status()

        if status == "complete":
            print "Status:", status, "  Not submitting."
            sys.stdout.flush()

            # ensure job marked as complete in db
            if self.auto:
                for j in id:
                    job = db.select_job(j)
                    if job["taskstatus"] == "Incomplete":
                        try:
                            pbs.complete_job(jobid=j)
                        except (pbs.PBSError, pbs.JobDBError, pbs.EligibilityError) as e:
                            print str(e)
                            sys.stdout.flush()

            # ensure results report written
            if not os.path.isfile(os.path.join(self.calcdir, "properties.calc.json")):
                self.finalize()

            return

        elif status == "not_converging":
            print "Status:", status, "  Not submitting."
            sys.stdout.flush()
            return

        elif status != "incomplete":
            raise vaspwrapper.VaspWrapperError("unexpected relaxation status: '" + status + "' and task: '" + task + "'")
            sys.stdout.flush()
            return


        print "Preparing to submit a VASP relaxation PBS job"
        sys.stdout.flush()

        # cd to configdir, submit jobs from configdir, then cd back to currdir
        currdir = os.getcwd()
        os.chdir(self.calcdir)

        # determine the number of atoms in the configuration
        #print "Counting atoms in the POSCAR"
        #sys.stdout.flush()
        #pos = vasp.io.Poscar(os.path.join(self.configdir,"POS"))
        #N = len(pos.basis)

        # determine the number of nodes from settings #TODO
        # still has to implement a option to use either atoms_per_proc or nodes_per_image or images_per_node
        nodes = int(self.n_images) * float(self.settings["nodes_per_image"])
        # construct command to be run
        cmd = ""
        if self.settings["preamble"] is not None:
        # Append any instructions given in the 'preamble' file, if given
            preamble = self.casm_directories.settings_path_crawl(self.settings["preamble"],
                                                                 self.configname, self.clex,
                                                                 self.calc_subdir)
            with open(preamble) as my_preamble:
                cmd += "".join(my_preamble)
        # Or just execute a single prerun line, if given
        if self.settings["prerun"] is not None:
            cmd += self.settings["prerun"] + "\n"
        cmd += "python -c \"import casm.vaspwrapper; casm.vaspwrapper.Neb('" + self.configdir + "').run()\"\n"
        if self.settings["postrun"] is not None:
            cmd += self.settings["postrun"] + "\n"

        print "Constructing a PBS job"
        sys.stdout.flush()
        # construct a pbs.Job
        job = pbs.Job(name=casm.jobname(self.configdir),\
                      account=self.settings["account"],\
                      nodes=int(nodes),
                      ppn=int(self.settings["ppn"]),\
                      walltime=self.settings["walltime"],\
                      pmem=self.settings["pmem"],\
                      qos=self.settings["qos"],\
                      queue=self.settings["queue"],\
                      message=self.settings["message"],\
                      email=self.settings["email"],\
                      priority=self.settings["priority"],\
                      command=cmd,\
                      auto=self.auto)

        print "Submitting"
        sys.stdout.flush()
        # submit the job
        job.submit()
        self.report_status("submitted")

        # return to current directory
        os.chdir(currdir)

        print "CASM VASPWrapper relaxation PBS job submission complete\n"
        sys.stdout.flush()


    def run_settings(self):
        """ Set default values based on runtime environment"""
        settings = dict(self.settings)

        # set default values

        if settings["npar"] == "CASM_DEFAULT":
            if "PBS_NUM_NODES" in os.environ:
                settings["npar"] = int(os.environ["PBS_NUM_NODES"])
            elif "SLURM_JOB_NUM_NODES" in os.environ:
                settings["npar"] = int(os.environ["SLURM_JOB_NUM_NODES"])
            else:
                settings["npar"] = None
        elif settings["npar"] == "VASP_DEFAULT":
            settings["npar"] = None

        if settings["npar"] is None:
            if settings["ncore"] == "CASM_DEFAULT":
                if "PBS_NUM_PPN" in os.environ:
                    settings["ncore"] = int(os.environ["PBS_NUM_PPN"])
                elif "SLURM_CPUS_ON_NODE" in os.environ:
                    settings["ncore"] = int(os.environ["SLURM_CPUS_ON_NODE"])
                else:
                    settings["ncore"] = None
            elif settings["ncore"] == "VASP_DEFAULT":
                settings["ncore"] = 1
        else:
            settings["ncore"] = None

        if settings["ncpus"] is None or settings["ncpus"] == "CASM_DEFAULT":
            if "PBS_NP" in os.environ:
                settings["ncpus"] = int(os.environ["PBS_NP"])
            elif "SLURM_NTASKS" in os.environ:
                settings["ncpus"] = int(os.environ["SLURM_NTASKS"])
            else:
                settings["ncpus"] = None

        if settings["run_limit"] is None or settings["run_limit"] == "CASM_DEFAULT":
            settings["run_limit"] = 10

        return settings

class run(object):
        """ Setup input files, run a vasp relaxation, and report results """
    def __init__(self, configdir):
        self.configdir = configdir
        # construct the Neb object
        calculation = vasp.Neb(self.calcdir, self.run_settings())

        # check the current status
        (status, task) = calculation.status()

        if status == "complete":
            print "Status:", status
            sys.stdout.flush()

            # mark job as complete in db
            if self.auto:
                try:
                    pbs.complete_job()
                except (pbs.PBSError, pbs.JobDBError, pbs.EligibilityError) as e:
                    print str(e)
                    sys.stdout.flush()

            # write results to properties.calc.json
            self.finalize()
            return

        elif status == "not_converging":
            print "Status:", status
            self.report_status("failed","run_limit")
            print "Returning"
            sys.stdout.flush()
            return

        elif status == "incomplete":

            if task == "setup":
                self.setup()

            self.report_status("started")
            (status, task) = calculation.run()

        else:
            self.report_status("failed", "unknown")
            raise vaspwrapper.VaspWrapperError("unexpected relaxation status: '" + status + "' and task: '" + task + "'")
            sys.stdout.flush()


        # once the run is done, update database records accordingly

        if status == "not_converging":

            # mark error
            if self.auto:
                try:
                    pbs.error_job("Not converging")
                except (pbs.PBSError, pbs.JobDBError) as e:
                    print str(e)
                    sys.stdout.flush()

            print "Not Converging!"
            sys.stdout.flush()
            self.report_status("failed", "run_limit")

            # print a local settings file, so that the run_limit can be extended if the
            #   convergence problems are fixed

            config_set_dir = self.casm_directories.configuration_calc_settings_dir(self.configname, self.clex, self.calc_subdir)

            try:
                os.makedirs(config_set_dir)
            except:
                pass
            settingsfile = os.path.join(config_set_dir, "neb.json")
            vaspwrapper.write_settings(self.settings, settingsfile)

            print "Writing:", settingsfile
            print "Edit the 'run_limit' property if you wish to continue."
            sys.stdout.flush()
            return

        elif status == "complete":

            # mark job as complete in db
            if self.auto:
                try:
                    pbs.complete_job()
                except (pbs.PBSError, pbs.JobDBError, pbs.EligibilityError) as e:
                    print str(e)
                    sys.stdout.flush()

            # write results to properties.calc.json
            self.finalize()

        else:
            self.report_status("failed", "unknown")
            raise vaspwrapper.VaspWrapperError("vasp relaxation complete with unexpected status: '" + status + "' and task: '" + task + "'")
            sys.stdout.flush()

    def report_status(self, status, failure_type=None):
        """Report calculation status to status.json file in configuration directory.

        Args:
            status: string describing calculation status. Currently used values are
                 not_submitted
                 submitted
                 complete
                 failed
             failure_type: optional string describing reason for failure. Currently used values are
                 unknown
                 electronic_convergence
                 run_limit"""

        output = dict()
        output["status"] = status
        if failure_type is not None:
            output["failure_type"] = failure_type

        outputfile = os.path.join(self.calcdir, "status.json")
        with open(outputfile, 'w') as file:
            file.write(json.dumps(output, file, cls=casm.NoIndentEncoder, indent=4, sort_keys=True))
        print "Wrote " + outputfile
        sys.stdout.flush()


    def finalize(self):
        if self.is_converged():
        # write properties.calc.json
            vaspdir = os.path.join(self.calcdir, "run.final")
            speciesfile = self.casm_directories.settings_path_crawl("SPECIES", self.configname, self.clex, self.calc_subdir)
            output = self.properties(vaspdir, True, speciesfile)
            outputfile = os.path.join(self.calcdir, "properties.calc.json")
            with open(outputfile, 'w') as file:
                file.write(json.dumps(output, file, cls=casm.NoIndentEncoder, indent=4, sort_keys=True))
            print "Wrote " + outputfile
            all_image_folders = [int(i.strip().split('_')[-1]) for i in os.listdir(self.calcdir) if "N_images" in i]
            num_images = [int(self.calcdir.strip().split('_')[-1])]
            if num_images == all_image_folders:
                shutil.copy(os.path.join(self.calcdir, "properties.calc.json"),
                            os.path.join(os.path.split(self.calcdir)[0], "properties.calc.json"))
                print "As the present run has highest number of images copied {0} to {1}".format(os.path.join(self.calcdir, "properties.calc.json"),
                                                                                             os.path.join(os.path.split(self.calcdir)[0],
                                                                                                          "properties.calc.json"))
            sys.stdout.flush()
            self.report_status('complete')

    def is_converged(self):
        # Check for electronic convergence in completed calculations. Returns True or False.

        # Verify that the last relaxation reached electronic convergence
        calculation = vasp.Neb(self.calcdir, self.run_settings())
        for i in range(len(calculation.rundir)):
            try:
                #vrun = vasp.io.Vasprun(os.path.join(self.calcdir, calculation.rundir[-i-1], "vasprun.xml"))
                vrun_oszicar = vasp.io.Oszicar(os.path.join(self.calcdir, calculation.rundir[-i-1], "01", "OSZICAR"))
                vrun_nelm = vasp.io.get_incar_tag("NELM", os.path.join(self.calcdir, calculation.rundir[-i-1]))
                if vrun_nelm == None: vrun_nelm = 60 ##pushing the default. may be write an addon to get it from outcar
                if len(vrun_oszicar.num_elm[-1]) >= vrun_nelm:
                    print('The last relaxation run (' +
                          os.path.basename(relaxation.rundir[-i-1]) +
                          ') failed to achieve electronic convergence; properties.calc.json will not be written.\n')
                    self.report_status('failed', 'electronic_convergence')
                    return False
                break
            except:
                pass

        # Verify that the final static run reached electronic convergence
        #vrun = vasp.io.Vasprun(os.path.join(self.calcdir, "run.final", "vasprun.xml"))
        vrun_oszicar = vasp.io.Oszicar(os.path.join(self.calcdir, "run.final", "01", "OSZICAR"))
        vrun_nelm = vasp.io.get_incar_tag("NELM", os.path.join(self.calcdir, "run.final"))
        if vrun_nelm == None: vrun_nelm = 60 ##pushing the default. may be write an addon to get it from outcar
        if vrun_oszicar.num_elm[-1] >= vrun_nelm:
            print('The final run failed to achieve electronic convergence; properties.calc.json will not be written.\n')
            self.report_status('failed', 'electronic_convergence')
            return False

        return True

    @staticmethod
    def properties(vaspdir, use_poscarfile=False, speciesfile=None):
        """Report results to properties.calc.json file in configuration directory, after checking for electronic convergence."""
        final_output = []
        num_images = vasp.io.get_incar_tag("IMAGES", vaspdir)
        for img in [str(j).zfill(2) for j in range(1, num_images+1)]:
            output = dict()
            #vrun = vasp.io.Vasprun( os.path.join(vaspdir, "vasprun.xml") )
            vrun_oszicar = vasp.io.Oszicar(os.path.join(vaspdir, img, "OSZICAR"))
            vrun_outcar = vasp.io.Outcar(os.path.join(vaspdir, img, "OUTCAR"))

            # the calculation is run on the 'sorted' POSCAR, need to report results 'unsorted'

            if (use_poscarfile is not None) and (speciesfile is not None):
                species_settings = vasp.io.species_settings(speciesfile)
                super_poscar = vasp.io.Poscar(os.path.join(vaspdir, img, "POSCAR"), species_settings)
                super_contcar = vasp.io.Poscar(os.path.join(vaspdir, img, "CONTCAR"), species_settings)
                unsort_dict = super_poscar.unsort_dict()
            else: # not implemented #TODO
                # fake unsort_dict (unsort_dict[i] == i)
                unsort_dict = dict(zip(range(0, len(vrun.basis)), range(0, len(vrun.basis))))
                super_poscar = vasp.io.Poscar(os.path.join(vaspdir, "POSCAR"))

            # unsort_dict:
            #   Returns 'unsort_dict', for which: unsorted_dict[orig_index] == sorted_index;
            #   unsorted_dict[sorted_index] == orig_index
            #   For example:
            #     'unsort_dict[0]' returns the index into the unsorted POSCAR of the first atom in the sorted POSCAR

            output["Image_number"] = img
            output["atom_type"] = super_poscar.type_atoms
            output["atoms_per_type"] = super_poscar.num_atoms
            output["coord_mode"] = super_poscar.coord_mode

            # as lists
            output["relaxed_forces"] = [ None for i in range(len(vrun_outcar.forces))]
            for i, v in enumerate(vrun_outcar.forces):
                output["relaxed_forces"][unsort_dict[i] ] = casm.NoIndent(vrun_outcar.forces[i])

            output["relaxed_lattice"] = [casm.NoIndent(list(v)) for v in super_contcar.lattice()]
            output["relaxed_basis"] = [None for i in range(len(super_contcar.basis))]
            for i, v in enumerate(super_contcar.basis):
                output["relaxed_basis"][unsort_dict[i]] = casm.NoIndent(list(super_contcar.basis[i].position))

            output["relaxed_energy"] = vrun_oszicar.E[-1]
            final_output.append(output)

        return final_output