"""Module containing a class to represent a group of simulations defined by one
or more simulation sequences.

"""

import os
import sys
import copy
import json
import shutil
import time
from datetime import datetime

import numpy as np
import jsonpickle
import jsonpickle.ext.numpy as jsonpickle_numpy
jsonpickle_numpy.register_handlers()

from matsim import (JS_TEMPLATE_DIR, utils, parse_opt, CONFIG, DB_CONFIG,
                    database, OPTSPEC)
from matsim.simulation import BaseUpdate, apply_base_update
from matsim.simulation import process
from matsim.simulation.sequence import SimSequence
from matsim.atomistic.simulation.castep import CastepSimulation
from matsim.atomistic.simulation.lammps import LammpsSimulation
from matsim.utils import (nest, merge, prt, dict_from_list, mut_exc_args,
                          set_nested_dict, get_recursive, update_dict)
from matsim.resources import Stage, Scratch, Archive, ResourceConnection
from matsim.readwrite import replace_in_file, delete_line, add_line, write_list_file

# Maps a given software to a Simulation class.
SOFTWARE_CLASS_MAP = {
    'castep': CastepSimulation,
    'lammps': LammpsSimulation
}


def write_jobscript(path, calc_paths, method, num_cores, is_sge, job_array, executable,
                    scratch_os=None, scratch_path=None, parallel_env=None,
                    selective_submission=False, module_load=None,
                    job_name=None, seedname=None):
    """
    Write a jobscript file whose execution runs calculation input files.

    Parameters
    ----------
    path : str
        Directory in which to save the generated jobscript file.
    calc_paths : list of str
        Directories in which calculations are to be run.
    method : str
        Either 'castep' or 'lammps'
    num_cores : int
        Number of processor cores to use for the calculations.
    is_sge : bool
        If True, jobscript is generated for the SGE batch scheduler. If False,
        jobscript is generated to immediately run the calculations.
    job_array : bool
        Only applicable if `sge` is True. If True, calculations are submitted
        as an SGE job array. If False, calculations are submitted in one go. If
        the number of calculations is one (i.e. len(`calc_paths`) == 1), this
        will be set to False. Setting this to False can be handy for many small
        calculations which won't take a long time to complete. If submitted as
        a job array, a significant fraction of the total completion time may be
        queuing.
    executable : str
        Executable command.
    scratch_os : str, optional
        Either 'nt' (Windows) or 'posix' (Unix-like, MacOS). The operating
        system on which the jobscript file will be executed. Default is to
        query to system on which this script is invoked i.e. `os.name`.
    scratch_path : str, optional
        Directory in which the jobscript is to be executed. Specify this path
        if the jobscript will be executed in a location different to the
        directory `path`, in which it is generated. By default, this is set to
        the same string as `path`.
    parallel_env : str, optional
        The SGE parallel environment on which to submit the calculations. Only
        applicable in `num_cores` > 1. Default is None.
    selective_submission : bool, optional
        Only applicable if `sge` is True. If True, the SGE task id flag `-t`
        [1] will be excluded from the jobscript file and instead this flag will
        be expected as a command line argument when executing the jobscript:
        e.g. "qsub jobscript.sh -t 1-10:2". Default is False.
    module_load : str, optional
        A string representing the path to a module to load within the
        jobscript. If specified, the statement "module load `module_load`" will
        be added to the top of the jobscript.
    job_name : str, optional
        Only applicable if `sge` is True. Default is None.
    seedname : str, optional
        Must be set if `method` is 'castep'.

    Returns
    -------
    None

    References
    ----------
    [1] http://gridscheduler.sourceforge.net/htmlman/htmlman1/qsub.html

    TODO:
    -   Add option for specifying parallel environment type for sge jobscripts

    """

    # General validation:

    if method == 'castep' and seedname is None:
        raise ValueError('`seedname` must be specified for CASTEP jobscripts.')

    if num_cores > 1 and parallel_env is None:
        raise ValueError('`parallel_env` must be set if `num_cores` > 1.')

    num_calcs = len(calc_paths)

    if num_cores <= 0:
        raise ValueError('Num cores not valid.')
    elif num_cores == 1:
        multi_type = 'serial'
    else:
        multi_type = 'parallel'

    if is_sge:
        sge_str = 'sge'
    else:
        sge_str = 'no_sge'

    if num_calcs == 1:
        job_array = False

    if job_array:
        job_arr_str = 'job_array'
    else:
        job_arr_str = 'single_job'

    if scratch_path is None:
        scratch_path = path

    if scratch_os is None:
        scratch_os = os.name

    if scratch_os == 'nt':
        scratch_path_sep = '\\'
    elif scratch_os == 'posix':
        scratch_path_sep = '/'

    # Get the template file name:
    tmp_fn = method + '_' + sge_str + '_' + \
        multi_type + '_' + scratch_os + '_' + job_arr_str + '.txt'

    # Get the template file path:
    tmp_path = os.path.join(JS_TEMPLATE_DIR, tmp_fn)

    # Write text file with all calc paths
    dirlist_fn = 'dir_list.txt'
    dir_list_path_stage = os.path.join(path, dirlist_fn)
    dir_list_path_scratch = scratch_path_sep.join([scratch_path, dirlist_fn])
    write_list_file(dir_list_path_stage, calc_paths)

    if scratch_os == 'posix':
        js_ext = 'sh'
    elif scratch_os == 'nt':
        js_ext = 'bat'

    js_name = 'jobscript.' + js_ext

    # Copy template file to path
    js_path = os.path.join(path, js_name)
    shutil.copy(tmp_path, js_path)

    # Add module load to jobscript:
    if module_load is not None:
        add_line(js_path, 1, '')
        add_line(js_path, 2, 'module load {}'.format(module_load))
        delete_line(js_path, '#$ -V')
        replace_in_file(js_path, '#!/bin/bash', '#!/bin/bash --login')

    # Make replacements in template file:
    replace_in_file(js_path, '<replace_with_dir_list>', dir_list_path_scratch)

    if multi_type == 'parallel':
        replace_in_file(js_path, '<replace_with_num_cores>', str(num_cores))
        replace_in_file(js_path, '<replace_with_pe>', parallel_env)

    if is_sge:
        if job_name is not None:
            replace_in_file(js_path, '<replace_with_job_name>', job_name)
        else:
            delete_line(js_path, '<replace_with_job_name>')

        if selective_submission:
            delete_line(js_path, '#$ -t')
        else:
            replace_in_file(js_path, '<replace_with_job_index_range>',
                            '1-' + str(num_calcs))

    if method == 'castep':
        replace_in_file(js_path, '<replace_with_seed_name>', seedname)

    replace_in_file(js_path, '<replace_with_executable>', executable)

    # For `method` == 'lammps', `sge` == True, `job_array` == False, we need
    # a helper jobscript, called by the SGE jobscript:
    if method == 'lammps' and is_sge and not job_array:

        if multi_type != 'serial' or scratch_os != 'posix':
            raise NotImplementedError('Jobscript parameters not supported.')

        help_tmp_path = os.path.join(
            JS_TEMPLATE_DIR, 'lammps_no_sge_serial_posix_single_job.txt')

        help_js_path = os.path.join(path, 'lammps_single_job.sh')
        shutil.copy(help_tmp_path, help_js_path)
        replace_in_file(help_js_path, '<replace_with_dir_list>',
                        dir_list_path_scratch)

    if is_sge:
        # Make a directory for SGE-related output. E.g. .o and .e files from CSF.
        os.makedirs(os.path.join(path, 'sge'))

    return job_array


def write_process_jobscript(path, job_name, dependency, num_calcs, human_id,
                            run_group_idx, job_array, selective_submission):
    """Write a job array dependency jobscript to auto-process a run group."""

    # Get the template file path
    tmp_path = os.path.join(JS_TEMPLATE_DIR, 'process_sge.txt')

    # Copy template file to path
    js_path = os.path.join(path, 'process_jobscript.sh')
    shutil.copy(tmp_path, js_path)

    # Make replacements in template file:
    replace_in_file(js_path, '<replace_with_process_job_name>', job_name)
    job_range = '1-' + str(num_calcs)
    replace_in_file(js_path, '<replace_with_job_index_range>', job_range)
    replace_in_file(js_path, '<replace_with_run_group_job_name>', dependency)
    replace_in_file(js_path, '<replace_with_sim_group_human_id>', human_id)
    replace_in_file(js_path, '<replace_with_run_group_idx>',
                    str(run_group_idx))

    job_dep_cmd = 'hold_jid_ad' if job_array else 'hold_jid'
    replace_in_file(js_path, '<replace_with_job_dependancy_cmd>', job_dep_cmd)

    if not job_array:
        delete_line(js_path, '#$ -t')
        replace_in_file(js_path, '$(($SGE_TASK_ID - 1))', '0')

    elif selective_submission:
        delete_line(js_path, '#$ -t')


class SimGroup(object):
    """Class to represent a group of related simulations.

    Things we'd want to do with a SimGroup:
    -   Set up directories on stage:
            - write sim input files
            - write jobscripts
            - write json serialisation of self
    -   Write input files for one or more simulation runs on stage/scratch
    -   Submit runs
    -   Add results to one or more runs
    -   Save as a JSON file
    -   Load from a JSON file (for adding results, or runs)


    """

    # Default path formatting options:
    path_options_default = {
        'parallel_sims_join': '_+_',
        'sim_numbers': True,
        'sim_numbers_join': '__',
        'sequence_names': False,
        'sub_dirs': [],
        'run_fmt': '{0}',
        'calcs_path': 'calcs',
        'human_id': '%Y-%m-%d-%H%M_%%r%%r%%r%%r%%r',
    }

    def __init__(self, base_sim_options, run_options, path_options, sim_updates,
                 sequences, human_id, name, sims=None, db_id=None,
                 vis_options=None):
        """Initialise a SimGroup object."""

        self.base_sim_options = base_sim_options
        self.run_options = run_options
        self.path_options = path_options
        self.vis_options = vis_options
        self.sim_updates = sim_updates
        self.sequences = sequences
        self.human_id = human_id
        self.name = name
        self.sims = sims

        self.software_name = run_options['software_name']
        self.stage = run_options.pop('stage')
        self.scratch = run_options.pop('scratch')
        self.archive = run_options.pop('archive')

        self.db_id = db_id

    @property
    def num_sims(self):
        """Return the number of simulations in this group."""
        return len(self.sim_updates)

    @property
    def sim_class(self):
        """Get the Simulation class associated with this group."""
        return SOFTWARE_CLASS_MAP[self.software_name]

    @property
    def sequence_lengths(self):
        # TODO remove this logic from SimGroup

        all_nest = [i.nest_idx for i in self.sequences]
        _, uniq_idx = np.unique(all_nest, return_index=True)
        ret = tuple([self.sequences[i].num_vals for i in uniq_idx])

        return ret

    def get_path_labels(self, sim_idx):
        # TODO remove this logic from SimGroup

        num_elems = self.sequence_lengths
        ind = [sim_idx]
        for j in range(len(num_elems) - 1):
            prod_range = num_elems[::-1][:j + 1]
            prod = np.product(prod_range)
            ind.append(prod * int(ind[-1] / prod))

        return ind[::-1]

    def save_state(self, resource_type, path=None):

        print('Saving SimGroup state -- PENDING')

        json_fn = 'sim_group.json'

        if resource_type == 'stage':
            json_path = self.stage.path.joinpath(json_fn)

        elif resource_type == 'scratch':
            json_path = self.scratch.path.joinpath(json_fn)

        if path:
            json_path = json_path.with_name(json_path.name + '.' + path)

        jsonpickle.set_encoder_options('json', sort_keys=True, indent=4)
        state = jsonpickle.encode(self)

        pick = jsonpickle.pickler.Pickler()

        state = {
            'base_sim_options': pick.flatten(self.base_sim_options),
            'sim_updates': pick.flatten(self.sim_updates),
            'sims': pick.flatten(self.sims),
            'sequences': pick.flatten(self.sequences),
        }

        with json_path.open('w') as sim_group_fp:
            json.dump(state, sim_group_fp, indent=2)

        print('Saving SimGroup state -- DONE')

    def generate_sim_group_sims(self):
        """Generate a list of simulations."""

        print('Generating Simulation objects -- PENDING')
        sims = []
        for sim_upds_idx, sim_upds in enumerate(self.sim_updates):

            print('GEN SIM: {}'.format(sim_upds_idx))
            sim_opt = copy.deepcopy(self.base_sim_options)

            for upd in sim_upds:
                sim_opt = apply_base_update(sim_opt, upd)

            sims.append(self.sim_class(sim_opt))

        print('Generating Simulation objects -- DONE')
        self.sims = sims

    @classmethod
    def load_state(cls, human_id, resource_type):
        """Load from database/JSON file representing a SimGroup object.

        Parameters
        ----------
        resource_type : str
            One of: "stage", "scratch".

        """

        # First get info from database
        sg_params = database.get_sim_group(human_id)

        # Instantiate resource objects:
        stage_name = sg_params['run_opt'].pop('stage')['name']
        scratch_name = sg_params['run_opt'].pop('scratch')['name']
        archive_name = sg_params['run_opt'].pop('archive')['name']

        add_path = sg_params['path_options']['sub_dirs']
        add_path += [sg_params['human_id']]
        stage = Stage(stage_name, add_path)
        scratch = Scratch(scratch_name, add_path)
        archive = Archive(archive_name, add_path)

        # Now want to open the JSON file:
        json_fn = 'sim_group.json'

        if resource_type == 'stage':
            json_path = stage.path.joinpath(json_fn)

        elif resource_type == 'scratch':
            json_path = scratch.path.joinpath(json_fn)

        elif resource_type == 'archive':
            json_path = archive.path.joinpath(json_fn)

        unpick = jsonpickle.unpickler.Unpickler()

        with open(json_path, 'r') as sim_group_fp:
            json_data = json.load(sim_group_fp)

        sg_params['run_opt'].update({
            'stage': stage,
            'scratch': scratch,
            'archive': archive,
        })

        state = {
            'base_sim_options': unpick.restore(json_data['base_sim_options']),
            'sim_updates': unpick.restore(json_data['sim_updates']),
            'sims': unpick.restore(json_data['sims']),
            'sequences': unpick.restore(json_data['sequences']),
            'path_options': sg_params['path_options'],
            'run_options': sg_params['run_opt'],
            'human_id': sg_params['human_id'],
            'name': sg_params['name'],
            'db_id': sg_params['db_id'],
        }

        return cls(**state)

    def make_visualisations(self):
        """Add plots to the simulation directories on stage."""

        # Get the deepest sequence nest index whose values affect the structure
        affect_struct = [i.affects_structure for i in self.sequences]
        deepest_idx = len(affect_struct) - affect_struct[::-1].index(True) - 1

        all_plot_paths = []
        vis_path_final = ''
        for sim_idx, sim in enumerate(self.sims):

            sim_path = self.get_sim_path(sim_idx)
            vis_path = sim_path[:(deepest_idx + 2)]

            if vis_path[-1] != vis_path_final:

                vis_path_final = vis_path[-1]
                vis_dir_path = self.stage.path.joinpath(*vis_path, 'plots')
                vis_dir_path.mkdir(parents=True)
                plot_path = str(vis_dir_path.joinpath('structure.html'))
                all_plot_paths.append(vis_path + ['plots'])
                vis_args = {
                    'save': True,
                    'save_args': {
                        'filename': plot_path
                    },
                    'plot_2d': 'xyz',
                    'group_atoms_by': self.vis_options.get('group_atoms_by'),
                    'style': self.vis_options.get('style'),
                }
                sim.structure.visualise(**vis_args)

        return all_plot_paths

    def write_initial_runs(self):
        """Write input files on stage."""

        if self.sims is not None:
            raise ValueError('Simulations have already been generated for this'
                             ' SimGroup object.')

        # Check this machine is the stage machine
        if self.stage.machine_name != CONFIG['machine_name']:
            raise ValueError('This machine does not have the same name as that '
                             'of the Stage associated with this SimGroup.')

        print('Writing simulation inputs -- PENDING')

        stage_path = self.stage.path
        scratch_path = self.scratch.path
        stage_path.mkdir(parents=True)

        # Copy makesims input file:
        shutil.copy2(CONFIG['option_paths']['makesims'], stage_path)

        sim_params = self.base_sim_options['params'][self.software_name]

        self.sim_class.copy_reference_data(
            sim_params, stage_path, scratch_path)

        self.generate_sim_group_sims()
        self._all_plot_paths = self.make_visualisations()

        for sim_idx, _ in enumerate(self.sims):
            self.sims[sim_idx].runs = [
                {'result': None, }
                for _ in range(len(self.run_options['groups']))]

        # Loop through each requested run group:
        for rg_idx, run_group in enumerate(self.run_options['groups']):

            prt(rg_idx, 'rg_idx')

            # Loop over each simulation in this run group
            sim_paths_stage = []
            sim_paths_scratch = []
            soft_inst = run_group['software_instance']
            for sim_idx in run_group['sim_idx']:

                run_path = self.get_run_path(sim_idx, rg_idx)
                sm_pth_stg = stage_path.joinpath(*run_path)
                sm_pth_sct = scratch_path.joinpath(*run_path)
                sim_paths_stage.append(sm_pth_stg)
                sim_paths_scratch.append(str(sm_pth_sct))

                print('write sim_idx: {}'.format(sim_idx))

                # Write simulation input files:
                self.sims[sim_idx].write_input_files(str(sm_pth_stg))

            # Write supporting files: jobscript, dirlist, options records
            rg_path = ['run_groups', str(rg_idx)]
            rg_path_stage = stage_path.joinpath(*rg_path)
            rg_path_scratch = scratch_path.joinpath(*rg_path)
            rg_path_stage.mkdir(parents=True)

            js_params = {
                'path': str(rg_path_stage),
                'calc_paths': sim_paths_scratch,
                'method': self.software_name,
                'num_cores': run_group['num_cores'],
                'is_sge': self.scratch.sge,
                'job_array': run_group['job_array'],
                'selective_submission': run_group['selective_submission'],
                'scratch_os': self.scratch.os_type,
                'scratch_path': str(rg_path_scratch),
                'parallel_env': soft_inst['parallel_env'],
                'module_load': soft_inst['module_load'],
                'job_name': '{}_{}'.format(self.name, rg_idx),
                'seedname': 'sim',
                'executable': soft_inst['executable'],
            }
            write_jobscript(**js_params)

            if run_group['auto_process'] and self.scratch.sge:

                # Add process jobscript:
                pjs_params = {
                    'path': str(rg_path_stage),
                    'job_name': 'p_' + self.name[2:],
                    'dependency': self.name,
                    'num_calcs': len(sim_paths_scratch),
                    'human_id': self.human_id,
                    'run_group_idx': rg_idx,
                    'job_array': run_group['job_array'],
                    'selective_submission': run_group['selective_submission'],
                }
                write_process_jobscript(**pjs_params)

        print('Writing simulation inputs -- DONE')

    def auto_submit_initial_runs(self):
        """Submit initial runs according the the run group flag `auto_submit`"""

        auto_sub_idx = []
        ask_sub_idx = []
        no_sub_idx = []

        for i_idx, i in enumerate(self.run_options['groups']):

            auto_sub = i.get('auto_submit', 'ask')
            if auto_sub == 'ask':
                ask_sub_idx.append(i_idx)
            elif auto_sub is True:
                auto_sub_idx.append(i_idx)
            else:
                no_sub_idx.append(i_idx)

        if ask_sub_idx:
            for ask_sub in ask_sub_idx:
                if utils.confirm('Submit run group #{}?'.format(ask_sub)):
                    self.submit_run_group(ask_sub)
                else:
                    print('Run group #{} was NOT submitted.'.format(ask_sub))

        if no_sub_idx:
            print('Run groups: {} will not be auto-submitted.'.format(no_sub_idx))

        if auto_sub_idx:
            print('Auto-submitting run groups: {}'.format(auto_sub_idx))
            for rg_idx in auto_sub_idx:
                self.submit_run_group(rg_idx)

    def submit_run_group(self, run_group_idx, run_idx=None):
        """Submit a run groups on Scratch. Can be invoked either on
        Stage (for submitting initial runs) or Scratch. ?? check this.

        Parameters
        ----------
        run_group_idx : int
            Run group index to submit (order of the run group within the sim
            group).
        run_idx : str
            Run indices (order of runs within run group) to submit in the case
            this run group has `selective_submission`. Must be in the SGE "-t"
            option format i.e. "n[-m[:s]]"

        """

        run_group = self.run_options['groups'][run_group_idx]

        # Validation
        if not run_group['selective_submission'] and run_idx is not None:
            msg = ('`run_idx` should only be supplied to `submit_run_group` '
                   'if `selective_submission` if True.')
            raise ValueError(msg)

        if run_group['selective_submission'] and run_idx is None:
            msg = ('`run_idx` must be supplied to `submit_run_group` '
                   'if `selective_submission` if True.')
            raise ValueError(msg)

        conn = ResourceConnection(self.get_machine_resource()[0], self.scratch)
        jobscript_ext = 'sh' if self.scratch.os_type == 'posix' else 'bat'
        jobscript_fn = 'jobscript.{}'.format(jobscript_ext)

        # Get ID in run_group table:
        rg_id = run_group['db_id']

        submit_time = database.get_run_group(rg_id)['submit_time']
        if submit_time:
            msg = 'Run group index: {} was already submitted (at {}).'
            raise ValueError(msg.format(run_group_idx, submit_time))

        sub_msg = ('Submitting run group: {}')
        print(sub_msg.format(run_group_idx))

        rg_path = '/'.join(['run_groups', str(run_group_idx)])

        if self.scratch.sge:
            cmd = ['qsub']
            if run_idx:
                cmd.append('-t {}'.format(run_idx))
            cmd.append('{}'.format(jobscript_fn))

        else:
            cmd = [jobscript_fn]

        # Execute command to get the hostname on Scratch:
        hostname = conn.run_command(['hostname'], block=True)

        # Execute command to submit the jobscript:
        submit_proc = conn.run_command(cmd, cwd=rg_path, block=False)

        # Get approximate time of execution, in a MySQL format:
        dt_fmt = '%Y-%m-%d %H:%M:%S'
        submit_time = datetime.strftime(datetime.now(), dt_fmt)

        # Update the database:
        # Set run state to 3 "in_queue" or 5 "running" (if not SGE)
        run_state_id = 3 if self.scratch.sge else 5
        database.set_run_group_submitted(
            rg_id, hostname, submit_time, run_state_id)

        # Check if submit process has ended:
        if self.scratch.sge:
            msg = 'Submitting to SGE - PENDING'
            msg_done = 'Submitting to SGE - DONE'
        else:
            msg = 'Running sims - PENDING'
            msg_done = 'Running sims - DONE'

        print(msg)

        while True:
            if submit_proc.poll() is not None:
                break
            time.sleep(0.5)

        # At this point either the simulations have been run, or
        # submitted (if SGE).

        with submit_proc.stdout as submit_out:
            submit_stdout = submit_out.read()

        with submit_proc.stderr as submit_err:
            submit_stderr = submit_err.read()

        if not submit_stderr:
            print(msg_done)

        else:
            msg = 'There was a problem with job submission: \n\n\t{}\n'
            raise ValueError(msg.format(submit_stderr))

        if self.scratch.sge:
            # Get the Job-ID from the submit return
            job_id_task_id = submit_stdout.split()[2]
            job_id = int(job_id_task_id.split('.')[0])
            database.set_run_group_sge_jobid(rg_id, job_id)

        else:
            # Set run_states to 6 "pending_process":
            database.set_all_run_states(rg_id, 6)

        if run_group['auto_process']:

            if self.scratch.sge:

                # Execute command to submit the process jobscript:

                cmd = ['qsub']
                if run_idx:
                    cmd.append('-t {}'.format(run_idx))
                cmd.append('process_jobscript.sh')

                submit_process_proc = conn.run_command(
                    cmd, cwd=rg_path, block=False)

                msg = 'Submitting process job to SGE - PENDING {}\r'
                msg_done = 'Submitting process job to SGE - COMPLETE '
                print(msg)
                while True:
                    if submit_process_proc.poll() is not None:
                        break
                    time.sleep(0.1)

                # At this point either the simulations have been run, or
                # submitted (if SGE).

                with submit_process_proc.stdout as submit_out:
                    submit_stdout = submit_out.read()

                with submit_process_proc.stderr as submit_err:
                    submit_stderr = submit_err.read()

                if not submit_stderr:
                    print(msg_done)

                else:
                    msg = ('There was a problem with process job '
                           'submission: \n\n\t{}\n')
                    raise ValueError(msg.format(submit_stderr))

            else:
                process.main(self, run_group_idx)

    def copy_to_scratch(self):
        """Copy group from Stage to Scratch."""

        print('Copying SimGroup to Scratch -- PENDING')
        self.check_is_stage_machine()

        # Copy directory to scratch:
        conn = ResourceConnection(self.stage, self.scratch)
        conn.copy_to_dest(ignore=['plots'])

        # Copy plots directly to archive:
        conn_arch = ResourceConnection(self.stage, self.archive)
        for plot_path in self._all_plot_paths:
            conn_arch.copy_to_dest(subpath=plot_path)

        # Change state of all runs to 2 ("pending_run")
        run_groups = database.get_run_groups(self.db_id)
        for rg_id in [i['id'] for i in run_groups]:
            database.set_all_run_states(rg_id, 2)

        print('Copying SimGroup to Scratch -- DONE')

    def get_machine_resource(self):
        """Returns a list of Resource (Stage and/or Scratch) this machine
        corresponds to.

        Returns
        -------
        list of Resource
            If this machine is the Stage, the Stage resource will be the first
            element in the returned list.

        """

        ret = []
        if CONFIG['machine_name'] == self.stage.machine_name:
            ret.append(self.stage)

        if CONFIG['machine_name'] == self.scratch.machine_name:
            ret.append(self.scratch)

        return ret

    def check_is_stage_machine(self):
        """Check the current machine is the Stage machine for this group."""

        if CONFIG['machine_name'] != self.stage.machine_name:
            raise ValueError('This machine does not have the same ID as that '
                             'of the Stage associated with this SimGroup.')

    def check_is_scratch_machine(self):
        """Check the current machine is the Scratch machine for this group."""

        if CONFIG['machine_name'] != self.scratch.machine_name:
            raise ValueError('This machine does not have the same ID as that '
                             'of the Scratch associated with this SimGroup.')

    def add_runs(self, opt):
        """Add additional run to one or more simulation on scratch."""

        # Need to check current machine is Scratch machine (config.yml)

    def get_sim_path(self, sim_idx):
        """Get the path to a simulation directory.

        Parameters
        ----------
        sim_idx : int
            Index of simulation within the group.

        """
        seq_id = self.sims[sim_idx].options['sequence_id']

        path = []
        nest_idx = seq_id['nest_idx'][0] - 1

        for sid_idx in range(len(seq_id['paths'])):

            add_path = seq_id['paths'][sid_idx]

            if self.path_options['sequence_names']:
                add_path = seq_id['names'][sid_idx] + '_' + add_path

            if seq_id['nest_idx'][sid_idx] == nest_idx:
                path[-1] += self.path_options['parallel_sims_join'] + add_path
            else:
                path.append(add_path)

            nest_idx = seq_id['nest_idx'][sid_idx]

        if self.path_options['sim_numbers']:
            path_labs = self.get_path_labels(sim_idx)
            path = [str(i) + self.path_options['sim_numbers_join'] + j
                    for i, j in zip(path_labs, path)]

        path = [self.path_options['calcs_path']] + path

        return path

    def get_run_path(self, sim_idx, run_idx):
        """Get the path to a simulation run directory."""

        sim_path = self.get_sim_path(sim_idx)
        run_path = sim_path + [self.path_options['run_fmt'].format(run_idx)]

        return run_path

    def check_run_success(self, sim_idx, run_idx):
        """Check a given run of a given sim has succeeded."""

        run_path = self.get_run_path(sim_idx, run_idx)
        run_path_full = self.scratch.path.joinpath(*run_path)
        success = self.sims[sim_idx].check_success(run_path_full)

        return success

    def parse_result(self, sim_idx, run_idx):
        """Parse results for a given run of a given sim and add to the sim
        results attribute."""

        run_path = self.get_run_path(sim_idx, run_idx)
        run_path_full = self.scratch.path.joinpath(*run_path)

        self.sims[sim_idx].parse_result(run_path_full, run_idx)
