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
                    database as dbs, OPTSPEC)
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
                 sequences, human_id, name, sims=None, dbid=None,
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

        self.software = run_options['software']
        self.stage = run_options.pop('stage')
        self.scratch = run_options.pop('scratch')
        self.archive = run_options.pop('archive')

        self.dbid = dbid

    @property
    def num_sims(self):
        """Return the number of simulations in this group."""
        return len(self.sim_updates)

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

    @classmethod
    def get_default_state_location(cls, human_id, stage, scratch, archive):
        """Return either stage or scratch str."""

        state_id = dbs.get_sim_group_state_id(human_id)

        if state_id in [1, 2]:
            resource_type = 'stage'
            if CONFIG['machine_name'] != stage.machine_name:
                raise ValueError('This is not the stage machine.')

        elif state_id in [3, 4, 5]:
            resource_type = 'scratch'
            if CONFIG['machine_name'] != scratch.machine_name:
                raise ValueError('This is not the scratch machine.')

        return resource_type

    def save_state(self, resource_type=None, path=None):

        print('Saving SimGroup state -- PENDING')

        json_fn = 'sim_group.json'

        if not resource_type:
            resource_type = SimGroup.get_default_state_location(
                self.human_id, self.stage, self.scratch, self.archive
            )

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

            sim = SOFTWARE_CLASS_MAP[self.software](sim_opt)
            sims.append(sim)

        print('Generating Simulation objects -- DONE')
        self.sims = sims

    @classmethod
    def load_state(cls, human_id, resource_type=None):
        """Load from database/JSON file representing a SimGroup object.

        Parameters
        ----------
        resource_type : str
            One of: "stage", "scratch", "archive", determines where to load the
            JSON file from. By default, the sim_group_state will be checked. If
            state is 1 "on_stage_initial" or 2 "on_stage_run_groups", 
            resource_type is set to "stage"; if state is 3, 4, or 5,
            resource_type is set to "scratch".

        """

        # First get info from database
        sg_params = dbs.get_sim_group_by_human_id(human_id)

        prt(sg_params, 'sg_params')

        # Get sim IDs from database as well:
        sims_db = dbs.get_sim_group_sims(sg_params['id'])

        # Instantiate resource objects:
        stage_name = sg_params['run_opt'].pop('stage')['name']
        scratch_name = sg_params['run_opt'].pop('scratch')['name']
        archive_name = sg_params['run_opt'].pop('archive')['name']

        add_path = sg_params['path_options']['sub_dirs']
        add_path += [sg_params['human_id']]
        stage = Stage(stage_name, add_path)
        scratch = Scratch(scratch_name, add_path)
        archive = Archive(archive_name, add_path)

        if not resource_type:
            resource_type = SimGroup.get_default_state_location(
                human_id, stage, scratch, archive
            )

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
            'dbid': sg_params['id'],
        }

        # Reset sim IDs from database, rather than relying on JSON file.
        for sim_idx, _ in enumerate(state['sims']):
            state['sims'][sim_idx].dbid = sims_db[sim_idx]['id']

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

    def initialise(self):
        """Write helper files on stage and add sim group to database."""

        # Check if this sim group is in database
        prt(dbs.get_sim_group_state_id(self.human_id),
            'dbs.get_sim_group_state_id(self.human_id)')
        if dbs.get_sim_group_state_id(self.human_id) is not False:
            raise ValueError('Simulations have already been generated for this'
                             ' SimGroup object.')

        # Check this machine is the stage machine
        if self.stage.machine_name != CONFIG['machine_name']:
            raise ValueError('This machine does not have the same name as that '
                             'of the Stage associated with this SimGroup.')

        print('Writing initial files -- PENDING')

        stage_path = self.stage.path
        scratch_path = self.scratch.path
        stage_path.mkdir(parents=True)

        # Copy makesims input file:
        shutil.copy2(CONFIG['option_paths']['makesims'], stage_path)

        sim_params = self.base_sim_options['params'][self.software]

        SOFTWARE_CLASS_MAP[self.software].copy_reference_data(
            sim_params, stage_path, scratch_path)

        self.generate_sim_group_sims()
        self._all_plot_paths = self.make_visualisations()

        print('Writing initial files -- DONE')

        print('Adding SimGroup to the database -- PENDING')
        dbs.add_sim_group(self)
        print('Adding SimGroup to the database -- DONE')
        prt(self.dbid, 'self.dbid')

        self.save_state('stage')

    def _validate_run_group_defn(self, rg_defn):
        """Validate a run group options dict."""

        msg_soft_1 = (
            'Software instance "{}" is not allowed on scratch "{}"'
        )
        msg_soft_2 = (
            'Invalid software instances "{}", since this does not use '
            'software "{}"'
        )
        msg_invalid_cores = (
            '{} core(s) is not supported on the specified software instance.'
        )
        sim_idx_msg_1 = (
            'Run group `sim_idx` must be either "all" or a list of '
            'integers that index the sims.'
        )
        sim_idx_msg_2 = (
            'Run group `sim_idx` cannot contain repeated indices.'
        )
        job_arr_msg = (
            '`job_array` cannot be `True` for a Scratch which is not SGE'
        )
        sel_sub_msg_1 = (
            '`selective_submission` cannot be `True` if `job_array` is `False`'
        )
        sel_sub_msg_2 = (
            '`selective_submission` cannot be `True` for a Scratch which is '
            'not SGE, and if `job_array` if `False`.'
        )
        sel_sub_msg_3 = (
            '`auto_submit` cannot be `True` or `ask` if `selective_submission`'
            ' is `True`.'
        )

        soft_inst_name = rg_defn['software_instance']
        soft_inst = dbs.get_software_instance_by_name(soft_inst_name)

        # Check software instance has correct software name:
        rg_software_name = soft_inst['software_name']
        if rg_software_name != self.software:
            raise ValueError(msg_soft_2.format(soft_inst_name, self.software))

        # Check software instance is allowed on this scratch:
        ok_scratch_ids = dbs.get_software_instance_ok_scratch(soft_inst_name)
        if self.scratch.scratch_id not in ok_scratch_ids:
            raise ValueError(msg_soft_1.format(
                soft_inst_name, self.scratch.name))

        # Replace the name of the software instance with the software instance
        # dict from the database:
        rg_defn['software_instance'] = soft_inst

        # Validate `num_cores`:
        ncores = rg_defn['num_cores']
        cores_good = soft_inst['min_cores'] <= ncores <= soft_inst['max_cores']
        if not cores_good:
            raise ValueError(msg_invalid_cores.format(ncores))

        # Validate `sim_idx`
        rg_sim_idx = rg_defn['sim_idx']

        if rg_sim_idx == 'all':
            rg_defn['sim_idx'] = list(range(self.num_sims))

        elif isinstance(rg_sim_idx, list):

            if min(rg_sim_idx) < 0 or max(rg_sim_idx) > (self.num_sims - 1):
                raise ValueError(sim_idx_msg_1)

            if len(set(rg_sim_idx)) != len(rg_sim_idx):
                raise ValueError(sim_idx_msg_2)

        else:
            raise ValueError(sim_idx_msg_1)

        # Validate `job_array` and `selective_submission`:
        if self.scratch.sge:

            job_array = rg_defn.get('job_array', len(rg_defn['sim_idx']) > 1)
            sel_sub = rg_defn.get('selective_submission', False)
            if sel_sub is True and not job_array:
                raise ValueError(sel_sub_msg_1)

        else:
            if rg_defn.get('job_array') is True:
                raise ValueError(job_arr_msg)
            if rg_defn.get('selective_submission') is True:
                raise ValueError(sel_sub_msg_2)

            job_array = False
            sel_sub = False

        if sel_sub and rg_defn['auto_submit'] in [True, 'ask']:
            raise NotImplementedError(sel_sub_msg_3)

        rg_defn['job_array'] = job_array
        rg_defn['selective_submission'] = sel_sub

    def submit_initial_runs(self):
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

    def add_run_group(self, run_group_defn):
        """Add a run group to the sim group.

        Parameters
        ----------
        run_group_defn : dict with keys:
            sim_idx : list of int or str
                Indexes simulations to include in this run group.
            software_instance : str
                Name of software instance.
            num_cores : int
            auto_submit: bool or str
                If True, submit this run group immediately. If False, do not submit.
                If "ask", provide a prompt to ask whether to submit.
            auto_process : bool
                If True, automatically process this run group after jobs complete.
            job_array : bool
                If True, submit runs as a job array.
            selective_submission : bool
                If True, allow submitting the runs in this run group selectively.

        """

        self._validate_run_group_defn(run_group_defn)

        # Find the new order for the new run group (within the sim group):
        order_in_sim_group = dbs.get_count_sim_group_run_groups(self.dbid)

        # Get resource type
        sg_state_id = dbs.get_sim_group_state_id(self.human_id)
        prt(sg_state_id, 'sg_state_id')
        if sg_state_id in [1, 2]:
            resource_type = 'stage'
        elif sg_state_id in [3, 4, 5]:
            resource_type = 'scratch'

        stage_path = self.stage.path
        scratch_path = self.scratch.path
        sim_paths_stage = []
        sim_paths_scratch = []
        sim_dbids = []
        run_order_in_sims = []

        # Add to the simulation runs list:
        for idx, sim_idx in enumerate(run_group_defn['sim_idx']):

            sim = self.sims[sim_idx]
            run_order_in_sim = len(sim.runs)
            run_order_in_sims.append(run_order_in_sim)
            sim_run_defn = {
                'run_group_order_in_sim_group': order_in_sim_group,
                'run_order_in_run_group': idx,
                'run_order_in_sim': run_order_in_sim,
                'run_params': run_group_defn['run_params'],
                'dbid': run_group_defn['dbid'],
                'result': None
            }
            sim.runs.append(sim_run_defn)

            run_path = self.get_run_path(sim_idx, order_in_sim_group)
            sm_pth_stg = stage_path.joinpath(*run_path)
            sm_pth_sct = scratch_path.joinpath(*run_path)
            sim_paths_stage.append(sm_pth_stg)
            sim_paths_scratch.append(str(sm_pth_sct))
            sim_dbids.append(sim.dbid)

            if resource_type == 'stage':
                input_path = sm_pth_stg
            elif resource_type == 'scratch':
                input_path = sm_pth_sct

            # Write simulation input files:
            self.sims[sim_idx].write_input_files(
                run_order_in_sim, str(input_path))

        # Write supporting files: jobscript, dirlist, options records
        rg_path = ['run_groups', str(order_in_sim_group)]
        rg_path_stage = stage_path.joinpath(*rg_path)
        rg_path_scratch = scratch_path.joinpath(*rg_path)

        if resource_type == 'stage':
            rg_path_resource = rg_path_stage
        elif resource_type == 'scratch':
            rg_path_resource = rg_path_scratch

        rg_path_resource.mkdir(parents=True)

        soft_inst = run_group_defn['software_instance']
        js_params = {
            'path': str(rg_path_resource),
            'calc_paths': sim_paths_scratch,
            'method': self.software,
            'num_cores': run_group_defn['num_cores'],
            'is_sge': self.scratch.sge,
            'job_array': run_group_defn['job_array'],
            'selective_submission': run_group_defn['selective_submission'],
            'scratch_os': self.scratch.os_type,
            'scratch_path': str(rg_path_scratch),
            'parallel_env': soft_inst['parallel_env'],
            'module_load': soft_inst['module_load'],
            'job_name': '{}_{}'.format(self.name, order_in_sim_group),
            'seedname': 'sim',
            'executable': soft_inst['executable'],
        }
        write_jobscript(**js_params)

        if run_group_defn['auto_process'] and self.scratch.sge:

            # Add process jobscript:
            pjs_params = {
                'path': str(rg_path_resource),
                'job_name': 'p_' + self.name[2:],
                'dependency': self.name,
                'num_calcs': len(sim_paths_scratch),
                'human_id': self.human_id,
                'run_group_idx': order_in_sim_group,
                'job_array': run_group_defn['job_array'],
                'selective_submission': run_group_defn['selective_submission'],
            }
            write_process_jobscript(**pjs_params)

        # Add run group to the database:
        if resource_type == 'stage':
            run_state = 1
        elif resource_type == 'scratch':
            run_state = 2

        rg_db = dbs.add_run_group(self.dbid, run_group_defn, sim_dbids,
                                  run_order_in_sims, run_state,
                                  self.scratch.sge)
        run_group_defn.update({
            'dbid': rg_db[0]
        })

        # Add to the run groups list:
        self.run_options['groups'].append(run_group_defn)

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

        # Need to identify if we're on Stage or Scratch
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

        state_id = dbs.get_sim_group_state_id(self.human_id)
        if state_id in [1, 2]:
            src_resource = self.stage
        elif state_id in [3, 4, 5]:
            if self.stage.machine_name == CONFIG['machine_name']:
                src_resource = self.stage
            else:
                src_resource = self.scratch

        conn = ResourceConnection(src_resource, self.scratch)

        if dbs.get_sim_group_state_id(self.human_id) < 4:
            dbs.set_sim_group_state_id(self.human_id, 4)

        jobscript_ext = 'sh' if self.scratch.os_type == 'posix' else 'bat'
        jobscript_fn = 'jobscript.{}'.format(jobscript_ext)

        # Get ID in run_group table:
        rg_id = run_group['dbid']

        submit_time = dbs.get_run_group(rg_id)['submit_time']
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
        dbs.set_run_group_submitted(
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
            dbs.set_run_group_sge_jobid(rg_id, job_id)

        else:
            # Set run_states to 6 "pending_process":
            dbs.set_all_run_states(rg_id, 6)

        if run_group['auto_process']:

            if self.scratch.sge:

                # Execute command to submit the process jobscript:

                cmd = ['qsub']
                if run_idx:
                    cmd.append('-t {}'.format(run_idx))
                cmd.append('process_jobscript.sh')

                submit_process_proc = conn.run_command(
                    cmd, cwd=rg_path, block=False)

                msg = 'Submitting process job to SGE - PENDING'
                msg_done = 'Submitting process job to SGE - COMPLETE'
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
        run_groups = dbs.get_sim_group_run_groups(self.dbid)
        for rg_id in [i['id'] for i in run_groups]:
            dbs.set_all_run_states(rg_id, 2)

        dbs.set_sim_group_state_id(self.human_id, 3)

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

    def get_run_path(self, sim_idx, run_group_idx):
        """Get the path to a simulation run directory."""

        sim_path = self.get_sim_path(sim_idx)
        run_path = sim_path
        run_path += [self.path_options['run_fmt'].format(run_group_idx)]

        return run_path

    def check_run_success(self, sim_idx, run_group_idx):
        """Check a given run of a given sim has succeeded."""

        run_path = self.get_run_path(sim_idx, run_group_idx)
        run_path_full = self.scratch.path.joinpath(*run_path)
        success = self.sims[sim_idx].check_success(run_path_full)

        return success

    def parse_result(self, sim_idx, run_group_idx):
        """Parse results for a given run of a given sim and add to the sim
        results attribute."""

        run_path = self.get_run_path(sim_idx, run_group_idx)
        run_path_full = self.scratch.path.joinpath(*run_path)

        sim = self.sims[sim_idx]

        run_order_in_sim = utils.dict_from_list(
            sim.runs, {'run_group_order_in_sim_group': run_group_idx}
        )['run_order_in_sim']

        # prt(run_path, 'run_path')
        # prt(run_path_full, 'run_path_full')
        # prt(run_order_in_sim, 'run_order_in_sim')

        run_res = sim.runs[run_order_in_sim]['result']
        if run_res:
            msg = 'Result has already been parsed for run_idx {}'
            raise ValueError(msg.format(run_order_in_sim))

        sim.parse_result(run_path_full, run_order_in_sim)
