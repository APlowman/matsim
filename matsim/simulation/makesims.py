"""matsim.simulation.makesims.py"""

import copy

from matsim import database as dbs
from matsim import utils
from matsim.utils import prt
from matsim.resources import Stage, Scratch, Archive
from matsim.simulation import apply_base_update
from matsim.simulation.sequence import get_sim_updates
from matsim.simulation.simgroup import SimGroup, SOFTWARE_CLASS_MAP


def validate_run_group_sim_idx(run_options, num_sims):
    """Validate the user input for the `sim_idx` key of each run group."""

    sim_idx_msg = ('Run group `sim_idx` must be either "all" or a list of '
                   'integers that index the sims.')

    for rg_idx in range(len(run_options['groups'])):

        run_group = run_options['groups'][rg_idx]
        rg_sim_idx = run_group['sim_idx']

        if rg_sim_idx == 'all':
            run_group['sim_idx'] = list(range(num_sims))

        elif isinstance(rg_sim_idx, list):

            if min(rg_sim_idx) < 0 or max(rg_sim_idx) > (num_sims - 1):
                raise ValueError(sim_idx_msg)

        else:
            raise ValueError(sim_idx_msg)


def validate_run_options(run_options):
    """Validate the run options for making a new SimGroup.

    The only missing validation check is checking the `sim_idx` of each run
    group, since we need to know the number of sims in the group to check this,
    which requires instantiating the SimSequenceGroup, which in turn requires
    knowing the software name, which is validated in this function.

    """

    # TODO: raise exception if sge.jobarray False and len(sim_idx) > 1

    scratch = run_options['scratch']
    software_name = None

    msg_soft_ins = 'Software instance "{}" is not allowed on scratch "{}"'
    msg_soft_name = 'All run groups must use the same software (name)!'
    msg_invalid_cores = ('{} core(s) is not supported on the specified '
                         'software instance.')

    for run_group in run_options['groups']:

        # Check software instance for this run group is allowed on this scratch
        soft_inst_name = run_group['software_instance']
        soft_inst = dbs.get_software_instance_by_name(soft_inst_name)
        run_group['software_instance'] = soft_inst
        rg_software_name = soft_inst['software_name']

        ok_scratch_ids = dbs.get_software_instance_ok_scratch(soft_inst_name)
        if scratch.scratch_id not in ok_scratch_ids:
            raise ValueError(msg_soft_ins.format(soft_inst_name, scratch.name))

        # Set software name for this sim group:
        if software_name is None:
            software_name = rg_software_name

        elif software_name != rg_software_name:
            raise ValueError(msg_soft_name)

        # Check valid num cores specified:
        ncores = run_group['num_cores']
        cores_good = soft_inst['min_cores'] <= ncores <= soft_inst['max_cores']
        if not cores_good:
            raise ValueError(msg_invalid_cores.format(ncores))

    run_options['software_name'] = software_name
    return software_name


def add_sim_group_to_db(sim_group):
    """Add the new sim group to the database."""

    print('Adding SimGroup to the database -- PENDING')
    db_ret = dbs.add_sim_group(sim_group)

    # Add sim grounp database ID:
    sim_group.db_id = db_ret['sim_group_id']

    # Add run grounp database IDs:
    run_groups = sim_group.run_options['groups']
    for rg_idx, _ in enumerate(run_groups):
        run_group = run_groups[rg_idx]
        rg_ids = db_ret['run_group_ids'][rg_idx]
        run_group.update({
            'db_id': rg_ids[0],
            'sge_db_id': rg_ids[1]
        })

    print('Adding SimGroup to the database -- DONE')


def make_new_simgroup(options):
    """Used to instantiate SimGroup when generating a new simulation
    group, as opposed to when loading an existing SimGroup from a JSON file.

    Parameters
    ----------
    options : dict
        Parsed and validated options dict.

    Returns
    -------
    sim_group : SimGroup

    """

    print('Making new SimGroup -- PENDING')

    path_options = options.pop('path_options')
    # TODO: deal with importing structure (and group):
    # import_options = options.pop('import_options')
    run_options = options.pop('run')
    seq_options = options.pop('sequences')

    # All remaning options are base simulation options:
    base_sim_options = options

    path_options = {**SimGroup.path_options_default, **path_options}
    sub_dirs = [utils.parse_times(i)[0] for i in path_options['sub_dirs']]
    path_options.update({'sub_dirs': sub_dirs})

    hid, hid_num = utils.parse_times(path_options['human_id'])
    human_id = hid
    name = 'j_' + hid_num
    add_path = path_options['sub_dirs'] + [human_id]

    run_options['stage'] = Stage(run_options['stage'], add_path)
    run_options['scratch'] = Scratch(run_options['scratch'], add_path)
    run_options['archive'] = Archive(run_options['archive'], add_path)

    software_name = validate_run_options(run_options)
    sim_class = SOFTWARE_CLASS_MAP[software_name]
    base_sim = sim_class(base_sim_options)
    sim_updates, sequences = get_sim_updates(seq_options, base_sim)
    validate_run_group_sim_idx(run_options, len(sim_updates))

    sim_group = SimGroup(base_sim_options, run_options, path_options,
                         sim_updates, sequences, human_id, name)

    print('Making new SimGroup -- DONE')

    return sim_group


def main(options):
    """Main function to generate a simulation group."""

    # Generate a new SimGroup object
    sim_group = make_new_simgroup(options)

    # Generate input files on stage (current machine):
    sim_group.write_initial_runs()

    # Add records to database:
    add_sim_group_to_db(sim_group)

    # Save current state of sim group as JSON file:
    sim_group.save_state('stage')

    # Copy to scratch
    do_copy = sim_group.copy_to_scratch()

    if do_copy:
        # Submit initial run groups on scratch
        sim_group.auto_submit_initial_runs()

    print('Finished setting up simulation group: {}'.format(sim_group.human_id))
