"""matsim.simulation.makesims.py"""

import copy

from matsim import database as dbs
from matsim import utils, CONFIG
from matsim.utils import prt
from matsim.resources import Stage, Scratch, Archive, ResourceConnection
from matsim.simulation import apply_base_update
from matsim.simulation.sequence import get_sim_updates
from matsim.simulation.simgroup import SimGroup, SOFTWARE_CLASS_MAP


def validate_run_options(run_options):
    """Validate the run options for making a new SimGroup.

    To check:
      1.) `software` exists on database
      2.) The machine of stage is the current machine
      3.) resource connections exists between:
          a)  stage -> scratch
          b)  stage -> archive
          c)  scratch -> archive

    """
    # Check software exists on database:
    _ = dbs.get_software(run_options['software'])

    stage = run_options['stage']
    scratch = run_options['scratch']
    archive = run_options['archive']

    # Check the machine of the stage is the current machine:
    if stage.machine_name != CONFIG['machine_name']:
        raise ValueError('Specified Stage must belong to the same machine as'
                         'that of this computer, specified in the config file.')

    # Check resource connections exist between resources:
    _ = ResourceConnection.check_exists(stage, scratch)
    _ = ResourceConnection.check_exists(stage, archive)
    _ = ResourceConnection.check_exists(scratch, archive)


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
    vis_options = options.pop('visualise')
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

    validate_run_options(run_options)

    sim_class = SOFTWARE_CLASS_MAP[run_options['software']]
    base_sim = sim_class(base_sim_options)
    sim_updates, sequences = get_sim_updates(seq_options, base_sim)

    sim_group = SimGroup(base_sim_options, run_options, path_options,
                         sim_updates, sequences, human_id, name,
                         vis_options=vis_options)

    print('Making new SimGroup -- DONE')

    return sim_group


def main(options):
    """Main function to generate a simulation group."""

    # Extract out run group options:
    run_group_options = options['run'].pop('groups', [])
    options['run']['groups'] = []
    prt(run_group_options, 'run_group_options')

    # Generate a new SimGroup object
    sim_group = make_new_simgroup(options)
    prt(sim_group.run_options, 'sim_group.run_options')

    # Generate helper files on stage (current machine), and add to database:
    sim_group.initialise()

    # Generate input files for initial run groups, add to database:
    for _, run_group_defn in enumerate(run_group_options):
        sim_group.add_run_group(run_group_defn)

    # Save current state of sim group as JSON file:
    sim_group.save_state()
    dbs.set_sim_group_state_id(sim_group.human_id, 2)

    # Copy to scratch
    copy_scratch = sim_group.run_options['copy_to_scratch']
    do_copy = False
    if copy_scratch == 'ask':
        if utils.confirm('Copy sim group to scratch?'):
            do_copy = True
    elif copy_scratch is True:
        do_copy = True

    if not do_copy:
        print('Sim group was NOT copied to scratch.')

    else:
        sim_group.copy_to_scratch()

        # Submit initial run groups on scratch
        sim_group.submit_initial_runs()

    print('Finished setting up simulation group: {}'.format(sim_group.human_id))
