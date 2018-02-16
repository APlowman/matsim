"""matsim.simulation.makesims.py"""

import copy

from matsim import database as dbs
from matsim import utils, CONFIG
from matsim.utils import prt, get_recursive
from matsim.resources import Stage, Scratch, Archive, ResourceConnection
from matsim.simulation import apply_base_update
from matsim.simulation.sequence import get_sim_updates
from matsim.simulation.simgroup import SimGroup, SOFTWARE_CLASS_MAP

IMPORT_MAP = [
    {
        'condition': (['params', 'castep', 'continuation'], ['*']),
        'import_files': ['*.check'],
        'import_args': {}
    },
    {
        'condition': (['params', 'castep', 'reuse'], ['*']),
        'import_files': ['*.check'],
        'import_args': {}
    },
    {
        'condition': (['params', 'optados', 'task'], ['*']),
        'import_files': ['*.bands'],
        'import_args': {}
    },
    {
        'condition': (['params', 'optados', 'task'], ['pdos']),
        'import_files': ['*.pdos_bin'],
        'import_args': {}
    },
    {
        'condition': (['params', 'optados', 'broadening'], ['adaptive', 'linear']),
        'import_files': ['*ome_bin'],
        'import_args': {}
    },
    {
        'condition': (['structure', 'source'], ['import']),
        'import_files': [],
        'import_args': {
            'opt_idx': ['structure', 'import', 'opt_idx'],
        }
    },
]


def process_imports(import_options, base_sim_options):
    """Sort out SimGroup imports."""

    # multiple import_options dict ==> one import_data dict
    # Add `import_file_paths` to import_data containing paths of copied files

    prt(base_sim_options, 'base_sim_options')

    matching_imp = {
        'import_files': [],
        'import_args': {},
    }

    # Find if any conditions in the `IMPORT_MAP` are matched in `base_sim_options`:
    for imp_map in IMPORT_MAP:

        # prt(imp_map['condition'][0], "imp_map['condition'][0]")
        # prt(imp_map['condition'][1], "imp_map['condition'][1]")

        try:
            cnd_val = get_recursive(base_sim_options, imp_map['condition'][0])
            # prt(cnd_val, 'cnd_val')

            if cnd_val in imp_map['condition'][1]:

                matching_imp['import_files'].extend(
                    copy.copy(imp_map['import_files'])
                )

                # Resolve any import args
                for key, val in imp_map['import_args'].items():

                    res_val = get_recursive(base_sim_options, val)
                    matching_imp['import_args'].update({
                        key: res_val
                    })

        except AttributeError as _:
            continue

    prt(matching_imp, 'matching_imp')

    # Instantiate import Simulations

    import_data = [
        copy.deepcopy(matching_imp)
        for _ in enumerate(import_options['runs'])
    ]

    # Assume import options within the base_sim_options are lists of the same
    # length as the import options runs list:
    for idx, imp_dat in enumerate(import_data):
        for key, val in imp_dat['import_args'].items():
            import_data[idx]['import_args'][key] = val[idx]

    res_type = import_options['source']

    # First load sim groups
    imp_sim_groups = []
    for imp_idx, imp_hid in enumerate(import_options['sim_groups']):
        imp_sg = SimGroup.load_state(imp_hid, resource_type=res_type)
        imp_sim_groups.append(imp_sg)

    prt(imp_sim_groups, 'imp_sim_groups')

    for imp_idx, imp_run in enumerate(import_options['runs']):

        sim_idx = imp_run['sim_idx']
        run_group_idx = imp_run['run_group_idx']
        sg_idx = imp_run['sim_group_idx']
        sim_group = imp_sim_groups[sg_idx]

        import_data[imp_idx].update({
            'sim': sim_group.sims[sim_idx],
            'run_group_idx': run_group_idx,
            'run_path': sim_group.get_run_path(sim_idx, run_group_idx),
        })

    prt(import_data, 'import_data')

    exit()


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
    run_options = options.pop('run')
    seq_options = options.pop('sequences')
    import_options = options.pop('import')

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
    if import_options:
        process_imports(import_options, base_sim_options)

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
