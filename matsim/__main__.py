"""matsim.__main__.py"""

import sys

import yaml

from matsim import OPTSPEC, MAKESIMS_FN, UPDATE_FN, PROCESS_FN, SEQ_DEFN, parse_opt
from matsim.utils import prt


def main(args):

    # Parse the update.yml options:
    with open(UPDATE_FN, 'r') as up_opts_fp:
        up_opts_raw = yaml.safe_load(up_opts_fp)
    up_opts = parse_opt(up_opts_raw, OPTSPEC['update'])

    if len(args) == 1 or args[1] == 'make':
        from matsim.simulation import makesims

        # Parse the makesims.yml options:
        with open(MAKESIMS_FN, 'r') as ms_opts_fp:
            ms_opts_raw = yaml.safe_load(ms_opts_fp)
        ms_opts = parse_opt(ms_opts_raw, OPTSPEC['makesims'])

        makesims.main(ms_opts, ms_opts_raw, SEQ_DEFN)

    elif args[1] == 'load':
        from matsim.simulation.simgroup import SimGroup

        sim_group = SimGroup.load_state(args[2], 'stage', SEQ_DEFN)
        prt(sim_group, 'sim_group')

        print('loaded, now saving...')

        sim_group.save_state('stage', path='testy_test')

    elif args[1] == 'process':
        from matsim.analysis import process

        # Parse the process.yml options:
        with open(PROCESS_FN, 'r') as pr_opts_fp:
            pr_opts_raw = yaml.safe_load(pr_opts_fp)
        pr_opts = parse_opt(pr_opts_raw, OPTSPEC['process'])

        process.main(pr_opts, SEQ_DEFN, up_opts)

    elif args[1] == 'update':
        from matsim import update
        update.main(up_opts)


if __name__ == '__main__':
    main(sys.argv)


# opt_lkup_fn = OPT_FILE_NAMES['lookup']
# opt_def_fn = OPT_FILE_NAMES['defaults']

# OPT_CONFIG_FN = os.path.join(SET_UP_PATH, OPT_FILE_NAMES['config'])
# OPT_RESOURCES_FN = os.path.join(SET_UP_PATH, OPT_FILE_NAMES['resources'])
# OPT_SOFTWARE_FN = os.path.join(SET_UP_PATH, OPT_FILE_NAMES['software'])

# if len(args) == 1 or args[1] == 'make':
#     from matsim.simulation import makesims

#     MS_OPT_FN = os.path.join(SET_UP_PATH, OPT_FILE_NAMES['makesims'])

#     with open(OPT_CONFIG_FN, 'r') as config_fs:
#         OPT_CONFIG = yaml.safe_load(config_fs)

#     with open(OPT_RESOURCES_FN, 'r') as resources_fs:
#         OPT_RESOURCES = yaml.safe_load(resources_fs)

#     with open(OPT_SOFTWARE_FN, 'r') as software_fs:
#         OPT_SOFTWARE = yaml.safe_load(software_fs)

#     makesims.main2(MS_OPT_FN, OPT_SPEC_FN, OPT_CONFIG,
#                    OPT_RESOURCES, OPT_SOFTWARE)

# elif args[1] == 'process':

#     from matsim.analysis import process

#     PR_OPT_FN = os.path.join(SET_UP_PATH, OPT_FILE_NAMES['process'])

#     with open(OPT_CONFIG_FN, 'r') as config_fs:
#         OPT_CONFIG = yaml.safe_load(config_fs)

#     process.main(PR_OPT_FN, OPT_SPEC_FN, OPT_CONFIG)

# elif args[1] == 'load':

#     from matsim.simgroup import SimGroup
#     sim_group = SimGroup.load_state(args[2], OPT_SPEC_FN)
#     prt(sim_group, 'sim_group')

#     with open(OPT_CONFIG_FN, 'r') as config_fs:
#         OPT_CONFIG = yaml.safe_load(config_fs)

#     # a = sim_group.sims[0].structure
#     # prt(a, 'a')

#     # cs = sim_group.sims[0].structure.crystal_structures[0]
#     # prt(cs, 'cs')
#     sim_group.save_state()
#     SimGroup.set_machine_id(OPT_CONFIG['machine_id'])
#     # sim_group.write_initial_runs()

# # elif args[1] == 'process':
# #     if len(args) != 3:
# #         raise ValueError('Specify SID to process.')
# #     from matsim.analysis import process
# #     opt_fn = OPT_FILE_NAMES['process']
# #     ps_opt = options_parser.validate_ps_opt(opt_fn, opt_lkup_fn)
# #     process.main(ps_opt, args[2])

# # elif args[1] == 'submit_process':
# #     if len(args) != 3:
# #         raise ValueError('Specify SID to process.')
# #     from matsim.analysis import submit_process
# #     submit_process.main(args[2])

# # elif args[1] == 'harvest':
# #     from matsim.analysis import harvest
# #     opt_fn = OPT_FILE_NAMES['harvest']
# #     hv_opt = options_parser.validate_hv_opt(opt_fn, opt_lkup_fn, opt_def_fn)
# #     harvest.main(hv_opt)

# # elif args[1] == 'plot':
# #     from matsim.analysis import makeplots
# #     opt_fn = OPT_FILE_NAMES['makeplots']
# #     mp_opt = options_parser.validate_mp_opt(opt_fn, opt_lkup_fn, opt_def_fn)
# #     makeplots.main(mp_opt)

# # elif args[1] == 'series_helper':
# #     from matsim.simulation import serieshelper
# #     opt_fn = OPT_FILE_NAMES['serieshelper']
# #     sh_opt = options_parser.validate_sh_opt(opt_fn, opt_lkup_fn, opt_def_fn)
# #     serieshelper.main(sh_opt)
# #     # series_helpers.refine_gamma_surface(args[2])

# else:
#     print('Invalid option.')
