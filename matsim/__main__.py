"""matsim.__main__.py"""

import sys
import argparse

import yaml

from matsim import OPTSPEC, MAKESIMS_FN, UPDATE_FN, PROCESS_FN, SEQ_DEFN, parse_opt
from matsim.utils import prt


def main(args=sys.argv):

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

    elif args[1] == 'submit':
        from matsim.simulation.simgroup import SimGroup

        sim_group = SimGroup.load_state(args[2], 'stage', SEQ_DEFN)
        run_group_idx = int(args[3])
        sim_group.submit_run_groups([run_group_idx])

    elif args[1] == 'process':
        from matsim.analysis import process

        # Parse the process.yml options:
        with open(PROCESS_FN, 'r') as pr_opts_fp:
            pr_opts_raw = yaml.safe_load(pr_opts_fp)
        pr_opts = parse_opt(pr_opts_raw, OPTSPEC['process'])

        process.main(pr_opts, args[2], SEQ_DEFN, up_opts)

    elif args[1] == 'update':
        from matsim import update
        update.main(up_opts)
