"""matsim.__main__.py"""

import sys
import argparse

import yaml

from matsim import OPTSPEC, MAKESIMS_FN, UPDATE_FN, PROCESS_FN, SEQ_DEFN, parse_opt
from matsim.utils import prt


def main(args=sys.argv):

    if len(args) == 1 or args[1] == 'make':

        from matsim.simulation import makesims

        # Parse the makesims.yml options:
        with open(MAKESIMS_FN, 'r') as ms_opts_fp:
            ms_opts_raw = yaml.safe_load(ms_opts_fp)
        ms_opts = parse_opt(ms_opts_raw, OPTSPEC['makesims'])

        makesims.main(ms_opts)

    elif args[1] == 'load':

        from matsim.simulation.simgroup import SimGroup

        sim_group = SimGroup.load_state(args[2], 'stage')
        prt(sim_group, 'sim_group')

        print('loaded, now saving...')

        sim_group.save_state('stage', path='testy_test')

    elif args[1] == 'submit':

        from matsim.simulation.simgroup import SimGroup

        sim_group = SimGroup.load_state(args[2], 'stage')

        prt(sim_group.run_options, 'run opt')

        run_group_idx = int(args[3])
        sim_group.submit_run_groups([run_group_idx])

    elif args[1] == 'process':

        from matsim.simulation import process
        from matsim.simulation.simgroup import SimGroup

        sim_group = SimGroup.load_state(args[2], 'scratch')
        proc_args = [sim_group]
        if len(args) == 4:
            run_group_idx = int(args[3])
            proc_args += [run_group_idx]
        process.main(*proc_args)

    elif args[1] == 'update':
        from matsim import update
        update.main()
