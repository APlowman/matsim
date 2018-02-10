"""matsim.__main__.py"""

import sys
import getopt

import yaml

from matsim import OPTSPEC, MAKESIMS_FN, UPDATE_FN, PROCESS_FN, SEQ_DEFN, parse_opt
from matsim.utils import prt


def main(args=sys.argv[1:]):

    try:
        long_args = [
            'make',
            'load=',
            'submit=',
            'process=',
            'force-process=',
            'run-group=',
            'update',
            'no-update'
        ]
        opts, args = getopt.getopt(args, 'ml:s:p:f:r:u', long_args)

    except getopt.GetoptError:
        print('error')
        sys.exit(2)

    if opts[0][0] in ('-m', '--make'):

        from matsim.simulation import makesims

        print('Make.')

        # Parse the makesims.yml options:
        with open(MAKESIMS_FN, 'r') as ms_opts_fp:
            ms_opts_raw = yaml.safe_load(ms_opts_fp)
        ms_opts = parse_opt(ms_opts_raw, OPTSPEC['makesims'])

        makesims.main(ms_opts)

    elif opts[0][0] in ('-l', '--load'):

        from matsim.simulation.simgroup import SimGroup

        print('Load.')

        hid = opts[0][1]

        sim_group = SimGroup.load_state(hid, 'stage')
        prt(sim_group, 'sim_group')

        print('loaded, now saving...')

        sim_group.save_state('stage', path='testy_test')

    elif opts[0][0] in ('-s', '--submit'):

        from matsim.simulation.simgroup import SimGroup

        print('Submit.')

        hid = opts[0][1]
        run_group_idx = None
        for opt, arg in opts:
            if opt in ('-r', '--run-group'):
                run_group_idx = arg
        sim_group = SimGroup.load_state(hid, 'stage')
        sim_group.submit_run_groups([run_group_idx])

    elif opts[0][0] in ('-p', '--process'):

        from matsim.simulation import process
        from matsim.simulation.simgroup import SimGroup

        hid = opts[0][1]
        do_update = True
        run_group_idx = None
        force_process_run_idx = None
        for opt, arg in opts:
            if opt == '--no-update':
                do_update = False
            if opt in ('-r', '--run-group'):
                run_group_idx = arg
            if opt in ('-f', '--force-process'):
                force_process_run_idx = [int(i) for i in arg.split(',')]

        print('Process: {}. rg_idx: {}. do_update: {}. force_process_run_idx: {}'.format(
            hid, run_group_idx, do_update, force_process_run_idx))

        sim_group = SimGroup.load_state(hid, 'scratch')
        proc_args = {
            'sim_group': sim_group,
            'run_group_idx': run_group_idx,
            'do_update': do_update,
            'force_process': force_process_run_idx,
        }
        process.main(**proc_args)

    elif opts[0][0] in ('-u', '--update'):

        from matsim import update
        update.main()
