"""matsim.__main__.py"""

import sys
import getopt
import time
import random

import yaml

from matsim import (OPTSPEC, MAKESIMS_FN, UPDATE_FN, PROCESS_FN, SEQ_DEFN,
                    parse_opt, ADD_RG_FN, database as dbs, utils)
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
            'run-idx=',
            'update',
            'no-update',
            'add-run-group'
        ]
        opts, args = getopt.getopt(args, 'ml:s:p:f:r:i:ua', long_args)

    except getopt.GetoptError:
        print('Error parsing arguments: {}'.format(args))
        sys.exit(2)

    if not opts:
        raise ValueError('Must supply arguments.')

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
        run_idx = None
        for opt, arg in opts:
            if opt in ('-r', '--run-group'):
                run_group_idx = int(arg)
            if opt in ('-i', '--run-idx'):
                run_idx = arg

        if run_group_idx is None:
            raise ValueError('Must specify `run_group_idx` using option: `-r`'
                             ' or `--run-group`.')

        sim_group = SimGroup.load_state(hid, 'stage')
        sim_group.submit_run_group(run_group_idx, run_idx=run_idx)

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
                run_group_idx = int(arg)
            if opt in ('-f', '--force-process'):
                force_process_run_idx = [int(i) for i in arg.split(',')]

        print('Process: {}. rg_idx: {}. do_update: {}. force_process_run_idx: {}'.format(
            hid, run_group_idx, do_update, force_process_run_idx))

        sim_group = SimGroup.load_state(hid, 'scratch')

        # First check if the sim_group is currently being processed.
        # Wait for a random time to split up process invokations which are
        # closely spaced:
        rand_wait_time = random.random() * 60
        print('Waiting for a random wait time of: {}'.format(rand_wait_time))
        time.sleep(rand_wait_time)

        msg = ('Another processing instance is currently working on this '
               'SimGroup. Waiting for {} seconds before re-polling..')
        poll_interval = 2
        while dbs.is_sim_group_processing(sim_group.dbid):
            print(msg.format(poll_interval))
            time.sleep(poll_interval)

        # Reload
        print('Reloading sim group.')
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

    elif opts[0][0] in ('-a', '--add-run-group'):

        from matsim.simulation.simgroup import SimGroup

        print('Add run group.')

        # Parse the add_run_group.yml options:
        with open(ADD_RG_FN, 'r') as add_rg_opts_fp:
            add_rg_opts_raw = yaml.safe_load(add_rg_opts_fp)
        add_rg_opts = parse_opt(add_rg_opts_raw, OPTSPEC['add_run_group'])

        hid = add_rg_opts.pop('id')
        sim_group = SimGroup.load_state(hid)
        sim_group.add_run_group(add_rg_opts)

        # Save current state of sim group as JSON file:
        sim_group.save_state()

        sg_state = dbs.get_sim_group_state_id(sim_group.human_id)
        if sg_state < 2:
            dbs.set_sim_group_state_id(sim_group.human_id, 2)

        elif sg_state > 2:

            run_groups = sim_group.run_options['groups']
            new_rg = run_groups[-1]
            rg_idx = len(run_groups) - 1

            if new_rg['auto_submit'] == 'ask':

                if utils.confirm('Submit run group #{}?'.format(rg_idx)):
                    sim_group.submit_run_group(rg_idx)
                else:
                    print('Run group #{} was NOT submitted.'.format(rg_idx))

            elif new_rg['auto_submit']:
                print('Auto-submitting run group: #{}'.format(rg_idx))
                sim_group.submit_run_group(rg_idx)

            else:
                print('Run group: {} will not be auto-submitted.'.format(rg_idx))
