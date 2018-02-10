"""matsim.analysis.process"""

from matsim import database
from matsim import update
from matsim.resources import ResourceConnection
from matsim.utils import prt


def main(sim_group, run_group_idx=None, do_update=True, force_process=None):
    """
    Process a given SimGroup:
    -   Run update to update run states in the database
    -   Find all runs in state 6 ("pending_process")
    -   Invoke check success method on each run
    -   If check success True, parse results and add to result attribute in
        SimGroup -> sim -> runs.
    -   Overwrite JSON file to Scratch (make backup of previous on Scratch)
    -   Copy new JSON file Archive (make backup of previous on Archive)

    Parameters
    ----------
    run_group_idx : int, optional, zero-indexed.
        If set, only process runs belonging to the this run group. Otherwise,
        process all run groups.
    force_process, list of int, optional, zero-indexed.
        List of run indices within specified run group to initially set to 
        `pending_process` regardless of checking `qstat`.

    """

    # prt(sim_group, 'sim_group')
    sg_id = sim_group.db_id
    # prt(sg_id, 'sg_id')

    if do_update:
        # Update (SGE) run states in database:
        update.main()

    sim_group.check_is_scratch_machine()

    if run_group_idx is None and force_process:
        raise ValueError('Must specify `run_group_idx`, if `force_process` is '
                         'set.')

    if force_process:
        rg_id = database.get_run_groups(sg_id)[run_group_idx]['id']
        prt(rg_id, 'rg_id')

        rg_runs = database.get_run_group_runs(rg_id)
        prt(rg_runs, 'rg_runs')

        force_run_ids = [rg_runs[i]['id'] for i in force_process]
        prt(force_run_ids, 'force_run_ids')

        database.set_many_run_states(force_run_ids, 6)

    # Find all runs belonging to this run group in states 6 "pending_process"
    # or state 8 "process_no_errors" (which normally should only be a transient
    # state):
    pending_process = database.get_sim_group_runs(sg_id, [6, 8])

    print('Found {} run(s) to process.'.format(len(pending_process)))
    # prt(pending_process, 'pending_process runs')
    # prt(run_group_idx, 'run_group_idx')

    if run_group_idx:

        # Only keep runs belonging to given run group
        for pen_run_idx in range(len(pending_process)):

            pen_run = pending_process[pen_run_idx]
            if pen_run['run_group_order_id'] == run_group_idx + 1:
                pending_process = [pen_run]
                break
            else:
                pending_process = []

    # Set state to 7 ("processing") for these runs
    run_ids = [i['id'] for i in pending_process]
    database.set_many_run_states(run_ids, 7)

    no_errs_pen_idx = []
    errs_pen_idx = []
    sim_run_idx = []
    run_group_sge_job_id = []
    run_order_id = []
    run_group_is_job_arr = []
    for pen_run_idx, pen_run in enumerate(pending_process):

        # Get path on scratch of run:
        sim_idx = pen_run['sim_order_id'] - 1
        run_idx = pen_run['run_group_order_id'] - 1
        sim_run_idx.append([sim_idx, run_idx])
        run_group_sge_job_id.append(pen_run['sge_job_id'])
        run_order_id.append(pen_run['order_id'])
        run_group_is_job_arr.append(bool(pen_run['sge_job_array']))
        run_success = sim_group.check_run_success(sim_idx, run_idx)

        if run_success:
            no_errs_pen_idx.append(pen_run_idx)
        else:
            errs_pen_idx.append(pen_run_idx)

    # Update states to 9 ("process_errors")
    err_ids = [pending_process[i]['id'] for i in errs_pen_idx]
    database.set_many_run_states(err_ids, 9)

    # Parse full output and add to sim.results[run_idx]
    for pen_run_idx in no_errs_pen_idx:
        sim_group.parse_result(*sim_run_idx[pen_run_idx])

    # Update states to 8 ("process_no_errors")
    no_err_ids = [pending_process[i]['id'] for i in no_errs_pen_idx]
    database.set_many_run_states(no_err_ids, 8)

    if no_errs_pen_idx:

        # Copy new sim_group.json to Archive location
        arch_conn = ResourceConnection(sim_group.scratch, sim_group.archive)

        # Overwrite sim_group.json with new results:
        sim_group.save_state('scratch')

        if not database.check_archive_started(sg_id):
            # Copy everything to archive apart from calcs directory and SGE
            # log files:
            arch_conn.copy_to_dest(ignore=['calcs', 'sge'])
            database.set_archive_started(sg_id)

        else:
            # Copy updated sim_group.json (and backup old one on archive)
            subpath = ['sim_group.json']
            arch_conn.copy_to_dest(subpath=subpath, file_backup=True)

        archived_ids = []
        for pen_run_idx in no_errs_pen_idx:
            # Copy relevent sim/run/ directories to Archive location
            subpath = sim_group.get_run_path(*sim_run_idx[pen_run_idx])
            arch_conn.copy_to_dest(subpath=subpath)

            archived_ids.append(pending_process[pen_run_idx]['id'])

            if sim_group.scratch.sge:
                # Copy SGE log files for these runs:

                # Need to form log stdout/err filenames:
                # Need SGE job_id from table `run_group_sge`
                sge_job_id_str = str(run_group_sge_job_id[pen_run_idx])
                prt(sge_job_id_str, 'sge_job_id_str')

                rg_idx_str = str(sim_run_idx[pen_run_idx][1])
                prt(rg_idx_str, 'rg_idx_str')

                stdoe_fn = sim_group.name + '_' + rg_idx_str + '.{}' + sge_job_id_str

                if run_group_is_job_arr[pen_run_idx]:
                    r_idx_str = str(run_order_id[pen_run_idx])
                    prt(r_idx_str, 'r_idx_str')

                    stdoe_fn += '.' + r_idx_str

                stdout_fn = stdoe_fn.format('o')
                stderr_fn = stdoe_fn.format('e')

                prt(stdoe_fn, 'stdoe_fn')
                prt(stdout_fn, 'stdout_fn')
                prt(stderr_fn, 'stderr_fn')

                stdeo_path = ['run_groups', rg_idx_str, 'sge']
                prt(stdeo_path, 'stdeo_path')

                stdout_path = stdeo_path + [stdout_fn]
                stderr_path = stdeo_path + [stderr_fn]

                prt(stdeo_path, 'stdeo_path')
                prt(stdout_path, 'stdout_path')
                prt(stderr_path, 'stderr_path')
                arch_conn.copy_to_dest(subpath=stdout_path)
                arch_conn.copy_to_dest(subpath=stderr_path)

        # Update states to 10 ("archived")
        database.set_many_run_states(archived_ids, 10)
