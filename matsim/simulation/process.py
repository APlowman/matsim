"""matsim.analysis.process"""

from matsim import database
from matsim import update
from matsim.resources import ResourceConnection
from matsim.utils import prt


def main(sim_group, run_group_idx=None):
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
    run_group_idx : int, optional
        If set, only process runs belonging to the this run group. Otherwise,
        process all run groups.

    """

    # Update (SGE) run states in database:
    update.main()

    sim_group.check_is_scratch_machine()
    sg_id = sim_group.db_id

    # prt(sg_id, 'sg_id')

    # Find all runs belonging to this run group in states 6 "pending_process"
    # or state 8 "process_no_errors" (which normally should only be a transient
    # state):
    pending_process = database.get_sim_group_runs(sg_id, [6, 8])
    prt(pending_process, 'pending_process runs')
    prt(run_group_idx, 'run_group_idx')

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

    for pen_run_idx, pen_run in enumerate(pending_process):

        # Get path on scratch of run:
        sim_idx = pen_run['sim_order_id'] - 1
        run_idx = pen_run['run_group_order_id'] - 1
        sim_run_idx.append([sim_idx, run_idx])

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
            # Copy everything to archive apart from calcs directory:
            arch_conn.copy_to_dest(ignore=['calcs'])
            database.set_archive_started(sg_id)

        else:
            # Copy updated sim_group.json
            subpath = ['sim_group.json']
            arch_conn.copy_to_dest(subpath=subpath, file_backup=True)

        archived_ids = []
        for pen_run_idx in no_errs_pen_idx:
            # Copy relevent sim/run/ directories to Archive location
            subpath = sim_group.get_run_path(*sim_run_idx[pen_run_idx])
            arch_conn.copy_to_dest(subpath=subpath)

            archived_ids.append(pending_process[pen_run_idx]['id'])

        # Update states to 10 ("archived")
        database.set_many_run_states(archived_ids, 10)
