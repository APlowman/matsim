"""`matsim.database_atm.py`"""

import pathlib
import json
import copy

from matsim import database as dbs, CONFIG
from matsim.utils import prt


def get_atm_sim(path, split_path=True, user_cred=None):

    user_cred = user_cred or CONFIG['user']
    user_id = dbs.get_user_id(user_cred)

    if split_path:
        path = json.dumps(pathlib.Path(path).parts)

    sql = (
        'select * from atm_sim '
        'where registered_path = %s '
        'and user_account_id = %s'
    )

    atm_sims = dbs.exec_select(sql, (path, user_id))

    return atm_sims


def get_atm_sim_by_sim_id(sim_id, user_cred=None):

    user_cred = user_cred or CONFIG['user']
    user_id = dbs.get_user_id(user_cred)

    sql = (
        'select ats.* from atm_sim ats '
        'inner join user_account ua on ua.id = ats.user_account_id '
        'where ats.sim_id = %s '
        'and ua.id = %s'
    )
    atm_sim = dbs.exec_select(sql, (sim_id, user_id))

    return atm_sim


def get_atm_sim_castep_by_atm_sim_id(atm_sim_id, user_cred=None):
    """Get an atm_sim_castep row from the atm_sim_id."""

    user_cred = user_cred or CONFIG['user']
    user_id = dbs.get_user_id(user_cred)

    sql = (
        'select atsc.* from atm_sim_castep atsc '
        'inner join atm_sim ats on ats.id = atsc.atm_sim_id '
        'inner join user_account ua on ua.id = ats.user_account_id '
        'where atsc.atm_sim_id = %s '
        'and ua.id = %s'
    )
    atm_sim_cst = dbs.exec_select(sql, (atm_sim_id, user_id))

    return atm_sim_cst


def validate_new_atm_sim(params):
    """Validate new atomistic sim before adding to database."""

    not_none = ['repeats', 'supercell_type', 'num_atoms', 'path']
    for key in not_none:
        if not params.get(key):
            raise ValueError('Key "{}" must not be `None`'.format(key))


def add_atm_sim(params, sim_id=None, run_id=None, user_cred=None):
    """Add atomistic simulation to the database.

    `sim_id` and `run_id` are passed if adding from matsim, and are not passed if
    adding from external source.

    """

    validate_new_atm_sim(params)

    user_cred = user_cred or CONFIG['user']
    user_id = dbs.get_user_id(user_cred)

    # Check it doesn't already exist:
    if get_atm_sim(params['path'], split_path=False, user_cred=user_cred):
        raise ValueError('Simulation already exists.')

    fields = [
        'supercell_type',
        'repeats',
        'num_atoms',
        'registered_path',
        'atom_constraints',
        'software',
        'user_account_id',
    ]
    args = [
        params['supercell_type'],
        params['repeats'],
        params['num_atoms'],
        params['path'],
        params['atom_constraints'],
        params['software'],
        user_id,
    ]
    if sim_id:
        fields.append('sim_id')
        args.append(sim_id)

    values = ', '.join(['%s'] * len(fields))
    fields_str = ', '.join(fields)

    sql = ('insert into atm_sim ({}) values ({})').format(fields_str, values)

    atm_sim_id = dbs.exec_insert(sql, tuple(args))

    return atm_sim_id
    # atm_sim_castep_id = None
    # atm_run_castep_id = None
    # atm_sim_lammps_id = None
    # atm_run_lammps_id = None

    # if params['software'] == 'castep':

    #     cst_params = params['castep_params']

    #     atm_sim_castep_id = add_atm_castep_sim(atm_sim_id, cst_params)

    #     atm_run_castep_id = add_atm_castep_run(
    #         atm_sim_castep_id, cst_params, run_id=run_id)

    # elif params['software'] == 'lammps':

    #     lmp_params = params['lammps_params']
    #     atm_sim_lammps_id = add_atm_lammps_sim(atm_sim_id, lmp_params)

    #     atm_run_lammps_id = add_atm_lammps_run(
    #         atm_sim_lammps_id, run_id=run_id)

    # atm_supercell_gb_csl_id = None
    # atm_supercell_gb = None
    # gb_params = params['gb_params']

    # if params['supercell_type'].startswith('gb'):

    #     atm_supercell_gb_csl_id = add_atm_supercell_gb_csl(
    #         gb_params, atm_sim_id)

    # if params['supercell_type'] == 'gb':

    #     atm_supercell_gb = add_atm_supercell_gb(gb_params, atm_sim_id)

    # vac_params = params['vacuum_expansions']
    # if vac_params:
    #     vac_exp_ids = add_supercell_expansions(vac_params, atm_sim_id)

    # ret = {
    #     'atm_sim_id': atm_sim_id,
    #     'atm_sim_castep_id': atm_sim_castep_id,
    #     'atm_run_castep_id': atm_run_castep_id,
    #     'atm_sim_lammps_id': atm_sim_lammps_id,
    #     'atm_run_lammps_id': atm_run_lammps_id,
    #     'atm_supercell_gb_csl_id': atm_supercell_gb_csl_id,
    #     'atm_supercell_gb': atm_supercell_gb,
    #     'vac_exp_ids': vac_exp_ids,
    # }

    # return ret


def add_supercell_expansions(vac_params, atm_sim_id):
    """Add to the `atm_supercell_expansion` table."""

    vac_exp_ids = []
    for vac_exp_idx, vac_exp in enumerate(vac_params):

        fields_run_vac = [
            'atm_sim_id',
            'magnitude',
            'direction',
            'func',
            'order_id',
        ]
        args_run_vac = [
            atm_sim_id,
            vac_exp['magnitude'],
            vac_exp.get('direction'),
            vac_exp['func'],
            vac_exp_idx,
        ]
        values_run_vac = ', '.join(['%s'] * len(fields_run_vac))
        fields_run_vac = ', '.join(fields_run_vac)

        sql_run_vac = (
            'insert into atm_supercell_expansion ({}) values ({})'
        ).format(fields_run_vac, values_run_vac)

        vac_exp_id = dbs.exec_insert(sql_run_vac, tuple(args_run_vac))
        vac_exp_ids.append(vac_exp_id)

    return vac_exp_ids


def add_atm_supercell_gb_csl(gb_params, atm_sim_id):
    """Add to the `atm_csl_supercell` table."""
    fields_csl = [
        'atm_sim_id',
        'sigma',
        'plane_type',
        'tilt_idx',
    ]
    values_csl = ', '.join(['%s'] * len(fields_csl))
    fields_csl = ', '.join(fields_csl)

    args_csl = [
        atm_sim_id,
        gb_params['sigma'],
        gb_params['plane_type'],
        gb_params.get('tilt_idx'),
    ]

    sql_csl = (
        'insert into atm_csl_supercell ({}) values ({})'
    ).format(fields_csl, values_csl)

    atm_supercell_gb_csl_id = dbs.exec_insert(sql_csl, tuple(args_csl))

    return atm_supercell_gb_csl_id


def add_atm_supercell_gb(gb_params, atm_sim_id):
    """Add to the `atm_grain_boundary_supercell` table."""

    fields_gb = [
        'atm_sim_id',
        'relative_shift_grain_idx',
        'relative_shift_magnitude',
    ]
    values_gb = ', '.join(['%s'] * len(fields_gb))
    fields_gb = ', '.join(fields_gb)
    args_gb = [
        atm_sim_id,
        gb_params.get('relative_shift_grain_idx'),
        gb_params.get('relative_shift_magnitude'),
    ]

    sql_gb = (
        'insert into atm_supercell_gb ({}) values ({})'
    ).format(fields_gb, values_gb)

    atm_supercell_gb_id = dbs.exec_insert(sql_gb, tuple(args_gb))

    return atm_supercell_gb_id


def add_atm_run_lammps(atm_sim_lammps_id, run_id=None):
    """Add to the `atm_run_castep` table."""

    fields_run_lmp = [
        'atm_sim_lammps_id',
    ]
    args_run_lmp = [
        atm_sim_lammps_id,
    ]
    if run_id:
        fields_run_lmp.append('run_id')
        args_run_lmp.append(run_id)

    values_run_lmp = ', '.join(['%s'] * len(fields_run_lmp))
    fields_run_lmp = ', '.join(fields_run_lmp)

    sql_run_lmp = (
        'insert into atm_run_lammps ({}) values ({})'
    ).format(fields_run_lmp, values_run_lmp)

    atm_run_lammps_id = dbs.exec_insert(sql_run_lmp, tuple(args_run_lmp))

    return atm_run_lammps_id


def add_atm_sim_lammps(atm_sim_id, lmp_params):
    """Add to the `atm_sim_castep` table."""

    fields_sim_lmp = [
        'atm_sim_id',
        'potential_name',
    ]
    values_sim_lmp = ', '.join(['%s'] * len(fields_sim_lmp))
    fields_sim_lmp = ', '.join(fields_sim_lmp)

    args_sim_lmp = [
        atm_sim_id,
        lmp_params['potential_name'],
    ]

    sql_sim_lmp = (
        'insert into atm_sim_lammps ({}) values ({})'
    ).format(fields_sim_lmp, values_sim_lmp)

    atm_sim_lammps_id = dbs.exec_insert(sql_sim_lmp, tuple(args_sim_lmp))

    return atm_sim_lammps_id


def add_atm_sim_castep(atm_sim_id, cst_params):
    """Add to the `atm_sim_castep` table."""

    fields_sim_cst = [
        'atm_sim_id',
        'cut_off_energy',
        'task',
        'kpoint_mp_grid',
        'kpoint_mp_offset',
        'kpoint_num',
        'mixing_scheme',
        'smearing_width',
        'metals_method',
        'xc_functional',
        'geom_method',
    ]
    values_sim_cst = ', '.join(['%s'] * len(fields_sim_cst))
    fields_sim_cst = ', '.join(fields_sim_cst)

    args_sim_cst = [
        atm_sim_id,
        cst_params['cut_off_energy'],
        cst_params['task'],
        cst_params['kpoint_mp_grid'],
        cst_params['kpoint_mp_offset'],
        cst_params['kpoint_num'],
        cst_params['mixing_scheme'],
        cst_params['smearing_width'],
        cst_params['metals_method'],
        cst_params['xc_functional'],
        cst_params['geom_method'],
    ]

    sql_sim_cst = (
        'insert into atm_sim_castep ({}) values ({})'
    ).format(fields_sim_cst, values_sim_cst)

    atm_sim_castep_id = dbs.exec_insert(sql_sim_cst, tuple(args_sim_cst))

    return atm_sim_castep_id


def add_atm_run_castep(atm_sim_castep_id, cst_params, run_id=None):
    """Add to the `atm_run_castep` table"""

    fields_run_cst = [
        'atm_sim_castep_id',
        'elec_energy_tol',
        'geom_energy_tol',
        'geom_force_tol',
        'geom_disp_tol',
        'geom_stress_tol',
        'opt_strategy',
    ]

    args_run_cst = [
        atm_sim_castep_id,
        cst_params['elec_energy_tol'],
        cst_params['geom_energy_tol'],
        cst_params['geom_force_tol'],
        cst_params['geom_disp_tol'],
        cst_params['geom_stress_tol'],
        cst_params['opt_strategy'],
    ]

    if run_id:
        fields_run_cst.append('run_id')
        args_run_cst.append(run_id)

    values_run_cst = ', '.join(['%s'] * len(fields_run_cst))
    fields_run_cst = ', '.join(fields_run_cst)

    sql_run_cst = (
        'insert into atm_run_castep ({}) values ({})'
    ).format(fields_run_cst, values_run_cst)

    atm_run_castep_id = dbs.exec_insert(sql_run_cst, tuple(args_run_cst))

    return atm_run_castep_id


def filter_sims(user_cred=None, **kwargs):

    user_cred = user_cred or CONFIG['user']
    user_id = dbs.get_user_id(user_cred)

    is_parse_json = [
        'repeats',
    ]
    is_parse_bool = [
        'atom_constraints'
    ]

    variables = {
        'repeats': 'sim',
        'sigma': 'csl',
        'kpoint_num': 'scst',
        'atom_constraints': 'sim',
        'cut_off_energy': 'scst',
        'elec_energy_tol': 'rcst',
        'kpoint_mp_grid': 'scst',
        'task': 'scst',
        'geom_method': 'scst',
        'geom_energy_tol': 'rcst',
        'geom_force_tol': 'rcst',
        'geom_disp_tol': 'rcst',
        'geom_stress_tol': 'rcst',
        'mixing_scheme': 'scst',
        'metals_method': 'scst',
        'smearing_width': 'scst',
        'opt_strategy': 'rcst',
        'num_atoms': 'sim',
        'supercell_type': 'sim',
        'plane_type': 'csl',
        'tilt_idx': 'csl',
        'software': 'sim',
        'relative_shift_grain_idx': 'gbs',
        'relative_shift_magnitude': 'gbs',
    }

    sql = (
        'select sim.id from atm_sim sim '
        'inner join user_account ua on ua.id = sim.user_account_id '
        'inner join atm_csl_supercell csl on csl.atm_sim_id = sim.id '
        'inner join atm_grain_boundary_supercell gbs on gbs.atm_sim_id = sim.id '
        'inner join atm_sim_castep scst on scst.atm_sim_id = sim.id '
        'inner join atm_run_castep rcst on rcst.atm_sim_castep_id = scst.id '
        'where ua.id = %s '
    )

    args = [user_id]

    for filter_name, filter_val in kwargs.items():

        if filter_name not in variables:
            raise ValueError('Not valid filter variable.')

        # prt(filter_name, 'filt name')
        # prt(type(filter_val), 'filt val')

        var_table = variables[filter_name]

        sql += (
            'and {}.{} = %s'.format(var_table, filter_name)
        )

        if filter_name in is_parse_json:
            filter_val = json.dumps(filter_val)

        if filter_name in is_parse_bool:
            filter_val = bool(filter_val)

        # prt(type(filter_val), 'filt val')

        args.append(filter_val)

    # TODO: one to many relationships e.g. supercell expansions

    # print('sql: {}'.format(sql))
    # print('args: {}'.format(args))

    atm_sim_ids = dbs.exec_select(sql, tuple(args), fetch_all=True)

    # print('atm_sim_ids: {}'.format(atm_sim_ids))
    atm_sim_ids = [i['id'] for i in atm_sim_ids]

    print('atm_sim_ids: {}'.format(atm_sim_ids))

    # exit(0)

    sim_id_wherein = ', '.join(['%s'] * len(atm_sim_ids))
    sql_distinct = (
        'select distinct {}.{} from atm_sim sim '
        'inner join user_account ua on ua.id = sim.user_account_id '
        'inner join atm_csl_supercell csl on csl.atm_sim_id = sim.id '
        'inner join atm_grain_boundary_supercell gbs on gbs.atm_sim_id = sim.id '
        'inner join atm_sim_castep scst on scst.atm_sim_id = sim.id '
        'inner join atm_run_castep rcst on rcst.atm_sim_castep_id = scst.id '
        'where ua.id = %s '
        'and sim.id in ({})'
    )

    distinct_args = [user_id] + atm_sim_ids

    # prt(distinct_args, 'distin arg')

    # Get distinct variables: use this instead (to do in one query):
    # select
    # group_concat(distinct repeats) as reps_grp,
    # group_concat(distinct num_atoms) as natoms_grp,
    # group_concat(distinct software) as software_grp,
    # ...
    # from atm_sim

    distinct_vars = {}
    common_vars = {}

    for var_name, var_tab in variables.items():

        sql_distinct_var = sql_distinct.format(
            var_tab, var_name, sim_id_wherein)

        # prt(sql_distinct_var, 'sql_distinct_var')

        dst_vars = dbs.exec_select(
            sql_distinct_var, tuple(distinct_args), fetch_all=True)

        dst_vars = [i[var_name] for i in dst_vars]
        # prt(dst_vars, 'dst_vars')

        if var_name in is_parse_json:
            dst_vars = [json.loads(i) for i in dst_vars]

        if var_name in is_parse_bool:
            dst_vars = [bool(i) for i in dst_vars]

        # prt(dst_vars, 'dst_vars')

        if len(dst_vars) == 1:
            common_vars.update({
                var_name: dst_vars[0]
            })
        else:
            distinct_vars.update({
                var_name: dst_vars
            })

    return atm_sim_ids, common_vars, distinct_vars
