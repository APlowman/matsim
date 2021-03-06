"""matsim.database.py"""
import json

import pymysql.cursors
from passlib.hash import pbkdf2_sha256

from matsim import DB_CONFIG, CONFIG
from matsim.utils import prt

_AUTO_LOOKUP = {
    True: 'true',
    False: 'false',
    'ask': 'ask',
    'true': True,
    'false': False,
}


def connect_db(dict_curs=True):
    """Connect to the database."""
    if dict_curs:
        connection = pymysql.connect(**DB_CONFIG, charset='utf8mb4',
                                     cursorclass=pymysql.cursors.DictCursor)
    else:
        connection = pymysql.connect(**DB_CONFIG, charset='utf8mb4')
    return connection


def exec_insert(sql, args):
    """Execute an insert SQL statement.

    Parameters
    ----------
    sql : str
    args : tuple

    Returns
    -------
    last_row_id : int

    """

    last_row_id = None
    conn = connect_db()
    try:
        with conn.cursor() as cursor:
            cursor.execute(sql, args)
            last_row_id = cursor.lastrowid

        conn.commit()

    finally:
        conn.close()

    return last_row_id


def exec_select(sql, args, fetch_all=False, dict_curs=True):
    """Execute a select SQL statement.

    Parameters
    ----------
    sql : str
    args : tuple
    fetch_all : bool
        If True, use `cursor.fetchall()` else use `cursor.fetchone()`.

    Returns
    -------
    result : dict

    """

    conn = connect_db(dict_curs)
    try:
        with conn.cursor() as cursor:

            cursor.execute(sql, args)

            if fetch_all:
                result = cursor.fetchall()
            else:
                result = cursor.fetchone()

    finally:
        conn.close()

    return result


def exec_update(sql, args):
    """Execute an update SQL statement.

    Parameters
    ----------
    sql : str
    args : tuple

    Returns
    -------


    """

    conn = connect_db()
    try:
        with conn.cursor() as cursor:
            cursor.execute(sql, args)

        conn.commit()

    finally:
        conn.close()


def add_user_account(user_name, password, email):
    """Add a user account."""

    sql = (
        'insert into user_account (user_name, password, email) '
        'values (%s, %s, %s)'
    )
    pass_hash = pbkdf2_sha256.hash(password)
    user_id = exec_select(sql, (user_name, pass_hash, email))

    return user_id


def get_user_id(user_cred):
    """Get user account id of given user."""

    password = user_cred['password']
    user_name = user_cred['name']

    sql = (
        'select id, password from user_account where user_name = %s'
    )
    msg = 'Incorrect user credentials.'

    user = exec_select(sql, (user_name,))
    if not user:
        raise ValueError(msg)

    pass_hash = user['password']
    user_id = user['id']

    is_verified = pbkdf2_sha256.verify(password, pass_hash)
    if is_verified:
        return user_id

    else:
        raise ValueError(msg)


def get_machines(user_cred=None):
    """Retrieve all machines belonging to a user."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select * from machine where user_account_id = %s'
    )
    result = exec_select(sql, (user_id,), fetch_all=True)

    return result


def add_machine(name, os_type, is_dropbox, user_cred=None):
    """Add a new machine associated with a given user."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    for mach in get_machines(user_cred):
        if mach['name'] == name:
            msg = 'Machine with name "{}" already exists.'.format(name)
            raise ValueError(msg)

    sql = (
        'insert into machine (user_account_id, name, os_type, is_dropbox) '
        'values (%s, %s, %s, %s)'
    )
    insert_id = exec_insert(sql, (user_id, name, os_type, is_dropbox))

    return insert_id


def add_stage(name, base_path, machine_name, user_cred=None):
    """Add a stage resource associated with a given user."""

    user_cred = user_cred or CONFIG['user']
    _ = get_user_id(user_cred)

    # Get machine ID:
    all_mach = get_machines(user_cred)
    machine_id = None
    for mach in all_mach:
        if mach['name'] == machine_name:
            machine_id = mach['id']

            # Check it's not dropbox
            if mach['is_dropbox']:
                msg = 'Cannot add a Stage resource to a dropbox machine.'
                raise ValueError(msg)

            break

    # Check machine id is owned by user:
    if not machine_id:
        raise ValueError('You do not have a machine with that name.')

    for stage in get_stages(user_cred):
        if stage['name'] == name:
            msg = 'Stage with name: "{}" already exists.'.format(name)
            raise ValueError(msg)

    sql_res = (
        'insert into resource (name, base_path, machine_id) '
        'values (%s, %s, %s)'
    )
    sql_stage = (
        'insert into stage (resource_id) '
        'values (%s)'
    )

    res_insert_id = exec_insert(sql_res, (name, base_path, machine_id))
    stage_insert_id = exec_insert(sql_stage, (res_insert_id,))

    return stage_insert_id


def add_scratch(name, base_path, machine_name, is_sge, user_cred=None):
    """Add a scratch resource associated with a given user."""

    user_cred = user_cred or CONFIG['user']
    _ = get_user_id(user_cred)

    # Get machine ID:
    all_mach = get_machines(user_cred)
    machine_id = None
    for mach in all_mach:
        if mach['name'] == machine_name:
            machine_id = mach['id']

            # Check it's not dropbox
            if mach['is_dropbox']:
                msg = 'Cannot add a Scratch resource to a dropbox machine.'
                raise ValueError(msg)

            break

    # Check machine id is owned by user:
    if not machine_id:
        msg = 'You do not have a machine called "{}".'
        raise ValueError(msg.format(machine_name))

    for scratch in get_scratches(user_cred):
        if scratch['name'] == name:
            msg = 'Scratch with name: "{}" already exists.'.format(name)
            raise ValueError(msg)

    sql_res = (
        'insert into resource (name, base_path, machine_id) '
        'values (%s, %s, %s)'
    )
    sql_scratch = (
        'insert into scratch (is_sge, resource_id) '
        'values (%s, %s)'
    )

    res_insert_id = exec_insert(sql_res, (name, base_path, machine_id))
    scratch_insert_id = exec_insert(sql_scratch, (is_sge, res_insert_id))

    return scratch_insert_id


def add_archive(name, base_path, machine_name, user_cred=None):
    """Add an archive resource associated with a given user."""

    user_cred = user_cred or CONFIG['user']
    _ = get_user_id(user_cred)

    # Get machine ID:
    all_mach = get_machines(user_cred)
    machine_id = None
    for mach in all_mach:
        if mach['name'] == machine_name:
            machine_id = mach['id']
            break

    # Check machine id is owned by user:
    if not machine_id:
        raise ValueError('You do not have a machine with that name.')

    for archive in get_archives(user_cred):
        if archive['name'] == name:
            msg = 'Archive with name: "{}" already exists.'.format(name)
            raise ValueError(msg)

    sql_res = (
        'insert into resource (name, base_path, machine_id) '
        'values (%s, %s, %s)'
    )
    sql_archive = (
        'insert into archive (resource_id) '
        'values (%s)'
    )

    res_insert_id = exec_insert(sql_res, (name, base_path, machine_id))
    archive_insert_id = exec_insert(sql_archive, (res_insert_id,))

    return archive_insert_id


def get_scratches(user_cred=None):
    """Retrieve all scratch resources belonging to a user."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select s.id "scratch_id", r.id "resource_id", m.id "machine_id", '
        'm.name "machine_name", r.name, r.base_path, '
        'm.is_dropbox, s.is_sge '
        'from scratch s '
        'inner join resource r on s.resource_id = r.id '
        'inner join machine m on r.machine_id = m.id '
        'where m.user_account_id = %s'
    )

    return exec_select(sql, (user_id,), fetch_all=True)


def get_stages(user_cred=None):
    """Retrieve all stage resources belonging to a user."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select s.id "stage_id", r.id "resource_id", m.id "machine_id", '
        'm.name "machine_name", r.name, r.base_path, '
        'm.is_dropbox '
        'from stage s '
        'inner join resource r on s.resource_id = r.id '
        'inner join machine m on r.machine_id = m.id '
        'where m.user_account_id = %s'
    )

    return exec_select(sql, (user_id,), fetch_all=True)


def get_machine_by_name(name, user_cred=None):
    """Get a machine by its name for a given user."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select * '
        'from machine '
        'where '
        'name = %s and '
        'user_account_id = %s'
    )
    machine = exec_select(sql, (name, user_id))

    if not machine:
        msg = ('No machine with name "{}" exists for given user.')
        raise ValueError(msg.format(name))

    return machine


def get_archive_by_name(name, user_cred=None):
    """Get an archive resource by name."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select a.id "archive_id", r.id "resource_id", r.name "name" '
        'from archive a '
        'inner join resource r on a.resource_id = r.id '
        'inner join machine m on r.machine_id = m.id '
        'where '
        'm.user_account_id = %s and '
        'r.name = %s'
    )
    archive = exec_select(sql, (user_id, name))

    if not archive:
        msg = ('No archive with name "{}" exists for given user.')
        raise ValueError(msg.format(name))

    return archive


def get_archive_by_id(archive_id, user_cred=None):
    """Get an archive resource by ID."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select a.id "archive_id", r.id "resource_id", r.name "name" '
        'from archive a '
        'inner join resource r on a.resource_id = r.id '
        'inner join machine m on r.machine_id = m.id '
        'where '
        'm.user_account_id = %s and '
        'a.id = %s'
    )
    archive = exec_select(sql, (user_id, archive_id))

    if not archive:
        msg = ('No archive with ID "{}" exists for given user.')
        raise ValueError(msg.format(archive_id))

    return archive


def get_scratch_by_name(name, user_cred=None):
    """Get a scratch resource by name."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select s.id "scratch_id", s.is_sge "scratch_is_sge", '
        'r.id "resource_id", r.name "name" '
        'from scratch s '
        'inner join resource r on s.resource_id = r.id '
        'inner join machine m on r.machine_id = m.id '
        'where '
        'm.user_account_id = %s and '
        'r.name = %s'
    )
    scratch = exec_select(sql, (user_id, name))

    if not scratch:
        msg = ('No scratch with name "{}" exists for given user.')
        raise ValueError(msg.format(name))

    return scratch


def get_scratch_by_id(scratch_id, user_cred=None):
    """Get a scratch resource by its ID."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select s.id "scratch_id", s.is_sge "scratch_is_sge", '
        'r.id "resource_id", r.name "name" '
        'from scratch s '
        'inner join resource r on s.resource_id = r.id '
        'inner join machine m on r.machine_id = m.id '
        'where '
        'm.user_account_id = %s and '
        's.id = %s'
    )
    scratch = exec_select(sql, (user_id, scratch_id))

    if not scratch:
        msg = ('No scratch with ID "{}" exists for given user.')
        raise ValueError(msg.format(scratch_id))

    return scratch


def get_scratch_by_resource_id(resource_id, user_cred=None):
    """Get a scratch resource by its resource ID."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select s.id "scratch_id", s.is_sge "scratch_is_sge", '
        'r.id "resource_id", r.name "name" '
        'from scratch s '
        'inner join resource r on s.resource_id = r.id '
        'inner join machine m on r.machine_id = m.id '
        'where '
        'm.user_account_id = %s and '
        'r.id = %s'
    )
    scratch = exec_select(sql, (user_id, resource_id))

    if not scratch:
        msg = ('No scratch with resource ID "{}" exists for given user.')
        raise ValueError(msg.format(resource_id))

    return scratch


def get_stage_by_resource_id(resource_id, user_cred=None):
    """Get a stage resource by its resource ID."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select s.id "stage_id" '
        'r.id "resource_id", r.name "name" '
        'from stage s '
        'inner join resource r on s.resource_id = r.id '
        'inner join machine m on r.machine_id = m.id '
        'where '
        'm.user_account_id = %s and '
        'r.id = %s'
    )
    stage = exec_select(sql, (user_id, resource_id))

    if not stage:
        msg = ('No stage with resource ID "{}" exists for given user.')
        raise ValueError(msg.format(resource_id))

    return stage


def get_stage_by_name(name, user_cred=None):
    """Get a stage resource by name."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select s.id "stage_id", r.id "resource_id", r.name "name"'
        'from stage s '
        'inner join resource r on s.resource_id = r.id '
        'inner join machine m on r.machine_id = m.id '
        'where '
        'm.user_account_id = %s and '
        'r.name = %s'
    )
    stage = exec_select(sql, (user_id, name))

    if not stage:
        msg = ('No stage with name "{}" exists for given user.')
        raise ValueError(msg.format(name))

    return stage


def get_stage_by_id(stage_id, user_cred=None):
    """Get a stage resource by its ID."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select s.id "stage_id", r.id "resource_id", r.name "name" '
        'from stage s '
        'inner join resource r on s.resource_id = r.id '
        'inner join machine m on r.machine_id = m.id '
        'where '
        'm.user_account_id = %s and '
        's.id = %s'
    )
    stage = exec_select(sql, (user_id, stage_id))

    if not stage:
        msg = ('No stage with ID "{}" exists for given user.')
        raise ValueError(msg.format(stage_id))

    return stage


def get_resource(resource_id, user_cred=None):
    """Get a resource by id."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select r.id "resource_id", r.name "resource_name", '
        'r.base_path "resource_base_path", m.id "machine_id", '
        'm.os_type "machine_os_type", m.is_dropbox "machine_is_dropbox", '
        'm.name "machine_name" '
        'from resource r '
        'inner join machine m on r.machine_id = m.id '
        'where '
        'm.user_account_id = %s and '
        'r.id = %s'
    )
    resource = exec_select(sql, (user_id, resource_id))

    if not resource:
        msg = ('No resource with ID "{}" exists for given user.')
        raise ValueError(msg.format(resource_id))

    return resource


def get_resource_connection(src_res_id, dst_res_id, user_cred=None):
    """Get a resource connection by source and destination resource ID."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select '
        'rc.host, rc.is_remote '
        'from resource_connection rc '
        'inner join resource rs on rs.id = rc.source_id '
        'inner join resource rd on rd.id = rc.destination_id '
        'inner join machine ms on rs.machine_id = ms.id '
        'inner join machine md on rd.machine_id = md.id '
        'where '
        'ms.user_account_id = %s and '
        'md.user_account_id = %s and '
        'rs.id = %s and '
        'rd.id = %s'
    )
    res_conn = exec_select(sql, (user_id, user_id, src_res_id, dst_res_id))
    if not res_conn:
        msg = ('No resource connection between source ID {} and destination '
               'ID {} is configured.')
        raise ValueError(msg.format(src_res_id, dst_res_id))

    return res_conn


def get_resource_connections_by_source(src_res_id, user_cred=None):
    """Get a all resource connection for a given source resource ID."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select '
        'rc.host, rc.is_remote, rc.destination_id '
        'from resource_connection rc '
        'inner join resource rs on rs.id = rc.source_id '
        'inner join resource rd on rd.id = rc.destination_id '
        'inner join machine ms on rs.machine_id = ms.id '
        'inner join machine md on rd.machine_id = md.id '
        'where '
        'ms.user_account_id = %s and '
        'md.user_account_id = %s and '
        'rs.id = %s'
    )
    res_conn = exec_select(sql, (user_id, user_id, src_res_id), fetch_all=True)
    if not res_conn:
        msg = ('No resource connection with source ID {} is configured.')
        raise ValueError(msg.format(src_res_id))

    return res_conn


def get_scratches_from_stage(stage_res_id, user_cred=None):
    """Get resource IDs of all scratches connectable from a given Stage
    resource ID.

    """

    check_resource_ownership(stage_res_id, user_cred)

    sql = (
        'select rc.destination_id from resource_connection rc '
        'inner join scratch sc on sc.resource_id = rc.destination_id '
        'where rc.source_id = %s'
    )

    return exec_select(sql, (stage_res_id,), fetch_all=True)


def get_archives(user_cred=None):
    """Retrieve all archive resources belonging to a user."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select a.id "archive_id", r.id "resource_id", m.id "machine_id", '
        'm.name "machine_name", r.name, r.base_path, '
        'm.is_dropbox '
        'from archive a '
        'inner join resource r on a.resource_id = r.id '
        'inner join machine m on r.machine_id = m.id '
        'where m.user_account_id = %s'
    )

    return exec_select(sql, (user_id,), fetch_all=True)


def get_resource_by_name_type(res_name, res_type, user_cred=None):
    """Find a resource by its name and type (one of stage, scratch or
    archive).

    """
    user_cred = user_cred or CONFIG['user']

    msg = ('Cannot find a resource with name: "{}" and type: '
           '"{}"'.format(res_name, res_type))

    if res_type == 'stage':
        all_res = get_stages(user_cred)

    elif res_type == 'scratch':
        all_res = get_scratches(user_cred)

    elif res_type == 'archive':
        all_res = get_archives(user_cred)

    else:
        raise ValueError(msg)

    for res in all_res:
        if res['name'] == res_name:
            return res

    raise ValueError(msg)


def add_resource_connection(source_name, source_type, dest_name, dest_type,
                            host=None, user_cred=None):
    """Add a resource connection between a source and destination."""

    user_cred = user_cred or CONFIG['user']
    _ = get_user_id(user_cred)

    # Find resource id of source and dest:
    src = get_resource_by_name_type(source_name, source_type, user_cred)
    dst = get_resource_by_name_type(dest_name, dest_type, user_cred)
    src_id = src['resource_id']
    dst_id = dst['resource_id']

    is_remote = False
    if src['machine_name'] != dst['machine_name']:
        is_remote = True

    sql = (
        'insert into resource_connection (host, is_remote, source_id, '
        'destination_id) '
        'values (%s, %s, %s, %s)'
    )
    res_conn_id = exec_insert(sql, (host, is_remote, src_id, dst_id))

    return res_conn_id


def get_resource_connections(user_cred=None):
    """Get resource connections for a given user."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    # Big sql statement since we want to identify what type the source and
    # destination resources are (i.e. stage, scratch or archive):
    sql = (
        'select '
        'rc.host, '
        'rc.is_remote, '
        'rc.source_id "src_id", '
        'rc.destination_id "dst_id", '
        'rs.name "src_name", '
        'rs.machine_id "src_machine_id", '
        'rd.name "dst_name", '
        'rd.machine_id "dst_machine_id", '
        'ms.user_account_id "machine_src_user_id", '
        'md.user_account_id "machine_dst_user_id", '
        ''
        'stgs.resource_id "stage_src_resource_id", '
        'stgd.resource_id "stage_dst_resource_id", '
        'scrs.resource_id "scratch_src_resource_id", '
        'scrd.resource_id "scratch_dst_resource_id", '
        'arcs.resource_id "archive_src_resource_id", '
        'arcd.resource_id "archive_dst_resource_id" '
        ''
        'from resource_connection rc '
        'inner join resource rs on rs.id = rc.source_id '
        'inner join resource rd on rd.id = rc.destination_id '
        'inner join machine ms on ms.id = rs.machine_id '
        'inner join machine md on md.id = rd.machine_id '
        ''
        'left outer join stage stgs on stgs.resource_id = rs.id '
        'left outer join stage stgd on stgd.resource_id = rd.id '
        'left outer join scratch scrs on scrs.resource_id = rs.id '
        'left outer join scratch scrd on scrd.resource_id = rd.id '
        'left outer join archive arcs on arcs.resource_id = rs.id '
        'left outer join archive arcd on arcd.resource_id = rd.id '
        ''
        'where ms.user_account_id = %s and md.user_account_id = %s'
    )

    res_conn_all = exec_select(sql, (user_id, user_id), fetch_all=True)

    bad_msg = ('Resource connection {} ID does not refer to a Stage, '
               'Scratch or Archive resource!.')

    # Now add a type key (stage, scratch, archive) with respective ids:
    for res_conn_idx in range(len(res_conn_all)):

        res_conn = res_conn_all[res_conn_idx]

        # Source:
        src_ids = [
            res_conn['stage_src_resource_id'],
            res_conn['scratch_src_resource_id'],
            res_conn['archive_src_resource_id']
        ]
        src_id_notnone = [i is not None for i in src_ids]

        if sum(src_id_notnone) != 1:
            raise ValueError(bad_msg.format('source'))

        src_id_type_idx = src_id_notnone.index(True)

        source = {}

        if src_id_type_idx == 0:
            src_type = 'stage'
            source.update({
                'stage_id': res_conn['stage_src_resource_id'],
            })
        elif src_id_type_idx == 1:
            src_type = 'scratch'
            source.update({
                'scratch_id': res_conn['scratch_src_resource_id'],
            })
        elif src_id_type_idx == 2:
            src_type = 'archive'
            source.update({
                'archive_id': res_conn['archive_src_resource_id'],
            })
        source.update({
            'name': res_conn['src_name'],
        })

        # Destination:
        dst_ids = [
            res_conn['stage_dst_resource_id'],
            res_conn['scratch_dst_resource_id'],
            res_conn['archive_dst_resource_id']
        ]
        dst_id_notnone = [i is not None for i in dst_ids]

        if sum(dst_id_notnone) != 1:
            raise ValueError(bad_msg.format('destination'))

        dst_id_type_idx = dst_id_notnone.index(True)

        destination = {}
        if dst_id_type_idx == 0:
            dst_type = 'stage'
            destination.update({
                'stage_id': res_conn['stage_dst_resource_id'],
            })
        elif dst_id_type_idx == 1:
            dst_type = 'scratch'
            destination.update({
                'scratch_id': res_conn['scratch_dst_resource_id'],
            })
        elif dst_id_type_idx == 2:
            dst_type = 'archive'
            destination.update({
                'archive_id': res_conn['archive_dst_resource_id'],
            })
        destination.update({
            'name': res_conn['dst_name'],
        })

        # Now overwrite the list with more useful keys:
        res_conn_all[res_conn_idx] = {
            'host': res_conn['host'],
            'is_remote': res_conn['is_remote'],
            'src_resource_type': src_type,
            'dst_resource_type': dst_type,
            'src': source,
            'dst': destination,
        }

    return res_conn_all


def add_software(name, user_cred=None):
    """Add a software for this user."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    for soft in get_all_softwares(user_cred):
        if soft['name'] == name:
            msg = 'Software with name "{}" already exists.'.format(name)
            raise ValueError(msg)

    sql = (
        'insert into software (user_account_id, name) '
        'values (%s, %s)'
    )
    insert_id = exec_insert(sql, (user_id, name))

    return insert_id


def get_software(name, user_cred=None):
    """Get a particular software by name, for this user."""

    user_cred = user_cred or CONFIG['user']
    _ = get_user_id(user_cred)

    sql = (
        'select * '
        'from software '
        'where name = %s'
    )
    software = exec_select(sql, (name,))

    if not software:
        msg = ('No software with name "{}" exists.')
        raise ValueError(msg.format(name))

    return software


def get_all_softwares(user_cred=None):
    """Get all softwares for this user."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select s.id, s.name '
        'from software s '
        'where s.user_account_id = %s'
    )

    return exec_select(sql, (user_id,), fetch_all=True)


def add_software_instance(name, module_load, core_range, executable, vers,
                          parallel_env, software_name, scratches, user_cred=None):
    """Add a software instance for this user."""

    user_cred = user_cred or CONFIG['user']
    _ = get_user_id(user_cred)

    for soft_ins in get_all_software_instances(user_cred):
        if soft_ins['name'] == name:
            msg = 'Software instance with name "{}" already exists.'
            raise ValueError(msg.format(name))

    # Get scratch IDs:
    scratch_ids = []
    for scratch_name in scratches:
        s_id = get_scratch_by_name(scratch_name)['scratch_id']
        scratch_ids.append(s_id)

    # Get software_id:
    software_id = get_software(software_name, user_cred)['id']
    min_cores = core_range[0]
    max_cores = core_range[1]

    sql = (
        'insert into software_instance (name, module_load, '
        'min_cores, max_cores, executable, version, parallel_env, '
        'software_id) '
        'values (%s, %s, %s, %s, %s, %s, %s, %s)'
    )
    si_insert_id = exec_insert(sql, (
        name, module_load, min_cores, max_cores, executable, vers,
        parallel_env, software_id
    ))

    # Set scratch IDs this is allowed on:
    scratch_software_insert_ids = []
    for s_id in scratch_ids:

        sql_scratch_soft_inst = (
            'insert into scratch_software_instance '
            '(scratch_id, software_instance_id) '
            'values (%s, %s)'
        )
        s_si_insert_id = exec_insert(
            sql_scratch_soft_inst, (s_id, si_insert_id))

        scratch_software_insert_ids.append(s_si_insert_id)

    return si_insert_id, scratch_software_insert_ids


def get_all_software_instances(user_cred=None):
    """Get all software instances for this user."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select si.* '
        'from software_instance si '
        'inner join software s on s.id = si.software_id '
        'where s.user_account_id = %s'
    )

    return exec_select(sql, (user_id,), fetch_all=True)


def get_software_instance_by_name(name, user_cred=None):
    """Get a particular software instance by name, for this user."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select si.*, s.name "software_name" '
        'from software_instance si '
        'inner join software s on s.id = si.software_id '
        'where si.user_account_id = %s and '
        'si.name = %s'
    )

    software_inst = exec_select(sql, (user_id, name))

    if not software_inst:
        msg = ('No software instance with name "{}" exists for given user.')
        raise ValueError(msg.format(name))

    return software_inst


def get_software_instance_by_id(software_instance_id, user_cred=None):
    """Get a particular software instance by ID."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select si.* , s.name "software_name" '
        'from software_instance si '
        'inner join software s on s.id = si.software_id '
        'where si.user_account_id = %s and '
        'si.id = %s'
    )
    software_inst = exec_select(sql, (user_id, software_instance_id))

    if not software_inst:
        msg = ('No software instance with ID "{}" exists for given user.')
        raise ValueError(msg.format(software_instance_id))

    return software_inst


def get_software_instance_ok_scratch(name, user_cred=None):
    """Get the scratch IDs that allow a given software instance."""

    user_cred = user_cred or CONFIG['user']
    _ = get_user_id(user_cred)

    soft_inst_id = get_software_instance_by_name(name, user_cred)['id']

    sql = (
        'select ssi.scratch_id '
        'from scratch_software_instance ssi '
        'where ssi.software_instance_id = %s'
    )
    scratch_ids = exec_select(sql, (soft_inst_id,), fetch_all=True)
    scratch_ids = [i['scratch_id'] for i in scratch_ids]

    if not scratch_ids:
        msg = (
            'No scratch IDs allow software instance with name "{}" '
            'for given user.'
        )
        raise ValueError(msg.format(name))

    return scratch_ids


def get_sim_group_state_id(human_id, user_cred=None):
    """Get the state of a given sim group, return False if does not exist."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select sim_group_state_id from sim_group '
        'where human_id = %s '
        'and user_account_id = %s'
    )
    state_id = exec_select(sql, (human_id, user_id))

    if not state_id:
        return False
    else:
        return state_id['sim_group_state_id']


def set_sim_group_state_id(human_id, state_id, user_cred=None):
    """Get sim group state ID."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'update sim_group '
        'set sim_group_state_id = %s '
        'where human_id = %s '
        'and user_account_id = %s'
    )
    args = (state_id, human_id, user_id)
    exec_update(sql, args)


def get_sim_group_sims(sim_group_id, user_cred=None):
    """Get all sims associated with a given sim group."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select s.* '
        'from sim s '
        'inner join sim_group sg on sg.id = s.sim_group_id '
        'where sg.id = %s '
        'and sg.user_account_id = %s '
        'order by s.order_in_sim_group'

    )
    sims = exec_select(sql, (sim_group_id, user_id), fetch_all=True)

    return sims


def get_sim_group_by_human_id(human_id, user_cred=None):
    """Get a SimGroup from the database by human_id"""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select * '
        'from sim_group sg '
        'where sg.human_id = %s and '
        'sg.user_account_id = %s'
    )
    sim_group = exec_select(sql, (human_id, user_id))

    if not sim_group:
        msg = ('No sim group with human_id "{}" exists for given user.')
        raise ValueError(msg.format(human_id))

    sg_id = sim_group['id']

    # Deserialise path options:
    sim_group['path_options'] = json.loads(sim_group['path_options'])

    # Get stage, scratch, archive by ID:
    run_opt = {
        'stage': get_stage_by_id(sim_group['stage_id']),
        'scratch': get_scratch_by_id(sim_group['scratch_id']),
        'archive': get_archive_by_id(sim_group['archive_id']),
        'groups': [],
    }

    # Get run groups:
    software_name = None
    all_run_groups = get_sim_group_run_groups(sg_id)
    for rg_idx, _ in enumerate(all_run_groups):

        run_group = all_run_groups[rg_idx]

        # get software instance name:
        soft_inst = get_software_instance_by_id(
            run_group['software_instance_id'])

        rg_soft_name = soft_inst['software_name']

        if software_name is None:
            software_name = rg_soft_name
            sim_group.update({
                'software_name': software_name
            })
        elif software_name != rg_soft_name:
            raise ValueError('Software name problemo!')

        run_group.update({
            'dbid': run_group.pop('id'),
            'auto_submit': _AUTO_LOOKUP[run_group['auto_submit']],
            'auto_process': _AUTO_LOOKUP[run_group['auto_process']],
            'job_array': bool(run_group['run_group_sge_job_array']),
            'selective_submission': bool(
                run_group['run_group_sge_selective_submission']),
        })

        run_opt['groups'].append({
            **run_group,
            'software_instance': soft_inst,
            'software': rg_soft_name,
        })

    run_opt.update({
        'software': rg_soft_name,
    })
    sim_group.update({
        'run_opt': run_opt
    })

    return sim_group


def is_sim_group_processing(sim_group_id, user_cred=None):
    """Check if a sim group is being processed."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select processing_now '
        'from sim_group sg '
        'inner join user_account ua on ua.id = sg.user_account_id '
        'where sg.id = %s '
        'and ua.id = %s '
    )
    args = (sim_group_id, user_id)
    result = exec_select(sql, args)
    if not result:
        msg = 'No sim group with id: {}, belonging to current user.'
        raise ValueError(msg.format(sim_group_id))

    return bool(result['processing_now'])


def get_machine_sync_client_path(dropbox_machine_id,
                                 current_machine_name=CONFIG['machine_name'],
                                 user_cred=None):
    """Get path of sync client folder on current machine."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    # Get machine ID of current machine
    cur_mach_id = get_machine_by_name(
        current_machine_name, user_cred=user_cred)['id']

    sql = (
        'select ms.* '
        'from machine_sync_client_path ms '
        'inner join machine m on m.id = ms.dropbox_machine_id '
        'where ms.dropbox_machine_id = %s '
        'and m.is_dropbox = %s '
        'and ms.sync_client_machine_id = %s '
        'and m.user_account_id = %s'
    )
    args = (dropbox_machine_id, 1, cur_mach_id, user_id)
    result = exec_select(sql, args)

    if not result:
        return False

    return result, current_machine_name


def set_sim_group_processing(sim_group_id, is_processing, user_cred=None):
    """Set a sim group to processing state."""

    user_cred = user_cred or CONFIG['user']
    _ = get_user_id(user_cred)

    # # First check ownership (MySQL warning raised when doing it in one query
    # # for some reason.)
    # sql = (
    #     'select sg.id from sim_group sg '
    #     'inner join user_account ua on ua.id = sg.user_account_id '
    #     'where sg.id = %s '
    #     'and ua.id = %s'
    # )
    # result = exec_select(sql, (sim_group_id, user_id))
    # if not result:
    #     raise ValueError('No sim group.')

    sql = (
        'update sim_group sg '
        'set sg.processing_now = %s '
        'where sg.id = %s '
    )
    args = (is_processing, sim_group_id)
    exec_update(sql, args)


def get_software_id_by_name(software_name):
    """Get the software ID by it's (unique) name."""

    sql = (
        'select * from software where name = %s'
    )
    software_id = exec_select(sql, (software_name, ))['id']
    return software_id


def add_sim_group(sim_group, user_cred=None):
    """Add a SimGroup object to the database."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    # Check the sim group is not already added
    if sim_group.dbid:
        raise ValueError('Sim group has already been added to the database.')

    # Add to the sim_group table:
    sql = (
        'insert into sim_group '
        '(human_id, user_account_id, stage_id, scratch_id, archive_id, '
        'path_options, name, sim_group_state_id, software_id, processing_now) '
        'values (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'
    )
    path_opt = json.dumps(sim_group.path_options)
    software_id = get_software_id_by_name(sim_group.software)
    args = (
        sim_group.human_id,
        user_id,
        sim_group.stage.stage_id,
        sim_group.scratch.scratch_id,
        sim_group.archive.archive_id,
        path_opt,
        sim_group.name,
        1,
        software_id,
        False
    )

    # prt(sql, 'sql')
    # prt(args, 'args')

    sg_id = exec_insert(sql, args)
    sim_group.dbid = sg_id

    # Add sims to sim table
    for sim_idx in range(sim_group.num_sims):

        sql_sim = (
            'insert into sim '
            '(sim_group_id, order_in_sim_group) '
            'values (%s, %s)'
        )
        sim_insert = exec_insert(sql_sim, (sg_id, sim_idx))
        sim_group.sims[sim_idx].dbid = sim_insert

    # # Add to the run_group table
    # rg_insert_ids = []
    # for run_group in sim_group.run_options['groups']:

    #     sge = sim_group.scratch.sge
    #     rg_insert = add_run_group(sg_id, run_group, sim_insert_ids, 1, sge)
    #     rg_insert_ids.append(rg_insert)


def get_sim_group_runs(sim_group_id, states=None, user_cred=None):
    """Get all runs associated with a given sim group, optionally in given
    states."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select r.*, s.order_in_sim_group "sim_order_in_sim_group", '
        'rg.order_in_sim_group "run_group_order_in_sim_group", '
        'rgsge.job_id "run_group_sge_job_id", '
        'rgsge.job_array "run_group_sge_job_array" '
        'from run_group rg '
        'inner join sim_group sg on sg.id = rg.sim_group_id '
        'inner join run r on r.run_group_id = rg.id '
        'inner join sim s on s.id = r.sim_id '
        'left outer join run_group_sge rgsge on rgsge.run_group_id = rg.id '
        'where sg.id = %s and '
        'sg.user_account_id = %s'
    )
    args = (sim_group_id, user_id)

    if states:
        sql += ' and r.run_state_id in %s'
        args = (sim_group_id, user_id, states)

    sql += ' order by rg.order_in_sim_group'

    run_groups = exec_select(sql, args, fetch_all=True)

    return run_groups


def get_run_group_runs(run_group_id, states=None, user_cred=None):
    """Get runs belonging to a given run group, optionally in given states."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select r.* '
        'from run r '
        'inner join run_group rg on rg.id = r.run_group_id '
        'inner join sim_group sg on sg.id = rg.sim_group_id '
        'where run_group_id = %s '
        'and sg.user_account_id = %s'
    )
    args = [run_group_id, user_id]

    if states:
        sql += ' and r.run_state_id in %s'
        args += [states]

    sql += ' order by r.order_in_run_group'

    runs = exec_select(sql, tuple(args), fetch_all=True)

    if not runs:
        msg = 'No runs in run group with ID {} is associated with this user.'
        raise ValueError(msg.format(run_group_id))

    return runs


def get_run_group(run_group_id, user_cred=None):
    """Get a run group."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select rg.* '
        'from run_group rg '
        'inner join sim_group sg on sg.id = rg.sim_group_id '
        'where rg.id = %s and '
        'sg.user_account_id = %s'
    )
    run_group = exec_select(sql, (run_group_id, user_id))

    if not run_group:
        msg = 'No run group with ID {} is associated with this user.'
        raise ValueError(msg.format(run_group_id))

    return run_group


def get_sim_group_run_groups(sim_group_id, user_cred=None):
    """Get all run groups associated with a given sim group."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select rg.*, rgsge.job_array "run_group_sge_job_array", '
        'rgsge.selective_submission "run_group_sge_selective_submission" '
        'from run_group rg '
        'inner join sim_group sg on sg.id = rg.sim_group_id '
        'left outer join run_group_sge rgsge on rgsge.run_group_id = rg.id '
        'where sg.id = %s and '
        'sg.user_account_id = %s '
        'order by rg.order_in_sim_group'
    )
    run_groups = exec_select(sql, (sim_group_id, user_id), fetch_all=True)

    return run_groups


def get_count_sim_group_run_groups(sim_group_id, user_cred=None):
    """Return number of run groups for a given sim group."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select count(*) from run_group rg '
        'inner join sim_group sg on sg.id = rg.sim_group_id '
        'where sg.id = %s '
        'and sg.user_account_id = %s'
    )
    count = exec_select(sql, (sim_group_id, user_id))['count(*)']

    return count


def add_run_group(sim_group_id, run_group_defn, sim_insert_ids, run_order_in_sims,
                  run_state, is_sge, user_cred=None):
    """Add a run group to a given sim group."""

    # First get existing run groups:
    pre_rgs = get_sim_group_run_groups(sim_group_id, user_cred)
    order_in_sim_group = len(pre_rgs)

    # Add to run_group table:
    sql = (
        'insert into run_group '
        '(sim_group_id, software_instance_id, num_cores, order_in_sim_group, '
        'auto_submit, auto_process) '
        'values (%s, %s, %s, %s, %s, %s)'
    )

    rg_insert_id = exec_insert(sql, (
        sim_group_id,
        run_group_defn['software_instance']['id'],
        run_group_defn['num_cores'],
        order_in_sim_group,
        _AUTO_LOOKUP[run_group_defn['auto_submit']],
        _AUTO_LOOKUP[run_group_defn['auto_process']]
    ))

    rg_sge_insert_id = None
    if is_sge:
        # Also add to run_group_sge table:
        sql_sge = (
            'insert into run_group_sge '
            '(run_group_id, job_array, selective_submission, '
            'resource_flag) '
            'values (%s, %s, %s, %s)'
        )
        rg_sge_insert_id = exec_insert(sql_sge, (
            rg_insert_id,
            run_group_defn['job_array'],
            run_group_defn['selective_submission'],
            run_group_defn.get('resource_flag')
        ))

    # Add to the run table
    run_insert_ids = []
    for order_in_run_group, _ in enumerate(run_group_defn['sim_idx']):

        sim_id = sim_insert_ids[order_in_run_group]
        order_in_sim = run_order_in_sims[order_in_run_group]

        sql_run = (
            'insert into run (sim_id, run_group_id, run_state_id, '
            'order_in_run_group, order_in_sim) '
            'values (%s, %s, %s, %s, %s)'
        )
        args = (
            sim_id, rg_insert_id, run_state, order_in_run_group, order_in_sim
        )
        rn_insert = exec_insert(sql_run, args)

        rn_sge_insert = None
        if is_sge:

            # Also add to run_sge table:
            sql_run_sge = ('insert into run_sge (run_id) values (%s)')
            rn_sge_insert = exec_insert(sql_run_sge, (rn_insert))

        run_insert_ids.append([rn_insert, rn_sge_insert])

    return rg_insert_id, rg_sge_insert_id, run_insert_ids


def check_run_group_ownership(run_group_id, user_cred=None):
    """Check given user owns given run group."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select * '
        'from run_group rg '
        'inner join sim_group sg on sg.id = rg.sim_group_id '
        'inner join user_account u on u.id = sg.user_account_id '
        'where rg.id = %s and '
        'u.id = %s'
    )
    check_res = exec_select(sql, (run_group_id, user_id))
    if not check_res:
        msg = 'No run group with ID {} is associated with this user.'
        raise ValueError(msg.format(run_group_id))


def check_many_runs_ownership(run_ids, user_cred=None):
    """Check given user owns all in a set of runs."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select r.id '
        'from run r '
        'inner join run_group rg on rg.id = r.run_group_id '
        'inner join sim_group sg on sg.id = rg.sim_group_id '
        'inner join user_account u on u.id = sg.user_account_id '
        'where r.id in %s and '
        'u.id = %s'
    )
    check_res = exec_select(sql, (run_ids, user_id), fetch_all=True)

    owned_run_ids = sorted([i['id'] for i in check_res])
    if sorted(run_ids) != owned_run_ids:
        msg = 'Not all runs with IDs {} are associated with this user.'
        raise ValueError(msg.format(run_ids))


def check_run_ownership(run_id, user_cred=None):
    """Check given user owns given run."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select * '
        'from run r '
        'inner join run_group rg on rg.id = r.run_group_id '
        'inner join sim_group sg on sg.id = rg.sim_group_id '
        'inner join user_account u on u.id = sg.user_account_id '
        'where r.id = %s and '
        'u.id = %s'
    )
    check_res = exec_select(sql, (run_id, user_id))
    if not check_res:
        msg = 'No run with ID {} is associated with this user.'
        raise ValueError(msg.format(run_id))


def set_run_state(run_id, state_id, user_cred=None):
    """Set the run state for a given run."""

    user_cred = user_cred or CONFIG['user']
    _ = get_user_id(user_cred)

    check_run_ownership(run_id)

    sql_runs = (
        'update run '
        'set run_state_id = %s '
        'where id = %s'
    )
    exec_update(sql_runs, (state_id, run_id))


def set_many_run_states(run_ids, state_id, user_cred=None):
    """Set many runs to a given state."""

    if not run_ids:
        return

    user_cred = user_cred or CONFIG['user']
    _ = get_user_id(user_cred)

    check_many_runs_ownership(run_ids)

    sql_runs = (
        'update run '
        'set run_state_id = %s '
        'where id in %s'
    )
    exec_update(sql_runs, (state_id, run_ids))


def set_all_run_states(run_group_id, state_id, user_cred=None):
    """Change the run state for all runs within a run group."""

    user_cred = user_cred or CONFIG['user']
    _ = get_user_id(user_cred)

    check_run_group_ownership(run_group_id)

    sql_runs = (
        'update run '
        'set run_state_id = %s '
        'where run_group_id = %s'
    )
    exec_update(sql_runs, (state_id, run_group_id))


def set_run_group_submitted(run_group_id, hostname, submit_time, run_state_id, user_cred=None):
    """Set a run group to the submitted state."""

    user_cred = user_cred or CONFIG['user']
    _ = get_user_id(user_cred)

    check_run_group_ownership(run_group_id)

    if run_state_id not in [3, 5]:
        msg = 'Run state should be set to 3 ("in_queue") or 5 ("running")'
        raise ValueError(msg)

    sql = (
        'update run_group '
        'set submit_time = %s, hostname = %s '
        'where id = %s'
    )
    exec_update(sql, (submit_time, hostname, run_group_id))
    set_all_run_states(run_group_id, run_state_id)


def set_run_group_sge_jobid(run_group_id, job_id, user_cred=None):
    """Add job ID to a run_group_sge row."""

    user_cred = user_cred or CONFIG['user']
    _ = get_user_id(user_cred)

    check_run_group_ownership(run_group_id)

    sql = (
        'update run_group_sge '
        'set job_id = %s '
        'where run_group_id = %s'
    )
    exec_update(sql, (job_id, run_group_id))


def check_scratch_ownership(scratch_id, user_cred=None):
    """Check given user owns given scratch."""

    get_scratch_by_id(scratch_id, user_cred)


def check_resource_ownership(resource_id, user_cred=None):
    """Check given resource ID owned by given user."""

    get_resource(resource_id, user_cred)


def get_runs_by_scratch(scratch_id, run_state_ids=None, user_cred=None):
    """Get all runs on a given scratch for this user, optionally filtered by
    run state.

    """

    user_cred = user_cred or CONFIG['user']
    _ = get_user_id(user_cred)

    check_scratch_ownership(scratch_id, user_cred)

    sql = (
        'select r.*, rgs.job_id "run_group_sge_job_id" '
        'from run r '
        'inner join run_group rg on rg.id = r.run_group_id '
        'inner join sim_group sg on sg.id = rg.sim_group_id '
        'left outer join run_group_sge rgs on rgs.run_group_id = rg.id '
        'where sg.scratch_id = %s'
    )

    if run_state_ids is None:
        run_state_ids = []

    if run_state_ids:
        sql += ' and r.run_state_id in %s'
        args = (scratch_id, run_state_ids)

    else:
        args = (scratch_id,)

    runs = exec_select(sql, args, fetch_all=True)
    return runs


def check_archive_started(sim_group_id, user_cred=None):
    """Check if archiving has started for given sim group and user."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'select archive_started '
        'from sim_group '
        'where user_account_id = %s and '
        'id = %s'
    )
    started = exec_select(sql, (user_id, sim_group_id))['archive_started']

    return bool(started)


def set_archive_started(sim_group_id, user_cred=None):
    """Check if archiving has started for given sim group and user."""

    user_cred = user_cred or CONFIG['user']
    user_id = get_user_id(user_cred)

    sql = (
        'update sim_group '
        'set archive_started = %s '
        'where user_account_id = %s and '
        'id = %s'
    )
    exec_update(sql, (1, user_id, sim_group_id))
