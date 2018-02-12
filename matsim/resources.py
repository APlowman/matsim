"""`matsim.resources.py`"""
import os
import pathlib
import subprocess
from datetime import datetime
import copy
import shutil

from matsim.copy_multi import copytree_multi
from matsim import utils, database, dbhelpers as dbh, CONFIG
from matsim.utils import prt


class Resource(object):
    """Class to represent a location on the local or a remote machine."""

    def __init__(self, resource_id, add_path=None):
        """Initialise resource object from ID and 'database' file."""

        res_defn = database.get_resource(resource_id)

        if add_path is None:
            add_path = []

        # Instantiate a "pure" path object in this case, which does not have
        # access to the file system.
        if res_defn['machine_os_type'] == 'nt':
            path_class = pathlib.PureWindowsPath

        elif res_defn['machine_os_type'] == 'posix':
            path_class = pathlib.PurePosixPath

        self.base_path = path_class(res_defn['resource_base_path'])
        self.path = self.base_path.joinpath(*add_path)

        # Check base_path is absolute
        if not self.base_path.is_absolute():
            msg = ('Resource `base_path` "{}" must be an absolute path.')
            raise ValueError(msg.format(self.base_path))

        self.resource_id = res_defn['resource_id']
        self.machine_id = res_defn['machine_id']
        self.machine_name = res_defn['machine_name']
        self.os_type = res_defn['machine_os_type']
        self.is_dropbox = res_defn['machine_is_dropbox']
        self.sync_client_paths = None

        if self.is_dropbox:

            sync_client_data = database.get_machine_sync_client_path(
                self.machine_id)

            if sync_client_data:

                # Use local files from sync client rather than Dropbox API:
                sync_client, mach_name = sync_client_data

                print('Using local sync client for Dropbox resource.')

                sync_client_path_parts = pathlib.Path(
                    sync_client['sync_path']).parts

                self.is_dropbox = False
                self.os_type = os.name
                self.machine_id = sync_client['sync_client_machine_id']
                self.machine_name = mach_name

                self.base_path = pathlib.Path(*(sync_client_path_parts +
                                                self.base_path.parts[1:]))

                self.path = pathlib.Path(*(sync_client_path_parts +
                                           self.path.parts[1:]))

        try:
            self.make_paths_concrete()
        except NotImplementedError:
            pass

    def make_paths_concrete(self):
        """Convert `base_path` and `path` to "concrete" path objects.

        Convert to "concrete" path objects (as opposed to "pure") to allow
        file operations.

        """

        if self.os_type == 'nt':
            path_class = pathlib.WindowsPath

        elif self.os_type == 'posix':
            path_class = pathlib.PosixPath

        self.base_path = path_class(self.base_path)
        self.path = path_class(self.path)


class Stage(Resource):
    """Class to represent the area on the local machine in which simulations
    input files are generated."""

    def __init__(self, name, add_path=None):

        stage_defn = database.get_stage_by_name(name)
        res_id = stage_defn['resource_id']
        super().__init__(res_id, add_path)
        self.name = name
        self.stage_id = stage_defn['stage_id']


class Scratch(Resource):
    """Class to represent the area on a machine in which simulations are to be
    run."""

    def __init__(self, name, add_path=None):

        scratch_defn = database.get_scratch_by_name(name)
        res_id = scratch_defn['resource_id']
        super().__init__(res_id, add_path)
        self.name = name
        self.scratch_id = scratch_defn['scratch_id']
        self.sge = scratch_defn['scratch_is_sge']


class Archive(Resource):
    """Class to represent the area on a machine in which completed simulations
    are archived."""

    def __init__(self, name, add_path=None):

        arch_defn = database.get_archive_by_name(name)
        res_id = arch_defn['resource_id']
        super().__init__(res_id, add_path)
        self.name = name
        self.archive_id = arch_defn['archive_id']


class ResourceConnection(object):
    """Class to represent a connection between a local resource and another."""

    @classmethod
    def check_exists(cls, src, dst):
        """Check a resource connection between src and dst exists on the 
        database, and return the result if it does. Otherwise raise ValueError.

        """

        if src == dst:
            res_conn = {
                'host': None,
                'is_remote': False,
            }

        else:
            res_conn = database.get_resource_connection(
                src.resource_id, dst.resource_id)

            if res_conn is None:
                msg = ('No resource connection information between source'
                       'and destination resources can be found.')
                raise ValueError(msg)

        return res_conn

    def __init__(self, src, dst):

        res_conn = ResourceConnection.check_exists(src, dst)

        self.src = src
        self.dst = dst
        self.host = res_conn['host']
        self.remote = src.machine_id != dst.machine_id
        self.os_types = (self.src.os_type, self.dst.os_type)

        if self.remote:

            # These (src, dst) OS types are allowed for remote dst:
            ok_remote_os_types = {
                ('nt', 'posix'),
                ('posix', 'posix'),
            }
            if self.os_types not in ok_remote_os_types:
                msg = ('This combination of source and destination '
                       '`os_type` {} is not supported for remote destination.')
                raise ValueError(msg.format(self.os_types))

        # Source is always local:
        self.src.make_paths_concrete()

        # Destination may be remote:
        if not self.remote:
            self.dst.make_paths_concrete()

    def check_conn(self):
        """Check a connection can be made between source and destination
        resources, and that the base paths of both source and destination exist.

        TODO: separate remote connection and destination base path check.

        """

        dst_msg = 'Destination `base_path` "" does not exist.'

        # Check source base path exists:
        if not self.src.base_path.exists():
            msg = 'Source `base_path` "{}" does not exist.'
            raise ValueError(msg.format(self.src.base_path))

        if self.remote:

            if self.dst.is_dropbox:
                dbx = dbh.get_dropbox()
                if not dbh.is_folder(dbx, str(self.dst.base_path)):
                    raise ValueError(dst_msg.format(self.dst.base_path))

            else:
                # Check remote connection can be made and destination base path exists:
                ssh_cm = 'ssh {} "[ -d {} ]"'
                ssh_cm = ssh_cm.format(self.host, self.dst.base_path)
                comp_proc = subprocess.run(['bash', '-c', ssh_cm])

                if comp_proc.returncode == 1:
                    msg = ('Remote connection to host "{}" could not be made, '
                           'or destination `base_path` "{}" does not exist')
                    raise ValueError(msg.format(self.host, self.dst.base_path))

        else:

            # Check destination base path exists:
            if not self.dst.base_path.exists():
                raise ValueError(dst_msg.format(self.dst.base_path))

    def file_exist_on_src(self, subpath=None):
        """TODO"""
        pass

    def file_exist_on_dst(self, subpath=None):
        """TODO"""
        pass

    def copy_to_dest(self, subpath=None, file_backup=False, ignore=None):
        """Copy content from the source resource to the destination resource.

        Contents is mirrored from source to destination.

        Parameters
        ----------
        subpath : list, optional
            Copy from a given subpath with source (to the same subpath on
            destination). Default is copy from the source `path`.
        file_backup : bool, optional
            Only applicable if `subpath` resolves to a file and the file path
            already exists on destination. If True, the name of the file on
            destination is prepended with a date time stamp.
        ignore : sequence of str
            Only applicable if `subpath` resolves to a directory. Patterns to
            ignore when copying (glob style).

        TODO: support Dropbox destination

        """

        self.check_conn()

        msg = ('\nCopying from resource "{}" to{}{} resource "{}".\n')
        is_rem_str = ' remote' if self.remote else ''
        is_db_str = ' Dropbox' if self.dst.is_dropbox else ''
        msg = msg.format(self.src.name, is_rem_str, is_db_str, self.dst.name)
        msg += '-' * (len(msg) - 2)
        print(msg)

        if subpath:
            src_path = self.src.path.joinpath(*subpath)
            dst_path = self.dst.path.joinpath(*subpath)
        else:
            src_path = self.src.path
            dst_path = self.dst.path

        if not src_path.exists():
            msg = 'Path {} does not exist.'
            raise ValueError(msg.format(src_path))

        msg = 'Source:      "{}"\nDestination: "{}"\n'
        msg += 'File backup: {}\n'
        msg += 'Is file:     {}\n'
        msg += 'Is dir:      {}\n'
        msg_fmt = [src_path, dst_path, file_backup,
                   src_path.is_file(), src_path.is_dir()]
        if ignore:
            msg += 'Excluding:   {}\n'
            msg_fmt += [','.join(['"{}"'.format(i) for i in ignore])]
        print(msg.format(*msg_fmt))

        if file_backup:
            dt_fmt = '%Y-%m-%d-%H%M%S.'
            bk_ts = datetime.strftime(datetime.now(), dt_fmt)
            bk_dst_path = str(dst_path.with_name(bk_ts + dst_path.name))

        if self.remote and self.dst.is_dropbox:

            dbx = dbh.get_dropbox()
            dst_path = str(dst_path)

            if src_path.is_file():

                src_path = str(src_path)

                if file_backup and dbh.is_file(dbx, dst_path):
                    # Rename destination file
                    dbh.rename_file(dbx, dst_path, bk_dst_path)

                dbh.upload_dropbox_file(dbx, src_path, dst_path)

            else:

                src_path = str(src_path)
                dbh.upload_dropbox_dir(dbx, src_path, dst_path, exclude=ignore)

        elif self.remote:

            orig_src_path = copy.copy(src_path)

            if self.src.os_type == 'nt':

                # Convert path to posix style for use within "Bash on Windows":
                path_args = ['/mnt', src_path.drive[0].lower(),
                             *src_path.parts[1:]]
                src_path = pathlib.PurePosixPath(*path_args)

            if orig_src_path.is_file():

                mkdirs = False

                # Check if file exists on destination:
                cmd = '[ -f {} ] && echo found'.format(dst_path)
                check_file_exist = self.run_command([cmd]).rstrip('\n')

                if check_file_exist == 'found' and file_backup:

                    # Rename on destination:
                    cmd = 'mv {} {}'.format(dst_path, bk_dst_path)
                    self.run_command([cmd])

            else:
                # Add trailing slash to path (for rsync to copy contents of src):
                src_path = str(src_path) + '/'
                mkdirs = True

            utils.rsync_remote(src_path, self.host, dst_path,
                               mkdirs=mkdirs, exclude=ignore)

        else:

            if src_path.is_file():

                if dst_path.is_file() and file_backup:

                    # Rename destination file
                    shutil.move(str(dst_path), bk_dst_path)

                # copy2 (and copy) will overwrite existing file:
                shutil.copy2(src_path, dst_path)

            else:
                if ignore:
                    ignore_func = shutil.ignore_patterns(*ignore)
                else:
                    ignore_func = None
                copytree_multi(src_path, dst_path, ignore=ignore_func)

    def run_command(self, cmd, cwd=None, block=True):
        """Execute a command on the destination resource.

        Parameters
        ----------
        cmd : list
            The first element is the command to execute and additional elements
            are optional arguments for that command.
        cwd : str, optional
            The directory in which to execute the command, relative to the
            `path` of the destination resource.

        """

        if self.dst.is_dropbox:
            msg = 'Cannot run a command on a Dropbox resource.'
            raise NotImplementedError(msg)

        self.check_conn()

        # no_dir_msg = 'Directory does not exist on scratch. Aborting.'

        if cwd is not None:
            cwd = str(self.dst.path.joinpath(cwd))
        else:
            cwd = str(self.dst.path)

        if self.remote:

            # Concatentate the command, since we're passing via SSH:
            cmd_str = ' '.join(cmd)

            run_args = ['bash', '-c']
            ssh_cmd = 'ssh {} "cd {} && {}"'
            run_args.append(ssh_cmd.format(self.host, cwd, cmd_str))
            cmd = run_args

        else:

            # Change to desired directory:
            os.chdir(cwd)

            if self.dst.os_type == 'nt':
                pass

            elif self.dst.os_type == 'posix':
                pass

        run_params = {
            'args': cmd,
            'encoding': 'utf-8',
            'stdout': subprocess.PIPE,
            'stderr': subprocess.PIPE,
        }

        if block:

            cmd_proc = subprocess.run(**run_params)
            if cmd_proc.stderr:
                msg = ('ResourceConnection could not run command '
                       'on destination: {}'.format(cmd_proc.stderr))
                raise ValueError(msg)
            else:
                return cmd_proc.stdout

        else:
            cmd_proc = subprocess.Popen(**run_params)

            return cmd_proc

        #     if self.os_name == 'nt' and scratch.os_name == 'nt':
        #         if not os.path.isdir(scratch.path):
        #             raise ValueError(no_dir_msg)

        #         js_path = os.path.join(scratch.path, 'jobscript.bat')
        #         # Run batch file in a new console window:
        #         subprocess.Popen(js_path,
        #                          creationflags=subprocess.CREATE_NEW_CONSOLE)

        #     elif self.os_name == 'posix' and scratch.os_name == 'posix':
        #         # Use rsync/scp
        #         js_path = os.path.join(scratch.path, 'jobscript.sh')
        #         os.chmod(js_path, 0o744)
        #         subprocess.Popen(js_path, shell=True)
