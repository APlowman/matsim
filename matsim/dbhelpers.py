"""`matsim.dbhelpers.py

Dropbox helper functions.

"""
import dropbox
import os
import posixpath
import fnmatch
from matsim import CONFIG


def get_dropbox(key=None):
    if not key:
        key = CONFIG['dropbox']
    return dropbox.Dropbox(key)


def is_file(dbx, path):
    """Check given path on dropbox is a file."""
    meta = dbx.files_get_metadata(path)
    return isinstance(meta, dropbox.files.FileMetadata)


def is_folder(dbx, path):
    """Check given path on dropbox is a folder."""
    try:
        meta = dbx.files_get_metadata(path)
        return isinstance(meta, dropbox.files.FolderMetadata)
    except dropbox.exceptions.ApiError as api_err:
        raise ValueError('Path not found on Dropbox: {}'.format(path))


def rename_file(dbx, src_path, dst_path):
    """Rename a file from one path to another."""
    if is_file(dbx, src_path):
        print('Renaming file: {} to {}'.format(src_path, dst_path))
        dbx.files_move_v2(src_path, dst_path)
    else:
        raise ValueError('Cannot rename a file that does not exist.')


def download_dropbox_file(dbx, dropbox_path, local_path):
    dbx.files_download_to_file(local_path, dropbox_path)


def upload_dropbox_file(dbx, local_path, dropbox_path, overwrite=False,
                        autorename=False):
    """
    Parameters
    ----------
    dbx: Dropbox
    local_path : str
        Path of file on local computer to upload to dropbox.
    dropbox_path : str
        Directory on dropbox to upload the file to.
    overwrite : bool
        If True, the file overwrites an existing file with the same name.
    autorename : bool
        If True, rename the file if there is a conflict.

    """

    if overwrite:
        mode = dropbox.dropbox.files.WriteMode('overwrite', None)
    else:
        mode = dropbox.dropbox.files.WriteMode('add', None)

    with open(local_path, mode='rb') as f:

        dbx.files_upload(f.read(), dropbox_path,
                         mode=mode, autorename=autorename)


def upload_dropbox_dir(dbx, local_path, dropbox_path, overwrite=False,
                       autorename=False, exclude=None):
    """
    Parameters
    ----------
    dbx: Dropbox
    local_path : str
        Path of file on local computer to upload to dropbox.
    dropbox_path : str
        Directory on dropbox to upload the file to.
    overwrite : bool
        If True, the file overwrites an existing file with the same name.
    autorename : bool
        If True, rename the file if there is a conflict.
    exclude : list, optional
        List of file or directory names to exclude, matched with `fnmatch` for
        files, or compared directly for directories.

    Notes
    -----
    Does not upload empty directories.

    """

    # Validation
    if not os.path.isdir(local_path):
        raise ValueError(
            'Specified `local_path` is not a directory: {}'.format(local_path))

    for root, dirs, files in os.walk(local_path):

        dirs[:] = [d for d in dirs if d not in exclude]

        print('Uploading from root directory: {}'.format(root))

        for file_name in files:

            up_file = False

            if exclude is not None:
                if not any([fnmatch.fnmatch(file_name, i) for i in exclude]):
                    up_file = True
            else:
                up_file = True

            if up_file:

                # Source path
                src_fn = os.path.join(root, file_name)

                # Destination path
                rel_path = os.path.relpath(src_fn, local_path)
                fn_db_path = os.path.join(dropbox_path, rel_path)
                dst_fn = posixpath.join(*fn_db_path.split(os.path.sep))
                if not dst_fn.startswith('/'):
                    dst_fn = '/' + dst_fn

                print('Uploading file: {}'.format(file_name))
                upload_dropbox_file(dbx, src_fn, dst_fn,
                                    overwrite=overwrite, autorename=autorename)
