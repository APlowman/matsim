"""matsim.simulation.__init__.py"""

from collections import namedtuple

from matsim.utils import get_recursive, set_nested_dict, update_dict

BaseUpdate = namedtuple(
    'BaseUpdate',
    ['address', 'val', 'val_seq_type', 'mode']
)


def apply_base_update(base_dict, update):
    """Apply a BaseUpdate namedtuple to a dict.

    Parameters
    ----------
    base_dict : dict
        dict to which an update is applied.
    update : BaseUpdate
        BaseUpdate namedtuple with fields: `address`, `val`, `val_seq_type`
        and `mode`.

    """
    if update.mode == 'replace':
        upd_val = update.val

    elif update.mode == 'append':
        upd_val = get_recursive(base_dict, update.address, [])
        upd_val.append(update.val)

    upd_dict = set_nested_dict(update.address, upd_val)
    new_base_dict = update_dict(base_dict, upd_dict)

    return new_base_dict
