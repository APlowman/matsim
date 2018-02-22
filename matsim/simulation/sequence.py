"""Module containing a class to represent a sequence of simulations."""

import warnings
import copy
import numpy as np
import yaml

from matsim import SEQ_DEFN
from matsim.utils import (set_nested_dict, get_recursive, update_dict, prt,
                          mut_exc_args, merge, nest)
from matsim import readwrite
from matsim.simulation import BaseUpdate
from matsim.simulation import sequence_funcs

SEQUENCE_FUNC_LOOKUP = {
    'kpoint_mp_spacing': sequence_funcs.process_kpoint_sequence,
}


def merge_parallel_seqences(sequences):
    """Get sequence updates, merging those from parallel sequences (those with
    the same `nest_idx`), and sort by `nest_idx`.

    TODO: check if i need to sort, probably not now since i've sorted sequences
    """

    # prt(sequences, 'sequences')

    seq_upds_mergd = {}
    for seq in sequences:

        # prt(seq, 'seq')

        nest_idx = seq.nest_idx
        upds = seq.updates

        if nest_idx in seq_upds_mergd:
            # print('nest_idx in seq_upds')
            seq_upds_mergd[nest_idx] = merge(seq_upds_mergd[nest_idx], upds)
        else:
            # print('nest_idx NOT in seq_upds')
            # prt(upds, 'upds')
            seq_upds_mergd.update({nest_idx: upds})

        prt(readwrite.format_dict(seq_upds_mergd), 'seq_upds_mergd')

    # exit()

    seq_upds_srtd = [val for _, val in sorted(seq_upds_mergd.items())]

    return seq_upds_srtd


def get_sim_updates(seq_options, base_sims):
    """Get a list of updates required to form each simulation in the group."""

    all_sim_updates = []

    for idx, base_sim in enumerate(base_sims):

        if idx == 0:
            sequences = [SimSequence(i, base_sim) for i in seq_options]
            sequences.sort(key=lambda x: x.nest_idx)

        # prt(sequences, 'sequences')
        # exit()

        seq_upds = merge_parallel_seqences(sequences)
        grp_upd = nest(*seq_upds)

        # prt(seq_upds, 'seq_upds')

        sim_updates = []
        for upd_lst in grp_upd:
            upd_lst_flat = [j for i in upd_lst for j in i]
            sim_updates.append(upd_lst_flat)

        all_sim_updates.append(sim_updates)

    return all_sim_updates, sequences


class SimSequence(object):
    """Options which parameterise a sequence of simulations."""

    def __init__(self, spec, base_sim):

        params, spec = self._validate_spec(spec, SEQ_DEFN)

        # Store this for easily saving/loading as JSON:
        self.spec = copy.deepcopy(spec)

        self.base_dict = params['base_dict']
        self.range_allowed = params['range_allowed']
        self.update_mode = params['update_mode']
        self.val_seq_type = params['val_seq_type']
        self.map_to_dict = params['map_to_dict']
        self.val_name = params['val_name']
        self.affects_structure = params['affects_structure']

        self.name = spec.pop('name')
        self.nest_idx = spec.pop('nest_idx')
        self.val_fmt = spec.pop('val_fmt')
        self.path_fmt = spec.pop('path_fmt')

        # Remove vals from spec to parse and to get remaining additional spec:
        vals_spec_keys = ['vals', 'start', 'step', 'stop']
        vals_spec = {i: spec.pop(i, None) for i in vals_spec_keys}

        self.vals = self._parse_vals(**vals_spec)
        self.additional_spec = spec

        self.updates = self._get_updates(base_sim)

    def to_jsonable(self):
        """Generate a dict representation that can be JSON serialised."""
        return {'spec': self.spec}

    @classmethod
    def from_jsonable(cls, state, seqn_defn):
        """Generate new instance from JSONable dict"""
        return cls(state['spec'], seqn_defn)

    def _validate_spec(self, spec, seq_defn):
        """
        TODO: if `map_to_dict`, check `val_name` is not None.

        """
        msg = 'SimSequence failed validation: '

        # Keys allowed in the sequence spec (from makesims options):
        req_specs = [
            'name',
            'nest_idx',
        ]
        ok_specs = req_specs + [
            'val_fmt',
            'path_fmt',
            'vals',
            'start',
            'step',
            'stop',
        ]

        # Keys allowed in the sequence definition (from sequences.yml):
        req_params = [
            'base_dict',
            'range_allowed',
            'val_seq_type',
            'update_mode',
            'defaults',
            'additional_spec',
            'val_name',
            'map_to_dict',
            'affects_structure',
        ]
        ok_params = req_params

        for i in req_specs:
            if i not in spec:
                req_spec_msg = msg + 'Sequence spec '
                if spec.get('name'):
                    req_spec_msg += '"{}" '.format(spec.get('name'))
                req_spec_msg += 'must have `{}` key.'.format(i)
                raise ValueError(req_spec_msg)

        params = seq_defn[spec['name']]

        for i in req_params:
            if i not in params:
                raise ValueError(msg + 'Sequence definition must have '
                                 '`{}` key.'.format(i))

        for param_key in params:
            if param_key not in ok_params:
                raise ValueError(msg + 'Sequence definition key "{}" is not '
                                 'allowed in sequence "{}".'.format(param_key, spec['name']))

        spec = {**params['defaults'], **spec}

        for spec_key in spec.keys():
            if spec_key not in ok_specs and spec_key not in params['additional_spec']:
                ok_spec_msg = (msg + 'Sequence spec key "{}" is not allowed in'
                               ' sequence "{}".'.format(spec_key, spec['name']))
                raise ValueError(ok_spec_msg)

        if spec['name'] not in seq_defn:
            raise ValueError(
                msg + 'Sequence name "{}" not known. Sequence name should be '
                'one of: {}'.format(spec['name'], list(
                    seq_defn.keys()))
            )

        if not params['range_allowed']:
            if any([i in spec for i in ['start', 'step', 'stop']]):
                raise ValueError(msg + 'Range is not allowed for '
                                 'sequence name: {}'.format(spec['name']))
        return params, spec

    def _parse_vals(self, vals=None, start=None, step=None, stop=None):
        """Parse sequence spec vals and """

        # print('in parse_vals')

        mut_exc_args({'vals': vals},
                     {'start': start, 'step': step, 'stop': stop})

        if vals is None:

            if not self.range_allowed:
                raise ValueError('Specifying a range for sequence "{}" is not '
                                 'allowed.'.format(self.name))

            step_wrn_msg = ('Spacing between series values will not be exactly'
                            ' as specified.')

            if not np.isclose((start - stop) % step, 0):
                warnings.warn(step_wrn_msg)

            diff = start - stop if start > stop else stop - start
            num = int(np.round((diff + step) / step))
            vals = np.linspace(start, stop, num=num)

        # Parse vals Numpy array or tuple if necessary (lists are represented
        # natively by YAML):
        if self.val_seq_type in ['array', 'tuple']:

            vals_prsd = []
            for val in vals:

                if self.val_seq_type == 'array':
                    vals_prsd.append(np.array(val))

                elif self.val_seq_type == 'tuple':
                    vals_prsd.append(tuple(val))

            vals = vals_prsd

        return vals

    @property
    def num_vals(self):
        """Get the number of values (simulations) in the sequence."""
        return len(self.vals)

    def _get_updates(self, base_sim):
        """
        Build a list of update dicts to be applied to the base options for each
        element in the group.

        """

        # print('in _get_updates')

        # Run any additional processing on the `vals`
        func = SEQUENCE_FUNC_LOOKUP.get(self.name)
        if func:
            func(self, base_sim)

        fmt_arr_opt_path = {
            'col_delim': '_',
            'row_delim': '__',
            'format_spec': self.path_fmt,
        }
        fmt_arr_opt_val = {
            'col_delim': '_',
            'row_delim': '__',
            'format_spec': self.val_fmt,
        }

        name_add = ['sequence_id', 'names']
        paths_add = ['sequence_id', 'paths']
        vals_add = ['sequence_id', 'vals']
        nest_add = ['sequence_id', 'nest_idx']

        updates = []
        for _, val in enumerate(self.vals):

            # prt(val, 'val')

            if self.val_seq_type == 'array':
                path_str = readwrite.format_arr(
                    val, **fmt_arr_opt_path)[:-2]

            elif self.val_seq_type == 'list':
                path_str = '_'.join([self.path_fmt.format(i) for i in val])

            else:
                path_str = self.path_fmt.format(val)

            val_formatted = copy.deepcopy(val)
            if self.val_fmt:
                # Format the value for the actual update but not the sequence
                # id update:
                if self.val_seq_type == 'array':
                    val_formatted = readwrite.format_arr(
                        val_formatted, **fmt_arr_opt_val)[:-2]
                else:
                    val_formatted = self.val_fmt.format(val_formatted)

            # If `map_to_dict` replace the val which updates the options with a
            # dict mapping `val_name`: val and all other `additional_spec`:
            upd_val = val_formatted
            if self.map_to_dict:
                # TODO move this check to _validate_spec
                if not self.val_name:
                    msg = ('`val_name` must be set if `map_to_dict` True.')
                    raise ValueError(msg)

                upd_val = {
                    self.val_name: val_formatted,
                    **self.additional_spec,
                }

            # Update that affects the generation of the Simulation:
            elem_upd = BaseUpdate(self.base_dict, upd_val,
                                  self.val_seq_type, self.update_mode)

            # Updates that parameterise this effect:
            seqid_upd = [
                BaseUpdate(['sequence_id'], {}, None, 'replace'),
                BaseUpdate(name_add, self.name, None, 'append'),
                BaseUpdate(paths_add, path_str, None, 'append'),
                BaseUpdate(vals_add, val, self.val_seq_type, 'append'),
                BaseUpdate(nest_add, self.nest_idx, None, 'append'),
            ]

            updates.append([elem_upd] + seqid_upd)

        return updates
