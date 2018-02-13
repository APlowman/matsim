"""matsim.atomistic.simulation.castep.py"""

import copy
import json

import numpy as np

from matsim import utils, database_atm as dbs_atm
from matsim.atomistic.software import castep as castepio
from matsim.atomistic.simulation import SUPERCELL_TYPE_LOOKUP
from matsim.atomistic.simulation.sim import AtomisticSimulation


class CastepSimulation(AtomisticSimulation):
    """Class to represent a CASTEP simulation."""

    def __init__(self, options):
        """Initialise a CastepSimulation."""
        super().__init__(options)
        self._process_options()

    def _process_options(self):
        """Additional processing on CASTEP options to prepare for writing input
        files.

        """
        super()._process_options()

        castep_opt = self.options['params']['castep']

        # Sort out checkpointing:
        if castep_opt.get('checkpoint') is True:
            if castep_opt.get('backup_interval') is not None:
                castep_opt['param'].update(
                    {'backup_interval': castep_opt['backup_interval']})

        else:
            castep_opt['param'].update({'write_checkpoint': 'none'})

        castep_opt.pop('backup_interval', None)
        castep_opt.pop('checkpoint', None)

        # Remove geometry optimisation parameters if not necessary:
        castep_task = castep_opt['param']['task'].upper()
        geom_opt_str = ['GEOMETRYOPTIMISATION', 'GEOMETRYOPTIMIZATION']
        if castep_task not in geom_opt_str:

            geom_keys = []
            for param_key in castep_opt['param']:
                if param_key.upper().startswith('GEOM_'):
                    geom_keys.append(param_key)

            for geom_k in geom_keys:
                castep_opt['param'].pop(geom_k, None)

        # Remove constraints if task is single point:
        if castep_task == 'SINGLEPOINT':

            constraints = self.options['structure']['constraints']
            cell_const = constraints['cell']
            atom_const = constraints['atoms']

            cell_const.pop('cell_angles_equal', None)
            cell_const.pop('cell_lengths_equal', None)
            cell_const.pop('fix_cell_angles', None)
            cell_const.pop('fix_cell_lengths', None)
            atom_const.pop('fix_xy_idx', None)
            atom_const.pop('fix_xz_idx', None)
            atom_const.pop('fix_yz_idx', None)
            atom_const.pop('fix_xyz_idx', None)

        # Add symmetry operations:
        if castep_opt['find_inv_sym']:

            if castep_opt['cell'].get('symmetry_generate') is True:
                msg = ('Cannot add inversion symmetry operation to CASTEP '
                       'CELL file if `symmetry_generate` is `True`.')
                raise ValueError(msg)

            sym_ops = self.structure.get_sym_ops()

            sym_rots = sym_ops['rotations']
            sym_trans = sym_ops['translations']
            inv_sym_rot = -np.eye(3, dtype=int)
            inv_sym_idx = np.where(
                np.all(sym_rots == inv_sym_rot, axis=(1, 2)))[0]

            if not inv_sym_idx:
                msg = 'The supercell does not have inversion symmetry.'
                raise ValueError(msg)

            inv_sym_trans = sym_trans[inv_sym_idx[0]]

            castep_opt['sym_ops'] = [
                np.vstack([np.eye(3), np.zeros((3,))]),
                np.vstack([inv_sym_rot, inv_sym_trans])
            ]

    def write_input_files(self, run_idx, path):
        """Write input files necessary to perform a CASTEP simulation. """

        run_parameters = super().get_run_parameters(run_idx)
        common_params = super().get_common_atomistic_parameters()
        cst_in_params = {
            'path': path,
            'seedname': run_parameters['castep']['seedname'],
            'cell': run_parameters['castep']['cell'],
            'param': run_parameters['castep']['param'],
            # 'sym_ops': run_parameters['castep']['sym_ops'],
            **common_params
        }
        castepio.write_castep_inputs(**cst_in_params)

    def to_jsonable(self):
        """Generate a dict representation that can be JSON serialised."""

        # `structure` and `options` are dealt with in the super-class:
        ret = super().to_jsonable()

        # Add `result` key from each run in `runs`:
        for i, j in zip(ret['runs'], self.runs):
            if j['result'] is not None:
                i['result'] = {
                    # TODO, check which CASTEP results vals need `tolist()`ing.
                    **j['result'],
                }

        return ret

    @classmethod
    def from_jsonable(cls, state):
        """Instantiate from a JSONable dict."""

        runs_native = state['runs']
        for idx, _ in enumerate(runs_native):
            if runs_native[idx]['result'] is not None:
                runs_native[idx]['result'] = {
                    # TODO, check which CASTEP results vals need `array`ing.
                    **runs_native[idx]['result'],
                }

        sup_types_str = state['structure']['meta'].get('supercell_type')
        sup_type_class = 'default'
        for sup_type in sup_types_str:
            if sup_type in SUPERCELL_TYPE_LOOKUP:
                sup_type_class = sup_type
                break
        struct_class = SUPERCELL_TYPE_LOOKUP[sup_type_class]

        state.update({
            'structure': struct_class.from_jsonable(state['structure']),
            'runs': runs_native,
        })
        return cls(state=state)

    def check_success(self, path):
        """Check a given run of this simulation has succeeded."""

        castep_opt = self.options['params']['castep']
        castep_task = castep_opt['param']['task'].upper()
        geom_opt_str = ['GEOMETRYOPTIMISATION', 'GEOMETRYOPTIMIZATION']

        if castep_task in geom_opt_str:
            success_func = castepio.check_success_geom_opt
        else:
            success_func = castepio.check_success_single_point

        return success_func(path)

    def parse_result(self, path, run_idx):
        """Parse results from path and add to runs[run_idx]['result']"""

        cst = castepio.read_castep_output(path)
        self.runs[run_idx]['result'] = cst

        cst_params = cst['params']

        cst_sim_param_names = [
            'cut_off_energy',
            'kpoint_num',
            'kpoint_mp_grid',
            'kpoint_mp_offset',
            'xc_functional',
            'mixing_scheme',
            'metals_method',
            'smearing_width',
            'geom_method',
        ]
        cst_run_param_names = [
            'elec_energy_tol',
            'opt_strategy',
            'geom_energy_tol',
            'geom_force_tol',
            'geom_disp_tol',
            'geom_stress_tol',
        ]
        cst_json_param_names = [
            'kpoint_mp_grid',
            'kpoint_mp_offset',
        ]

        # Need to add to table `atm_run_castep`, and if not already
        # added, also `atm_sim`, `atm_sim_castep`

        atm_sim = dbs_atm.get_atm_sim_by_sim_id(self.dbid)

        if atm_sim:

            atm_sim_id = atm_sim['id']
            atm_sim_castep_id = dbs_atm.get_atm_sim_castep_by_atm_sim_id(
                atm_sim_id
            )['id']

        else:

            # Supercell type
            sup_types = self.structure.meta['supercell_type']
            if 'bicrystal' in sup_types:
                supercell_type = 'gb'
            elif 'bulk_bicrystal' in sup_types:
                supercell_type = 'gb_bulk'
            elif 'surface_bicrystal' in sup_types:
                supercell_type = 'gb_surface'
            elif 'bulk' in sup_types:
                supercell_type = 'bulk'
            elif 'surface' in sup_types:
                supercell_type = 'surface'

            # Unit repeats in the supercell
            repeats = None
            if supercell_type == 'bulk':
                repeats = self.options['structure']['box_lat']
            elif supercell_type.startswith('gb'):
                # For a CSL GB, this is `gb_size` * identity matrix
                repeats = self.options['structure']['gb_size']

            # Atom constraints. Not sure if this is right, but we're saying 3
            # are for CoM:
            atm_constraints = (cst['params']['ion_constraints_num'] > 3)

            # Add to `atm_sim` table:
            atm_sim_params = {
                'supercell_type': supercell_type,
                'repeats': json.dumps(repeats),
                'num_atoms': cst['num_ions'],
                'path': path,
                'atom_constraints': atm_constraints,
                'software': 'castep'
            }
            atm_sim_id = dbs_atm.add_atm_sim(atm_sim_params, self.dbid)

            # Add to `atm_sim_castep` table:
            castep_sim_params = {
                'task': cst_params['calc_type']
            }
            for par_name in cst_sim_param_names:

                par = cst_params[par_name]
                if par_name in cst_json_param_names:
                    par = json.dumps(par)
                castep_sim_params.update({par_name: par})

            atm_sim_castep_id = dbs_atm.add_atm_sim_castep(
                atm_sim_id, castep_sim_params)

        # Add to the `atm_run_castep` table
        castep_run_params = {}
        for par_name in cst_run_param_names:

            par = cst_params[par_name]
            if par_name in cst_json_param_names:
                par = json.dumps(par)

            castep_run_params.update({par_name: par})

        dbs_atm.add_atm_run_castep(
            atm_sim_castep_id, castep_run_params, self.runs[run_idx]['dbid'])
