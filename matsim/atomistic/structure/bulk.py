"""`matsim.atomistic.structure.bulk.py`"""
import copy

import numpy as np
from matsim import vectors

from matsim.atomistic.structure.atomistic import AtomisticStructure
from matsim.atomistic.structure.crystal import CrystalBox


class BulkCrystal(AtomisticStructure):
    """Class to represent a bulk crystal."""

    def __init__(self, as_params, repeats):
        """Constructor method for BulkCrystal object."""

        super().__init__(**as_params)
        self.repeats = repeats
        self.meta.update({'supercell_type': ['bulk']})


def bulk_crystal_from_parameters(crystal_structure, repeats, overlap_tol=1,
                                 tile=None):
    """Generate a BulkCrystal object given a `CrystalStructure` object and
    an integer array of column vectors representing the multiplicity of each new
    edge vector.

    Parameters
    ----------
    crystal_structure : CrystalStructure
    repeats : ndarray of int of shape (3, 3)
    tiles : sequence of int of length 3.

    """

    # Validation
    if any([i in vectors.num_equal_cols(repeats) for i in [2, 3]]):
        raise ValueError(
            'Identical columns found in repeats: \n{}\n'.format(repeats))

    supercell = np.dot(crystal_structure.bravais_lattice.vecs, repeats)

    crystal_box = CrystalBox(crystal_structure, supercell)
    atom_sites = crystal_box.atom_sites
    lattice_sites = crystal_box.lattice_sites

    crystal_idx_lab_atm = {
        'crystal_idx': (
            np.array([0]),
            np.zeros(atom_sites.shape[1], dtype=int)
        ),
    }
    crystal_idx_lab_lat = {
        'crystal_idx': (
            np.array([0]),
            np.zeros(lattice_sites.shape[1], dtype=int)
        ),
    }

    atom_labels = copy.deepcopy(crystal_box.atom_labels)
    atom_labels.update({**crystal_idx_lab_atm})

    lattice_labels = copy.deepcopy(crystal_box.lattice_labels)
    lattice_labels.update({**crystal_idx_lab_lat})

    int_sites, int_labels = None, None
    if crystal_box.interstice_sites is not None:

        crystal_idx_lab_int = {
            'crystal_idx': (
                np.array([0]),
                np.zeros(crystal_box.interstice_sites.shape[1], dtype=int)
            ),
        }

        int_labels = copy.deepcopy(crystal_box.interstice_labels)
        int_labels.update({**crystal_idx_lab_int})

    crystals = [{
        'crystal': np.copy(supercell),
        'origin': np.zeros((3, 1)),
        'cs_idx': 0,
        'cs_orientation': np.eye(3),
        'cs_origin': [0, 0, 0],
    }]

    as_params = {
        'supercell': supercell,
        'atom_sites': atom_sites,
        'atom_labels': atom_labels,
        'lattice_sites': lattice_sites,
        'lattice_labels': lattice_labels,
        'interstice_sites': int_sites,
        'interstice_labels': int_labels,
        'crystals': crystals,
        'crystal_structures': [crystal_structure],
        'overlap_tol': overlap_tol,
        'tile': tile,
    }
    bulk_crystal = BulkCrystal(as_params, repeats)

    return bulk_crystal
