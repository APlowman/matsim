"""matsim.atomistic.structure.atomistic.py"""

import copy
import warnings

import numpy as np
import spglib
from mendeleev import element
from vecmaths import rotation

from matsim import geometry, vectors, utils
from matsim.utils import RestrictedDict, mut_exc_args, prt
from matsim.atomistic.structure import site_labs_from_jsonable, site_labs_to_jsonable
from matsim.atomistic.structure.crystal import CrystalBox, CrystalStructure
from matsim.atomistic.structure.visualise import visualise as struct_visualise


def modify_crystal_structure(cs, vol_change, ca_change):
    """
    Regenerate a CrystalStructure with a modified bravais lattice.

    Parameters
    ----------
    cs : CrystalStructure object
    vol_change : float
        Percentage change in volume
    ca_change : float
        Percentage change in c/a ratio

    Returns
    -------
    CrystalStructure

    """
    bl = cs.bravais_lattice

    # Modify hexagonal CrystalStructure
    if bl.lattice_system != 'hexagonal':
        raise NotImplementedError('Cannot modify non-hexagonal crystal '
                                  'structure.')

    # Generate new a and c lattice parameters based on originals and volume and
    # c/a ratio changes:
    v = get_hex_vol(bl.a, bl.c)
    v_new = v * (1 + vol_change / 100)

    ca_new = (bl.c / bl.a) * (1 + ca_change / 100)
    a_new = get_hex_a(ca_new, v_new)
    c_new = ca_new * a_new

    bl_new = BravaisLattice('hexagonal', a=a_new, c=c_new)
    cs_new = CrystalStructure(bl_new, copy.deepcopy(cs.motif))
    return cs_new


class AtomisticStructureException(Exception):
    pass


class AtomisticStructure(object):
    """
    Class to represent crystals of atoms

    Attributes
    ----------
    atom_sites : ndarray of shape (3, N)
        Array of column vectors representing the atom positions.
    supercell : ndarray of shape (3, 3)
        Array of column vectors representing supercell edge vectors.
    lattice_sites : ndarray of shape (3, M), optional
        Array of column vectors representing lattice site positions.
    crystals : list of dict of (str : ndarray or int), optional
        Each dict contains at least these keys:
            `crystal` : ndarray of shape (3, 3)
                Array of column vectors representing the crystal edge vectors.
            `origin` : ndarray of shape (3, 1)
                Column vector specifying the origin of this crystal.
        Additional keys are:
            'cs_idx': int
                Index of `crystal_structures`, defining to which
                CrystalStructure this crystal belongs.
            `cs_orientation`: ndarray of shape (3, 3)
                Rotation matrix which rotates the CrystalStructure lattice
                unit cell from the initialised BravaisLattice object to some
                other desired orientation.
            'cs_origin': list of float or int
                Origin of the CrystalStructure unit cell in multiples of the
                CrystalStructure unit cell vectors. For integer values, this
                will not affect the atomic structure of the crystal. To
                produce a rigid translation of the atoms within the crystal,
                non-integer values can be used.

    crystal_structures : list of CrystalStructure, optional
    crystal_idx : ndarray of shape (N,), optional
        Defines to which crystal each atom belongs.
    lat_crystal_idx : ndarray of shape (M,), optional
        Defines to which crystal each lattice site belongs
    species_idx : ndarray of shape (N,), optional
        Defines to which species each atom belongs, indexed within the atom's
        crystal_structure. For atom index `i`, this indexes
        `crystal_structures[
            crystals[crystal_idx[i]]['cs_idx']]['species_set']`
        Either specify (`all_species` and `all_species_idx`) or (`species_idx`
        and `motif_idx`), but not both.
    motif_idx : ndarray of shape (N,), optional
        Defines to which motif atom each atom belongs, indexed within the
        atom's crystal_structure. For atom index `i`, this indexes
        `crystal_structures[
            crystals[crystal_idx[i]]['cs_idx']]['species_motif']`
        Either specify (`all_species` and `all_species_idx`) or (`species_idx`
        and `motif_idx`), but not both.
    all_species : ndarray of str, optional
        1D array of strings representing the distinct species. Either specify
        (`all_species` and `all_species_idx`) or (`species_idx` and
        `motif_idx`), but not both.
    all_species_idx : ndarray of shape (N, ), optional
        Defines to which species each atom belongs, indexed over the whole
        AtomisticStructure. This indexes `all_species`. Either specify
        (`all_species` and `all_species_idx`) or (`species_idx` and
        `motif_idx`), but not both.

    atom_sites_frac
    num_atoms_per_crystal
    num_atoms
    num_crystals
    reciprocal_supercell

    Methods
    -------
    todo

    TODO:
    -   Re-write docstrings.
    -   Consolidate atom/lattice/interstice into a list of Sites objects.

    """

    # Keys allowed in the meta attribute:
    ok_meta = ['supercell_type']

    def __init__(self, supercell=None, atom_sites=None, atom_labels=None, origin=None,
                 lattice_sites=None, lattice_labels=None, interstice_sites=None,
                 interstice_labels=None, crystals=None, crystal_structures=None,
                 overlap_tol=1, state=None):
        """Constructor method for AtomisticStructure object."""

        mut_exc_args(
            {
                'supercell': supercell,
                'atom_sites': atom_sites,
                'atom_labels': atom_labels,
            },
            {'state': state}
        )

        if state:

            self.origin = state['origin']

            self.atom_sites = state['atom_sites']
            self.atom_labels = state['atom_labels']
            self.supercell = state['supercell']
            self.meta = state['meta']

            self.lattice_sites = state.get('lattice_sites')
            self.lattice_labels = state.get('lattice_labels')
            self.interstice_sites = state.get('interstice_sites')
            self.interstice_labels = state.get('interstice_labels')

            self.crystals = state['crystals']

            # print('loaded AtomisticStructure from jsonable.')
            # prt(self.crystals, 'self.crystals')

            self.crystal_structures = state['crystal_structures']
            self._overlap_tol = state['_overlap_tol']

        else:

            if origin is None:
                origin = np.zeros((3, 1))

            self.origin = utils.to_col_vec(origin)

            self.atom_sites = atom_sites
            self.atom_labels = atom_labels
            self.supercell = supercell
            self.meta = RestrictedDict(AtomisticStructure.ok_meta)

            self.lattice_sites = lattice_sites
            self.lattice_labels = lattice_labels
            self.interstice_sites = interstice_sites
            self.interstice_labels = interstice_labels
            self.crystals = crystals
            self.crystal_structures = crystal_structures
            self._overlap_tol = overlap_tol

            self.check_overlapping_atoms(overlap_tol)

            # Check handedness:
            if self.volume < 0:
                raise ValueError('Supercell does not form a right - handed '
                                 'coordinate system.')

    def to_jsonable(self):
        """Generate a dict representation that can be JSON serialised."""

        # Only add meta keys allowed in this class, since a subclass may need
        # to jsonify meta keys allowed in it's own class.
        meta_js = {key: val for key, val in self.meta.items()
                   if key in AtomisticStructure.ok_meta}

        # prt(self.crystals, 'self.crystals')

        crystals_js = []
        for crys in self.crystals:
            crystals_js.append({
                'crystal': crys['crystal'].tolist(),
                'origin': crys['origin'].tolist(),
                'cs_idx': crys['cs_idx'],
                'cs_orientation': crys['cs_orientation'].tolist(),
                'cs_origin': crys['cs_origin'],
            })

        crys_struct_js = [i.to_jsonable() for i in self.crystal_structures]

        state = {
            'origin': self.origin.tolist(),
            'supercell': self.supercell.tolist(),
            'atom_sites': self.atom_sites.tolist(),
            'atom_labels': site_labs_to_jsonable(self.atom_labels),
            'lattice_sites': None,
            'lattice_labels': None,
            'interstice_sites': None,
            'interstice_labels': None,
            'meta': meta_js,
            'crystals': crystals_js,
            'crystal_structures': crys_struct_js,
            '_overlap_tol': self._overlap_tol,
        }

        if self.lattice_sites is not None:
            state.update({
                'lattice_sites': self.lattice_sites.tolist(),
                'lattice_labels': site_labs_to_jsonable(self.lattice_labels),
            })

        if self.interstice_sites is not None:
            state.update({
                'interstice_sites': self.interstice_sites.tolist(),
                'interstice_labels': site_labs_to_jsonable(self.interstice_labels),
            })

        return state

    @classmethod
    def from_jsonable(cls, state):
        """Instantiate from a JSONable dict."""

        crystals_native = []
        for crys in state['crystals']:
            crystals_native.append({
                'crystal': np.array(crys['crystal']),
                'origin': np.array(crys['origin']),
                'cs_idx': crys['cs_idx'],
                'cs_orientation': np.array(crys['cs_orientation']),
                'cs_origin': crys['cs_origin'],
            })

        crys_struct_native = [CrystalStructure.from_jsonable(i)
                              for i in state['crystal_structures']]

        meta_native = RestrictedDict(
            AtomisticStructure.ok_meta,
            {key: val for key, val in state['meta'].items()
             if key in AtomisticStructure.ok_meta},
        )

        state.update({
            'origin': np.array(state['origin']),
            'supercell': np.array(state['supercell']),
            'atom_sites': np.array(state['atom_sites']),
            'atom_labels': site_labs_from_jsonable(state['atom_labels']),
            'meta': meta_native,
            'crystals': crystals_native,
            'crystal_structures': crys_struct_native,
            '_overlap_tol': state['_overlap_tol'],
        })

        if state.get('lattice_sites') is not None:
            state.update({
                'lattice_sites': np.array(state['lattice_sites']),
                'lattice_labels': site_labs_from_jsonable(state['lattice_labels']),
            })

        if state.get('interstice_sites') is not None:
            state.update({
                'interstice_sites': np.array(state['interstice_sites']),
                'interstice_labels': site_labs_from_jsonable(state['interstice_labels']),
            })

        return state

    def translate(self, shift):
        """
        Translate the AtomisticStructure.

        Parameters
        ----------
        shift : list or ndarray of size 3

        """

        shift = utils.to_col_vec(shift)
        self.origin += shift
        self.atom_sites += shift

        if self.lattice_sites is not None:
            self.lattice_sites += shift

        if self.interstice_sites is not None:
            self.interstice_sites += shift

        if self.crystals is not None:
            for c_idx in range(len(self.crystals)):
                self.crystals[c_idx]['origin'] += shift

    def rotate(self, rot_mat):
        """
        Rotate the AtomisticStructure about its origin according to a rotation
        matrix.

        Parameters
        ----------
        rot_mat : ndarray of shape (3, 3)
            Rotation matrix that pre-multiplies column vectors in order to 
            rotate them about a particular axis and angle.

        """

        origin = np.copy(self.origin)
        self.translate(-origin)

        self.supercell = np.dot(rot_mat, self.supercell)
        self.atom_sites = np.dot(rot_mat, self.atom_sites)

        if self.lattice_sites is not None:
            self.lattice_sites = np.dot(rot_mat, self.lattice_sites)

        if self.interstice_sites is not None:
            self.interstice_sites = np.dot(rot_mat, self.interstice_sites)

        if self.crystals is not None:

            for c_idx in range(len(self.crystals)):

                c = self.crystals[c_idx]

                c['crystal'] = np.dot(rot_mat, c['crystal'])
                c['origin'] = np.dot(rot_mat, c['origin'])

                if 'cs_orientation' in c.keys():
                    c['cs_orientation'] = np.dot(rot_mat, c['cs_orientation'])

        self.translate(origin)

    def visualise(self, **kwargs):
        struct_visualise(self, **kwargs)

    def reorient_to_lammps(self):
        """
        Reorient the supercell and its contents to a LAMMPS-compatible
        orientation. Also translate the origin to (0,0,0).

        Returns
        -------
        ndarray of shape (3, 3)
            Rotation matrix used to reorient the supercell and its contents

        """

        # Find rotation matrix which rotates to a LAMMPS compatible orientation
        sup_lmps = rotation.align_xy(self.supercell)
        R = np.dot(sup_lmps, self.supercell_inv)

        # Move the origin to (0,0,0): (I presume this is necessary for LAMMPS?)
        self.translate(-self.origin)

        # Rotate the supercell and its contents by R
        self.rotate(R)

        return R

    def wrap_sites_to_supercell(self, sites='all', dirs=None):
        """
        Wrap sites to within the supercell.

        Parameters
        ----------
        sites : str
            One of "atom", "lattice", "interstice" or "all".
        dirs : list of int, optional
            Supercell direction indices to apply wrapping. Default is None, in
            which case atoms are wrapped in all directions.            

        """

        # Validation
        if dirs is not None:
            if len(set(dirs)) != len(dirs):
                raise ValueError('Indices in `dirs` must not be repeated.')

            if len(dirs) not in [1, 2, 3]:
                raise ValueError('`dirs` must be a list of length 1, 2 or 3.')

            for d in dirs:
                if d not in [0, 1, 2]:
                    raise ValueError('`dirs` must be a list whose elements are'
                                     '0, 1 or 2.')

        allowed_sites_str = [
            'atom',
            'lattice',
            'interstice',
            'all',
        ]

        if not isinstance(sites, str) or sites not in allowed_sites_str:
            raise ValueError('`sites` must be a string and one of: "atom", '
                             '"lattice", "interstice" or "all".')

        if sites == 'all':
            sites_arr = [
                self.atom_sites,
                self.lattice_sites,
                self.interstice_sites
            ]
        elif sites == 'atom':
            sites_arr = [self.atom_sites]
        elif sites == 'lattice':
            sites_arr = [self.lattice_sites]
        elif sites == 'interstice':
            sites_arr = [self.interstice_sites]

        for s_idx in range(len(sites_arr)):

            s = sites_arr[s_idx]
            if s is None:
                continue

            # Get sites in supercell basis:
            s_sup = np.dot(self.supercell_inv, s)

            # Wrap atoms:
            s_sup_wrp = np.copy(s_sup)
            s_sup_wrp[dirs] -= np.floor(s_sup_wrp[dirs])

            # Snap to 0:
            s_sup_wrp = vectors.snap_arr_to_val(s_sup_wrp, 0, 1e-12)

            # Convert back to Cartesian basis
            s_std_wrp = np.dot(self.supercell, s_sup_wrp)

            # Update attributes:
            sites_arr[s_idx][:] = s_std_wrp

    def add_point_defects(self, point_defects):
        """
        Add point defects to the structure.

        Parameters
        ----------
        point_defects : list of PointDefect objects

        """
        pass

    def add_atom(self, coords, species, crystal_idx, is_frac_coords=False):
        """Add an atom to the structure."""
        pass

    def remove_atom(self, atom_idx):
        """Remove an atom from the structure."""
        pass

    @property
    def supercell_inv(self):
        return np.linalg.inv(self.supercell)

    @property
    def atom_sites_frac(self):
        return np.dot(self.supercell_inv, self.atom_sites)

    @property
    def lattice_sites_frac(self):
        if self.lattice_sites is not None:
            return np.dot(self.supercell_inv, self.lattice_sites)
        else:
            return None

    @property
    def interstice_sites_frac(self):
        if self.interstice_sites is not None:
            return np.dot(self.supercell_inv, self.interstice_sites)
        else:
            return None

    @property
    def species(self):
        return self.atom_labels['species'][0]

    @property
    def species_idx(self):
        return self.atom_labels['species'][1]

    @property
    def all_species(self):
        """Get the species of each atom as a string array."""
        return self.species[self.species_idx]

    @property
    def spglib_cell(self):
        """Returns a tuple representing valid input for the spglib library."""

        cell = (self.supercell.T,
                self.atom_sites_frac.T,
                [element(i).atomic_number for i in self.all_species])
        return cell

    @property
    def num_atoms_per_crystal(self):
        """Computes number of atoms in each crystal, returns a list."""

        if self.crystals is None:
            return None

        na = []
        for c_idx in range(len(self.crystals)):
            crystal_idx_tup = self.atom_labels['crystal_idx']
            crystal_idx = crystal_idx_tup[0][crystal_idx_tup[1]]
            na.append(np.where(crystal_idx == c_idx)[0].shape[0])

        return na

    @property
    def num_atoms(self):
        """Computes total number of atoms."""
        return self.atom_sites.shape[1]

    @property
    def num_crystals(self):
        """Returns number of crystals."""
        return len(self.crystals)

    @property
    def reciprocal_supercell(self):
        """Returns the reciprocal supercell as array of column vectors."""

        v = self.supercell
        cross_1 = np.cross(v[:, 1], v[:, 2])
        cross_2 = np.cross(v[:, 0], v[:, 2])
        cross_3 = np.cross(v[:, 0], v[:, 1])

        B = np.zeros((3, 3))
        B[:, 0] = 2 * np.pi * cross_1 / (np.dot(v[:, 0], cross_1))
        B[:, 1] = 2 * np.pi * cross_2 / (np.dot(v[:, 1], cross_2))
        B[:, 2] = 2 * np.pi * cross_3 / (np.dot(v[:, 2], cross_3))

        return B

    def get_kpoint_grid(self, separation):
        """
        Get the MP kpoint grid size for a given kpoint separation.

        Parameters
        ----------
        separation : float or int or ndarray of shape (3, )
            Maximum separation between kpoints, in units of inverse Angstroms.
            If an array, this is the separations in each reciprocal supercell
            direction.

        Returns
        -------
        ndarray of int of shape (3, )
            MP kpoint grid dimensions along each reciprocal supercell
            direction.

        """

        recip = self.reciprocal_supercell
        grid = np.ceil(np.round(
            np.linalg.norm(recip, axis=0) / (separation * 2 * np.pi),
            decimals=8)
        ).astype(int)

        return grid

    def get_kpoint_spacing(self, grid):
        """
        Get the kpoint spacing given an MP kpoint grid size.

        Parameters
        ----------
        grid : list of length 3
            Grid size in each of the reciprocal supercell directions.

        Returns
        -------
        ndarray of shape (3, )
            Separation between kpoints in each of the reciprocal supercell
            directions.

        """

        recip = self.reciprocal_supercell
        seps = np.linalg.norm(recip, axis=0) / (np.array(grid) * 2 * np.pi)

        return seps

    @property
    def crystal_centres(self):
        """Get the midpoints of each crystal in the structure."""

        return [geometry.get_box_centre(c['crystal'], origin=c['origin'])
                for c in self.crystals]

    def tile_supercell(self, tiles):
        """
        Tile supercell and its sites by some integer factors in each supercell 
        direction.

        Parameters
        ----------
        tiles : tuple or list of length 3 or ndarray of size 3
            Number of repeats in each supercell direction.

        """
        invalid_msg = ('`tiles` must be a tuple or list of three integers '
                       'greater than 0.')

        if isinstance(tiles, np.ndarray):
            tiles = np.squeeze(tiles).tolist()

        if len(tiles) != 3:
            raise ValueError(invalid_msg)

        for t in tiles:
            if not isinstance(t, int) or t < 1:
                raise ValueError(invalid_msg)

        tl_atm, tl_atm_lb = self.get_tiled_sites(
            self.atom_sites, self.atom_labels, tiles)

        self.atom_sites = tl_atm
        self.atom_labels = tl_atm_lb

        if self.lattice_sites is not None:
            tl_lat, tl_lat_lb = self.get_tiled_sites(
                self.lattice_sites, self.lattice_labels, tiles)

            self.lattice_sites = tl_lat
            self.lattice_labels = tl_lat_lb

        if self.interstice_sites is not None:
            tl_int, tl_int_lb = self.get_tiled_sites(
                self.interstice_sites, self.interstice_labels, tiles)

            self.interstice_sites = tl_int
            self.interstice_labels = tl_int_lb

        self.supercell *= tiles

    def get_tiled_sites(self, sites, site_labels, tiles):
        """
        Get sites (atoms, lattice, interstice) and site labels tiled by some
        integer factors in each supercell direction.

        Sites are tiled in the positive supercell directions.

        Parameters
        ----------
        tiles : tuple or list of length 3 or ndarray of size 3
            Number of repeats in each supercell direction.

        Returns
        -------
        sites_tiled : ndarray
        labels_tiled : dict

        """

        invalid_msg = ('`tiles` must be a tuple or list of three integers '
                       'greater than 0.')

        if isinstance(tiles, np.ndarray):
            tiles = np.squeeze(tiles).tolist()

        if len(tiles) != 3:
            raise ValueError(invalid_msg)

        sites_tiled = np.copy(sites)
        labels_tiled = {k: tuple(np.copy(i) for i in v)
                        for k, v in site_labels.items()}

        for t_idx, t in enumerate(tiles):

            if t == 1:
                continue

            if not isinstance(t, int) or t < 1:
                raise ValueError(invalid_msg)

            v = self.supercell[:, t_idx:t_idx + 1]
            v_range = v * np.arange(1, t)

            all_t = v_range.T[:, :, np.newaxis]

            sites_stack = all_t + sites_tiled
            add_sites = np.hstack(sites_stack)
            sites_tiled = np.hstack([sites_tiled, add_sites])

            labels_tiled_new = {}
            for k, v in labels_tiled.items():

                add_label_idx = np.tile(v[1], t - 1)
                new_label_idx = np.concatenate((v[1], add_label_idx))
                labels_tiled_new.update({
                    k: (v[0], new_label_idx)
                })

            labels_tiled = labels_tiled_new

        return sites_tiled, labels_tiled

    def get_interatomic_dist(self, periodic=True):
        """
        Find the distances between unique atom pairs across the whole
        structure.

        Parameters
        ----------
        periodic : bool
            If True, the atom sites are first tiled in each supercell direction
            to ensure that distances between periodic cells are considered.
            Currently, this is crude, and so produces interatomic distances
            between like atoms (i.e. of one supercell vector length).

        Returns
        ------
        ndarray of shape (N,)

        TODO:
        -   Improve consideration of periodicity. Maybe instead have a function
            `get_min_interatomic_dist` which gets the minimum distances of each
            atom and every other atom, given periodicity.

        """
        if periodic:
            atms = self.get_tiled_sites(
                self.atom_sites, self.atom_labels, [2, 2, 2])[0]
        else:
            atms = self.atom_sites

        return vectors.get_vec_distances(atms)

    def check_overlapping_atoms(self, tol=None):
        """
        Checks if any atoms are overlapping within a tolerance.

        Parameters
        ----------
        tol : float, optional
            Distance below which atoms are considered to be overlapping. By,
            default uses the value assigned on object initialisation as
            `_overlap_tol`.

        Raises
        ------
        AtomisticStructureException
            If any atoms are found to overlap within `tol`.

        """
        if tol is None:
            tol = self._overlap_tol

        dist = self.get_interatomic_dist()
        if np.any(dist < tol):
            raise AtomisticStructureException('Found overlapping atoms. '
                                              'Minimum separation: '
                                              '{:.3f}'.format(np.min(dist)))

    def get_sym_ops(self):
        return spglib.get_symmetry(self.spglib_cell)

    def shift_atoms(self, shift, wrap=False):
        """
        Perform a rigid shift on all atoms, in fractional supercell coordinates.

        Parameters
        ----------
        shift : list or tuple of length three or ndarry of shape (3,) of float
            Fractional supercell coordinates to translate all atoms by.
        wrap : bool
            If True, wrap atoms to within the supercell edges after shift.
        """

        shift = np.array(shift)[:, np.newaxis]
        shift_std = np.dot(self.supercell, shift)
        self.atom_sites += shift_std

        if wrap:
            self.wrap_sites_to_supercell(sites='atom')

    def add_vac(self, thickness, dir_idx, position=1):
        """
        Extend the supercell in a given direction.

        Supercell vector given by direction index `dir_idx` is extended such
        that it's component in the direction normal to the other two supercell
        vectors is a particular `thickness`.

        Parameters
        ----------
        thickness : float
            Thickness of vacuum to add
        dir_idx : int 0, 1 or 2
            Supercell direction in which to add vacuum
        position : float
            Fractional coordinate along supercell vector given by `dir_idx` at
            which to add the vacuum. By default, adds vacuum to the far face of
            the supercell, such that atom Cartesian coordinates are not
            affected. Must be between 0 (inclusive) and 1 (inclusive).
        """

        # TODO: validate it does what we want. Maybe revert back to calling it
        # `add_surface_vac`.

        warnings.warn('!! Untested function... !!')

        if dir_idx not in [0, 1, 2]:
            raise ValueError('`dir_idx` must be 0, 1 or 2.')

        if position < 0 or position > 1:
            raise ValueError('`position` must be between 0 (inclusive) and 1 '
                             '(inclusive).')

        non_dir_idx = [i for i in [0, 1, 2] if i != dir_idx]
        v1v2 = self.supercell[:, non_dir_idx]
        v3 = self.supercell[:, dir_idx]

        n = np.cross(v1v2[:, 0], v1v2[:, 1])
        n_unit = n / np.linalg.norm(n)
        v3_mag = np.linalg.norm(v3)
        v3_unit = v3 / v3_mag
        d = thickness / np.dot(n_unit, v3_unit)

        v3_mag_new = v3_mag + d
        v3_new = v3_unit * v3_mag_new

        self.supercell[:, dir_idx] = v3_new

        asf = self.atom_sites_frac
        shift_idx = np.where(asf[dir_idx] > position)[0]

        self.atom_sites[:, shift_idx] += (n_unit * thickness)

    def check_atomic_environment(self, checks_list):
        """Invoke checks of the atomic environment."""

        allowed_checks = {
            'atoms_overlap': self.check_overlapping_atoms,
        }

        for chk, func in allowed_checks.items():
            if chk in checks_list:
                func()

    @property
    def volume(self):
        """Get the volume of the supercell."""
        sup = self.supercell
        return np.dot(np.cross(sup[:, 0], sup[:, 1]), sup[:, 2])


class BulkCrystal(AtomisticStructure):
    """

    Attributes
    ----------
    crystal_structure : CrystalStructure

    TODO:
    -   Add proper support for cs_orientation and cs_origin. Maybe allow one of
        `box_lat` or `box_std` for the more general case.

    """

    def __init__(self, crystal_structure, box_lat, overlap_tol=1):
        """Constructor method for BulkCrystal object."""

        # Validation
        if any([i in vectors.num_equal_cols(box_lat) for i in [2, 3]]):
            raise ValueError(
                'Identical columns found in box_lat: \n{}\n'.format(box_lat))

        supercell = np.dot(crystal_structure.bravais_lattice.vecs, box_lat)

        cb = CrystalBox(crystal_structure, supercell)
        atom_sites = cb.atom_sites
        lattice_sites = cb.lattice_sites

        crystal_idx_lab = {
            'crystal_idx': (
                np.array([0]),
                np.zeros(atom_sites.shape[1], dtype=int)
            ),
        }

        atom_labels = copy.deepcopy(cb.atom_labels)
        atom_labels.update({**crystal_idx_lab})

        lattice_labels = copy.deepcopy(cb.lattice_labels)
        lattice_labels.update({**crystal_idx_lab})

        int_sites, int_labels = None, None
        if cb.interstice_sites is not None:
            int_labels = copy.deepcopy(cb.interstice_labels)
            int_labels.update({**crystal_idx_lab})

        crystals = [{
            'crystal': supercell,
            'origin': np.zeros((3, 1)),
            'cs_idx': 0,
            'cs_orientation': np.eye(3),
            'cs_origin': [0, 0, 0]
        }]

        super().__init__(supercell,
                         atom_sites,
                         atom_labels,
                         lattice_sites=lattice_sites,
                         lattice_labels=lattice_labels,
                         interstice_sites=int_sites,
                         interstice_labels=int_labels,
                         crystals=crystals,
                         crystal_structures=[crystal_structure],
                         overlap_tol=overlap_tol)

        self.meta.update({'supercell_type': ['bulk']})


class PointDefect(object):
    """
    Class to represent a point defect embedded within an AtomisticStructure

    Attributes
    ----------
    defect_species : str
        Chemical symbol of the defect species or "v" for vacancy
    host_species : str
        Chemical symbol of the species which this defect replaces or "i" for
        interstitial.
    index : int
        The atom or interstitial site index within the AtomisticStructure.
    charge : float
        The defect's electronic charge relative to that of the site it occupies.
    interstice_type : str
        Set to "tetrahedral" or "octahedral" if `host_species` is "i".

    """

    def __init__(self, defect_species, host_species, index=None, charge=0,
                 interstice_type=None):

        # Validation
        if interstice_type not in [None, 'tetrahedral', 'octahedral']:
            raise ValueError('Interstice type "{}" not understood.'.format(
                interstice_type))

        if host_species != 'i' and interstice_type is not None:
            raise ValueError('Non-interstitial defect specified but '
                             '`interstice_type` also specified.')

        if defect_species == 'v' and host_species == 'i':
            raise ValueError('Cannot add a vacancy defect to an '
                             'interstitial site!')

        if host_species == 'i' and interstice_type is None:
            raise ValueError('`interstice_type` must be specified for '
                             'interstitial point defect.')

        self.defect_species = defect_species
        self.host_species = host_species
        self.index = index
        self.charge = charge
        self.interstice_type = interstice_type

    def __str__(self):
        """
        References
        ----------
        https://en.wikipedia.org/wiki/Kr%C3%B6ger%E2%80%93Vink_notation

        """
        # String representation of the charge in Kroger-Vink notation
        if self.charge == 0:
            charge_str = 'x'
        elif self.charge > 0:
            charge_str = '•' * abs(self.charge)
        elif self.charge < 0:
            charge_str = '′' * abs(self.charge)

        out = '{}_{}^{}'.format(self.defect_species, self.host_species, charge_str,
                                self.index)

        if self.index is not None:
            idx_str_int = 'interstitial' if self.host_species == 'i' else 'atom'
            idx_str = 'at {} index {}'.format(idx_str_int, self.index)
        else:
            idx_str = ''

        if self.interstice_type is not None:
            out += ' ({}'.format(self.interstice_type)
            out += ' ' + idx_str + ')'
        else:
            out += ' (' + idx_str + ')'

        return out

    def __repr__(self):
        return ('PointDefect({!r}, {!r}, index={!r}, charge={!r}, '
                'interstice_type={!r})').format(
                    self.defect_species, self.atom_site, self.index, self.charge,
                    self.interstice_type)
