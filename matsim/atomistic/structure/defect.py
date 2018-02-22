"""`matsim.atomistic.structure.defect.py`"""


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

        out = '{}_{}^{}'.format(self.defect_species,
                                self.host_species, charge_str)

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
                    self.defect_species, self.host_species, self.index, self.charge,
                    self.interstice_type)
