"""matsim.simulation.sim.py"""

import copy

from matsim import utils


class Simulation(object):
    """Class to represent a simulation."""

    def __init__(self):

        self.dbid = None
        self.runs = []
        self.options = {}

    def get_run_parameters(self, run_idx):
        """Merge the simulation parameters with those of a particular run to
        get parameters for writing the input files for that run."""

        sim_params = copy.deepcopy(self.options['params'])
        run_params = self.runs[run_idx]['run_params']
        utils.update_dict(sim_params, run_params)

        return sim_params
