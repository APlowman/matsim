"""`matsim.simulation.sequence_funcs.py"""


def process_kpoint_sequence(sequence, base_sim):
    """Process values in a kpoint spacing sequence to avoid """

    # Get the MP grid from the supercell and modify vals to remove those
    # which generate duplicate MP grids.
    unique_grids = {}
    for v in sequence.vals:

        kpt_gd = tuple(base_sim.structure.get_kpoint_grid(v))

        if unique_grids.get(kpt_gd) is None:
            unique_grids.update({kpt_gd: [v]})
        else:
            unique_grids[kpt_gd].append(v)

    unique_vals = []
    for k, v in unique_grids.items():
        unique_vals.append(sorted(v)[0])

    print('Found unique kpoint spacings: {}'.format(unique_vals))
    sequence.vals = sorted(unique_vals)
