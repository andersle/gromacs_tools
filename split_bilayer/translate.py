"""Translate bilayers."""
import pathlib
import sys
import numpy as np
from matplotlib import pyplot as plt
from gromacs import (
    read_gromacs_file,
    write_gromacs_gro_file,
)

plt.style.use('seaborn-talk')


def get_positions(frame):
    """Get positions given indices."""
    xpos = np.array([i for i in frame['x']])
    ypos = np.array([i for i in frame['y']])
    zpos = np.array([i for i in frame['z']])
    xyz = np.column_stack((xpos, ypos, zpos))
    return xyz


def merge_snapshots(frames):
    """Extract a subset of atoms from a given frame."""
    snapshot = {
        'header': 'Merged.',
        'box': frames[0]['box'],
        'residunr': [],
        'residuname': [],
        'atomname': [],
        'atomnr': [],
    }
    for key in ('residunr', 'residuname', 'atomname', 'atomnr'):
        for frame in frames:
            for item in frame[key]:
                snapshot[key].append(item)
    return snapshot


def main(upper_file, lower_file, delta_z):
    """Read frames and extract upper/lower bilayer."""
    upper = [i for i in read_gromacs_file(upper_file)][0]
    upper_xyz = get_positions(upper)
    lower = [i for i in read_gromacs_file(lower_file)][0]
    lower_xyz = get_positions(lower)
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    ax1.scatter(upper_xyz[:, 0], upper_xyz[:, 2], marker='o')
    ax1.scatter(lower_xyz[:, 0], lower_xyz[:, 2], marker='s')

    upper_xyz[:, 2] += delta_z
    lower_xyz[:, 2] -= delta_z

    ax2.scatter(upper_xyz[:, 0], upper_xyz[:, 2], marker='o')
    ax2.scatter(lower_xyz[:, 0], lower_xyz[:, 2], marker='s')
    write_gromacs_gro_file(
        'translated_{}.gro'.format(pathlib.Path(upper_file).stem),
        upper,
        upper_xyz,
        np.zeros_like(upper_xyz)
    )
    write_gromacs_gro_file(
        'translated_{}.gro'.format(pathlib.Path(lower_file).stem),
        lower,
        lower_xyz,
        np.zeros_like(lower_xyz)
    )
    merged = merge_snapshots((upper, lower))
    merged_xyz = np.vstack((upper_xyz, lower_xyz))
    write_gromacs_gro_file(
        'merged.gro',
        merged,
        merged_xyz,
        np.zeros_like(merged_xyz)
    )
    fig.tight_layout()
    plt.show()


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], float(sys.argv[3]))
