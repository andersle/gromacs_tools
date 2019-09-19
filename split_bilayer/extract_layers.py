"""Extract layers from a bilayer."""
import pathlib
import sys
import numpy as np
from matplotlib import pyplot as plt
from gromacs import (
    read_gromacs_file,
    write_gromacs_gro_file,
)

plt.style.use('seaborn-talk')


def select_residues(frame, residue_name):
    """Get indices for residues based on a name."""
    index = []
    for i, residue in enumerate(frame['residuname']):
        if residue == residue_name:
            index.append(i)
    return index


def select_residues_number(frame, residunr):
    """Select residues based on a residue numbers."""
    index = []
    for i, residue in enumerate(frame['residunr']):
        if residue in residunr:
            index.append(i)
    return index


def select_atom_in_residues(frame, residue_index, atom_name):
    """Select atoms in residues."""
    atoms = []
    for i in residue_index:
        if frame['atomname'][i] == atom_name:
            atoms.append(i)
    return atoms


def get_positions(frame, index):
    """Get positions given indices."""
    xpos = np.array([frame['x'][i] for i in index])
    ypos = np.array([frame['y'][i] for i in index])
    zpos = np.array([frame['z'][i] for i in index])
    xyz = np.column_stack((xpos, ypos, zpos))
    return xyz


def plot_heads(xyz, upper, lower):
    """Plot positions of heads."""
    fig = plt.figure()
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)

    pairs = ((0, 1), (0, 2), (1, 2))
    labels = ('x', 'y', 'z')

    for (idx1, idx2), axi in zip(pairs, (ax1, ax2, ax3)):
        scat1 = axi.scatter(xyz[upper, idx1], xyz[upper, idx2])
        scat2 = axi.scatter(xyz[lower, idx1], xyz[lower, idx2])
        axi.set_xlabel(labels[idx1])
        axi.set_ylabel(labels[idx2])
        axi.axvline(x=np.average(xyz[:, idx1]), color='#262626', ls=':')
        axi.axhline(y=np.average(xyz[:, idx2]), color='#262626', ls=':')
        axi.axhline(y=np.average(xyz[upper, idx2]),
                    color=scat1.get_facecolor()[0], ls='--')
        axi.axhline(y=np.average(xyz[lower, idx2]),
                    color=scat2.get_facecolor()[0], ls='--')
    fig.tight_layout()


def split_bilayer(frame, heads, split_along):
    """Split a bilayer along a given axis."""
    xyz = get_positions(frame, heads)
    split = ('x', 'y', 'z').index(split_along)
    pos = xyz[:, split]
    avg_pos = np.average(pos)
    upper = np.where(pos >= avg_pos)[0]
    lower = np.where(pos < avg_pos)[0]
    plot_heads(xyz, upper, lower)
    upper_idx = [heads[i] for i in upper]
    lower_idx = [heads[i] for i in lower]
    upper_resnum = [frame['residunr'][i] for i in upper_idx]
    lower_resnum = [frame['residunr'][i] for i in lower_idx]
    print('Residues in upper: {}'.format(len(upper_resnum)))
    print('Residues in lower: {}'.format(len(lower_resnum)))
    idx1 = select_residues_number(frame, upper_resnum)
    idx2 = select_residues_number(frame, lower_resnum)
    return idx1, idx2


def extract_snapshot(frame, index):
    """Extract a subset of atoms from a given frame."""
    snapshot = {
        'header': 'Extracted from bilayer.',
        'box': frame['box'],
        'residunr': [],
        'residuname': [],
        'atomname': [],
        'atomnr': [],
        'x': [],
        'y': [],
        'z': [],
    }
    for key in ('residunr', 'residuname', 'atomname', 'atomnr', 'x', 'y', 'z'):
        for i in index:
            snapshot[key].append(frame[key][i])
    xpos = np.array([i for i in snapshot['x']])
    ypos = np.array([i for i in snapshot['y']])
    zpos = np.array([i for i in snapshot['z']])
    xyz = np.column_stack((xpos, ypos, zpos))
    return snapshot, xyz


def main(grofile, lipid_type, atom_head='P', split_along='z'):
    """Read frames and extract upper/lower bilayer."""
    frames = read_gromacs_file(grofile)
    for i, frame in enumerate(frames):
        index_lipid = select_residues(frame, lipid_type)
        heads = select_atom_in_residues(frame, index_lipid, atom_head)
        upper, lower = split_bilayer(frame, heads, split_along)
        upper_gro, upper_xyz = extract_snapshot(frame, upper)
        lower_gro, lower_xyz = extract_snapshot(frame, lower)
        write_gromacs_gro_file(
            'upper_{}_{}.gro'.format(pathlib.Path(grofile).stem, i),
            upper_gro,
            upper_xyz,
            np.zeros_like(upper_xyz)
        )
        write_gromacs_gro_file(
            'lower_{}_{}.gro'.format(pathlib.Path(grofile).stem, i),
            lower_gro,
            lower_xyz,
            np.zeros_like(lower_xyz)
        )
    plt.show()


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
