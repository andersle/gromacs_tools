"""Remove water from a bilayer."""
import sys
import numpy as np
from matplotlib import pyplot as plt
from pyretis.inout.formats.gromacs import (
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


def select_atoms(frame, atom_name):
    """Select atoms based on a name."""
    index = []
    for i, atom in enumerate(frame['atomname']):
        if atom == atom_name:
            index.append(i)
    return index


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


def make_snapshots(frame, ignore):
    """Extract a subset of atoms from a given frame."""
    snapshot = {
        'header': 'Removed.',
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
        for i, item in enumerate(frame[key]):
            if i not in ignore:
                snapshot[key].append(item)
    return snapshot


def get_bilayer_boundary(frame, axis, bilayer_atom):
    """Get extent of a bilayer along a axis."""
    layer = select_atoms(frame, bilayer_atom)
    xyz = get_positions(frame)
    layer_xyz = xyz[layer, :]
    avg_layer = np.average(layer_xyz[:, axis])
    upper = np.where(layer_xyz[:, axis] >= avg_layer)[0]
    lower = np.where(layer_xyz[:, axis] < avg_layer)[0]
    upper_pos = np.average(layer_xyz[upper, axis])
    lower_pos = np.average(layer_xyz[lower, axis])
    bounds = [min(upper_pos, lower_pos), max(upper_pos, lower_pos)]
    return bounds


def select_waters(frame, bounds, axis, water_oxygen):
    """Select water molecules within a boundary."""
    xyz = get_positions(frame)
    waters = select_atoms(frame, water_oxygen)
    remove = np.where(
        (
            (bounds[0] < xyz[waters, axis]) &
            (bounds[1] > xyz[waters, axis])
        )
    )[0]
    idx_water = [waters[i] for i in remove]
    print('Will remove {} water molecule(s)'.format(len(idx_water)))
    water_resnum = {frame['residunr'][i] for i in idx_water}
    water_index = set(select_residues_number(frame, water_resnum))
    print('Atoms to remove: {}'.format(len(water_index)))
    return water_index


def main(gro_file, split_along='z'):
    """Read frames and extract upper/lower bilayer."""
    gro = [i for i in read_gromacs_file(gro_file)][0]

    split = ('x', 'y', 'z').index(split_along)
    bounds = get_bilayer_boundary(gro, split, 'P')
    water_index = select_waters(gro, bounds, split, 'OW')
    snapshot = make_snapshots(gro, water_index)
    xyz = get_positions(snapshot)
    write_gromacs_gro_file(
        'removed.gro',
        snapshot,
        xyz,
        np.zeros_like(xyz),
    )


if __name__ == '__main__':
    main(sys.argv[1])
