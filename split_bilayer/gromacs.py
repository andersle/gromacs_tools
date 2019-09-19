"""Methods for reading GROMACS .gro files."""
import numpy as np


_GRO_FMT = '{0:5d}{1:5s}{2:5s}{3:5d}{4:8.3f}{5:8.3f}{6:8.3f}'
_GRO_VEL_FMT = _GRO_FMT + '{7:8.4f}{8:8.4f}{9:8.4f}'
_GRO_BOX_FMT = '{:15.9f}'


def read_gromacs_lines(lines):
    """Read and parse GROMACS GRO data.

    This method will read a GROMACS file and yield the different
    snapshots found in the file.

    Parameters
    ----------
    lines : iterable
        Some lines of text data representing a GROMACS GRO file.

    Yields
    ------
    out : dict
        This dict contains the snapshot.

    """
    lines_to_read = 0
    snapshot = {}
    read_natoms = False
    gro = (5, 5, 5, 5, 8, 8, 8, 8, 8, 8)
    gro_keys = ('residunr', 'residuname', 'atomname', 'atomnr',
                'x', 'y', 'z', 'vx', 'vy', 'vz')
    gro_type = (int, str, str, int, float, float, float, float, float, float)
    for line in lines:
        if read_natoms:
            read_natoms = False
            lines_to_read = int(line.strip()) + 1
            continue  # just skip to next line
        if lines_to_read == 0:  # new snapshot
            if snapshot:
                yield snapshot
            snapshot = {'header': line.strip()}
            read_natoms = True
        elif lines_to_read == 1:  # read box
            snapshot['box'] = np.array(
                [float(i) for i in line.strip().split()]
            )
            lines_to_read -= 1
        else:  # read atoms
            lines_to_read -= 1
            current = 0
            for i, key, gtype in zip(gro, gro_keys, gro_type):
                val = line[current:current+i].strip()
                if not val:
                    # This typically happens if we try to read velocities
                    # and they are not present in the file.
                    break
                value = gtype(val)
                current += i
                try:
                    snapshot[key].append(value)
                except KeyError:
                    snapshot[key] = [value]
    if snapshot:
        yield snapshot


def read_gromacs_file(filename):
    """Read a GROMACS GRO file."""
    with open(filename, 'r') as fileh:
        for snapshot in read_gromacs_lines(fileh):
            yield snapshot


def write_gromacs_gro_file(outfile, txt, xyz, vel):
    """Write configuration in GROMACS GRO format.

    Parameters
    ----------
    outfile : string
        The name of the file to create.
    txt : dict of lists of strings
        This dict contains the information on residue-numbers, names,
        etc. required to write the GRO file.
    xyz : numpy.array
        The positions to write.
    vel : numpy.array
        The velocities to write.

    """
    print('Writing GROMACS gro file "{}"'.format(outfile))
    resnum = txt['residunr']
    resname = txt['residuname']
    atomname = txt['atomname']
    atomnr = txt['atomnr']
    npart = len(xyz)
    with open(outfile, 'w') as output:
        output.write('{}\n'.format(txt['header']))
        output.write('{}\n'.format(npart))
        for i in range(npart):
            buff = _GRO_VEL_FMT.format(
                resnum[i],
                resname[i],
                atomname[i],
                atomnr[i],
                xyz[i, 0],
                xyz[i, 1],
                xyz[i, 2],
                vel[i, 0],
                vel[i, 1],
                vel[i, 2])
            output.write('{}\n'.format(buff))
        box = ' '.join([_GRO_BOX_FMT.format(i) for i in txt['box']])
        output.write('{}\n'.format(box))
