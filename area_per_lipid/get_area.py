"""Calculate areas from xvg-files."""
import sys
import numpy as np
from matplotlib import pyplot as plt


plt.style.use('seaborn-talk')


LIPIDS = 64


def read_xvg_file(filename):
    """Return data in xvg file as numpy array."""
    data = []
    legends = []
    with open(filename, 'r') as fileh:
        for lines in fileh:
            if lines.startswith('@ s') and lines.find('legend') != -1:
                legend = lines.split('legend')[-1].strip()
                legend = legend.replace('"', '')
                legends.append(legend.lower())
            else:
                if lines.startswith('#') or lines.startswith('@'):
                    pass
                else:
                    data.append([float(i) for i in lines.split()])
    data = np.array(data)
    data_dict = {'step': data[:, 0]}
    for i, key in enumerate(legends):
        data_dict[key] = data[:, i+1]
    return data_dict


def plot_raw_data(all_data, match=None):
    """Make a plot of all the given data."""
    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    step_list = [data['step'] for data in all_data]
    step = np.concatenate(step_list).ravel()
    data_list = {}
    for data in all_data:
        if match is not None:
            keys = [i for i in data if i in match]
            for key in keys:
                if key not in data_list:
                    data_list[key] = []
                data_list[key].append(data[key])

        else:
            for key, val in data.items():
                if key != 'step':
                    if key not in data_list:
                        data_list[key] = []
                    data_list[key].append(val)
    for key, val in data_list.items():
        values = np.concatenate(val).ravel()
        line, = ax1.plot(step, values, label=key)
        ax1.axhline(y=np.average(values), ls=':', color=line.get_color())
    ax1.set_xlabel('Step')
    ax1.legend()
    fig.tight_layout()


def plot_xy(xdata, ydata, xlabel, ylabel):
    """Plot xy data in a new figure."""
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    line, = ax1.plot(xdata, ydata)
    ax1.axhline(y=np.average(ydata), color=line.get_color(), ls=':')
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(ylabel)
    fig.tight_layout()


def store_area(outfile, step, area):
    """Save the calculated area."""
    mat = np.column_stack((step, area))
    np.savetxt(outfile, mat, header='Step Area (nm^2)')


def main(xvgfiles):
    """Read xvg file and calculate area."""
    all_data = [read_xvg_file(i) for i in xvgfiles]
    plot_raw_data(all_data, match=['box-x', 'box-y', 'box-z'])
    area_list = [data['box-x'] * data['box-y'] for data in all_data]
    step_list = [data['step'] for data in all_data]
    area = np.concatenate(area_list).ravel()
    step = np.concatenate(step_list).ravel()
    plot_xy(step, area, 'Step', 'Area (xy) nm$^2$')
    plot_xy(step, 100.0 * area / LIPIDS, 'Step', 'Area per lipid (xy) Å$^2$')
    store_area('area.txt', step, area)
    plt.show()


if __name__ == '__main__':
    main(sys.argv[1:])
