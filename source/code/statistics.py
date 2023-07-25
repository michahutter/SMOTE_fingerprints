import matplotlib.pyplot as plt
from rdkit.DataStructs.cDataStructs import ExplicitBitVect, SparseBitVect, LongSparseIntVect, IntSparseIntVect

from source.code import utils


def create_statistics_output_name(out_file):
    """ Creates the name for the statistics files.

    Uses the given name and adds the suffix 'statistic.asc'

    Parameters
    ----------
    out_file : str
        Name of the output file.

    Returns
    -------
    str
        The new out_file name
    """

    try:
        if type(out_file != str):
            out_file = str(out_file)

        if len(out_file) <= 0:
            raise NameError

        if out_file[-4:] == '.asc':
            res = out_file[:len(out_file)-4] + '_statistic.asc'
        else:
            res = out_file + '_statistics.asc'

    except NameError:
        print('\nERROR: Name for output file needs to be at least one symbol long.\n')
        res = None

    return res


def create_statistics_output(data_raw, out_file):
    """ Creates a file with a detailed overview of set bits in the given list of fingerprints.

    Parameters
    ----------
    data_raw : list of list
        First list has length four and contains all compound SMILES in the first list, all compound names in the
        second list, all byte representation sizes in the third and all byte representations in the fourth list.
        Elements sharing an index belong together. Data_raw is provided by fingerprints.read_byte_fingerprint_file()
        function.
    out_file : str
        Name of the output file.
    """

    data = list()
    for i in range(len(data_raw[1])):
        data.append([data_raw[1][i], data_raw[3][i]])

    if type(data_raw[3][0]) == ExplicitBitVect or (data_raw[3][0]) == SparseBitVect:
        length = len(data_raw[3][0])
    elif type(data_raw[3][0]) == IntSparseIntVect or (data_raw[3][0]) == LongSparseIntVect:
        length = data_raw[3][0].GetLength()
    else:
        length = None

    write_statistics_output(out_file, length, data)


# Creates a file containing the info about what components have been used,
# As well as how many times certain bits have been set.
# Can be used as a basis to create nice graphs representing this info.
def write_statistics_output(out_file, length, data, threshold=None):
    """  Writes the statistics files.

    Helper function for create_statistics_output()

    Parameters
    ----------
    out_file : str
        Name of the output file.
    length : int
        length of the given fingerprints
    data : list of list
        The list contains two list. The first list contains the compound names and the second contains the
        corresponding fingerprints.
    threshold : float
        Used threshold for combined fingerprint
    """

    if length is None:
        return
    names = [names[0] for names in data]
    fps = [fps[1] for fps in data]
    out_file = utils.check_extension_output(out_file, '.asc')
    if out_file is None:
        return

    header1 = 'Used compounds:\n'
    header2 = 'Information about how many times particular bits are set:\n'
    body1 = '\n'
    body2 = '\n'
    count = 0
    probs = [0] * length

    for i in range(len(names)):
        body1 = body1 + names[i] + '\n'
    body1 = body1 + '\n'

    for i in range(length):
        for j in range(len(fps)):
            if fps[j][i] != 0:
                count += 1

        body2 = body2 + 'Bit' + str(i) + ' ' + str(count/len(fps)) + '\n'
        probs[i] = count/len(fps)
        count = 0

    with open(out_file, 'wb') as f:
        f.write(header1.encode('ascii'))
        f.write(body1.encode('ascii'))
        f.write(header2.encode('ascii'))
        f.write(body2.encode('ascii'))

    f.close()

    create_bar_chart([x for x in range(length)], probs, 'Bits', 'Rate',
                     (out_file[:len(out_file)-4] + '.png'), threshold)


def create_bar_chart(x, y, x_label, y_label, out_name, threshold=None):
    """ Creates a bar chart for visualization in png format.

    Parameters
    ----------
    x : list
        Contains all values that should be used for x-axis
    y : list
        Contains all values that should be used for y-axis
    x_label : str
        Name of x-axis label
    y_label : str
        Name of y-axis label
    out_name : str.
         Name of the output file.
    threshold : float
        Used threshold for combined fingerprint
    """

    plt.bar(x, y, color='red')
    plt.title('Rate of set bits')
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    if threshold is not None:
        threshold = [threshold] * len(y)
        plt.plot(x, threshold, color='black')

    plt.savefig(out_name)
