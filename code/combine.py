from rdkit.DataStructs.cDataStructs import TanimotoSimilarity

from source.code import statistics
from source.code.fingerprints import *


# noinspection PyArgumentList
def create_combine_output(ranks, fingerprints, fp_type, out_file, mode, amount, threshold=0.5):
    """ Creates files containing Information about the used compounds, used fingerprint type and combining result.

    The first file that gets created is a general overview of the combination process and the results. It uses
    the provided out_file name. The Second file that gets created, contains detailed information about how many
    times every single bit was set in the fingerprints list. It adds the suffix 'statistic' to the out_file name.
    The third file that gets created is a visual representation of the data in the second file. Representing every
    bit as a bar in a bar chart. Uses the same name as the second file but with a .png extension.

    Parameters
    ----------
    ranks : list of list
        Each list contains a list with the data found in one row of a ranks file. One list contains the name of an
        active, and it's summed up Manhatten distances to all inactives, in this order. Ranks is provided by
        ranks.read_ranking_output() function.
    fingerprints : list of list
        First list has length four and contains all compound SMILES in the first list, all compound names in the
        second list, all byte representation sizes in the third and all byte representations in the fourth list.
        Elements sharing an index belong together. Fingerprints is provided by
        fingerprints.read_byte_fingerprint_file() function.
    fp_type : str
        fingerprint type to be used provided by fingerprints.read_byte_header_type() function.
    out_file : str
        Name of the output file
    mode : str
        Choose between top k compounds or k-nearest to top compound.
    amount : int
        Number of compounds to be used in specified mode.
    threshold : float
        Determines how many times bit must be set in fingerprints list for this bit to be set in combined fingerprint.
        Values must be between 0 and 1.

    Returns
    -------
    Binary
        The binary representation of the resulting combined fingerprint.
    """

    try:
        if amount > len(fingerprints[3]):
            raise AssertionError

        out_file = utils.check_extension_output(out_file, '.asc')
        if out_file is None:
            return None

        if amount > len(ranks):
            raise ValueError

        # Check for mode
        if mode.lower() == "top".lower():

            data = get_top_fingerprints(ranks, fingerprints, amount)
            names = [n[0] for n in data]
            res_full = combine_fingerprints([d[1] for d in data], threshold)
            res = res_full[0].ToBinary()

            fingerprints = list()
            for i in range(len(data)):
                fingerprints.append(data[i][1].ToBinary())

        elif mode.lower() == "nearest".lower():

            data = get_k_nearest_fingerprints(ranks, fingerprints, amount)
            names = [n[0] for n in data]
            res_full = combine_fingerprints([d[1] for d in data], threshold)
            res = res_full[0].ToBinary()

            fingerprints = list()
            for i in range(len(data)):
                fingerprints.append(data[i][1].ToBinary())

        else:

            raise NameError

        write_combination_output(out_file, amount, mode, names, res, fp_type, threshold)
        out_stat_file = statistics.create_statistics_output_name(out_file)
        statistics.write_statistics_output(out_stat_file, res_full[1], data, threshold)

        return res

    except ValueError:
        print('\nERROR: Given amount of: ' + str(amount) + ' exceeds the entries in the ranking.'
              '\nAmount has to be smaller or equal to the entries in the ranking.\n'),
        return None

    except NameError:
        print('\nERROR: No Valid mode chosen.\n'
              'Valid modes are: top, nearest.\n')
        return None

    except AssertionError:
        print('\nERROR: Given amount of fingerprints cannot be bigger than actual amount of provided '
              'fingerprints: ' + str(len(fingerprints[3])) + '.\n')
        return None


def get_top_fingerprints(ranks, fingerprints, amount):
    """ Used to determine the wanted amount of fingerprints with top ranking.

    gets called if mode is set to 'top'

    Parameters
    ----------
    ranks : list of list
        Each list contains a list with the data found in one row of a ranks file. One list contains the name of an
        active, and it's summed up Manhatten distances to all inactives, in this order. Ranks is passed by
        combine.create_combine_output() function.
    fingerprints : list of list
        Each list contains a list with the data found in one row of a fingerprints file. One list contains the compound
        SMILES, name, fingerprint byte representation length and the fingerprint byte representation in this order.
        fingerprints is passed by combine.create_combine_output() function.
    amount : int
        Number of to compounds that should be chosen.
    Returns
    -------
    list of list
        Every list contains a list with the name of the top n compound and the corresponding byte representation
        of the fingerprint, in this order.
    """

    try:
        data = list()
        for i in range(amount):
            for j in range(len(fingerprints[1])):
                if fingerprints[1][j] == ranks[i][0]:
                    data.append([ranks[i][0], fingerprints[3][j]])
                    break

        if len(data) != amount:
            raise ValueError

    except ValueError:
        print('\nERROR: Fingerprints of top ' + str(amount) + ' ranked compounds are not\n'
                                                              'contained in the given fingerprint file.\n')
        data = None

    return data


# k-nearest are NOT determined via Manhatten distance but with Tanimoto index
def get_k_nearest_fingerprints(ranks, fingerprints, amount):
    """ Used to determine the top ranking compound and its k-nearest compounds.

    To determine the k-nearest compounds to the top the Tanimoto index is used.

    Parameters
    ----------
    ranks : list of list
        Each list contains a list with the data found in one row of a ranks file. One list contains the name of an
        active, and it's summed up Manhatten distances to all inactives, in this order. Ranks is passed by
        combine.create_combine_output() function.
    fingerprints : list of list
        Each list contains a list with the data found in one row of a fingerprints file. One list contains the compound
        SMILES, name, fingerprint byte representation length and the fingerprint byte representation in this order.
        fingerprints is passed by combine.create_combine_output() function.
    amount : int
        how many k-nearest compounds to the top should be respected

    Returns
    -------
    list of list
        Every list contains a list with the name of the top compound/k-nearest compound and the corresponding byte
        representation of the fingerprint, in this order.

    """
    try:
        data = list()
        fp = None
        # keeping the ranking order
        # getting the top one
        for i in range(len(ranks)):
            if fingerprints[1][i] == ranks[0][0]:
                fp = ([ranks[0][0], fingerprints[3][i]])
                break

        if fp is None:
            raise ValueError

        # looking for the k-nearest
        for i in range(len(fingerprints[1])):
            ti = TanimotoSimilarity(fp[1], fingerprints[3][i])
            data.append([fingerprints[1][i], fingerprints[3][i], ti])

        data.sort(key=lambda x: x[2], reverse=True)

    except ValueError:
        print('\nERROR: Fingerprint of top ranked compound: ' + str(ranks[0][0]) + ' is not\n'
              'contained in the given fingerprint file.\n')

        data = None

    return data[:amount + 1]


def combine_fingerprints(fp_list, threshold=0.5):
    """ Used to create a combined fingerprint from a list of fingerprints.

    Generally the combined fingerprint is created like a logical or operation on the fingerprint list.
    However, a bit is only set to 1, if the feature has a higher appearance rate in the fingerprint list
    than the given threshold.

    Parameters
    ----------
    fp_list : list
        list of fingerprint objects that should be combined
    threshold : float
        The minimum appearance rate bits needs to have in the fingerprint list, for it to be set to 1.

    Returns
    -------
    list
        First element is the combined fingerprint, second element is its length.

    """

    try:

        if threshold < 0 or threshold > 1:
            raise ValueError

        # Depending on fingerprint type different datatypes are used. Choose correct type and length.
        if type(fp_list[0]) == ExplicitBitVect:
            res_fp = ExplicitBitVect(len(fp_list[0]))
            length = len(fp_list[0])
        elif type(fp_list[0]) == IntSparseIntVect:
            res_fp = IntSparseIntVect(fp_list[0].GetLength())
            length = fp_list[0].GetLength()
        elif type(fp_list[0]) == LongSparseIntVect:
            res_fp = LongSparseIntVect(fp_list[0].GetLength())
            length = fp_list[0].GetLength()
        elif type(fp_list[0]) == SparseBitVect:
            res_fp = SparseBitVect(len(fp_list[0]))
            length = len(fp_list[0])
        else:
            raise NameError

        # Initial fingerprint has only 0s set 1s where amount of set bits in fingerprint list is >= threshold
        for j in range(length):
            count = 0
            for i in range(len(fp_list)):
                if fp_list[i][j] != 0:
                    count += 1
            if count / len(fp_list) >= threshold:
                res_fp[j] = 1

    except ValueError:
        print('\nERROR: Given value of threshold: ' + str(threshold) + ' is out of range'
              '\nValue of threshold needs to be between 0 and 1.\n')

        res_fp = None
        length = None

    except NameError:
        print('\nERROR: Invalid fingerprint type chosen.\n')
        res_fp = None
        length = None

    return [res_fp, length]


def write_combination_output(out_file, amount, mode, names, fp, fp_type, threshold):
    """ writes the combine output file.

    Helper function for create_combine_output() for writing the standard results file.

    Parameters
    ----------
    out_file : str
         Name of the output file
    amount : int
        Number of compounds to be used in specified mode.
    mode : str
        Choose between top k compounds or k-nearest to top compound.
    names : list
        Contains used compound names. Passed from get_top_fingerprints/get_k_nearest_fingerprints depending
        on chosen mode.
    fp : fingerprint
        the created combined fingerprint with the combine_fingerprints().
    fp_type : str
        Chosen fingerprint type. Retrieved from fingerprints.read_byte_header_type().
    threshold : float
        Chosen threshold in create_combine_output().
    """

    out_file = utils.check_extension_output(out_file, '.asc')
    if out_file is None:
        return

    header = 'Information about the combination process:\n\nUsed compounds:\n'
    body1 = ''
    fp_intro = 'Fingerprint Type:\n' + fp_type + '\n\nUsed threshold:\n' + str(threshold) +\
        '\n\nFingerprint byte representation with according length:\n' + str(len(fp)) + ' '

    if mode.lower() == 'top'.lower():
        for i in range(amount):
            body1 = body1 + 'Top' + str(i + 1) + '   ' + names[i] + '\n'
        body1 = body1 + '\n'

    elif mode.lower() == 'nearest'.lower():
        body1 = body1 + 'Top1' + '        ' + names[0] + '\n'
        for i in range(1, amount+1):
            body1 = body1 + str(i) + '-nearest' + '   ' + names[i] + '\n'
        body1 = body1 + '\n'

    else:
        print('\nNo Valid mode chosen.\n'
              'Valid modes are: top, nearest.\n')
        return

    with open(out_file, 'wb') as f:
        f.write(header.encode('ascii'))
        f.write(body1.encode('ascii'))
        f.write(fp_intro.encode('ascii'))
        f.write(fp)

    f.close()


def read_combine(combine_file):
    """ Reads generated combine files.

    Files generated with the combine command/create_combine_output function can be read with this function.


    Parameters
    ----------
    combine_file : str
         Path to the combine file to be read.

    Returns
    -------
    list
        First element is the used fingerprint type for the generated combined fingerprint, second element in
        the list is the byte representation of the generated combined fingerprint.
    """

    try:
        if type(combine_file) != str:
            raise TypeError

        data = list()
        eof = False
        with open(combine_file, 'rb') as in_file:
            while True:

                line = in_file.readline()
                if not line:
                    break

                if eof:
                    break

                if line.decode() == 'Used compounds:\n':
                    names = list()
                    line = next(in_file)
                    while line.decode() != '\n':
                        names.append(line.decode().split()[1])
                        line = next(in_file)

                if line.decode() == 'Fingerprint Type:\n':
                    fp_type = in_file.readline()[:-1].decode()
                    data.append(fp_type)

                if line.decode() == 'Fingerprint byte representation with according length:\n':
                    line = in_file.readline()
                    llist = line.split(maxsplit=1)

                    llist[0] = llist[0].decode()

                    if llist[1][-1] == 10:
                        llist[1] = llist[1][:-1]

                    if (len(llist[1])) != int(llist[0]):
                        while (len(llist[1])) != int(llist[0]):
                            try:
                                next_line = next(in_file)
                                if next_line[-1] == 10:
                                    llist[1] = llist[1] + b'\x0a' + next_line[:-1]
                                else:
                                    llist[1] = llist[1] + b'\x0a' + next_line

                            except:
                                llist[1] = llist[1] + b'\x0a'
                                eof = True

                    data.append(llist[1])

    except FileNotFoundError:
        print('\nERROR: Input combine_file: ' + str(combine_file) + ' could not be found.\n')
        data = None

    except TypeError:
        print('\nERROR: The function read_ranking_output got non str input.\n')
        data = None

    return [data, names]
