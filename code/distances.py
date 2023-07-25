from rdkit.DataStructs.cDataStructs import ExplicitBitVect, SparseBitVect, LongSparseIntVect, IntSparseIntVect

from source.code import utils

# ***************************************** Creation of distance files *****************************************


def create_distances_output(fps1, fps2, names1, names2, out_file):
    """ Creates a file containing the names of the fingerprints being compared and their fingerprint distance.

    Names1 and fps1 correspond to the first infile, names2 and fps2 correspond to the second infile.
    Infiles have to follow the format of files created by the fingerprints command/create_fingerprints_output()
    function.
    The output file will have three columns separated by whitespaces. The first column will contain names from
    names1, the second column will contain the names from names2, the third column contains the Manhatten distance of
    the fingerprints in fps1 and fps2.

    Parameters
    ----------
    fps1 : list
        Contains fingerprints regenerated from the byte representation in the first infile.
    fps2 : list
        Contains fingerprints regenerated from the byte representation in the second infile.
    names1 : list
        Contains names from the first infile.
    names2 : list
        Contains names from the second infile.
    out_file : str
        Name of the output file.

    Returns
    -------
    list
        Containing the information found in the outfile in a list where each entry represents one row. None if error.
    """

    out_file = utils.check_extension_output(out_file, '.asc')
    if out_file is None:
        return None

    data = list()
    # First line of output file
    s = 'active inactive distance\n'
    s = s.encode('ascii')
    with open(out_file, 'wb') as f:
        f.write(s)

        # Generated output respecting first line
        for i in range(len(fps1)):
            for j in range(len(fps2)):
                res = manhatten_distance(fps1[i], fps2[j])
                if res is None:
                    return None
                line = str(names1[i]) + ' ' + str(names2[j]) + ' ' + str(res) + '\n'
                line = line.encode('ascii')
                data.append([names1[i], names2[j], res])
                f.write(line)

    f.close()
    return data


# Calculates Manhatten distance between two bitvector representations
# For fingerprints using counts, counts will simply be seen as set bits
# since only the presence of feature is important here.
def manhatten_distance(m1, m2):
    """ Calculates Manhatten distance between two bitvector representations.

    For fingerprints using counts, counts will simply be seen as set bits when calculating the Manhatten distance.
    Bitvectors have to be of the same type.

    Parameters
    ----------
    m1 : bitvector
    m2 : bitvector

    Returns
    -------
    int
        Manhatten distance of the bitvectors. None if errors occur.
    """

    try:
        # if they do not have the same type they should not be compared.
        if type(m1) != type(m2):
            raise TypeError

        # Used by MACCS, Morgan, Pharmacophore
        if (type(m1) == ExplicitBitVect) or (type(m1) == SparseBitVect):

            # fingerprints should have the same length.
            if len(m1) != len(m2):
                raise TypeError

            res = (m1 ^ m2).GetNumOnBits()

        elif (type(m1) == LongSparseIntVect) or (type(m1) == IntSparseIntVect):
            # fingerprints should have same length.
            if m1.GetLength() != m2.GetLength():
                raise TypeError

            # Used by Atom Pairs, Topological Torsion
            # GetNonzeroElements returns the index of the bits that are not set to 0.
            # If an index only appears in one of the two fingerprints we get 1,
            # however, if the same index is in both fingerprints we need 0.
            bits_dict1 = m1.GetNonzeroElements()
            bits_dict2 = m2.GetNonzeroElements()

            bits1 = [key for key, val in bits_dict1.items()]
            bits2 = [key for key, val in bits_dict2.items()]

            lst = bits1 + bits2

            no_dupe = [x for x in lst if lst.count(x) == 1]
            res = len(no_dupe)

        else:
            raise ValueError

    except TypeError:
        print('\nError: Fingerprints trying to be compared are not of the same kind.\n')
        res = None
    except ValueError:
        print('\nError: No valid fingerprint type used as input.\n')
        res = None

    return res


def read_distances(in_file):
    """ Reads distance files.

    Files generated with the distance command/create_distances_output function can be read with this function.
    Output of this function will reflect the format of the in_file.

    Parameters
    ----------
    in_file : str
        Path to the input file

    Returns
    -------
    list of list
        Every list of list contains names of actives, names of inactives, Manhatten distance of the two
    """

    try:
        if type(in_file) != str:
            raise TypeError

        data = list()
        with open(in_file, 'r') as in_file:
            next(in_file)
            while True:

                line = in_file.readline()
                if not line:
                    break

                llist = line.split()
                if len(llist) != 3:
                    raise AssertionError

                data.append([llist[0], llist[1], int(llist[2])])

    except FileNotFoundError:
        print('\nERROR: Input distances_file: ' + str(in_file) + ' could not be found.\n')
        data = None

    except TypeError:
        print('\nERROR: Argument for function read_distances was not of type str.\n')
        data = None

    except AssertionError:
        print('\nERROR: Input file: ' + str(in_file) + ' does not have required format\n')
        data = None

    return data
