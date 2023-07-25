from rdkit.Chem import AllChem as Chem, AllChem
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Torsions
from rdkit.Chem.AtomPairs.Pairs import pyScorePair
from rdkit.Chem.Pharm2D import Generate
from rdkit.Chem.Pharm2D.SigFactory import SigFactory
from rdkit.DataStructs.cDataStructs import ExplicitBitVect, SparseBitVect, LongSparseIntVect, IntSparseIntVect

from source.code.globals import PHARMACOPHORE_FEATURES
from source.code import utils


# ***************************************** Creation of fingerprint files *****************************************

def create_fingerprints_output(in_file, out_file, fp_type, arguments=None):
    """ Creates a file containing the byte representations of a chosen fingerprint type for compounds in in_file.

    Output file will contain 4 columns separated by whitespaces. First column contains the smiles of the compound,
    second column contains the name, third column contains the length of the byte representation and the fourth
    column contains the byte representation. Used fingerprint type has to be one of the following: 'maccs', 'morgan'
    'Atom Pairs', 'Topological Torsion' or '2D Pharmacophore'.

    Parameters
    ----------
    in_file : str
        Name of the input file.
    out_file : str
        Name of the output file.
    fp_type : str
        name of the chosen fingerprint type.
    arguments : optional str
        Optional arguments to change details of fingerprint calculation.

    Returns
    -------
    list of list
        A multidimensional list, each entry contains a list representing one row in the output data.
    """

    # Reading the file and getting contents
    try:
        out_file = utils.check_extension_output(out_file, '.asc')
        if out_file is None:
            return None
        data = read_compound_file(in_file)
        if data is None:
            return None

        # Extract SMILES
        smiles = list()

        indexes = list()
        for i in range(len(data)):
            res = Chem.MolFromSmiles(data[i][0])
            # Ignore the ones with invalid SMILES
            if res is not None:
                smiles.append(res)
            else:
                indexes.append(i)
                print('\nSMILES: ' + data[i][0] + ' is not valid.\nContinuing without that entry.')

        # Remove not used entries
        if len(indexes) != 0:
            for i in range(len(indexes)):
                data.pop(indexes[i])

        # Create fingerprints
        fps = generate_fingerprints(mols=smiles, fp_type=fp_type, args=arguments)
        if fps is None:
            return None

        # Get Binary pickle of the fingerprints
        fps_bin = list()
        for i in range(len(smiles)):
            fps_bin.append(fps[i].ToBinary())

        # Write data into the output file.
        write_fingerprint_output(out_file, fp_type, data, fps_bin)

        res = [list(list(zip(*data))[0]), list(list(zip(*data))[1]), fps]

    except TypeError:
        print('\nERROR: Data contained in input file does not match required format.\n'
              'The file must consist of separate lines, with compounds in SMILES\n'
              'format in first position and the corresponding names in second\n'
              'position, separated by a whitespace.\n')
        res = None

    return res


def read_compound_file(in_file):
    """ Used to read input files.

    Input files need to have two columns separated by tabs. First column needs to contain smiles second column needs
    to contain the corresponding name.

    Parameters
    ----------
    in_file : str
        Path to the input file.

    Returns
    -------
    list
        A list containing the data read the file. None if errors occurred.
    """

    try:
        # Argument is of wrong type
        if type(in_file) != str:
            raise TypeError

        # Actual reading of file
        data = list()
        with open(in_file, 'r') as in_file:
            while True:

                line = in_file.readline()
                if not line:
                    break

                llist = line.split()
                data.append(llist)

    # In case file could not be found
    except FileNotFoundError:
        print('\nERROR: Input file: ' + str(in_file) + ' does not exist.\n')
        data = None

    except TypeError:
        print('\nERROR: Function read_compound_file received an argument of the wrong type.\n'
              ' Only strings representing files are allowed.\n')
        data = None
    return data


# Creates fingerprint depending on chosen fp_type
# Has options for MACCS, Morgan/circular (ECFP4), Atom pairs, Topological Torsion
# 2D Pharmacophore
# need to take another look at how they all work specifically 2d pharmacophore needs signature factory
# takes a dict of arguments to feed the fingerprint methods that can be personalized
def generate_fingerprints(mols, fp_type, args):
    """ Creates fingerprints depending on chosen fingerprint type.

    Possible fp_types are: 'MACCS', 'Morgan'(circular (ECFP4)), 'Atom pairs', 'Topological Torsion'
    '2D Pharmacophore'

    Parameters
    ----------
    mols : list
        list of molecule objects created with rdkid MolFromSmiles().
    fp_type : str
        type of fingerprint that should be generated.
    args : str
        Optional argument for using different Pharmacophore features file.

    Returns
    -------
    list of fingerprints
        Containing the generated fingerprints.
    """

    try:
        if fp_type.lower() == 'MACCS'.lower():
            fps = [MACCSkeys.GenMACCSKeys(x) for x in mols]

        elif fp_type.lower() == 'Morgan'.lower():

            fps = [AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=2048) for x in mols]

        elif fp_type.lower() == 'Atom Pairs'.lower():
            fps = [Pairs.GetHashedAtomPairFingerprint(x, nBits=2048) for x in mols]

        elif fp_type.lower() == 'Topological Torsion'.lower():
            fps = [Torsions.GetHashedTopologicalTorsionFingerprint(x, nBits=2048) for x in mols]

        elif fp_type.lower() == '2D Pharmacophore'.lower():
            if args:
                fdef_name = args
                factory = ChemicalFeatures.BuildFeatureFactory(fdef_name)
                sig_factory = SigFactory(factory, trianglePruneBins=False)
                sig_factory.SetBins([(0, 2), (2, 5), (5, 8)])
                sig_factory.Init()
                sig_factory.GetSigSize()
                fps = list()
                for i in range(len(mols)):
                    fps.append(Generate.Gen2DFingerprint(mols[i], sig_factory))
            else:

                factory = ChemicalFeatures.BuildFeatureFactory(PHARMACOPHORE_FEATURES)
                sig_factory = SigFactory(factory, trianglePruneBins=False)
                sig_factory.SetBins([(0, 2), (2, 5), (5, 8)])
                sig_factory.Init()
                sig_factory.GetSigSize()
                fps = list()
                for i in range(len(mols)):
                    fps.append(Generate.Gen2DFingerprint(mols[i], sig_factory))

        else:
            raise NameError

    except NameError:
        print('\nError: No valid fingerprint type chosen!\n'
              'Please choose one of the following valid fingerprint types: MACCS,\n'
              'Morgan, Atom Pairs, Topological Torsion, 2D Pharmacophore.\n')
        fps = None

    return fps


def write_fingerprint_output(out_file, fp_type, data, fps_bin):
    """ Helper function, takes care of writing the output file for create_fingerprints_output().

    Arguments of this file are all passed from the create_fingerprints_output() function.

    Parameters
    ----------
    out_file : str
        output file specified in create_fingerprints_output().
    fp_type : str
        Fingerprint type used in create_fingerprints_output().
    data : list of list
        Containing data generated in create_fingerprints_output().
    fps_bin : list
        Containing the byte representations of fingerprints created in create_fingerprints_output().
    """

    # First line of output file: Header
    s = 'compound name size ' + fp_type
    with open(out_file, 'wb+') as f:
        f.write(s.encode('ascii'))

        # Generate output respecting Header
        for i in range(len(data)):
            line = '\n' + str(data[i][0]).strip('\n') + ' ' + str(data[i][1]).strip('\n') \
                   + ' ' + str(len(fps_bin[i])) + ' '
            line = line.encode('ascii')
            f.write(line)
            f.write(fps_bin[i])

    f.close()


# ***************************************** Reading fingerprint files *****************************************

def regenerate_fingerprint(fp, fp_type):
    """ Uses byte representation of fingerprints to recreate the molecule object using original datatype.

    Parameters
    ----------
    fp : bytearray
        Byte representation of the fingerprint that should be restored.
    fp_type : str
        The fingerprint type it originally was.

    Returns
    -------
    ExplicitBitVect
        If original fingerprint was of type MACCS or Morgan.
    IntSparseIntVect
        If original fingerprint was of type Atom Pairs.
    LongSparseIntVect
        If original fingerprint was of type Topological Torsion.
    SparseBitVect
        If original fingerprint was of type 2D Pharmacophore.
    """

    try:
        if fp_type.lower() == 'MACCS'.lower():
            res = ExplicitBitVect(fp)
        elif fp_type.lower() == 'Morgan'.lower():
            res = ExplicitBitVect(fp)
        elif fp_type.lower() == 'Atom Pairs'.lower():
            res = IntSparseIntVect(fp)
        elif fp_type.lower() == 'Topological Torsion'.lower():
            res = LongSparseIntVect(fp)
        elif fp_type.lower() == '2D Pharmacophore'.lower():
            res = SparseBitVect(fp)
        else:
            raise NameError

        return res

    except NameError:
        print('\nERROR: Function regenerate_fingerprint got an invalid fingerprint type as argument.\n')
        return None

    except:
        return None


# noinspection PyTypeChecker
def read_byte_fingerprint_file(in_file):
    """ Reads generated fingerprint files.

    Files generated with the fingerprints command/create_fingerprints_output function can be read with this function.
    Since the fingerprints are stored as a pickled byte representation, it makes sure that the right amount of bytes
    is read, by using the given length info. Only stops reading bytes for a fingerprint when the given length is
    reached.

        Parameters
        ----------
        in_file : str
            Path to the input file.

        Returns
        -------
        list of list
            First list has length four and contains all compound smiles in the first list, all compound names in the
            second list, all byte representation sizes in the third and all byte representations in the fourth list.
            Elements sharing an index belong together.
        """

    try:
        if type(in_file) != str:
            raise TypeError

        in_file = utils.check_extension_input(in_file, '.asc')
        if in_file is None:
            return None

        data = list()
        eof = False
        with open(in_file, 'rb') as bin_file:

            # Reading of header and checking its format.
            line = next(bin_file)

            res = check_fingerprint_byte_file(line)
            if not res:
                raise AssertionError

            line = line.decode().split(' ', 3)
            fp_type = line[3][:-1]

            # Reading of data contained in file
            for line in bin_file:

                if not line:
                    break

                if eof:
                    break

                llist = line.split(maxsplit=3)
                if len(llist) != 4:
                    raise AssertionError

                # Turning the first three elements back from bytes
                llist[0] = llist[0].decode()
                llist[1] = llist[1].decode()
                llist[2] = llist[2].decode()

                # Do not want to cut of the last byte of the last element that is not followed by new line
                if llist[3][-1] == 10:
                    llist[3] = llist[3][:-1]

                # In case the byte representations contain a newline, adding those lines until it has the right length
                if (len(llist[3])) != int(llist[2]):
                    while (len(llist[3])) != int(llist[2]):
                        if eof:
                            break
                        try:
                            next_line = next(bin_file)
                            if next_line[-1] == 10:
                                llist[3] = llist[3] + b'\x0a' + next_line[:-1]
                            else:
                                llist[3] = llist[3] + b'\x0a' + next_line

                        except:
                            llist[3] = llist[3] + b'\x0a'
                            eof = True

                data.append(llist)
            bin_file.close()

            for i in range(len(data)):

                data[i][3] = regenerate_fingerprint(data[i][3], fp_type)
                if data[i][3] is None:
                    return None

        # Transpose list of lists to make pass to distances.py easier.
        data = list(map(list, zip(*data)))
        return data

    except NameError:
        print('\nERROR: Input file: ' + str(in_file) + ' does not use a valid file extension\n')
        return None
    except AssertionError:
        print('\nERROR: Input file: ' + str(in_file) + ' does not have the correct format.\n')
        return None
    except StopIteration:
        print('\nERROR: For input file: ' + str(in_file) + ' end of file was reached\n'
                                                           'without retrieving reading data. File might be empty\n')
        return None
    except FileNotFoundError:
        print('\nERROR: Input file: ' + str(in_file) + ' does not exist.\n')
        return None
    except TypeError:
        print('\nERROR: Input argument for function read_byte_fingerprint_file was not of type string.\n')
        return None


def read_byte_header_type(in_file):
    """ Gets the fingerprint type used in fingerprint file.

    Files generated with the fingerprints command/create_fingerprints_output function contain the used fingerprint
    type in their header. This function reads that value from the header.

    Parameters
    ----------
    in_file : str
        Path to the input file

    Returns
    -------
    str
        Returns the fingerprint type specified in the header

    """

    try:
        if type(in_file) != str:
            raise TypeError

        with open(in_file, 'rb') as bin_file:

            # Reading of header and checking its format.
            line = next(bin_file)
            line = line.decode().split(' ', 3)
            if len(line) != 4:
                raise AssertionError

            fp_type = line[3][:-1]
            if fp_type.lower() not in ["maccs", "morgan", "atom pairs", "topological torsion", "2d pharmacophore"]:
                raise NameError

            bin_file.close()

        return fp_type

    except NameError:
        print('\nERROR: Input file: ' + str(in_file) + ' does not use a valid fingerprint type.\n')
        return None
    except FileNotFoundError:
        return None
    except TypeError:
        print('\nERROR: Input argument for function read_byte_header_type was not of type string.\n')
        return None


def check_fingerprint_byte_file(line):
    """ This function is used to check if a fingerprint file has the correct format.

    Parameters
    ----------
    line : bytearray
        Header of file to be checked

    Returns
    -------
    bool
        True if header has correct form, False otherwise
    """

    # Reading of header and checking its format.
    try:
        line = line.decode().split(' ', 3)
        if len(line) != 4:
            return False

        if line[0] != 'compound' or line[1] != 'name' or line[2] != 'size':
            return False

        fp_type = line[3][:-1]
        if fp_type.lower() not in ["maccs", "morgan", "atom pairs", "topological torsion", "2d pharmacophore"]:
            return False

        return True

    except:
        return False
