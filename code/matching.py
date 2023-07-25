import contextlib
import json
import math
import mysql.connector
import requests
from mysql.connector import Error
from rdkit import DataStructs
from rdkit.DataStructs import ExplicitBitVect

from source.code import fingerprints
from source.code import utils


# noinspection PyBroadException
def get_matches(fp, fp_type, out_file, names, filtered=True):
    """ Used to find matches for a combined fingerprint in the SureChEMBL data.

    For the matching process a MySQL database is used. It contains pre-generated fingerprints to all SureChEMBL
    compounds for all available fingerprint types in this project: 'maccs', 'morgan' 'Atom Pairs', 'Topological Torsion'
    and '2D Pharmacophore'. The ten compounds with the highest similarity scores will be written to an outfile,
    together with their SID, SureChEMBL ID similarity score and all patent IDs associated with it. Patent Ids are
    retrieved using the PUG-REST service of the Pubchem website.

    Parameters
    ----------
    fp : fingerprint
        Combined fingerprint
    fp_type : str
        Chosen fingerprint type
    out_file : str
        Name of the output file
    names : list
        Names of compounds used for creating combined fingerprint
    filtered : bool
        True if filtering should be used, false otherwise
    """

    try:
        # Needed to establish connection to the local database.
        conn = mysql.connector.connect(host='localhost',
                                       database='surechembl_data',
                                       user='root',
                                       password='NPassWord_1411')

        # Table names cannot start with numbers
        if fp_type.lower() == '2d pharmacophore':
            fp_type = 'pharmacophore'

        # Queries used to get the fingerprints fitting to fp_type
        get_fps = """   SELECT surechembl_id, """ + fp_type.lower().replace(' ', '_') + """
                        FROM """ + fp_type.lower().replace(' ', '_') + """_fps;"""

        get_fps_filtered = """  SELECT schembl_smiles.surechembl_id, """ + fp_type.lower().replace(' ', '_') + """
                                FROM """ + fp_type.lower().replace(' ', '_') + """_fps, schembl_smiles
                                WHERE length <= 200
                                AND mw <= 600
                                AND one_molecule = 1
                                AND schembl_smiles.surechembl_id = """ + fp_type.lower().replace(' ', '_') + """_fps
                                .surechembl_id; """

        out_file = utils.check_extension_output(out_file, '.asc')
        if out_file is None:
            return False

        # Create cursor object to be able to pass queries.
        if conn.is_connected():
            cursor = conn.cursor()

            # Executing the fetching, depending on whether filtering is true or false.
            print('Start fetching SureChEMBL entries.\n')

            if filtered:
                cursor.execute(get_fps_filtered)
                data = cursor.fetchall()
            else:
                cursor.execute(get_fps)
                data = cursor.fetchall()

            # Readjust name for regeneration of fingerprints
            if fp_type.lower() == 'pharmacophore':
                fp_type = '2d pharmacophore'

            print('Finished fetching SureChEMBL entries.\n')

            # Regenerate the fingerprints and store the ones with top 10 highest scores
            print('Start comparing combine fingerprint with SureChEMBL entries.\n')
            top10 = [[0, None, None]] * 10
            all_scores = list()
            length = len(data)

            # getting list of SCHEMBL IDs of used compounds to filter them out
            # continue even if fails ??
            try:
                dupes = get_schembls(names)
            except:
                print('\nWARNING: Unable to finish process of finding duplicate entries.\n'
                      'Results may contain input compounds.\n')
                dupes = list()

            # Adjust fp only once if necessary
            if fp_type.lower() == 'morgan' or fp_type.lower() == '2d pharmacophore':
                fp = DataStructs.FoldFingerprint(fp, 8)
            elif fp_type.lower() == 'atom pairs' or fp_type.lower() == 'topological torsion':
                help1 = ExplicitBitVect(fp.GetLength())
                help1.SetBitsFromList(list(fp.GetNonzeroElements().keys()))
                fp = help1

            for i in range(length):
                curr_fp = fingerprints.regenerate_fingerprint(data[i][1], fp_type)
                if fp_type.lower() == 'maccs':
                    score = DataStructs.FingerprintSimilarity(fp, curr_fp, metric=DataStructs.TanimotoSimilarity)
                elif fp_type.lower() == 'morgan' or fp_type.lower() == '2d pharmacophore':
                    f_curr_fp = DataStructs.FoldFingerprint(curr_fp, 8)
                    score = DataStructs.FingerprintSimilarity(fp, f_curr_fp, metric=DataStructs.TanimotoSimilarity)

                else:
                    # Before folding they need to be turned into ExplicitBitVects
                    help2 = ExplicitBitVect(curr_fp.GetLength())
                    help2.SetBitsFromList(list(curr_fp.GetNonzeroElements().keys()))

                    # folding was uneccessary but is possible here
                    curr_fp = help2
                    score = DataStructs.FingerprintSimilarity(fp, curr_fp, metric=DataStructs.TanimotoSimilarity)

                # Used to generate a file just containing all scores calculated
                all_scores.append(score)

                if data[i][0] not in dupes:

                    if score > top10[9][0]:
                        top10.append([score, curr_fp, data[i][0]])
                        top10.sort(reverse=True, key=lambda x: x[0])  # Only consider first element when sorting
                        top10 = top10[:-1]

                    if i % (math.floor(length / 10)) == 0:
                        print(str(math.ceil((i / length) * 100)) + '%')

            print('\nStart fetching patent IDs and writing outfile.\n')

            # Ignore if no internet connection is possible

            url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sourceid/SureChEMBL/"
            for i in range(10):
                url = url + top10[i][2] + ','
            url = url[:-1]
            url = url + '/xrefs/PatentID/JSON'

            try:
                page = requests.get(url)
                patents = json.loads(page.text)
                patents = patents["InformationList"]["Information"]
            except:
                patents = list()
                print('\nWARNING: Unable to retrieve Patent IDs from PUBCHEM.\n'
                      'Continuing without them.\n')

            write_match_file_surechembl(out_file, patents, top10, fp_type)
            write_all_scores_file(out_file, all_scores)

    except Error as e:
        # Error occurred when connecting
        conn = None
        print("Error from server:", e)
        return False

    finally:
        # closing database connection.
        if conn.is_connected():
            cursor.close()
            conn.close()

    return True


def get_matches_external_base(fp, fp_type, out_file, data_file):
    # check if fingerprint file, if so use it. If no generate it and use the resulting data
    print('\nStart fetching fingerprints from custom file.')

    # Do not need error in this case when it fails.
    with contextlib.redirect_stdout(None):
        data = fingerprints.read_byte_fingerprint_file(data_file)
    if data is None:
        print('\nFile did not contain fingerprints, try generating fingerprints of ' + fp_type + ' type.')
        outfile = utils.cut_extension(data_file)
        outfile = outfile + '_fingerprints.asc'
        data = fingerprints.create_fingerprints_output(data_file, outfile, fp_type)

        if data is None:
            print('\nERROR: Provided input file does not contain valid compounds.\n'
                  'Either use a file containing SMILES and corresponding names.\n'
                  'or use a fingerprint file as input file.\n')
            return False

        else:
            print('\nCreation of fingerprint successful.')

    # Matching process
    print('\nStart comparing combine fingerprint with custom compounds.\n')
    top10 = [[0, None, None]] * 10
    all_scores = list()
    length = len(data[0])

    # Adjust fp only once if necessary
    if fp_type.lower() == 'morgan' or fp_type.lower() == '2d pharmacophore':
        fp = DataStructs.FoldFingerprint(fp, 8)
    elif fp_type.lower() == 'atom pairs' or fp_type.lower() == 'topological torsion':
        help1 = ExplicitBitVect(fp.GetLength())
        help1.SetBitsFromList(list(fp.GetNonzeroElements().keys()))
        fp = help1

    for i in range(length):
        curr_fp = data[len(data)-1][i]
        if fp_type.lower() == 'maccs':
            score = DataStructs.FingerprintSimilarity(fp, curr_fp, metric=DataStructs.TanimotoSimilarity)
        elif fp_type.lower() == 'morgan' or fp_type.lower() == '2d pharmacophore':
            f_curr_fp = DataStructs.FoldFingerprint(curr_fp, 8)
            score = DataStructs.FingerprintSimilarity(fp, f_curr_fp, metric=DataStructs.TanimotoSimilarity)

        else:
            # Before folding they need to be turned into ExplicitBitVects
            help2 = ExplicitBitVect(curr_fp.GetLength())
            help2.SetBitsFromList(list(curr_fp.GetNonzeroElements().keys()))

            # folding was not necessary but is possible here
            curr_fp = help2
            score = DataStructs.FingerprintSimilarity(fp, curr_fp, metric=DataStructs.TanimotoSimilarity)

        # Used to generate a file just containing all scores calculated
        all_scores.append(score)

        if score > top10[9][0]:
            top10.append([score, curr_fp, data[1][i]])
            top10.sort(reverse=True, key=lambda x: x[0])  # Only consider first element when sorting
            top10 = top10[:-1]

        if i % (math.floor(length / 10)) == 0:
            print(str(math.ceil((i / length) * 100)) + '%')

    print('\nFinished finding substances with highest similarities.')
    print('\nStart writing results.')

    write_match_file_external(out_file, top10, fp_type, data_file)
    write_all_scores_file(out_file, all_scores)

    return True


def write_match_file_surechembl(out_file, patents, top10, fp_type):
    """ Creates the file containing the results of the matching process performed on the SureChEMBL database.

    Helper function for get_matches()

    Parameters
    ----------
    out_file : str
        Name of the output file.
    patents : dict
        Contains information retrieved from the PUG-REST service of pubchem
    top10 : list of list
        Every list contains a list with the similarity score, the corresponding fingerprint and their SureChEMBL ID,
        in this order. The outer list is sorted by the similarity scores in descending order.
    fp_type : str
        Used fingerprint type.
    """

    with open(out_file, 'w+') as f:

        p = True
        if len(patents) < 1:
            p = False

        # Adding header and additional information.
        f.write('List of top 10 SureChEMBL compounds with highest similarity scores:\n\n')
        f.write('Fingerprint type:\n' + fp_type + '\n\n\n')

        for i in range(len(top10)):
            rank = '--------------- Top ' + str(i + 1) + ' ---------------\n\n'
            f.write(rank)

            if p:
                substance_data = 'SID:' + str(patents[i]["SID"]) + ',    SureChEMBL_ID:' + str(top10[i][2])
                substance_data = substance_data + ',   Similarity_score:' + str(top10[i][0]) + '\n'
                f.write(substance_data)

                patent_listing = 'PatentID: '
                if len(patents[i]) > 1:
                    for element in patents[i]["PatentID"]:
                        patent_listing = patent_listing + str(element) + ', '
                        # Formatting output file, don't let lines get  too long
                        if len(patent_listing) >= 100:
                            patent_listing = patent_listing + '\n'
                            f.write(patent_listing)
                            patent_listing = ''
                    # Cut of last unnecessary comma
                    patent_listing = patent_listing[:-2]
                    patent_listing = patent_listing + '\n\n\n'
                    f.write(patent_listing)
            else:
                substance_data = 'SureChEMBL_ID:' + str(top10[i][2])
                substance_data = substance_data + ',   Similarity_score:' + str(top10[i][0]) + '\n'
                f.write(substance_data)

    f.close()


def write_match_file_external(out_file, top10, fp_type, comp_name):
    """ Creates the file containing the results of the matching process performed on an external database.

        Helper function for get_matches_external_base()

        Parameters
        ----------
        comp_name : str
            Name of file containing compounds that matching is performed on.
        out_file : str
            Name of the output file.
        top10 : list of list
            Every list contains a list with the similarity score, the corresponding fingerprint and their SureChEMBL ID,
            in this order. The outer list is sorted by the similarity scores in descending order.
        fp_type : str
            Used fingerprint type.
        """

    with open(out_file, 'w+') as f:
        # Adding header and additional information.
        f.write('List of top 10 ' + comp_name + ' compounds with highest similarity scores:\n\n')
        f.write('Fingerprint type:\n' + fp_type + '\n\n\n')

        for i in range(len(top10)):
            rank = '--------------- Top ' + str(i + 1) + ' ---------------\n\n'
            f.write(rank)
            substance_data = 'Identifier:' + str(top10[i][2]) + ',   Similarity_score:' + str(top10[i][0]) + '\n\n\n'
            f.write(substance_data)

    f.close()


def write_all_scores_file(out_file, all_scores):
    """ Function used to create an overview of all calculated similarity scores during a matching process.

    Generated file will contain an enumeration for the similarity scores in the first column and the corresponding
    similarity scores in the second column.

    Parameters
    ----------
    out_file : str
        Name that was chosen for the initial matching out_file, new name will be based on this.
    all_scores : list
        Contains all generated similarity scores
    """

    # outfile = utils.cut_extension(outfile)
    for i in range(len(out_file)):
        if out_file[len(out_file) - 1] != '.':
            out_file = out_file[:-1]
        else:
            out_file = out_file[:-1]
            break

    out_file = out_file + '_all_scores.asc'

    with open(out_file, 'w+') as f:
        # Adding header and additional information.
        f.write('List of all similarity scores appearing while matching to SureChEMBL database:\n\n')
        for i in range(len(all_scores)):
            line = str(i + 1) + ' ' + str(all_scores[i]) + '\n'
            f.write(line)

        f.close()


def get_schembls(names):
    """ Get SureChEMBL Ids for a List of compounds using PubChems PUG-REST service.

    Helper function to filter out compounds that were already used in creation of the combined fingerprint.

    Parameters
    ----------
    names : list of strings
        Names of compounds

    Returns
    -------
    List
        All SureChEMBL IDs of the given compound names.
    """

    url1 = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'
    url2 = '/synonyms/JSON'

    schembls = list()
    for i in range(len(names)):
        url = url1 + names[i] + url2
        page = requests.get(url)
        if page.status_code == 200:
            synonyms = json.loads(page.text)
            synonyms = synonyms['InformationList']['Information'][0]['Synonym']
            schembl = [s for s in synonyms if 'SCHEMBL' in s]
            if len(schembl) == 1:
                schembls.append(schembl[0])

    return schembls


def read_matching_file(in_file):
    """ Reads generated matching files.

    Files generated with the matching command/get_matches function can be read with this function.

    Parameters
    ----------
    in_file : str
        Name of matching file to be read.

    Returns
    -------
    list of list
        First element of list contains original fingerprint type. Second element contains a list of the SureChEMBL Ids
        belonging to the compounds listed in the file.
    """

    try:
        if type(in_file) != str:
            raise TypeError

        schembls = list()
        with open(in_file, 'r') as file:
            while True:

                line = file.readline()
                if not line:
                    break

                if line == 'Fingerprint type:\n':
                    line = next(file)
                    fp_type = line[:-1]

                if line.find('SureChEMBL_ID:') >= 0:
                    split = line.split('SureChEMBL_ID:', 1)
                    split = split[1]
                    split = split.split(',', 1)
                    schembl = split[0]
                    schembls.append(schembl)

                elif line.find('Identifier:') >= 0:
                    split = line.split('Identifier:', 1)
                    split = split[1]
                    split = split.split(',', 1)
                    schembl = split[0]
                    schembls.append(schembl)

        return [fp_type, schembls]

    except FileNotFoundError:
        print('\nERROR: Input matching_file: ' + str(in_file) + ' does not exist.\n')
        return None
    except TypeError:
        print('\nERROR: Function read_matching_file received an argument of the wrong type.\n'
              ' Only strings representing files are allowed.\n')
        return None


def read_all_scores(in_file):
    """ Reads all_scores files generated during get_matches process.

    Files generated during the execution of the matching command/get_matches function can be read with this function.

    Parameters
    ----------
    in_file : str
        Name of all_scores file to be read.

    Returns
    -------
    list of list
        Every list of the list contains a line of the file being read. The first element in a list is the enumeration
        the second element is the similarity score belonging to the number in the first position.
    """

    try:
        if type(in_file) != str:
            raise TypeError

        data = list()
        with open(in_file, 'r') as f:
            while True:

                line = f.readline()
                if not line:
                    break

                llist = line.split()
                data.append(llist)

        return data[2:]

    except FileNotFoundError:
        print('\nERROR: Input all_scores_file: ' + str(in_file) + ' does not exist.\n')
        return None
    except TypeError:
        print('\nERROR: Function read_all_scores received an argument of the wrong type.\n'
              ' Only strings representing files are allowed.\n')
        return None
