import mysql.connector
from mysql.connector import Error
from rdkit import DataStructs
from rdkit.DataStructs import ExplicitBitVect

from source.code import utils
from source.code import fingerprints


def compare(fp_type, schembls, fps, names, out_file):
    """ Used to get average similarity scores for a set of compounds to the best hits of a combined fingerprint.

    Parameters
    ----------
    fp_type : str
        Used fingerprint type
    schembls : list
        SureChEMBL IDs of the compounds with the highest similarity scores to a combines fingerprint
    fps : list
        fingerprints of the compounds that get compared to the top compounds
    names : list
        Names of the compounds that get compared to the top compounds
    out_file : str
        Name of the output file
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

        out_file = utils.check_extension_output(out_file, '.asc')
        if out_file is None:
            return None

        # Create cursor object to be able to pass queries.
        if conn.is_connected():
            cursor = conn.cursor()

            schembl_entries = list()
            for i in range(len(schembls)):
                get_top = """   SELECT surechembl_id, """ + fp_type.lower().replace(' ', '_') + """
                                FROM """ + fp_type.lower().replace(' ', '_') + """_fps
                                WHERE surechembl_id LIKE '""" + schembls[i] + """';"""

                cursor.execute(get_top)
                curr_schembl = cursor.fetchall()
                schembl_entries.append([curr_schembl[0][0], curr_schembl[0][1]])

            # Readjust name for regeneration of fingerprints
            if fp_type.lower() == 'pharmacophore':
                fp_type = '2d pharmacophore'

            # Regenerating the actual fingerprints of the top 10 compounds
            for i in range(len(schembl_entries)):
                schembl_entries[i][1] = fingerprints.regenerate_fingerprint(schembl_entries[i][1], fp_type)

            scores_list = list()
            for i in range(len(fps)):
                score = 0
                for j in range(len(schembl_entries)):
                    if fp_type.lower() == 'maccs' or fp_type.lower() == 'morgan' or \
                            fp_type.lower() == '2d pharmacophore':

                        score = score + DataStructs.FingerprintSimilarity(fps[i], schembl_entries[j][1],
                                                                          metric=DataStructs.TanimotoSimilarity)

                    else:

                        help1 = ExplicitBitVect(fps[i].GetLength())
                        help1.SetBitsFromList(list(fps[i].GetNonzeroElements().keys()))
                        help2 = ExplicitBitVect(schembl_entries[j][1].GetLength())
                        help2.SetBitsFromList(list(schembl_entries[j][1].GetNonzeroElements().keys()))

                        score = score + DataStructs.FingerprintSimilarity(help1, help2,
                                                                          metric=DataStructs.TanimotoSimilarity)

                # Getting average score
                score = score / 10
                scores_list.append([names[i], score])

            # Sort by score so output is nicely ordered
            scores_list.sort(reverse=True, key=lambda x: x[1])
            avg = 0
            for element in scores_list:
                avg = avg + element[1]
            avg = avg / len(scores_list)

            write_compare_file(scores_list, avg, out_file)

    except Error as e:
        # Error occurred when connecting
        print("Error from server:", e)

    finally:
        # closing database connection.
        if conn.is_connected():
            cursor.close()
            conn.close()


def compare_external(hits, fps, names, out_file, database_fp, database_identifiers):
    out_file = utils.check_extension_output(out_file, '.asc')
    if out_file is None:
        return None

    database_entries = list()
    for i in range(len(hits)):
        idx = database_identifiers.index(hits[i])
        database_entries.append((database_identifiers[idx], database_fp[idx]))

    scores_list = list()
    for i in range(len(fps)):
        score = 0
        for j in range(len(database_entries)):
            score = score + DataStructs.FingerprintSimilarity(fps[i], database_entries[j][1],
                                                              metric=DataStructs.TanimotoSimilarity)

        # Getting average score
        score = score / 10
        scores_list.append([names[i], score])

    # Sort by score so output is nicely ordered
    scores_list.sort(reverse=True, key=lambda x: x[1])
    avg = 0
    for element in scores_list:
        avg = avg + element[1]
    avg = avg / len(scores_list)

    write_compare_file(scores_list, avg, out_file)


def write_compare_file(scores_list, avg, out_file):
    """ Used to write the output file for the compare function.

    Parameters
    ----------
    scores_list : list of list
        Contains names and corresponding similarity scores of compounds being compared to top ranking ones
    avg : float
        Average similarity score across all compounds. Provided by compare() function.
    out_file : str
        Name of the output file.
    """

    with open(out_file, 'w+') as f:
        # Header of out file
        f.write('Average Similarity scores of all compounds:\n\n')
        for i in range(len(scores_list)):
            line = scores_list[i][0] + ' ' + str(scores_list[i][1]) + '\n'
            f.write(line)

        f.write('\nAverage of all compounds:\n' + str(avg))

    f.close()


def read_compare_file(in_file):
    """ Reads generated compare files.

    Does not read the line at the end of the containing an overall average of all the found scores in the file.
    Only reads all the singular entries.

    Parameters
    ----------
    in_file : str
        Name of compare to be read.

    Returns
    -------
    list of list
        Every list contains a line from the file to be read. The first element in every list contains the name of the
        compound, the second element contains the average of the similarity scores of this compound to the top hits of a
        combined fingerprint.
    """

    try:
        # Argument is of wrong type
        if type(in_file) != str:
            raise TypeError

        # Actual reading of file
        data = list()
        with open(in_file, 'r') as f:
            while True:

                line = f.readline()
                if not line:
                    break

                llist = line.split()
                data.append(llist)

        return data[2:-3]

    except FileNotFoundError:
        print('\nERROR: Input compare_file: ' + str(in_file) + ' does not exist.\n')
        return None
    except TypeError:
        print('\nERROR: Function read_compare_file received an argument of the wrong type.\n'
              ' Only strings representing files are allowed.\n')
        return None
