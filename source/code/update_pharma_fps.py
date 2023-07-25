import mysql.connector
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import SanitizeFlags
from mysql.connector import Error

from source.code.fingerprints import generate_fingerprints


def update(filtered):
    """ Used for initializing the MySQL database.

    Uses the data in the table 'schembl_smiles' to generate '2D Pharmacophore' fingerprints. This process can take very
    long for longer SMILES. Which is why first the batch with properties desirable for possible matches is calculated.
    Which usually are compounds with shorter SMILES. The second batch might not be all that necessary to set up but can
    always be run in the background to slowly fill the whole database with all available data.
    """

    try:
        # Needed to establish connection to the local database.
        conn = mysql.connector.connect(host='localhost',
                                       database='surechembl_data',
                                       user='root',
                                       password='NPassWord_1411')

        if conn.is_connected():
            cursor = conn.cursor()

            get_filtered = """          SELECT surechembl_id, smiles
                                        FROM schembl_smiles
                                        WHERE length <= 200
                                        AND mw <= 600
                                        AND one_molecule = 1
                                        ORDER BY length ASC;"""

            get_rest = """              SELECT surechembl_id, smiles
                                        FROM
                                            (SELECT surechembl_id, smiles, length
                                            FROM schembl_smiles
                                            WHERE length > 200
                                            UNION
                                            SELECT surechembl_id, smiles, length
                                            FROM schembl_smiles
                                            WHERE one_molecule = 0
                                            UNION
                                            SELECT surechembl_id, smiles, length
                                            FROM schembl_smiles
                                            WHERE mw > 600
                                                    
                                        ORDER BY length ASC) AS sub;"""

            if filtered:
                load_fingerprints(get_filtered, cursor, conn, 1000)
                print('\n---------- Finished chunk 1/1 -----------\n')

            else:
                load_fingerprints(get_filtered, cursor, conn, 1000)
                print('\n---------- Finished chunk 1/2 -----------\n')
                load_fingerprints(get_rest, cursor, conn, 1)
                print('\n---------- Finished chunk 2/2 -----------\n')

    except Error as e:
        # Error occurred when connecting
        print("Error from server:", e)

    finally:
        # closing database connection.
        if conn.is_connected():
            cursor.close()
            conn.close()


def load_fingerprints(query, cursor, conn, count):
    """ Generated fingerprints and inserts them into 'pharmacophore_fps' table.

    Helper function for update_pharma_fps.update().

    Parameters
    ----------
    query : str
        MySQL query which data should be taken from 'schembl_smiles'.
    cursor : cursor
        For MySQL connection.
    conn : connector
        For MySQL connection.
    count : int
        Step size printed to console to follow progress of filling process.
    """

    cursor.execute(query)
    current = cursor.fetchall()

    cursor.execute("""SELECT Count(*) 
                      FROM pharmacophore_fps;""")
    i = cursor.fetchall()[0][0]

    cursor.execute("""SELECT Count(*)
                              FROM schembl_smiles
                              WHERE length <= 100
                              AND mw <= 600
                              AND one_molecule = 1;""")
    offset = cursor.fetchall()[0][0]

    if count == 1:
        i = i - offset
    if i >= offset:
        return


    for j in range(i, len(current)):
        temp = list()

        mol = Chem.MolFromSmiles(current[j][1], sanitize=False)
        mol.UpdatePropertyCache(strict=False)
        Chem.SanitizeMol(mol, SanitizeFlags.SANITIZE_FINDRADICALS | SanitizeFlags.SANITIZE_KEKULIZE |
                         SanitizeFlags.SANITIZE_SETAROMATICITY | SanitizeFlags.SANITIZE_SETCONJUGATION |
                         SanitizeFlags.SANITIZE_SETHYBRIDIZATION | SanitizeFlags.SANITIZE_SYMMRINGS,
                         catchErrors=True)

        temp.append(generate_fingerprints([mol], '2D Pharmacophore', None)[0].ToBinary())

        add_new_pharma = """    INSERT INTO pharmacophore_fps(surechembl_id, pharmacophore)
                                 VALUES
                                ('""" + str(current[j][0]) + """', %s);"""

        cursor.execute(add_new_pharma, temp)
        if j % count == 0:
            print(j)
            conn.commit()

    conn.commit()
