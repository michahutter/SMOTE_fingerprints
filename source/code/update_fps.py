import mysql.connector
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import SanitizeFlags
from mysql.connector import Error

from source.code.fingerprints import generate_fingerprints


def update():
    """ Used for initializing the MySQL database.

    Uses the data in the table 'schembl_smiles' to generate 'maccs', 'morgan', 'Atom Pairs' and 'Topological Torsion'
    fingerprints for all SureChEBML compounds. Since 2D Pharmacophore fingerprints appear to be taking longer than the
    other fingerprint types, those will be generated separated from those. For the four fingerprint types handled here
    four tables will be filled with the SureChEMBL ID of a compound and the byte representation of the generated
    fingerprint. Every table of those four represents a different fingerprint type.
    """

    try:
        # Needed to establish connection to the local database.
        conn = mysql.connector.connect(host='localhost',
                                       database='surechembl_data',
                                       user='root',
                                       password='NPassWord_1411')

        # Queries used to get the basic data: last updated timestamp and upload directory
        get_sid_smiles = """   SELECT surechembl_id, smiles
                                FROM   schembl_smiles
                                ORDER BY surechembl_id;"""

        # Create cursor object to be able to pass queries
        if conn.is_connected():
            cursor = conn.cursor()

            print('Retrieve data.')
            cursor.execute(get_sid_smiles)
            data = cursor.fetchall()

            cursor.execute("""SELECT Count(*) 
                              FROM topological_torsion_fps;""")

            j = cursor.fetchall()[0]

            print('Start generating fingerprints')
            for i in range(j, len(data)):

                temp1 = list()
                temp2 = list()
                temp3 = list()
                temp4 = list()

                mol = Chem.MolFromSmiles(data[i][1], sanitize=False)
                mol.UpdatePropertyCache(strict=False)
                Chem.SanitizeMol(mol, SanitizeFlags.SANITIZE_FINDRADICALS | SanitizeFlags.SANITIZE_KEKULIZE |
                                 SanitizeFlags.SANITIZE_SETAROMATICITY | SanitizeFlags.SANITIZE_SETCONJUGATION |
                                 SanitizeFlags.SANITIZE_SETHYBRIDIZATION | SanitizeFlags.SANITIZE_SYMMRINGS,
                                 catchErrors=True)

                temp1.append(generate_fingerprints([mol], 'maccs', None)[0].ToBinary())
                temp2.append(generate_fingerprints([mol], 'morgan', None)[0].ToBinary())
                temp3.append(generate_fingerprints([mol], 'atom pairs', None)[0].ToBinary())
                temp4.append(generate_fingerprints([mol], 'topological torsion', None)[0].ToBinary())

                add_new_maccs = """     INSERT INTO maccs_fps(surechembl_id, maccs)
                                        VALUES
                                            ('""" + str(data[i][0]) + """', %s);"""

                add_new_morgan = """    INSERT INTO morgan_fps(surechembl_id, morgan)
                                        VALUES
                                            ('""" + str(data[i][0]) + """', %s);"""

                add_new_atom_pairs = """INSERT INTO atom_pairs_fps(surechembl_id, atom_pairs)
                                        VALUES
                                            ('""" + str(data[i][0]) + """', %s);"""

                add_new_topo = """      INSERT INTO topological_torsion_fps(surechembl_id, topological_torsion)
                                        VALUES
                                            ('""" + str(data[i][0]) + """', %s);"""

                # In the case of this program taking up from where it left, it could have happend in between
                # fingerprint insertions

                try:
                    cursor.execute(add_new_maccs, temp1)
                except:
                    print('MACCS: ' + str(i) + 'already inserted')
                try:
                    cursor.execute(add_new_morgan, temp2)
                except:
                    print('Morgan: ' + str(i) + 'already inserted')
                try:
                    cursor.execute(add_new_atom_pairs, temp3)
                except:
                    print('Atom Pairs: ' + str(i) + 'already inserted')
                cursor.execute(add_new_topo, temp4)

                if i % 1000 == 0:
                    print(i)
                    conn.commit()

            conn.commit()

    except Error as e:
        # Error occurred when connecting
        print("Error from server:", e)

    finally:
        # closing database connection.
        if conn.is_connected():
            cursor.close()
            conn.close()


if __name__ == "__main__":
    update()
