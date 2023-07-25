import mysql.connector
from mysql.connector import Error
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import SanitizeFlags
from rdkit.Chem import Descriptors


def update():
    """ Used for initializing the MySQL database

    Uses the data in the tables 'cid_smiles' and 'sid_map' to extract the SureChEMBL data. This is achieved by
    joining the tables together on the pubchem cids and filtering for compounds with the origin name 'SureChEMBL'.
    Joining is necessary, since cid_smiles contains the needed SMILES for fingerprint generation and sid_map contains
    the origin names. The results will then be stored in a table named 'schembl_smiles'. For later filtering purposes
    more information about the compounds will then be generated and stored in the schembl_smiles table as well. This
    extra data involves: amount of non-H-atoms in the SMILES, molecular weight of the compounds and the info whether
    they only consist of one molecule (SMILES contains no '.').
    """

    try:
        # Needed to establish connection to the local database.
        conn = mysql.connector.connect(host='localhost',
                                       database='surechembl_data',
                                       user='root',
                                       password='NPassWord_1411')

        # Create cursor object to be able to pass queries
        if conn.is_connected():
            cursor = conn.cursor()

            # Other queries for updating database
            clear_schembl_smiles = """          TRUNCATE TABLE schembl_smiles"""

            fill_schembl_smiles = """           INSERT INTO schembl_smiles (surechembl_id, smiles)
                                                SELECT sid_map.origin_id, cid_smiles.smiles
                                                FROM cid_smiles, sid_map
                                                WHERE cid_smiles.cid = sid_map.cid
                                                AND origin_name = 'SureChEMBL'"""

            cursor.execute(clear_schembl_smiles)
            conn.commit()
            print('\nTable containing smiles and SureChEMBL data cleaned')

            cursor.execute(fill_schembl_smiles)
            conn.commit()
            print('\nTable containing smiles and SureChEMBL data repopulated')

            # Retrieving freshly created data
            get_sid_smiles = """    SELECT surechembl_id, smiles
                                    FROM schembl_smiles
                                    ORDER BY surechembl_id;"""

            cursor.execute(get_sid_smiles)
            raw = cursor.fetchall()

            for i in range(len(raw)):
                smiles = raw[i][1]
                length = get_smiles_length(smiles)
                mw = get_mw(smiles)
                one_molecule = get_one_molecule(smiles)

                update_entry = """UPDATE schembl_smiles
                                  SET length = """ + str(length) + """,
                                      mw = """ + str(mw) + """,
                                      one_molecule = """ + str(one_molecule) + """
                                  WHERE surechembl_id = '""" + raw[i][0] + """';"""

                cursor.execute(update_entry)

                if i % 1000 == 0:
                    print(i)
                    conn.commit()

            conn.commit()

    except Error as e:
        # Error occured when connecting
        print("Error from server:", e)

    finally:
        # closing database connection.
        if conn.is_connected():
            cursor.close()
            conn.close()


def get_smiles_length(smiles):
    """ Determines the number of non-H-atoms in a SMILES.

    Parameters
    ----------
    smiles : str
        SMILES to count non-H-atoms on.

    Returns
    -------
    int
        Number of non-H-atoms in given SMILES.
    """

    count = 0
    for i in range(len(smiles)):
        if smiles[i].isalpha():
            count = count + 1

    h_amount = smiles.count('H')
    h_amount = h_amount + smiles.count('h')

    return count - h_amount


def get_mw(smiles):
    """ Determines the molecular weight of a compound.

    Parameters
    ----------
    smiles : str
        Compound to calculate molecular weight of.

    Returns
    -------
    float
        Exact molecular weight of given compound.
    """

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol, SanitizeFlags.SANITIZE_FINDRADICALS | SanitizeFlags.SANITIZE_KEKULIZE |
                     SanitizeFlags.SANITIZE_SETAROMATICITY | SanitizeFlags.SANITIZE_SETCONJUGATION |
                     SanitizeFlags.SANITIZE_SETHYBRIDIZATION | SanitizeFlags.SANITIZE_SYMMRINGS,
                     catchErrors=True)

    mw = Descriptors.ExactMolWt(mol)

    return mw


def get_one_molecule(smiles):
    """ Determines if given SMILES represents only one molecule.

    Parameters
    ----------
    smiles : str
        SMILES to determine whether it only represents one molecule.

    Returns
    -------
    True if SMILES only represents one molecule, False if not.
    """

    if '.' in smiles:
        return 0
    else:
        return 1
