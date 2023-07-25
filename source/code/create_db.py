import mysql.connector
from mysql.connector import Error
from mysql.connector.constants import ClientFlag


def create_db():
    """ Used for creating the MYSQL database structure.

    This command will create the database and the needed tables that will contain the SureChEMBL data."""

    try:
        conn = mysql.connector.connect(host='localhost',
                                       user='root',
                                       password='NPassWord_1411',
                                       client_flags=[ClientFlag.LOCAL_FILES])

        if conn.is_connected():
            cursor = conn.cursor()

            cursor.execute('CREATE DATABASE IF NOT EXISTS surechembl_data')

            # Tables
            cursor.execute("""CREATE TABLE IF NOT EXISTS surechembl_data.cid_smiles
                           (
                                cid int,
                                smiles varchar(10000),
                                PRIMARY KEY(cid)

                           );""")

            cursor.execute("""CREATE TABLE IF NOT EXISTS surechembl_data.sid_map
                                       (
                                            sid int,
                                            origin_name varchar(200),
                                            origin_id varchar(500),
                                            cid int,
                                            PRIMARY KEY(sid)

                                       );""")

            cursor.execute("""CREATE TABLE IF NOT EXISTS surechembl_data.atom_pairs_fps
                                       (
                                            surechembl_id varchar(200),
                                            atom_pairs blob,
                                            PRIMARY KEY(surechembl_id)

                                       );""")

            cursor.execute("""CREATE TABLE IF NOT EXISTS surechembl_data.maccs_fps
                                       (
                                            surechembl_id varchar(200),
                                            maccs tinyblob,
                                            PRIMARY KEY(surechembl_id)

                                       );""")

            cursor.execute("""CREATE TABLE IF NOT EXISTS surechembl_data.morgan_fps
                                                   (
                                                        surechembl_id varchar(200),
                                                        maccs blob,
                                                        PRIMARY KEY(surechembl_id)

                                                   );""")

            cursor.execute("""CREATE TABLE IF NOT EXISTS surechembl_data.pharmacophore_fps
                                                   (
                                                        surechembl_id varchar(200),
                                                        maccs blob,
                                                        PRIMARY KEY(surechembl_id)

                                                   );""")

            cursor.execute("""CREATE TABLE IF NOT EXISTS surechembl_data.schembl_smiles
                                                   (
                                                        surechembl_id varchar(200),
                                                        smiles varchar(10000),
                                                        length int,
                                                        mw float,
                                                        one_molecule tinyint,
                                                        PRIMARY KEY(surechembl_id)
                                                        
                                                   );""")

            cursor.execute("""CREATE TABLE IF NOT EXISTS surechembl_data.topological_torsion_fps
                                                   (
                                                        surechembl_id varchar(200),
                                                        topological_torsion blob,
                                                        PRIMARY KEY(surechembl_id)

                                                   );""")

    except Error as e:
        # Error occurred when connecting
        print("Error from server:", e)

    finally:
        # closing database connection.
        if conn.is_connected():
            cursor.close()
            conn.close()
