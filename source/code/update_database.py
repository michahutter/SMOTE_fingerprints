import ftplib
import gzip
import os
import mysql.connector
from mysql.connector import Error
from mysql.connector.constants import ClientFlag


def update():
    """ Used for initializing the MySQL database

    This will update the data into the tables 'cid_smiles' and 'sid_map'. Those tables are needed to extract
    the data belonging to the SureChEMBL database. Using the FTP services of pubchem this will download the
    latest versions of the files 'CID-SMILES.csv' and 'SID-Map.csv', into the folder MySQL allows loading data
    from. The data contained in these files will then be inserted into MySQL tables with fitting columns. After
    this process is finished, the files will be deleted since they take quite a lot of space unzipped ~25 GB.
    """
    try:
        # Needed to establish connection to the local database.
        conn = mysql.connector.connect(host='localhost',
                                       database='surechembl_data',
                                       user='root',
                                       password='NPassWord_1411',
                                       client_flags=[ClientFlag.LOCAL_FILES])

        get_upload_path = """SHOW GLOBAL VARIABLES LIKE 'secure_file_priv';"""

        # Create cursor object to be able to pass queries
        if conn.is_connected():
            cursor = conn.cursor()

            # Getting path of directory with allowance to upload
            cursor.execute(get_upload_path)
            upload_status = cursor.fetchall()
            upload_path = upload_status[0][1].replace('\\', '/')

            # Other queries for updating database
            clear_cid_smiles = """TRUNCATE TABLE cid_smiles"""

            clear_sid_map = """TRUNCATE TABLE sid_map"""

            load_cid_smiles = """load data infile '""" + upload_path + """CID-SMILES.csv' 
                                 into table cid_smiles
                                 fields terminated by '\t'
                                 lines terminated by '\n'"""

            load_sid_map = """load data infile '""" + upload_path + """SID-Map.csv' 
                              into table sid_map
                              fields terminated by '\t'
                              lines terminated by '\n'
                              (sid, origin_name, origin_id, @cid)
                              SET 
                                  cid = NULLIF(@cid,'');"""

            # Get raw data
            print('\nStart downloading and unzipping substance data from PubChem.')
            get_file('/pubchem/Substance/Extras/', 'SID-Map.gz', 'SID-Map.csv', upload_path)
            print('\nFinished downloading and unzipping substance data from PubChem.')
            print('\nStart downloading and unzipping smiles data from PubChem.')
            get_file('/pubchem/Compound/Extras/', 'CID-SMILES.gz', 'CID-SMILES.csv', upload_path)
            print('\nFinished downloading and unzipping smiles data from PubChem.')

            # Make the tables clean
            cursor.execute(clear_sid_map)
            conn.commit()
            print('\nTable containing substance data cleaned')
            cursor.execute(clear_cid_smiles)
            conn.commit()
            print('\nTable containing smiles data cleaned')

            # Load data into tables
            cursor.execute(load_sid_map)
            conn.commit()
            print('\nTable containing substance data repopulated')
            cursor.execute(load_cid_smiles)
            conn.commit()
            print('\nTable containing smiles data repopulated')

            # Removing the raw files
            if os.path.exists(upload_path + 'CID-SMILES.gz'):
                os.remove(upload_path + 'CID-SMILES.gz')
                
            if os.path.exists(upload_path + 'CID-SMILES.csv'):
                os.remove(upload_path + 'CID-SMILES.csv')
                
            if os.path.exists(upload_path + 'SID-Map.gz'):
                os.remove(upload_path + 'SID-Map.gz')
                
            if os.path.exists(upload_path + 'SID-Map.csv'):
                os.remove(upload_path + 'SID-Map.csv')

    except Error as e:
        # Error occurred when connecting
        print("Error from server:", e)

    finally:
        # closing database connection.
        if conn.is_connected():
            cursor.close()
            conn.close()


def get_file(path, filename, end_file, upload_path):
    """ Uses the FTP-service of pubchem to download files and unzips them.

    Parameters
    ----------
    path : str
        Path of file to be downloaded
    filename : str
        Name of file in directory given by path that should be downloaded
    end_file : str
        Name of unzipped file
    upload_path : str
        Path to directory in which file should be downloaded and unzipped in.
    """

    ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
    ftp.login()
    ftp.cwd(path)
    ftp.retrbinary("RETR " + filename, open(upload_path + filename, 'wb').write)
    ftp.quit()

    with gzip.open(upload_path + filename, 'rb') as in_file, \
            open(upload_path + end_file, 'wb') as out_file:
        while True:
            # Read in chunks or it breaks
            block = in_file.read(65536)
            if not block:
                break
            else:
                out_file.write(block)
