#!/usr/env/python
from Bio.PyEMBOSS import Usa
from Bio.Seq import Seq
import os


def seqret(usa_str, database):
    '''
    Retrieve a sequence from an EMBOSS database. A valid USA is of
    the form database:ID([start:stop]) where parentheses denote an
    optional component.

    Paramters
    ---------
    usa_str : An EMBOSS format USA string
    database : An instance of PyEmboss.Database.Db()

    Returns
    -------
    A Usa object with the sequence from the EMBOSS database
    '''

    #Create empty USA from the input string
    usa = Usa.parse_usa(usa_str)
    initial_dir = os.getcwd()
    os.chdir(database.db_dir)
    #EMBOSS databases are in folders of the same name
    if not usa.db in os.listdir("."):
        raise LookupError("Could not locate database: %s" % usa.db)
    os.chdir(usa.db)

    #Has this database been processed previously?
    if usa.db in database.index_dict:
        #Use previous entries and records from entryname.idx
        index_fasta = database.index_dict[usa.db]
    else:
        #Get name of idx file
        idx_file = os.path.join(
            os.getcwd(),
            [f for f in os.listdir(".") if ".idx" in f][0]
        )
        #Read idx file to get record information
        index_fasta = database.process_index(usa.db, idx_file)

    #Has this database been processed prevously?
    if usa.db in database.fasta_dict:
        #Use prevous FASTA file location
        fasta = database.fasta_dict[usa.db]
    else:
        #Get name of lkp file
        lkp_file = os.path.join(
            os.getcwd(),
            [f for f in os.listdir(".") if ".lkp" in f][0]
        )
        #Read lkp file to get FASTA location
        fasta = database.get_fasta_file(usa.db, lkp_file)

    #Read the sequence from the FASTA file
    sequence = database.read_fasta(
        os.path.join(os.getcwd(), fasta),
        index_fasta, usa.id
    )

    #Slice the sequence if specified
    if usa.start or usa.stop:
        sequence = sequence[usa.start:usa.stop]

    #Return to original directory
    os.chdir(initial_dir)

    #Biopython Seq object of sequence
    usa.sequence = Seq(sequence)

    return usa
