from Bio.Seq import Seq
import os
import binascii


class Db(object):
    '''
    Allows the user to interface with an EMBOSS database. Includes tools to
    read idx, lkp, and FASTA files. Each instance represents a database
    directory. Processed idx and lkp files are stored in memory.

    Parameters
    ----------
    db_dir : an EMBOSS database directory

    Methods
    -------
    process_index : read an idx file
    get_fasta_file : read a lkp file
    read_fasta : read the specified entry from a FASTA file
    '''
    def __init__(self, db_dir):
        self.db_dir = db_dir
        self.index_dict = {}
        self.fasta_dict = {}

    def __str__(self):
        return self.db_dir

    def process_index(self, db_name, index_file):
        '''
        Reads an idx file to get the start locations of each sequence
        in the FASTA file.

        Parameters
        ----------
        db_name : the name of the database
        index_file : Path of the idx file

        Returns
        -------
        A dictionary of entries
            key : entry id
            value : start position in FASTA file
        '''
        entries = {}
        with open(index_file, "rb") as f:
            #Size of the idx file. Bytes belonging numbers are written in
            #reverse to this file
            f_size = int(binascii.hexlify(f.read(4)[::-1]), 16)
            f.seek(4, 1)
            #Each record is the same length
            rec_length = int(binascii.hexlify(f.read(2)[::-1]), 16)
            #The file has a 300 byte header
            f.seek(300)
            #Read each record
            for n in range((f_size-300)/rec_length):
                #The entry id
                entry_id = f.read(rec_length-10)
                #Start position in the FASTA file
                offset = int(binascii.hexlify(f.read(4)[::-1]), 16)
                entries[entry_id] = offset
                f.seek(6, 1)
        #Save processed file information to memory for future use
        self.index_dict[db_name] = entries
        return entries

    def get_fasta_file(self, db_name, lkp_file):
        '''
        Read a lkp file to get the name of the FASTA file.

        Paramters
        ---------
        db_name : database name
        lkp_file : path to the lkp file

        Returns
        -------
        The path to the FASTA file
        '''
        with open(lkp_file, "rb") as f:
            #FASTA path information is in the last 24 bytes of the file
            f.seek(-24, 2)
            fasta = f.read(24).strip('\0')
            #Save in memory for future retrievals
            self.fasta_dict[db_name] = fasta
        return fasta

    def read_fasta(self, fasta_file, index, entry):
        '''
        Retrieve the sequence of an entry from a FASTA file

        Paramters
        ---------
        fasta_file : path of the FASTA file
        index : index dictionary (returned from process index)
        entry : ID of the sequence to retrieve

        Returns
        -------
        Sequence of the requested entry as a string
        '''
        with open(fasta_file, "r") as f:
            if not entry in index:
                raise LookupError("Cannot find entry: " + entry)
            #Start location
            loc = index[entry]
            f.seek(loc)
            first_line = f.readline()
            #FASTA files begin with >
            if first_line[0] != ">":
                print ("Error processing file."
                       "FASTA entry does not start with '>'")
                exit()
            sequence = ""
            end_of_entry = False
            #Add to the sequence line by line. The end of the entry is found
            #when the next line begins with >
            while not end_of_entry:
                line = f.readline()
                if line and line[0] != ">":
                    sequence += line.strip()
                else:
                    end_of_entry = True
        return sequence
