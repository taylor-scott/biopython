from Bio.Seq import Seq
import os
import binascii


class Usa(object):
    def __init__(self, db, seq_id, start, stop, sequence=None):
        self.db = db
        self.id = seq_id
        self.start = start
        self.stop = stop
        self.sequence = sequence
        if self.start > 0:
            self.start -= 1

    def __str__(self):
        if not self.stop and not self.start:
            return "{}:{}".format(self.db, self.id)
        return "{}:{}[{}:{}]".format(
            self.db,
            self.id,
            *["" if n is None else n for n in (self.start, self.stop)]
        )

    def __len__(self):
        if self.sequence:
            return len(self.sequence)
        else:
            return self.stop - self.start + 1


class Database(object):
    def __init__(self, db_dir):
        self.db_dir = db_dir
        self.index_dict = {}
        self.fasta_dict = {}

    def __str__(self):
        return self.db_dir

    def _process_index(self, db_name, index_file):
        entries = {}
        with open(index_file, "rb") as f:
            f_size = int(
                binascii.hexlify(''.join([b for b in f.read(4)][::-1])),
                16
            )
            f.seek(4, 1)
            rec_length = int(
                binascii.hexlify(''.join([b for b in f.read(2)][::-1])),
                16
            )
            f.seek(300)
            for n in range((f_size-300)/rec_length):
                entry_id = ''.join([b for b in f.read(rec_length-10)])
                offset = int(
                    binascii.hexlify(''.join([b for b in f.read(4)][::-1])),
                    16
                )
                entries[entry_id] = offset
                f.seek(6, 1)
        self.index_dict[db_name] = entries
        return entries

    def _get_fasta_file(self, db_name, lkp_file):
        with open(lkp_file, "rb") as f:
            f.seek(-24, 2)
            fasta = f.read(24).strip('\0')
            self.fasta_dict[db_name] = fasta
        return fasta

    def _read_fasta(self, fasta_file, index, entry):
        with open(fasta_file, "r") as f:
            if not entry in index:
                raise LookupError("Cannot find entry: " + entry)
            loc = index[entry]
            f.seek(loc)
            first_line = f.readline()
            if first_line[0] != ">":
                print ("Error processing file."
                       "FASTA entry does not start with '>'")
                exit()
            sequence = ""
            end_of_entry = False
            while not end_of_entry:
                line = f.readline()
                if line and line[0] != ">":
                    sequence += line.strip()
                else:
                    end_of_entry = True
        return sequence

    def get_sequence(self, usa):
        initial_dir = os.getcwd()
        os.chdir(self.db_dir)
        if not usa.db in os.listdir("."):
            raise LookupError("Could not locate database: %s" % usa.db)
        os.chdir(usa.db)
        if usa.db in self.index_dict:
            index_fasta = self.index_dict[usa.db]
        else:
            idx_file = os.path.join(
                os.getcwd(),
                [f for f in os.listdir(".") if ".idx" in f][0]
            )
            index_fasta = self._process_index(usa.db, idx_file)

        if usa.db in self.fasta_dict:
            fasta = self.fasta_dict[usa.db]
        else:
            lkp_file = os.path.join(
                os.getcwd(),
                [f for f in os.listdir(".") if ".lkp" in f][0]
            )
            fasta = self._get_fasta_file(usa.db, lkp_file)

        if len(index_fasta) == 0:
            raise ValueError("Error processing database: No entries found")
        sequence = self._read_fasta(
            os.path.join(os.getcwd(), fasta),
            index_fasta, usa.id
        )
        if usa.start or usa.stop:
            sequence = sequence[usa.start:usa.stop]
        os.chdir(initial_dir)
        usa.sequence = Seq(sequence)
        return usa
