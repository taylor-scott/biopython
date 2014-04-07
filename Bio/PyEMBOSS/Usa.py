import re


class Usa(object):

    '''
    A container for an EMBOSS uniform sequence address (USA). Supports USAs of
    the form database:id[start:stop]. The components of the USA are stored as
    separate variable within the object. An additional variable, the sequence
    is stored as a Biopython Seq object.

    WARNING : EMBOSS subsequences are indexed differently than Python strings.
    EMBOSS indices begin at 1 are are inclusive. Python indices begin at 0 and
    the start index is inclusive while the stop index is exclusive. Always pass
    start and stop indices as expected by EMBOSS (i.e. add 1 to the indices).

    Parameters
    ----------
    db : The database
    seq_id : The ID of the sequence (in the database)
    start : The start position of the desired subsequence in EMBOSS form
    stop : The stop position of the desired subsequence
    utr_length : The sequence/subsequence

    Methods
    -------
    None

    Notes
    -----
    Pass start and stop indices as expected by EMBOSS. See the warning above
    for more information.
    '''

    def __init__(self, db, seq_id, start, stop, sequence=None):
        self.db = db
        self.id = seq_id
        self.start = start
        self.stop = stop
        self.sequence = sequence
        #Store the start index internally as a python index to ensure accurate
        #subsequence slicing
        if self.start > 0:
            self.start -= 1

    def __str__(self):
        if not self.stop and not self.start:
            return "{}:{}\n{}".format(
                self.db,
                self.id,
                self.sequence.__repr__()
            )
        else:
            #Display the start index as an EMBOSS index to ensure compatibility
            #with other EMBOSS programs
            adj_start = self.start
            if adj_start is not None:
                adj_start += 1
            return "{}:{}[{}:{}]\n{seq}".format(
                self.db,
                self.id,
                *["" if n is None else n for n in (adj_start, self.stop)],
                seq=self.sequence.__repr__()
            )

    def __len__(self):
        #More accurate in case the sequence changes
        if self.sequence:
            return len(self.sequence)
        #If the nucleotide sequence is not provided
        else:
            return self.stop - self.start


def parse_usa(usa_str):
    '''
    Create a Usa object from a USA string. The string is indexed with EMBOSS
    indices (see warning in Usa() docstring) to ensure combatibility with
    other EMBOSS programs.

    Parameters
    ----------
    usa_str : An EMBOSS format USA string

    Returns
    -------
    A Usa object
    '''
    usa_pattern = re.compile(
        "([A-Za-z_0-9]+):([A-Za-z0-9_]+)(?:\[(\d+)?:(\d+)?])?"
    )
    valid_usa = usa_pattern.search(usa_str)
    if valid_usa:
        usa_parts = list(valid_usa.groups())
        usa_parts[-2:] = [None if n is None else int(n)
                          for n in usa_parts[-2:]]
        usa = Usa(*usa_parts)
        return usa
    else:
        print "ERROR: {} is not a valid USA".format(usa_str)
        raise ValueError
