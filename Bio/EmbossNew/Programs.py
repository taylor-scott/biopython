from Bio.Emboss import Db
import os
import re

database = None


def set_db_dir(directory):
    try:
        global database
        curr = os.getcwd()
        os.chdir(os.path.abspath(directory))
        os.chdir(curr)
        if database is None or database.db_dir != os.path.abspath(directory):
            db = Db.Database(os.path.abspath(directory))
            database = db
    except OSError:
        raise OSError("Could not find the directory {}".format(directory))


def _parse_usa(usa_str):
    usa_pattern = re.compile(
        "([A-Za-z_0-9]+):([A-Za-z0-9_]+)(?:\[(\d+)?:(\d+)?])?"
    )
    valid_usa = usa_pattern.search(usa_str)
    if valid_usa:
        usa_parts = list(valid_usa.groups())
        usa_parts[-2:] = [None if n is None else int(n)
                          for n in usa_parts[-2:]]
        usa = Db.Usa(*usa_parts)
        return usa
    else:
        print "ERROR: {} is not a valid USA".format(usa_str)
        raise ValueError


def seqret(usa_str):
    global database

    parsed_usa = _parse_usa(usa_str)
    ret_usa = database.get_sequence(parsed_usa)
    return ret_usa
