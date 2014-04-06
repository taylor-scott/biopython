#!/usr/env/python
import os
import re
from Bio.EmbossNew.Database import Usa


def _parse_usa(usa_str):
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


def seqret(usa_str, database):
    parsed_usa = _parse_usa(usa_str)
    ret_usa = database.get_sequence(parsed_usa)
    return ret_usa
