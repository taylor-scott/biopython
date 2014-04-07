from Bio.Seq import Seq
from Bio.PyEMBOSS.Usa import Usa


def getorf(utr, min=None, max=None):
    frames = (0, 1, 2)
    START = "ATG"
    STOP = ("TGA", "TAA", "TAG")
    orf_list = []
    sequence = str(utr.sequence)
    for frame in frames:
        i = 0
        nuc_seq = sequence[frame:]
        orf_codons = []
        start_pos = None
        stop_pos = None
        while i < len(nuc_seq)/3 + 1:
            codon = nuc_seq[i*3:i*3+3]
            if not orf_codons and codon == START:
                orf_codons.append(codon)
                start_pos = 1 + i*3 + frame
            elif orf_codons and not codon in STOP:
                orf_codons.append(codon)
            elif orf_codons and codon in STOP:
                orf = "".join(orf_codons)
                if (len(orf) >= min or min is None) and\
                   (len(orf) <= max or max is None):
                    stop_pos = i*3 + frame
                    orf_usa = Usa(
                        utr.db,
                        "{}_{}".format(utr.id, len(orf_list)+1),
                        start_pos, stop_pos
                    )
                    seq = Seq(orf)
                    orf_usa.sequence = seq
                    orf_list.append(orf_usa)
                    orf_codons = []
                else:
                    orf_codons = []
            i += 1
        if orf_codons:
            orf = "".join(orf_codons)
            if (len(orf) >= min or min is None) and\
               (len(orf) <= max or max is None):
                stop_pos = len(sequence)
                orf_usa = Usa(
                    utr.db,
                    "{}_{}".format(utr.id, len(orf_list)+1),
                    start_pos,
                    stop_pos
                )
                seq = Seq(orf)
                orf_usa.sequence = seq
                orf_list.append(orf_usa)
                orf_codons = []
    return orf_list
