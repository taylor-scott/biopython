from Bio.Seq import Seq
from Bio.pyEmboss.Usa import Usa


def getorf(utr, min=None, max=None):
    '''
    Find all ORFs of the appropriate size in a Usa object sequence. Only finds
    ORFs on the sense strand. This mimics the output of
        getorf utr outfile -minsize min -maxsize max -find 3 -noreverse
    Does NOT include the stop codon

    Paramters
    ---------
    utr : A Usa object with a valid sequence
    min : minimum ORF size (nt)
    max : maximum ORF size (nt)

    Returns
    -------
    A list of ORFs as Usa objects with sequence information.

    Notes
    -----
    min and max sizes do not include a stop codon
    '''
    #3 reading frames
    frames = (0, 1, 2)
    START = "ATG"
    STOP = ("TGA", "TAA", "TAG")
    orf_list = []
    sequence = str(utr.sequence)
    #Iterate through each frame sequentially and find all ORFs
    for frame in frames:
        i = 0
        #Change frames
        nuc_seq = sequence[frame:]
        orf_codons = []
        start_pos = None
        stop_pos = None
        #Consider every codon in the frame
        while i < len(nuc_seq)/3 + 1:
            codon = nuc_seq[i*3:i*3+3]
            if not orf_codons and codon == START:
                orf_codons.append(codon)
                #Add 1 because EMBOSS indices start at 1 not 0
                start_pos = 1 + i*3 + frame
            elif orf_codons and not codon in STOP:
                orf_codons.append(codon)
            elif orf_codons and codon in STOP:
                orf = "".join(orf_codons)
                #Correct size?
                if (len(orf) >= min or min is None) and\
                   (len(orf) <= max or max is None):
                    #Do not add 1 because EMBOSS stop indices are inclusive
                    stop_pos = i*3 + frame
                    #Create Usa object
                    #Add integer to the end of the ID to distinguish ORFs
                    orf_usa = Usa(
                        utr.db,
                        "{}_{}".format(utr.id, len(orf_list)+1),
                        start_pos, stop_pos
                    )
                    #A Bipython Seq object
                    seq = Seq(orf)
                    #Add sequence to Usa
                    orf_usa.sequence = seq
                    orf_list.append(orf_usa)
                    orf_codons = []
                else:
                    orf_codons = []
            i += 1
        #End of sequence, no stop codon --> Add rest of sequence as ORF
        #This is the same behavior as getorf
        if orf_codons:
            orf = "".join(orf_codons)
            #Right size?
            if (len(orf) >= min or min is None) and\
               (len(orf) <= max or max is None):
                stop_pos = len(sequence)
                #Create Usa object
                #Add integer to the end of ID to distinguish ORFs
                orf_usa = Usa(
                    utr.db,
                    "{}_{}".format(utr.id, len(orf_list)+1),
                    start_pos,
                    stop_pos
                )
                #A Biopython Seq object
                seq = Seq(orf)
                #Add sequence to Usa
                orf_usa.sequence = seq
                orf_list.append(orf_usa)
                orf_codons = []

    return orf_list
