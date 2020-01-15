from Bio.Align.Applications import MafftCommandline, ClustalOmegaCommandline
from Bio import AlignIO
import io
import numpy

def align_with_mafft(filepath, localpair=False, maxiterate=1000):
    """
    Align a file with the given filepath using MAFFT
    :param filepath: The file to align
    :param localpair: Should we use the l-insi method
    :return: The MAFFT alignment
    """
    mafft_cline = MafftCommandline(input=filepath, localpair=localpair, maxiterate=maxiterate)
    stdout, stderr = mafft_cline()
    align = AlignIO.read(io.StringIO(stdout), "fasta")
    return align


def align_with_clustal_omega(filepath):
    """
    Align a file with the given filepath using Clustal Omega
    :param filepath: The file to align
    :return: The Clustal Omega alignment
    """
    clustal_omega_cline = ClustalOmegaCommandline(infile=filepath)

    stdout, stderr = clustal_omega_cline()
    align = AlignIO.read(io.StringIO(stdout), "fasta")
    return align


def read_alignment(filepath, filetype):
    """
    Read an alignment from file
    :param filepath: Filepath to read from
    :param filetype: The type of file we're reading
    :return: The alignment from file
    """

    return AlignIO.read(filepath, filetype)


def write_alignment(file, filepath, filetype):
    """
    Write an alignment to file
    
    :param file: The alignment
    :param filepath: The filepath to write to
    :param filetype: The type of file to create
    """
    AlignIO.write(file, filepath, filetype)


def get_percent_identity(seq1, seq2, count_gaps=False):
    """
    Calculate the percent identity between two aligned sequences
    :param seq1: First sequence
    :param seq2: Second sequence
    :param count_gaps: Whether we should count a gap as a mismatch or no
    :return:
    """

    # Make sure the sequence content is a string
    seq1 = str(seq1)
    seq2 = str(seq2)

    # print (seq1)
    # print (seq2)

    matches = sum(aa1 == aa2 for aa1, aa2 in zip(seq1, seq2) if aa1 != "-" and aa2 != "-")

    # Set the length based on whether we want identity to count gaps or not
    # length = len(seq1) if count_gaps else min(len(seq1.replace("-", ""))- seq2.count("-"), len(seq2.replace("-", "")) - seq1.count("-"))
    if count_gaps:
        length = len(seq1)
    else:
        length = sum ([1 for (aa1, aa2) in zip(seq1, seq2) if aa1 != "-" and aa2 != "-"])

    # print ('matches ', matches)
    # print ('length ', length)

    pct_identity = 100.0 * matches / length

    return pct_identity

def get_percent_identity_of_alignment(records, realign=False, count_gaps=False):
    """
    Return a numpy array of all the pairwise percent identities in an alignment
    :param records:
    :param realign:
    :param count_gaps:
    :return:
    """
    count = 0
    percent_identity = 0
    percent_identities = []
    for num, record in enumerate(records):
        record = records[num]
        for other_record in records[num:]:
            if record.name != other_record.name:
                # print(record.name, other_record.name)
                if realign:
                    print ("Realigning sequences isn't implemented yet")

                else:
                    # percent_identity += get_percent_identity(record.seq, other_record.seq, count_gaps)
                    percent_identities.append(get_percent_identity(record.seq, other_record.seq, count_gaps))
                    # print("Percent identity of %s and %s is %d%%" % (
                    # record.name, other_record.name, get_percent_identity(record.seq, other_record.seq, count_gaps)))
                    count += 1

    percent_identities_array = numpy.array(percent_identities)
    return percent_identities_array

def score_match(pair, matrix):
    if pair not in matrix:
        return matrix[(tuple(reversed(pair)))]
    else:
        return matrix[pair]

def score_pairwise(seq1, seq2, matrix, gap_s, gap_e):
    score = 0
    gap = False
    for i in range(len(seq1)):
        pair = (seq1[i], seq2[i])
        if not gap:
            if '-' in pair:
                gap = True
                score += gap_s
            else:
                score += score_match(pair, matrix)
        else:
            if '-' not in pair:
                gap = False
                score += score_match(pair, matrix)
            else:
                score += gap_e
    return score