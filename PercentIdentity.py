from __future__ import print_function, division

from operator import itemgetter
import os
import sys
import tempfile
import warnings

try:
    from Bio import pairwise2
    from Bio.SubsMat import MatrixInfo as matlist
except ImportError as exception:
    print("[!] Could not import Biopython modules", file=sys.stderr)
    raise exception

#
def align_sequences(sequence_A, sequence_B, **kwargs):
    """
    Performs a global pairwise alignment between two sequences
    using the BLOSUM62 matrix and the Needleman-Wunsch algorithm
    as implemented in Biopython. Returns the alignment, the sequence
    identity and the residue mapping between both original sequences.
    """

    def _calculate_identity(sequenceA, sequenceB):
        """
        Returns the percentage of identical characters between two sequences.
        Assumes the sequences are aligned.
        """

        sa, sb, sl = sequenceA, sequenceB, len(sequenceA)
        matches = [sa[i] == sb[i] for i in range(sl)]
        seq_id = (100 * sum(matches)) / sl

        gapless_sl = sum([1 for i in range(sl) if (sa[i] != '-' and sb[i] != '-')])
        gap_id = (100 * sum(matches)) / gapless_sl
        return (seq_id, gap_id)

    #
    matrix = kwargs.get('matrix', matlist.blosum62)
    gap_open = kwargs.get('gap_open', -10.0)
    gap_extend = kwargs.get('gap_extend', -0.2)

    alns = pairwise2.align.globalds(sequence_A, sequence_B,
                                    matrix, gap_open, gap_extend,
                                    penalize_end_gaps=(False, False) )

    best_aln = alns[0]
    aligned_A, aligned_B, score, begin, end = best_aln

    # Calculate sequence identity
    seq_id, g_seq_id = _calculate_identity(aligned_A, aligned_B)
    print(seq_id, g_seq_id)
    return ((aligned_A, aligned_B), seq_id, g_seq_id)


sequenceA = "MELDHRTSGGLHAYPGPRGGQVAKPNVILQIGKCRAEMLEHVRRTHRHLLAEVSKQVERELKGLHRSVGKLESNLDGYVPTSDSQRWKKSIKACLCRCQETIANLERWVKREMHVWREVFYRLERWADRLESTGGKYPVGSESARHTVSVGVGGPESYCHEADGYDYTVSPYAITPPPAAGELPGQEPAEAQQYQPWVPGEDGQPSPGVDTQIFEDPREFLSHLEEYLRQVGGSEEYWLSQIQNHMNGPAKKWWEFKQGSVKNWVEFKKEFLQYSEGTLSREAIQRELDLPQKQGEPLDQFLWRKRDLYQTLYVDADEEEIIQYVVGTLQPKLKRFLRHPLPKTLEQLIQRGMEVQDDLEQAAEPAGPHLPVEDEAETLTPAPNSESVASDRTQPE"
sequenceB = "MGAGARRGGRPRPRGRIKSRRRLQGAELGAAALCSSPLRRTLPRSLAAAAAAGPKGPSWLHGTAGPRGLCSLPPLASPVAAATCAPMELDHMTSGGLHAYPGPRGGPAAKPNVILQIGKCRAEMLEHVRRTHRHLLTEVSKQVERELKGLHRSVGKLESNLDGYVPTGDSQRWKKSIKACLCRCQETIANLERWVKREMHVWREVFYRLEKWADRLESMGGKYPVGNDPARHTVSVGVGGPEGYCQEADAYDYTVSPYAITPPPASGELPGQEPVDAQQYPPWGAGEDGQPSPGVDTQIFEDPREFLSHLEEYLRQVGGSEEYWLSQIQNHMNGPAKKWWEFKQGSVKNWVEFKKEFLQYSEGTLSREAIQRELDLPQKQGEPLDQFLWRKRDLYQTLYVDAEEEEIIQYVVGTLQPKLKRFLRHPLPKTLEQLIQRGMEVQDGLEQTAEPTGPRLPTEDEAEALTPALTSESVTSDRTQPE"

align_sequences(sequenceA,sequenceB)
