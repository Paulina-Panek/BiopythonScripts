from __future__ import print_function, division

import sys

try:
    from Bio import pairwise2
    from Bio.SubsMat import MatrixInfo as matlist
except ImportError as exception:
    print("[!] Could not import Biopython modules", file=sys.stderr)
    raise exception

def align_sequences(sequenceA, sequenceB, sl, **kwargs):
# Performs a global pairwise alignment between two sequences. Returns the alignment, the sequence identity and the residue mapping between both original sequences.

    def _calculate_identity(sequenceA, sequenceB, sl):
        """
        Returns the percentage of identical characters between two sequences.
        Assumes the sequences are aligned.
        """

        matches = [sequenceA[i] == sequenceB[i] for i in range(sl)]

        seq_id = (100 * sum(matches)) /sl
        return (seq_id)

    matrix = kwargs.get('matrix', matlist.blosum62)
    gap_open = kwargs.get('gap_open', -10.0)
    gap_extend = kwargs.get('gap_extend', -0.2)  #set as default in Clustal

    alns = pairwise2.align.globalds(sequenceA, sequenceB,
                                    matrix, gap_open, gap_extend,
                                    penalize_end_gaps=(True, True) )

    best_aln = alns[0]
    aligned_A, aligned_B, score, begin, end = best_aln
    #print(aligned_A)
    #print(aligned_B)

    # Calculate sequence identity
    seq_id = _calculate_identity(aligned_A, aligned_B, sl)
    #print(seq_id)
    #print(len(aligned_A))
    return ((aligned_A, aligned_B), seq_id)


#sequenceA = "MESVWIMKPWLRHLQQSSGPVNLTADGVARKGWVAGRPPLQEASWLLENFYMWLLLWPEVCSVTPAPRSQRAELCSDRQRPLSPPLSKGVCRPSHLRALTVLSPPPPFKNKEVRCRSAVVPRRGAAQPSSPQFHVSQRQRDEPLEQGQELSPLPDTGLEDTAPPVPQRALACSSLRLGFTAMACIFPSGGLRAEDAEGFGARVASAPSTFPARPPDTGCPRFSQKDALAGRVSSRQTRSTAPAPQPTQATGLAARQVPSSPARPPKSQQPSPAPATCAGARRGRSYRASGGLHAYPGPVAKPSVILQIGKCRAEMLEHVRRPHRHLLTEVSKQVERELKGLHRSVGKLEHNLDGYVPTGDSQRWKKSIKACLCRCQETIASLERWVKREMHVWREVFYRLERWADRLESMGGKYPVGNEPARHTVSVGVGGPESYCQETDGYDYTVSPYAITPPPAAAELPGQEPAEAQQYQPWVPGEDGQPSPGVDTQIFEDPREFLSHLEEYLRQVGGSEEFWLSQIQNHMNGPAKKWWEFKQGSVKNWVQFKKEFLQYSEGTLSREAIQRELELPQKQGEPLDQFLWRKRDLYQTLYVDADEEEVIQYVVGTLQPKLKRFLRHPLPKTLEQLIQRGKEGQDGLEQAAEPAGPHPPSEPEHGAITPAPTSESVASDRTQPE"
#sequenceB = "MELDHRTSGGLHAYPGPRGGQVAKPNVILQIGKCRAEMLEHVRRTHRHLLAEVSKQVERELKGLHRSVGKLESNLDGYVPTSDSQRWKKSIKACLCRCQETIANLERWVKREMHVWREVFYRLERWADRLESTGGKYPVGSESARHTVSVGVGGPESYCHEADGYDYTVSPYAITPPPAAGELPGQEPAEAQQYQPWVPGEDGQPSPGVDTQIFEDPREFLSHLEEYLRQVGGSEEYWLSQIQNHMNGPAKKWWEFKQGSVKNWVEFKKEFLQYSEGTLSREAIQRELDLPQKQGEPLDQFLWRKRDLYQTLYVDADEEEIIQYVVGTLQPKLKRFLRHPLPKTLEQLIQRGMEVQDDLEQAAEPAGPHLPVEDEAETLTPAPNSESVASDRTQPE"

#sl = len(sequenceA)

#seq_A_length = (len(sequenceA))
#align_sequences(sequenceA,sequenceB, sl)
