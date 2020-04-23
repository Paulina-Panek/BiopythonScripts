# Paulina Panek
# April 2020
# Script parsing result GenPept file (.gp) to get name of the sequence, classification, sequence, length

from Bio import Entrez
Entrez.email = "ppanek@hpu.edu"
from Bio import SeqIO
from Bio import Align
from Bio.SubsMat.MatrixInfo import blosum62
from PercentIdentity import *  #imports all functions from PercentIdentity.py

def numberRecords(ListRecords):
    #Function prints number of records
    records = list(SeqIO.parse(ListRecords, "genbank"))
    print("Found %i records in initial file " % len(records))

def CheckIfDuplicate(first_sequence_name, second_sequence_name, first_sequence, second_sequence):
    # returns 0 (same sequences), 1 (not same sequences, or 3 (something went wrong, function didn't work

    return_value = 3

    # if same species AND length of sequence is the same, check if the sequence is the same
    if (first_sequence_name == second_sequence_name):
        if first_sequence == second_sequence:
            return_value = 0  #same sequences
        else:
            return_value = 1
    else:
        return_value = 1

    return(return_value)

def RemoveLike(protein_name):
#if protein has word like in it's name, returns 1
    ret_val = 0

    if "like" in protein_name:
        ret_val = 1
    return ret_val

def unknown_aas(sequence):
    #returns number of unknown amino acids in sequence
    X_in_sequence = 0

    if 'X' in sequence:
        X_in_sequence = X_in_sequence + 1
    return X_in_sequence

numberRecords("arc_sequences_04202020.gp")

file = open("AllSpecies.fasta", "w")

aligner = Align.PairwiseAligner()

human_sequence = "MELDHRTSGGLHAYPGPRGGQVAKPNVILQIGKCRAEMLEHVRRTHRHLLAEVSKQVERELKGLHRSVGKLESNLDGYVPTSDSQRWKKSIKACLCRCQETIANLERWVKREMHVWREVFYRLERWADRLESTGGKYPVGSESARHTVSVGVGGPESYCHEADGYDYTVSPYAITPPPAAGELPGQEPAEAQQYQPWVPGEDGQPSPGVDTQIFEDPREFLSHLEEYLRQVGGSEEYWLSQIQNHMNGPAKKWWEFKQGSVKNWVEFKKEFLQYSEGTLSREAIQRELDLPQKQGEPLDQFLWRKRDLYQTLYVDADEEEIIQYVVGTLQPKLKRFLRHPLPKTLEQLIQRGMEVQDDLEQAAEPAGPHLPVEDEAETLTPAPNSESVASDRTQPE"

def Classify(ListRecords):
    #assigns group, write with sequence to a file, (in progress) remove duplicate sequences or unknown XXXX

    counter = 0
    counterRecs = 0
    duplicates = 0
    old_sequence_name = "empty"
    old_sequence = "no sequence yet"
    new_sequence_name = "empty2"
    sequence_title = "error! check what happened here"

    for seq_record in SeqIO.parse(ListRecords, "gb"):  #for every record in the list

        duplicates = duplicates + 1

        # setting up initial vatiables
        new_sequence_name = str(seq_record.seq) + "\n"
        new_sequence_length = len(seq_record)
        new_sequence = str(seq_record.seq)
        assignment = "UNASSIGNED FIX ME"
        Number_of_X = unknown_aas(new_sequence)
        prot_name = seq_record.description

        if (CheckIfDuplicate(new_sequence_name, old_sequence_name, new_sequence, old_sequence) == 1) and (Number_of_X == 0) and RemoveLike(prot_name) == 0:  # if not the same and no unknown aas (X) and no "like" in protein name, continue

            duplicates = duplicates - 1

#Classification block begins~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            if seq_record.annotations["taxonomy"][2] == "Ecdysozoa":  #classify as invertebrate
                assignment = "(I)"

            elif seq_record.annotations["taxonomy"][6] == "Amphibia":  # classify as amphibia
                assignment = "(A)"

            elif seq_record.annotations["taxonomy"][6] == "Actinopterygii":  # classify as fish
                assignment = "(F)"

            elif seq_record.annotations["taxonomy"][6] == "Archelosauria":  # classify as reptile or bird
                if seq_record.annotations["taxonomy"][11] == "Coelurosauria" or seq_record.annotations["taxonomy"][11] == "Aves": #bird
                    assignment = "(B)"
                else:
                    assignment = "(R)"

            elif seq_record.annotations["taxonomy"][6] == "Archosauria":  # classify as bird
                if seq_record.annotations["taxonomy"][11] == "Aves":  # bird
                    assignment = "(B)"
                else:
                    counter = counter + 1

            elif seq_record.annotations["taxonomy"][6] == "Lepidosauria" or  seq_record.annotations["taxonomy"][6] == "Testudines + Archosauria group":
                assignment = "(R)"

            elif seq_record.annotations["taxonomy"][6] == "Mammalia":
                if seq_record.annotations["taxonomy"][9] == "Primates":
                    assignment = "(P)"
                else:
                    assignment = "(M)"
            else:
                assignment = "UNCLASSIFIED FIX ME\n"
                counter = counter + 1
#end of classification block~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            counterRecs = counterRecs + 1  #counter of records that made it to file

# Begin pairwise alignment with human
            aligner.open_gap_score = -10
            aligner.extend_gap_score = -0.2
            aligner.substitution_matrix = blosum62
            alignments = aligner.align(human_sequence, new_sequence)

            sequence_title = (">" + str(counterRecs) + ". " + seq_record.annotations["source"] + ", " + assignment + "\n")
            file.write(sequence_title)
            file.write(new_sequence + "\n")


            old_sequence_length = new_sequence_length
            old_sequence_name = new_sequence_name
            old_sequence = new_sequence
            prot_name = seq_record.description
            print(prot_name)

    print("Number of unclassified species:", counter)
    print("Number of removed duplicates or skipped sequences with unknown aas:", duplicates)
    print("Number of records written to file: ", counterRecs)
    print(alignments[0], alignments[1])

    file.close()

Classify("arc_sequences_04202020.gp")


#for seq_record in SeqIO.parse("arc_sequences_04202020.gp","gb"): #uses GenPept file
#print(seq_record.description) #protein name [organism]
#print(seq_record.seq) # sequence
#print(seq_record.annotations["source"]) #name (common name)
#print(seq_record.annotations["taxonomy"][0])
#print(seq_record.annotations)
#print(len(seq_record)) #length of sequence
