# Paulina Panek
# April 2020
# Script parsing result xml file to get name of the sequence, classification, sequence, length

import Bio
from Bio import Entrez
Entrez.email = "ppanek@hpu.edu"
from Bio import SeqIO
from io import StringIO
from Bio import GenBank

def numberRecords(ListRecords):
    #Function prints number of records
    records = list(SeqIO.parse(ListRecords, "genbank"))
    print("Found %i records in initial file " % len(records))

def CheckIfDuplicate(first_name, second_name):
    ## returns 0 (same sequences), 1 (not same sequences, or 3 (something went wrong, function didn't work
    return_value = 3

    # if same species AND length of sequence is the same, check if the sequence is the same
    if (first_name == second_name):
        return_value = 0  #same sequences
    else:
        return_value = 1

    return(return_value)

def unknown_aas(sequence):
    #returns number of unknown amino acids in sequence
    X_in_sequence = 0

    if 'X' in sequence:
        X_in_sequence = X_in_sequence + 1
    return X_in_sequence

numberRecords("arc_sequences_04202020.gp")

file = open("UniqueNames.txt", "w")

def UniqueNames(ListRecords):
    #assigns group, write with sequence to a file, (in progress) remove duplicate sequences or unknown XXXX

    counter = 0
    counterRecs = 0
    duplicates = 0
    old_sequence_name = "empty"
    old_sequence_length = 0
    old_sequence = "no sequence yet"
    new_sequence_name = "empty2"
    new_sequence_length = 0
    sequence_title = "error! check what happened here"

    for seq_record in SeqIO.parse(ListRecords, "gb"):  #for every record in the list

        duplicates = duplicates + 1

        # setting up initial vatiables
        new_sequence_name = str(seq_record.annotations["organism"]) + "\n"
        new_sequence = str(seq_record.seq)

        Number_of_X = unknown_aas(new_sequence)

        if (CheckIfDuplicate(old_sequence_name, new_sequence_name) == 1) and (Number_of_X == 0):  # if not the same and no unknown aas (X), continue

            duplicates = duplicates - 1

            counterRecs = counterRecs + 1  #counter of records that made it to file


            sequence_title = (seq_record.annotations["organism"] + "\n")
            file.write(new_sequence_name)

            old_sequence_length = new_sequence_length
            old_sequence_name = new_sequence_name
            old_sequence = new_sequence

    print("Number of unclassified species:", counter)
    print("Number of removed duplicates or skipped sequences with unknown aas:", duplicates)
    print("Number of records written to file: ", counterRecs)
    file.close()
    print(seq_record.annotations["organism"])

UniqueNames("arc_sequences_04202020.gp")


#for seq_record in SeqIO.parse("arc_sequences_04202020.gp","gb"): #uses GenPept file
#print(seq_record.description) #protein name [organism]
#print(seq_record.seq) # sequence
#print(seq_record.annotations["source"]) #name (common name)
#print(seq_record.annotations["taxonomy"][0])
#print(seq_record.annotations)
#print(len(seq_record)) #length of sequence
