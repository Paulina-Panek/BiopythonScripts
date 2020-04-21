# Paulina Panek
# April 2020
# Script parsing result xml file to get name of the sequence, classification, sequence, length

import Bio
from Bio import Entrez
Entrez.email = "ppanek@hpu.edu"
from Bio import SeqIO
from io import StringIO
from Bio import GenBank
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.CheckSum import seguid


def numberRecords(ListRecords):
    #Function prints number of records
    records = list(SeqIO.parse(ListRecords, "genbank"))
    print("Found %i records in initial file: " % len(records))

numberRecords("arc_sequences_04202020.gp")

file = open("AllSpecies.fasta", "w")

def Classify(ListRecords):
    #assigns group, prints with sequence

    counter = 0
    counterRecs = 0

    for seq_record in SeqIO.parse(ListRecords, "gb"):

        sequence = str(seq_record.seq)
        assignment = "UNASSIGNED FIX ME"

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
            print(seq_record.annotations["source"] + ", UNCLASSIFIED FIX ME\n")
            counter = counter + 1

        file.write(">" + seq_record.annotations["source"] + ", " + assignment + "\n")
        file.write(sequence + "\n")
        counterRecs = counterRecs + 1

    print("Number of unclassified species:", counter)
    print("Number of records written to file: ", counterRecs)
    file.close()

Classify("arc_sequences_04202020.gp")


#for seq_record in SeqIO.parse("arc_sequences_04202020.gp","gb"): #uses GenPept file
    #print(seq_record.description) #protein name [organism]
    #print(seq_record.seq) # sequence
    #print(seq_record.annotations["source"]) #name (common name)
    #print(seq_record.annotations["taxonomy"][0])
    #print(seq_record.annotations)
    #print(len(seq_record)) #length of sequence


