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
    #Function rints number of records
    records = list(SeqIO.parse(ListRecords, "genbank"))
    print("Found %i records" % len(records))

#numberRecords("arc_sequences_04202020.gp")


def Classify(ListRecords):
    #assigns group, prints with sequence
    counter = 0
    for seq_record in SeqIO.parse(ListRecords, "gb"):
        if seq_record.annotations["taxonomy"][2] == "Ecdysozoa":  #classify as invertebrate
            print(seq_record.annotations["source"] + ", (I)")

        elif seq_record.annotations["taxonomy"][6] == "Amphibia":  # classify as amphibia
            print(seq_record.annotations["source"] + ", (A)")

        elif seq_record.annotations["taxonomy"][6] == "Actinopterygii":  # classify as fish
            print(seq_record.annotations["source"] + ", (F)")

        elif seq_record.annotations["taxonomy"][6] == "Archelosauria":  # classify as reptile or bird
            if seq_record.annotations["taxonomy"][11] == "Coelurosauria" or seq_record.annotations["taxonomy"][11] == "Aves": #bird
                print(seq_record.annotations["source"] + ", (B)")
            else:
                print(seq_record.annotations["source"] + ", (R)")

        elif seq_record.annotations["taxonomy"][6] == "Archosauria":  # classify as bird
            if seq_record.annotations["taxonomy"][11] == "Aves":  # bird
                print(seq_record.annotations["source"] + ", (B)")
            else:
                print(seq_record.annotations["source"] + ", UNCLASSIFIED FIX ME")
                counter = counter + 1

        elif seq_record.annotations["taxonomy"][6] == "Lepidosauria" or  seq_record.annotations["taxonomy"][6] == "Testudines + Archosauria group":
            print(seq_record.annotations["source"] + ", (R)")

        elif seq_record.annotations["taxonomy"][6] == "Mammalia":
            if seq_record.annotations["taxonomy"][9] == "Primates":
                print(seq_record.annotations["source"] + ", (P)")
            else:
                print(seq_record.annotations["source"] + ", (M)")
        else:
            print(seq_record.annotations["source"] + ", UNCLASSIFIED FIX ME")
            counter = counter + 1
        print(seq_record.seq)


    print("Number of unclassified species:", counter)

Classify("arc_sequences_04202020.gp")


#for seq_record in SeqIO.parse("arc_sequences_04202020.gp","gb"): #uses GenPept file
    #print(seq_record.description) #protein name [organism]
    #print(seq_record.seq) # sequence
    #print(seq_record.annotations["source"]) #name (common name)
    #print(seq_record.annotations["taxonomy"][0])
    #print(seq_record.annotations)
    #print(len(seq_record)) #length of sequence


