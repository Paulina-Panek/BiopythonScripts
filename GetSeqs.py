# Paulina Panek
# April 2020
# Script parsing result xml file to get name of the sequence, classification, sequence, length

import Bio
from Bio import Entrez
Entrez.email = "ppanek@hpu.edu"
from Bio import SeqIO


for seq_record in SeqIO.parse("arc_sequences_04202020.gp","gb"): #uses GenPept file
    #print(seq_record.description) #protein name [organism]
    #print(seq_record.seq) # sequence
    #print(seq_record.annotations["source"]) #name (common name)
    print(seq_record.annotations["taxonomy"]) #name (common name)
    #print(seq_record.annotations)