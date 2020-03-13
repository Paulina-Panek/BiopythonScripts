# Paulina Panek
# March 2020
# Script that searches NCBI protein database (db) by Accession Number (id)
# and prints full record in xml format

import Bio
from Bio import Entrez
from Bio import SeqIO
import os

Entrez.email = "ppanek@hpu.edu"


handle = Entrez.efetch(db="protein", id="GER55599", rettype = "gp", retmode = "text")
record = SeqIO.read(handle, "gb")
handle.close()
print("*************************")
print(record.id)
print(record.name)
