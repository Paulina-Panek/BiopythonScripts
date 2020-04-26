# Writes file with records that have multiple entries for a single organism

from Bio import SeqIO

def WriteRepeats(ListRecords):
    counter_old = 0
    old_species = "empty"
    new_species = "empty2"
    file = open("RepeatedSequences.txt", "w")
    old_record = "empty3"

    for seq_record in SeqIO.parse(ListRecords, "gb"):
        new_organism = str(seq_record.annotations["organism"])
        new_record = seq_record

        if (new_species == old_species):
            file.write(str(counter_old) + "\n\n" + str(old_record) +"\n\n" + str(new_record))
        else:
            old_record = new_record
            old_species = new_species
    file.close()

WriteRepeats("arc_sequences_04202020.gp")