from Bio import SeqIO
genome_record = SeqIO.read("ecoli.gb", "genbank")

fhand = open('data2.txt','w')

for feature in genome_record.features:
    if feature.type == "CDS":
        gene_sequence = feature.extract(genome_record.seq)
        if 'join' in str(feature.location):
            fhand.write(str(feature.location)+'\n')
            fhand.write(str(gene_sequence)+ '\n\n')
