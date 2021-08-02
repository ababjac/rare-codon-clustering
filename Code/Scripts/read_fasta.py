from Bio import SeqIO
import pandas as pd

records = []
with open("../Data/sequence.txt/sequence.txt", "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        nc = record.description.split()[0]
        gene = record.description.split()[1].split('=')[1].replace(']', '')
        seq = str(record.seq)
        l = [nc, gene, seq]
        records.append(l)

df = pd.DataFrame(records, columns = ['NC_ID', 'Gene', 'mRNA'])
df.to_csv('../Data/full2.csv')
print(df)
