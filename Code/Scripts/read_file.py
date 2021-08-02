from Bio import SeqIO
import pandas as pd

records = []
with open("../Data/full.txt", "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        gene = record.description.split()[1].split('=')[1].replace(']', '')
        uni = record.description.split()[3].split(':')[1].replace(']', '')
        l = [uni, gene]
        records.append(l)

df = pd.DataFrame(records, columns = ['Uniprot_ID', 'Gene'])
df.to_csv('../Data/uniprot_mapping.csv')
print(df)
