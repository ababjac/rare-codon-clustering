import random
import pandas as pd

codon_table = {'TCA': 'S', 'AAT': 'N', 'TGG': 'W', 'GAT': 'D', 'GAA': 'E', 'TTC': 'F', 'CCG': 'P',
           'ACT': 'T', 'GGG': 'G', 'ACG': 'T', 'AGA': 'R', 'TTG': 'L', 'GTC': 'V', 'GCA': 'A',
           'TGA': '*', 'CGT': 'R', 'CAC': 'H', 'CTC': 'L', 'CGA': 'R', 'GCT': 'A', 'ATC': 'I',
           'ATA': 'I', 'TTT': 'F', 'TAA': '*', 'GTG': 'V', 'GCC': 'A', 'GAG': 'E', 'CAT': 'H',
           'AAG': 'K', 'AAA': 'K', 'GCG': 'A', 'TCC': 'S', 'GGC': 'G', 'TCT': 'S', 'CCT': 'P',
           'GTA': 'V', 'AGG': 'R', 'CCA': 'P', 'TAT': 'Y', 'ACC': 'T', 'TCG': 'S', 'ATG': 'M',
           'TTA': 'L', 'TGC': 'C', 'GTT': 'V', 'CTT': 'L', 'CAG': 'Q', 'CCC': 'P', 'ATT': 'I',
           'ACA': 'T', 'AAC': 'N', 'GGT': 'G', 'AGC': 'S', 'CGG': 'R', 'TAG': '*', 'CGC': 'R',
           'AGT': 'S', 'CTA': 'L', 'CAA': 'Q', 'CTG': 'L', 'GGA': 'G', 'TGT': 'C', 'TAC': 'Y',
           'GAC': 'D'}

#amino acid table that maps amino acids to a list of codons that creates them
amino_acid_table = {'S': ['TCA', 'TCC', 'TCT', 'TCG', 'AGC', 'AGT'], 'N': ['AAT', 'AAC'], 'W': ['TGG'],
          'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'], 'F': ['TTC', 'TTT'], 'P': ['CCG', 'CCT', 'CCA', 'CCC'],
          'T': ['ACT', 'ACG', 'ACC', 'ACA'], 'G': ['GGG', 'GGC', 'GGT', 'GGA'],
          'R': ['AGA', 'CGT', 'CGA', 'AGG', 'CGG', 'CGC'], 'L': ['TTG', 'CTC', 'TTA', 'CTT', 'CTA', 'CTG'],
          'V': ['GTC', 'GTG', 'GTA', 'GTT'], 'A': ['GCA', 'GCT', 'GCC', 'GCG'], '*': ['TGA', 'TAA', 'TAG'],
          'H': ['CAC', 'CAT'], 'I': ['ATC', 'ATA', 'ATT'], 'K': ['AAG', 'AAA'], 'Y': ['TAT', 'TAC'],
          'M': ['ATG'], 'C': ['TGC', 'TGT'], 'Q': ['CAG', 'CAA']}



df = pd.read_csv('../Data/Ecoli/full.csv')
# rand_seqs = []
#
# for seq in df['seq']:
#     codon_list = [seq[s:s+3] for s in range(0, len(seq), 3)]
#     new_seq = ''
#
#     for codon in codon_list:
#         if len(codon) != 3:
#             new_codon = codon
#         else:
#             new_codon = random.choice(amino_acid_table[codon_table[codon]])
#         new_seq = new_seq + new_codon
#
#     rand_seqs.append(new_seq)

shuf_seqs = []

for seq in df['mRNA']:
    codon_list = [seq[s:s+3] for s in range(0, len(seq), 3)]

    bool = False
    if len(codon_list[-1]) != 3:
        end = codon_list.pop(-1)
        bool = True

    shuf_list = random.shuffle(codon_list)

    new_seq = ''.join(codon_list)
    if bool:
        new_seq = new_seq + end
        bool = False

    shuf_seqs.append(new_seq)

df['shuf_seqs'] = shuf_seqs
df.to_csv('../Data/Ecoli/shuf_seqs_ecoli.csv')
