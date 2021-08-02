import pandas as pd
import numpy as np
import statistics


WINDOW_SIZE = 9

#codon table that maps codon to its amino acid
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

codon_freq_ecoli = {'TTT': 22.38, 'TCT': 8.61, 'TAT': 16.36, 'TGT': 5.19, 'TTC': 16.21,
                      'TCC': 8.81, 'TAC': 12.15, 'TGC': 6.34, 'TTA': 13.83, 'TCA': 7.57,
                      'TAA': 2.03, 'TGA': 1.04, 'TTG': 13.37, 'TCG': 8.79, 'TAG': 0.3,
                      'TGG': 15.21, 'CTT': 11.44, 'CCT': 7.22, 'CAT': 12.84, 'CGT': 20.7,
                      'CTC': 10.92, 'CCC': 5.56, 'CAC': 9.44, 'CGC': 21.48, 'CTA': 3.93,
                      'CCA': 8.44, 'CAA': 15.1, 'CGA': 3.67, 'CTG': 52.1, 'CCG': 22.65,
                      'CAG': 29.21, 'CGG': 5.72, 'ATT': 30.21, 'ACT': 9.02, 'AAT': 18.26,
                      'AGT': 9.08, 'ATC': 24.6, 'ACC': 22.88, 'AAC': 21.47, 'AGC': 15.89,
                      'ATA': 4.88, 'ACA': 7.63, 'AAA': 33.94, 'AGA': 2.43, 'ATG': 27.59,
                      'ACG': 14.47, 'AAG': 10.7, 'AGG': 1.48, 'GTT': 18.39, 'GCT': 15.54,
                      'GAT': 32.43, 'GGT': 24.45, 'GTC': 15.07, 'GCC': 25.45, 'GAC': 19.14,
                      'GGC': 28.65, 'GTA': 10.97, 'GCA': 20.61, 'GAA': 39.55, 'GGA': 8.44,
                      'GTG': 25.9, 'GCG': 32.79, 'GAG': 18.24, 'GGG': 11.29}

codon_freq_scerevisiae = {  'TTT': 26.18, 'TCT': 23.35, 'TAT': 19.05, 'TGT': 7.82, 'TTC': 17.88,
                            'TCC': 14.07, 'TAC': 14.6, 'TGC': 4.75, 'TTA': 26.33, 'TCA': 19.05,
                            'TAA': 0.95, 'TGA': 0.6, 'TTG': 26.5, 'TCG': 8.71, 'TAG': 0.46,
                            'TGG': 10.35, 'CTT': 12.27, 'CCT': 13.57, 'CAT': 13.89, 'CGT': 6.26,
                            'CTC': 5.52, 'CCC': 6.91, 'CAC': 7.74, 'CGC': 2.63, 'CTA': 13.52,
                            'CCA': 17.81, 'CAA': 27.1, 'CGA': 3.1, 'CTG': 10.65, 'CCG': 5.42,
                            'CAG': 12.42, 'CGG': 1.82, 'ATT': 30.1, 'ACT': 20.24, 'AAT': 36.61,
                            'AGT': 14.6, 'ATC': 16.99, 'ACC': 12.48, 'AAC': 24.8, 'AGC': 9.96,
                            'ATA': 18.29, 'ACA': 18.18, 'AAA': 42.83, 'AGA': 21.05, 'ATG': 20.68,
                            'ACG': 8.15, 'AAG': 30.52, 'AGG': 9.45, 'GTT': 21.47, 'GCT': 20.28,
                            'GAT': 38.09, 'GGT': 22.59, 'GTC': 11.23, 'GCC': 12.14, 'GAC': 20.39,
                            'GGC': 9.78, 'GTA': 12.07, 'GCA': 16.26, 'GAA': 45.81, 'GGA': 11.19,
                            'GTG': 10.72, 'GCG': 6.17, 'GAG': 19.55, 'GGG': 6.06}


#---------------------------------------------------------------------------------------------------------------#

def get_aa_frequencies_list(codon_freqs, amino_acid_table, codon_table):
    """
    Parameters: codon_freqs(dictionary), amino_acid_table(dictionary), codon_table(dictionary)
    Description: Creates a list for each amino acid of all associated codon frequencies
    Return: dictionary of lists ({key=amino_acid, val=list(codon_freq))
    Todo: NONE
    """
    #n = len(codon_freqs)
    aa_freqs = {}
    aa_map = {}

    # for i in range(n):
    #     aa_freqs[i] = {}
    for aa, l in amino_acid_table.items():
        aa_freqs[aa] = []
        aa_map[aa] = []
        for codon in l:
            if codon in codon_freqs.keys():
                aa_freqs[aa].append(codon_freqs[codon])
                aa_map[aa].append(codon + ' ' + str(codon_freqs[codon]))

    return aa_freqs, aa_map

#---------------------------------------------------------------------------------------------------------------#

def calc_minMax_percent(seq, aa_avg_freq, codon_freq, codon_table, window_size):
    """
    Parameters: seq(String), aa_avg_freq(dictionary), codon_freq(dictionary), codon_table(dictionary), window_size(int)
    Description: Generates a list of %minMax values based on window_size for a given sequence
    Return: list
    Todo: NONE
    """
    codon_seq = [seq[s:s+3] for s in range(0, len(seq), 3)]
    n = len(codon_seq)
    skip = int(window_size/2)
    minMax_by_row = [0.0] * n

    actual = 0.0
    maximum = 0.0
    minimum = 0.0
    average = 0.0
    percent_max = 0.0
    percent_min = 0.0

    for i in range(skip-1, skip+n-window_size+1):
        codons_in_window = codon_seq[i-skip:i+window_size-skip]
        actual = maximum = minimum = average = percent_max = percent_min = 0.0

        for codon in codons_in_window:
            if len(codon) != 3:
                continue
            freq_for_window = aa_avg_freq[codon_table[codon]]

            if not freq_for_window:
                freq_for_window.append(0)

            m = len(freq_for_window)

            actual += codon_freq[codon]
            maximum += max(freq_for_window)
            minimum += min(freq_for_window)
            average += sum(freq_for_window) / m

        actual /= window_size
        maximum /= window_size
        minimum /= window_size
        average /= window_size

        if (maximum - average) != 0:
            percent_max = ((actual-average)/(maximum-average))*100
        if (average - minimum) != 0:
            percent_min = ((average-actual)/(average-minimum))*100

        if percent_max >= 0:
            minMax_by_row[i] = percent_max
        else:
            minMax_by_row[i] = -1 * percent_min

    return minMax_by_row

#---------------------------------------------------------------------------------------------------------------#

def add_col_to_df(MM_list, df):
    MM_min_list = []
    MM_max_list = []
    MM_avg_list = []
    MM_median_list = []

    for l in MM_list:
        if not l:
            MM_min_list.append(None)
            MM_avg_list.append(None)
            MM_max_list.append(None)
            MM_median_list.append(None)
        else:
            MM_min_list.append(min(l))
            MM_avg_list.append(sum(l) / len(l))
            MM_max_list.append(max(l))
            MM_median_list.append(statistics.median(l))

    #print(MM_min_list)
    df['MM_min'] = MM_min_list
    df['MM_avg'] = MM_avg_list
    df['MM_max'] = MM_max_list
    df['MM_median'] = MM_median_list

    return df

#---------------------------------------------------------------------------------------------------------------#

def reduce_matrix(df_path, vals_path, codon_freq, colname_seq, colname_id, out_path, perc):

    df = pd.read_csv(df_path)
    vals = pd.read_csv(vals_path)
    vals.drop(vals.columns[vals.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)

    vals.columns = df[colname_id]
    vals.index = df[colname_id]

    aa_freq, aa_map = get_aa_frequencies_list(codon_freq, amino_acid_table, codon_table)

    #get minMax for full data set
    minMax_by_row_full = []
    for seq in df[colname_seq]:
        minMax_by_row_full.append(calc_minMax_percent(seq, aa_freq, codon_freq, codon_table, WINDOW_SIZE))

    df = add_col_to_df(minMax_by_row_full, df)
    df_middle = df.sort_values(by='MM_avg', ascending=False)[round(perc*len(df)):round((1-perc)*len(df))]

    genes = df_middle[colname_id].values.tolist()
    remove_list = [r for r in df[colname_id] if r not in genes]

    vals.drop(index=remove_list, axis=0, inplace=True)
    vals.drop(labels=remove_list, axis=1, inplace=True)

    vals.to_csv(out_path)

#script
# reduce_matrix('../Data/Ecoli/full.csv', '../Files/MM/Full/Plain//HD_ecoli.csv', codon_freq_ecoli, 'mRNA', 'Gene', '../Files/MM/Omit10th/HD_ecoli.csv', 0.1)
# reduce_matrix('../Data/Ecoli/rand_seqs_ecoli.csv', '../Files/MM/Full/Plain//HD_ecoli_rand.csv', codon_freq_ecoli, 'rand_seqs', 'Gene', '../Files/MM/Omit10th/HD_ecoli_rand.csv', 0.1)
# reduce_matrix('../Data/Ecoli/shuf_seqs_ecoli.csv', '../Files/MM/Full/Plain//HD_ecoli_shuf.csv', codon_freq_ecoli, 'shuf_seqs', 'Gene', '../Files/MM/Omit10th/HD_ecoli_shuf.csv', 0.1)
#
# reduce_matrix('../Data/Ecoli/full.csv', '../Files/MM/Full/Plain//KS_stats_ecoli.csv', codon_freq_ecoli, 'mRNA', 'Gene', '../Files/MM/Omit10th/KS_stats_ecoli.csv', 0.1)
# reduce_matrix('../Data/Ecoli/rand_seqs_ecoli.csv', '../Files/MM/Full/Plain//KS_stats_ecoli_rand.csv', codon_freq_ecoli, 'rand_seqs', 'Gene', '../Files/MM/Omit10th/KS_stats_ecoli_rand.csv', 0.1)
# reduce_matrix('../Data/Ecoli/shuf_seqs_ecoli.csv', '../Files/MM/Full/Plain//KS_stats_ecoli_shuf.csv', codon_freq_ecoli, 'shuf_seqs', 'Gene', '../Files/MM/Omit10th/KS_stats_ecoli_shuf.csv', 0.1)

# reduce_matrix('../Data/Yeast/YeastSeqs.csv', '../Files/MM/Full/Plain//HD_yeast.csv', codon_freq_scerevisiae, 'seq', 'locus_tag', '../Files/MM/Omit10th/HD_yeast.csv', 0.1)
# reduce_matrix('../Data/Yeast/rand_seqs_yeast.csv', '../Files/MM/Full/Plain//HD_yeast_rand.csv', codon_freq_scerevisiae, 'rand_seqs', 'locus_tag', '../Files/MM/Omit10th/HD_yeast_rand.csv', 0.1)
# reduce_matrix('../Data/Yeast/shuf_seqs_yeast.csv', '../Files/MM/Full/Plain//HD_yeast_shuf.csv', codon_freq_scerevisiae, 'shuf_seqs', 'locus_tag', '../Files/MM/Omit10th/HD_yeast_shuf.csv', 0.1)
#
# reduce_matrix('../Data/Yeast/YeastSeqs.csv', '../Files/MM/Full/Plain//KS_stats_yeast.csv', codon_freq_scerevisiae, 'seq', 'locus_tag', '../Files/MM/Omit10th/KS_stats_yeast.csv', 0.1)
# reduce_matrix('../Data/Yeast/rand_seqs_yeast.csv', '../Files/MM/Full/Plain//KS_stats_yeast_rand.csv', codon_freq_scerevisiae, 'rand_seqs', 'locus_tag', '../Files/MM/Omit10th/KS_stats_yeast_rand.csv', 0.1)
# reduce_matrix('../Data/Yeast/shuf_seqs_yeast.csv', '../Files/MM/Full/Plain//KS_stats_yeast_shuf.csv', codon_freq_scerevisiae, 'shuf_seqs', 'locus_tag', '../Files/MM/Omit10th/KS_stats_yeast_shuf.csv', 0.1)


#script
reduce_matrix('../Data/Ecoli/full.csv', '../Files/MM/Full/Plain//HD_ecoli.csv', codon_freq_ecoli, 'mRNA', 'Gene', '../Files/MM/Omit30/HD_ecoli.csv', 0.3)
reduce_matrix('../Data/Ecoli/rand_seqs_ecoli.csv', '../Files/MM/Full/Plain//HD_ecoli_rand.csv', codon_freq_ecoli, 'rand_seqs', 'Gene', '../Files/MM/Omit30/HD_ecoli_rand.csv', 0.3)
reduce_matrix('../Data/Ecoli/shuf_seqs_ecoli.csv', '../Files/MM/Full/Plain//HD_ecoli_shuf.csv', codon_freq_ecoli, 'shuf_seqs', 'Gene', '../Files/MM/Omit30/HD_ecoli_shuf.csv', 0.3)

reduce_matrix('../Data/Ecoli/full.csv', '../Files/MM/Full/Plain//KS_stats_ecoli.csv', codon_freq_ecoli, 'mRNA', 'Gene', '../Files/MM/Omit30/KS_stats_ecoli.csv', 0.3)
reduce_matrix('../Data/Ecoli/rand_seqs_ecoli.csv', '../Files/MM/Full/Plain//KS_stats_ecoli_rand.csv', codon_freq_ecoli, 'rand_seqs', 'Gene', '../Files/MM/Omit30/KS_stats_ecoli_rand.csv', 0.3)
reduce_matrix('../Data/Ecoli/shuf_seqs_ecoli.csv', '../Files/MM/Full/Plain//KS_stats_ecoli_shuf.csv', codon_freq_ecoli, 'shuf_seqs', 'Gene', '../Files/MM/Omit30/KS_stats_ecoli_shuf.csv', 0.3)

reduce_matrix('../Data/Yeast/YeastSeqs.csv', '../Files/MM/Full/Plain//HD_yeast.csv', codon_freq_scerevisiae, 'seq', 'locus_tag', '../Files/MM/Omit30/HD_yeast.csv', 0.3)
reduce_matrix('../Data/Yeast/rand_seqs_yeast.csv', '../Files/MM/Full/Plain//HD_yeast_rand.csv', codon_freq_scerevisiae, 'rand_seqs', 'locus_tag', '../Files/MM/Omit30/HD_yeast_rand.csv', 0.3)
reduce_matrix('../Data/Yeast/shuf_seqs_yeast.csv', '../Files/MM/Full/Plain//HD_yeast_shuf.csv', codon_freq_scerevisiae, 'shuf_seqs', 'locus_tag', '../Files/MM/Omit30/HD_yeast_shuf.csv', 0.3)

reduce_matrix('../Data/Yeast/YeastSeqs.csv', '../Files/MM/Full/Plain//KS_stats_yeast.csv', codon_freq_scerevisiae, 'seq', 'locus_tag', '../Files/MM/Omit30/KS_stats_yeast.csv', 0.3)
reduce_matrix('../Data/Yeast/rand_seqs_yeast.csv', '../Files/MM/Full/Plain//KS_stats_yeast_rand.csv', codon_freq_scerevisiae, 'rand_seqs', 'locus_tag', '../Files/MM/Omit30/KS_stats_yeast_rand.csv', 0.3)
reduce_matrix('../Data/Yeast/shuf_seqs_yeast.csv', '../Files/MM/Full/Plain//KS_stats_yeast_shuf.csv', codon_freq_scerevisiae, 'shuf_seqs', 'locus_tag', '../Files/MM/Omit30/KS_stats_yeast_shuf.csv', 0.3)
