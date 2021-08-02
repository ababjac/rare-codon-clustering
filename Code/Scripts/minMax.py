import pandas as pd
import numpy as np
import statistics
import math
from sklearn import linear_model
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy import stats
import random

#---------------------------------------------------------------------------------------------------------------#
#constant window size -- change to be passed in as parameter or possibly try different values during regressions
#window sizes of 9 are very good based on recent paper but 17 was used historcally
WINDOW_SIZE = 9
BINS = np.arange(-100, 100, 5)

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
                      'TAA': 2.03, 'TGA': 1.04, 'TTG': 13.37, 'TCG': 8.79, 'TAG': 0.25,
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

def make_linear_regression(X, Y):
    #Y = Y.values.reshape(-1, 1)
    #X = X.values.reshape(-1, 1)

    n = int(-.9*len(X))
    X_train = X.head(n)
    Y_train = Y[:n]
    #X_test = X[n:]
    #Y_test = Y[n:]

    #plt.scatter(X_test, Y_test)
    #plt.title(plt_title)
    #plt.xlabel('%MinMax Min')
    #plt.ylabel('log10(PHI)')
    #plt.show()

    regr = linear_model.LinearRegression()
    regr.fit(X_train, Y_train)
    print(regr.score(X_train, Y_train))
    #print(regr.score(X_test, Y_test))
    #plt.plot(X_test, regr.predict(X_test))

    return regr

#---------------------------------------------------------------------------------------------------------------#

def make_scatterplot_data(X, Y, xlabel, ylabel, plt_title):
    plt.scatter(X, Y)
    plt.title(plt_title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()

def make_histogram_data(X, xlabel, plt_title, bins):
    plt.hist(X, bins)
    plt.title(plt_title)
    plt.xlabel(xlabel)
    plt.ylabel('Count')
    plt.show()

#---------------------------------------------------------------------------------------------------------------#

def run_regressions(df):
    X1 = sm.add_constant(df[['MM_avg']])
    #X1 = sm.add_constant(df[['MM_avg']])
    Y = df['log10.PHI']
    #print(Y)
    var = sm.OLS(Y, X1)
    model = var.fit()

    print(model.summary())

    X2 = sm.add_constant(df['MM_min'])
    #Y = np.log(df['abundance (molecules/cell)'])
    var = sm.OLS(Y, X2)
    model = var.fit()

    print(model.summary())
    # #Y = df['abundance (molecules/cell)']
    # X1 = df['MM_min']
    # X2 = df['MM_avg']
    # #X3 = df['MM_max']
    # X4 = df[['MM_min', 'MM_avg']]
    # #print(X4)
    #
    # #print(X1, X2, X3)
    #
    # plot_data(X1, Y, 'MM_min')
    # plot_data(X2, Y, 'MM_avg')
    # #plot_data(X3, Y, 'MM_max')
    #
    # #regr_X1 = make_linear_regression(X1, Y)
    # #regr_X2 = make_linear_regression(X2, Y)
    # #regr_X3 = make_linear_regression(X3, Y)
    # regr_X4 = make_linear_regression(X4, Y)

#---------------------------------------------------------------------------------------------------------------#

def get_negative_window_counts(nl, percent=False):
    counts = []

    for l in nl:
        count = 0
        for elem in l:
            if elem < 0:
                count += 1

        if(percent):
            counts.append(float(count/len(l)))
        else:
            counts.append(count)

    return counts

#---------------------------------------------------------------------------------------------------------------#

def geo_mean_overflow(iterable):
    a = np.log(iterable)
    return np.exp(a.mean())

#---------------------------------------------------------------------------------------------------------------#

def generate_minmax_plot(minmax_list, title, avgRRT=None, mean=False):
    indexes = list(range(1, len(minmax_list)+1))
    all_pos_mm = [elem+101 for elem in minmax_list]
    geo_mean = geo_mean_overflow(all_pos_mm) - 101

    minmax_list = [None if i==0.0 else i for i in minmax_list]

    plt.plot(indexes, minmax_list, label = "%MinMax", color = '#0000CC')

    if avgRRT != None:
        plt.plot(indexes, avgRRT, label = "Average RRT", color = '#808080')

    if mean:
        plt.plot([0,len(minmax_list)],[geo_mean,geo_mean],'r--')
        print(geo_mean)
    else:
        plt.plot([0,len(minmax_list)],[0,0],'r--')

    plt.legend(loc = 4)
    plt.xlim(0,len(minmax_list))
    plt.ylim(-100,100)
    #plt.fill_between(num, stdevbel0,stdevabv0, color = '#808080', alpha = .1)
    plt.title(title)
    plt.xlabel("Center Of Codon Window")
    plt.ylabel("%MinMax")

    plt.show()

#---------------------------------------------------------------------------------------------------------------#

def convert_codon_to_AA_seq(codon_seq):
    return [codon_table[codon] for codon in codon_seq]

#---------------------------------------------------------------------------------------------------------------#
#courtesy of Gabe
def calc_avgRRT(seq, aaMapDict, aaFreqDict, freqDict, mapDict, numRRT=1000):
    codon_seq = [seq[s:s+3] for s in range(0, len(seq), 3)]
    aaSeq = convert_codon_to_AA_seq(codon_seq)

    rrtCodonList = []
    for i in range(numRRT):
        rrtCodonList.append([])
        for j in range(len(aaSeq)):
            aa = aaSeq[j]
            data = aaMapDict[aa]
            codons = []
            freqs = []
            for duo in data:
                cAndF = duo.split()
                codons.append(cAndF[0])
                freqs.append(float(cAndF[1]))
            freqSum = sum(freqs)
            randomNum = random.uniform(0,freqSum)
            total = 0
            flag = 0
            for k in range(len(freqs)):
                if flag == 0:
                    total += freqs[k]
                    if total > randomNum:
                        rrtCodonList[i].append(codons[k])
                        flag = 1
    rrtList = []
    for i in range(numRRT):
        codonSequence = ''.join(rrtCodonList[i])
        mmValues = calc_minMax_percent(codonSequence, aaFreqDict, freqDict, mapDict, WINDOW_SIZE)
        rrtList.append(mmValues)

    #Average the RRTs
    avgRRT = []
    stdRRT = []
    for i in range(int(WINDOW_SIZE/2)):
        avgRRT.append(None)
        stdRRT.append(None)
    for i in range(int(WINDOW_SIZE/2),len(rrtList[0])-int(WINDOW_SIZE/2)):
        avg = 0
        stdlist = []
        for j in range(len(rrtList)):
            avg += rrtList[j][i]
            stdlist.append(rrtList[j][i])
        avg = avg/numRRT
        stdRRT.append(statistics.stdev(stdlist))
        avgRRT.append(round(avg,2))
    for i in range(int(WINDOW_SIZE/2)):
        avgRRT.append(None)
        stdRRT.append(None)

    return avgRRT

#---------------------------------------------------------------------------------------------------------------#

def fix_sequence(seq):
    START = 'ATG'
    STOP = ['TAG', 'TAA', 'TGA']

    if len(seq) % 3 != 0:
        pos = seq.find(START)
        seq = seq[pos:]

    codon_seq =  [seq[s:s+3] for s in range(0, len(seq), 3)]

    actual_seq = []
    bool = False
    for codon in codon_seq:
        if codon == START:
            bool = True
        if bool:
            if len(codon) == 3:
                actual_seq.append(codon)
        if codon in STOP:
            bool = False

    return ''.join(actual_seq)

#----------------------------------#---------------------------------------------------------------------------------------------------------------#

def get_bin_counts(list, bins):
    hist = np.histogram(list, bins)

    return hist[0].tolist()

#---------------------------------------------------------------------------------------------------------------#

def hellinger(p, q):
    return np.sqrt(np.sum((np.sqrt(p) - np.sqrt(q)) ** 2)) / np.sqrt(2)

#---------------------------------------------------------------------------------------------------------------#


#MAIN#
#shuf1 = 'ATGCATCACCATCACCATCACCATAACTATACAAAATTTGATGTAAAAAATTGGGTTCGCCGTGAGCATTTTGAGTTTTATCGGCATCGTTTACCATGTGGTTTTAGCTTAACAAGCAAAATTGATATCACGACGTTAAAAAAGTCATTGGATGATTCAGCGTATAAGTTTTACCCCGTGATGATATACTTAATTGCCCAAGCAGTTAACCAATTTGACGAGCTCCGAATGGCGATTAAAGATGATGAACTCATTGTGTGGGACTCCGTCGATCCGCAATTTACTGTATTCCATCAGGAGACTGAAACTTTTAGTGCGCTATCCTGCCCCTATAGCTCGGACATCGATCAATTTATGGTGAATTACCTGTCGGTCATGGAAAGGTATAAGTCCGACACAAAATTATTTCCCCAGGGCGTGACGCCAGAAAACCATCTAAACATATCCGCGCTGCCATGGGTTAACTTTGACTCGTTCAACTTAAATGTAGCTAACTTTACTGACTATTTTGCGCCGATCATAACCATGGCCAAATACCAACAGGAAGGGGATAGGCTCTTGTTACCCTTGAGCGTCCAGGTGCATCACGCCGTGTGTGACGGCTTTCATGTTGCACGCTTCATTAACAGATTACAGGAGTTATGTAATTCTAAATTGAAGTAA'
#
# aa_freq_ecoli, aa_map_ecoli = get_aa_frequencies_list(codon_freq_ecoli, amino_acid_table, codon_table)
# MM = calc_minMax_percent(fix_sequence(shuf1), aa_freq_ecoli, codon_freq_ecoli, codon_table, WINDOW_SIZE)
# generate_minmax_plot(MM, '', mean=True)



df = pd.read_csv('../Data/full.csv')
seq1 = df['mRNA'][df['Gene'] == 'mdtC'].values[0]
seq2 = df['mRNA'][df['Gene'] == 'cusA'].values[0]

aa_freq_ecoli, aa_map_ecoli = get_aa_frequencies_list(codon_freq_ecoli, amino_acid_table, codon_table)


MM_1 = calc_minMax_percent(fix_sequence(seq1), aa_freq_ecoli, codon_freq_ecoli, codon_table, WINDOW_SIZE)
generate_minmax_plot(MM_1, 'mdtC', mean=True)
MM_2 = calc_minMax_percent(fix_sequence(seq2), aa_freq_ecoli, codon_freq_ecoli, codon_table, WINDOW_SIZE)
generate_minmax_plot(MM_2, 'cusA', mean=True)

make_histogram_data(MM_1, '%MM', 'mdtC', BINS)
make_histogram_data(MM_2, '%MM', 'cusA', BINS)

H_1 = np.divide(get_bin_counts(MM_1, BINS), len(MM_1))
H_2 = np.divide(get_bin_counts(MM_2, BINS), len(MM_2))

HD = hellinger(H_1, H_2)
print(HD)


#
# new_seqs = []
# for seq in df['mRNA']:
#     corrected_seq = fix_sequence(seq)
#     new_seqs.append(corrected_seq)
#
# df['mRNA_adjusted'] = new_seqs
# #get codon counts for top 10%
# #df_top_10 = df.sort_values(by='disorder_ratio', ascending=False).head(round(len(df)*.1))
# #df_bottom_10 = df.sort_values(by='disorder_ratio', ascending=True).head(round(len(df)*.1))
#
# aa_freq_ecoli, aa_map_ecoli = get_aa_frequencies_list(codon_freq_ecoli, amino_acid_table, codon_table)
#
# #get minMax for full data set
# minMax_by_row_full = []
# for seq in df['mRNA_adjusted']:
#     minMax_by_row_full.append(calc_minMax_percent(seq, aa_freq_ecoli, codon_freq_ecoli, codon_table, WINDOW_SIZE))
#
# # #get minMax values by row for top 10% of log10(PHI) col from data set
# # minMax_by_row_top_10 = []
# # for seq in df_top_10['mRNA_adjusted']:
# #     minMax_by_row_top_10.append(calc_minMax_percent(seq, aa_freq_ecoli, codon_freq_ecoli, codon_table, WINDOW_SIZE))
# # #
# # minMax_by_row_bottom_10 = []
# # for seq in df_bottom_10['mRNA_adjusted']:
# #     minMax_by_row_bottom_10.append(calc_minMax_percent(seq, aa_freq_ecoli, codon_freq_ecoli, codon_table, WINDOW_SIZE))
# #
# # #run regressions
# df = add_col_to_df(minMax_by_row_full, df)
# df['MinMax'] = minMax_by_row_full


#print('Spearman Correlation for AVG %MM vs Disorder = ', stats.spearmanr(df['MM_avg'].values.tolist(), df['disorder_ratio'].values.tolist(), nan_policy='omit'))
#
# df_top_10 = add_col_to_df(minMax_by_row_top_10, df_top_10)
# df_bottom_10 = add_col_to_df(minMax_by_row_bottom_10, df_bottom_10)
#
#counts_top_10 = get_negative_window_counts(minMax_by_row_top_10)
# counts_bottom_10 = get_negative_window_counts(minMax_by_row_bottom_10)
#
#perc_top_10 = get_negative_window_counts(minMax_by_row_top_10, percent=True)
# perc_bottom_10 = get_negative_window_counts(minMax_by_row_bottom_10, percent=True)
#
#norm_counts_top_10 = [(float(elem)-min(counts_top_10))/(max(counts_top_10)-min(counts_top_10)) for elem in counts_top_10]
# norm_counts_bottom_10 = [(float(elem)-min(counts_bottom_10))/(max(counts_bottom_10)-min(counts_bottom_10)) for elem in counts_bottom_10]
#
# df_top_10['negCounts'] = counts_top_10
# df_top_10['negPerc'] = perc_top_10
# df_top_10['negNorm'] = norm_counts_top_10
#
# df_top_10['MinMax'] = minMax_by_row_top_10
#
# df_bottom_10['negCounts'] = counts_bottom_10
# df_bottom_10['negPerc'] = perc_bottom_10
# df_bottom_10['negNorm'] = norm_counts_bottom_10

# avg = sum(df_top_10['negNorm']) / len(df_top_10['negNorm'])
# std = df_top_10.loc[:,'negNorm'].std()
#
# print(avg, std)

#indexes = [print(row) for row in df_top_10.index if df_top_10['negPerc'][row] >= 0.13]
# print(df_top_10.loc[[83]])
# print(df_top_10.loc[[172]])
# print(df_top_10.loc[[187]])

#indexes2 = [print(row) for row in df_top_10.index if df_top_10['negNorm'][row] >= 0.8] #based on avg and 2 std deviations
#
# #correct protein ids
# # print(df_top_10.loc[[1128]])
# # print(df_top_10.loc[[508]])
# # print(df_top_10.loc[[1396]])
# print(df_top_10.loc[[24]])
# print(df_top_10.loc[[735]])
# print(df_top_10.loc[[783]])
# print(df_top_10.loc[[724]])
# print(df_top_10.loc[[1003]])
# print(df_top_10.loc[[1532]])
# print(df_top_10.loc[[773]])
# print(df_top_10.loc[[490]])
# print(df_top_10.loc[[1726]])

# df_counts = pd.DataFrame()
# df_counts['negNorm_bottom10'] = norm_counts_bottom_10
# df_counts['negNorm_top10'] = norm_counts_top_10
#
# df_counts.plot.hist(bins=10, alpha=0.7)
# plt.title("Comparison of bottom and top 10% log10(PHI) Negative Normalized %MinMax Frequency")
# plt.xlabel("%MinMax Negative Values Normalized Count")
# plt.show()

# make_histogram_data(perc_top_10, '%MinMax Negative Percentage per Window for top 10% log10(PHI)', 'Histogram of Negative %MinMax Percentage per Window for top 10% log10(PHI)', 10)
# make_histogram_data(perc_bottom_10, '%MinMax Negative Count Percentage per Window for bottom 10% log10(PHI)', 'Histogram of Negative %MinMax Percentage per Window for bottom 10% log10(PHI)', 10)
#
# make_histogram_data(counts_top_10, '%MinMax Negative Count per Window for top 10% log10(PHI)', 'Histogram of Negative %MinMax Count per Window for top 10% log10(PHI)', 10)
# make_histogram_data(counts_bottom_10, '%MinMax Negative Count per Window for bottom 10% log10(PHI)', 'Histogram of Negative %MinMax Count per Window for bottom 10% log10(PHI)', 10)
#
#make_histogram_data(norm_counts_top_10, '%MinMax Negative Count Normalized per Window for top 10% log10(PHI)', 'Histogram of Negative %MinMax Normalized Count per Window for top 10% log10(PHI)', 10)
# make_histogram_data(norm_counts_bottom_10, '%MinMax Negative Count Normalized per Window for bottom 10% log10(PHI)', 'Histogram of Negative %MinMax Normalized Count per Window for bottom 10% log10(PHI)', 10)

# make_histogram_data(df_top_10['MM_avg'], '%MinMax Avg for top 10% Disorder', 'Histogram of Average for top 10% Disorder', 10)
# make_histogram_data(df_bottom_10['MM_avg'], '%MinMax Avg for bottom 10% Disorder', 'Histogram of Average for bottom 10% Disorder', 10)
# make_histogram_data(df_top_10['MM_min'], '%MinMax Min for top 10% Disorder', 'Histogram of Minimum for top 10% Disorder', 10)
# make_histogram_data(df_bottom_10['MM_min'], '%MinMax Min for bottom 10% Disorder', 'Histogram of Minimum for bottom 10% Disorder', 10)

# l1 = np.log(df['abundance..molecules.cell.'])
# l2 = np.log(df['MM_median'])
#
# corr, pval = stats.spearmanr(l1, l2, nan_policy='omit')
# print(corr, pval)
#make_scatterplot_data(df['MM_avg'], l, "%MinMax Average", "Log10 of Abundance", "Plot of Average %MinMax vs Log10(Abundance) for Full Dataset")
# make_scatterplot_data(df_top_10['MM_min'], df_top_10['disorder_size'], "%MinMax Minimum", "Disorder Count", "Plot of Minimum %MinMax vs Disorder for Top 10% Dataset")
# make_scatterplot_data(df_bottom_10['MM_avg'], df_bottom_10['disorder_size'], "%MinMax Average", "Disorder Count", "Plot of Average %MinMax vs Disorder for Bottom 10% Dataset")
# make_scatterplot_data(df_bottom_10['MM_min'], df_bottom_10['disorder_size'], "%MinMax Minimum", "Disorder Count", "Plot of Minimum %MinMax vs Disorder for Bottom 10% Dataset")
#
# run_regressions(df)
# run_regressions(df_top_10)
# run_regressions(df_bottom_10)

#generate historic minMax plots
# row = df.loc[df['Gene'] == 'leuS']
# avgRRT = calc_avgRRT(list(row['mRNA_adjusted'])[0], aa_map_ecoli, aa_freq_ecoli, codon_freq_ecoli, codon_table)
# generate_minmax_plot(list(row['MinMax'])[0], 'leuS', avgRRT=avgRRT)
#
# row = df.loc[df['Gene'] == 'infB']
# avgRRT = calc_avgRRT(list(row['mRNA_adjusted'])[0], aa_map_ecoli, aa_freq_ecoli, codon_freq_ecoli, codon_table)
# generate_minmax_plot(list(row['MinMax'])[0], 'infB', avgRRT=avgRRT)

# row = df.loc[df['Gene'] == 'rpsT']
# avgRRT = calc_avgRRT(list(row['mRNA_adjusted'])[0], aa_map_ecoli, aa_freq_ecoli, codon_freq_ecoli, codon_table)
# generate_minmax_plot(list(row['MinMax'])[0], 'rpsT', mean=True, avgRRT=avgRRT)
# #
# row = df.loc[df['Gene'] == 'rplS']
# avgRRT = calc_avgRRT(list(row['mRNA_adjusted'])[0], aa_map_ecoli, aa_freq_ecoli, codon_freq_ecoli, codon_table)
# generate_minmax_plot(list(row['MinMax'])[0], 'rplS', mean=True, avgRRT=avgRRT)
# # #
# row = df.loc[df['Gene'] == 'rplD']
# avgRRT = calc_avgRRT(list(row['mRNA_adjusted'])[0], aa_map_ecoli, aa_freq_ecoli, codon_freq_ecoli, codon_table)
# generate_minmax_plot(list(row['MinMax'])[0], 'rplD', mean=True, avgRRT=avgRRT)
# # #
# row = df.loc[df['Gene'] == 'rplL']
# avgRRT = calc_avgRRT(list(row['mRNA_adjusted'])[0], aa_map_ecoli, aa_freq_ecoli, codon_freq_ecoli, codon_table)
# generate_minmax_plot(list(row['MinMax'])[0], 'rplL', mean=True, avgRRT=avgRRT)

#T-test code
#t, p = stats.ttest_ind(df_top_10['MM_avg'], df_bottom_10['MM_avg'], equal_var=False, nan_policy='omit')
#print("ttest_ind:            t = %g  p = %g" % (t, p))

# l = [elem for elem in df_top_10['MM_avg'] if not math.isnan(elem)]
# #print('avg =', statistics.mean(l), 'stdev =', statistics.stdev(l))
# abar, asd, na = statistics.mean(l), statistics.stdev(l), len(l)
#
# l = [elem for elem in df_bottom_10['MM_avg'] if not math.isnan(elem)]
# #print('avg =', statistics.mean(l), 'stdev =', statistics.stdev(l))
# bbar, bsd, nb = statistics.mean(l), statistics.stdev(l), len(l)
#
# t2, p2 = stats.ttest_ind_from_stats(abar, asd, na,
#                               bbar, bsd, nb,
#                               equal_var=False)
# print("ttest_ind_from_stats: t = %g  p = %g" % (t2, p2))
