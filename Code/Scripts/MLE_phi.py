#courtesy of Lu
#Github: https://github.com/luzhixiu1996/MLE-Phi/blob/master/Code/MLE.py

import math
import numpy as np
from scipy import optimize
from scipy.optimize import minimize_scalar
from sklearn.preprocessing import minmax_scale
from statistics import mean, median
import pandas as pd

deltaEtaFile="../Data/Scer_Selection.csv"
deltaMFile="../Data/Scer_Mutation.csv"


#----------------------------UNMODIFIED----------------------------#
deltaEtaFile=open(deltaEtaFile)
lines=deltaEtaFile.readlines()
etaDict=dict()
for line in lines[1:]:
    splitList=line.split(",")
    if len(splitList[0])>0 and len(splitList[1])>0:
        etaDict[splitList[1]]=splitList[2]

deltaMFile=open(deltaMFile)
lines=deltaMFile.readlines()
mDict=dict()
for line in lines[1:]:
    splitList=line.split(",")
    if len(splitList[0])>0 and len(splitList[1])>0:
        mDict[splitList[1]]=splitList[2]


codontable = {
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

inverseTable=  {'S': ['TCA', 'TCC', 'TCT', 'TCG', 'AGC', 'AGT'], 'N': ['AAT', 'AAC'], 'W': ['TGG'],
                'D': ['GAT', 'GAC'], 'E': ['GAA', 'GAG'], 'F': ['TTC', 'TTT'], 'P': ['CCG', 'CCT', 'CCA', 'CCC'],
              'T': ['ACT', 'ACG', 'ACC', 'ACA'], 'G': ['GGG', 'GGC', 'GGT', 'GGA'],
              'R': ['AGA', 'CGT', 'CGA', 'AGG', 'CGG', 'CGC'], 'L': ['TTG', 'CTC', 'TTA', 'CTT', 'CTA', 'CTG'],
              'V': ['GTC', 'GTG', 'GTA', 'GTT'], 'A': ['GCA', 'GCT', 'GCC', 'GCG'], '*': ['TGA', 'TAA', 'TAG'],
              'H': ['CAC', 'CAT'], 'I': ['ATC', 'ATA', 'ATT'], 'K': ['AAG', 'AAA'], 'Y': ['TAT', 'TAC'],
              'M': ['ATG'], 'C': ['TGC', 'TGT'], 'Q': ['CAG', 'CAA']}

synonymonDict={}
for key in inverseTable:
    valueList=inverseTable[key]
    for value in valueList:
        synonymonDict[value]=valueList


def loadSequence(sequence):
    startCodon="ATG"
    stopCodonList=["TAG","TAA","TGA"]
    codonList=[]
    i=0
    while(i<len(sequence)):
        codon=sequence[i:i+3]
        if len(codon)==3:
            codonList.append(codon)
        i+=3
    actualCodonList=[]
    started=False
    for codon in codonList:
        if codon in stopCodonList:
            break
        if started:
            actualCodonList.append(codon)
        if codon==startCodon:
            started=True
    codonList=actualCodonList
   # print "codon readed successful, the number of codon in this sequence is %d"%(len(codonList))
    return codonList

#this method removes sequence that cant be handled by mDict, etaDict and "TGG",
def parseSequence(sequence):
    i=0
    startCodon="ATG"
    stopCodonList=["TAG","TAA","TGA"]
    parsedSequence=""
    while(i<len(sequence)):
        codon=sequence[i:i+3]
        if len(codon)==3:
            if (codon in mDict and codon in etaDict) or (codon in startCodon) or (codon in stopCodonList) and ("TGG" not in codon):
                parsedSequence+=codon
        i+=3
    return parsedSequence


def cutSequence(seq):
    sequence=seq
    startCodon="ATG"
    stopCodonList=["TAG","TAA","TGA"]
    codonList=[]
    i=0
    while(i<len(sequence)):
        codon=sequence[i:i+3]
        if len(codon)==3:
            codonList.append(codon)
        i+=3
    actualCodonList=[]
    started=FalsewindowSize
    for codon in codonList:
        if codon in stopCodonList:
            break
        if started:
            actualCodonList.append(codon)
        if codon==startCodon:
            started=True
    codonList=actualCodonList
   # print "codon readed successful, the number of codon in this sequence is %d"%(len(codonList))
    return codonList


def roundList(lst):
    roundedList=[]
    decimalPlaces=6
    for i in lst:
        roundedList.append(round(i,decimalPlaces))
    return roundedList

phiDict=dict()
def method4(codonList):
    global mDict
    global etaDict
    global phiDict
    phiList=[]

    for codon in codonList:
        if codon in phiDict:
             phiList.append(phiDict[codon])
        else:
            maxprob=0.0
            selectedPhi=0.0
            rangeList=[]
            for i in range(1,101):
                rangeList.append(i/100.0)
            for phi in rangeList:
                if codon in mDict:
                    deltaM=float(mDict[codon])
                else:
                    deltaM=1.0
                if codon in etaDict:
                    deltaEta=float(etaDict[codon])
                else:
                    deltaEta=1.0
                global synonymonDict
                synoList=synonymonDict[codon]
                divisor=np.exp(-1.00*deltaM-(deltaEta*phi))
                dividant=0
                for syno in synoList:
                    if syno in mDict:
                        deltaM=float(mDict[syno])
                    else:
                        deltaM=1.0
                    if syno in etaDict:
                        deltaEta=float(etaDict[syno])
                    else:
                        deltaEta=1.0
                    tmp=np.exp((-1.00*deltaM)-(deltaEta*phi))
    #                print((-1.00*deltaM)-(deltaEta*phi))
    #                print (tmp)
                    dividant+=tmp
                if not dividant==0:
                    prob=divisor/dividant
                else:
                    prob=0

                if prob>maxprob:
                    maxprob=prob
                    selectedPhi=phi
            if not selectedPhi==0:
                phiDict[codon]=selectedPhi
                phiList.append(selectedPhi)
            else:
                phiList.append(0.001)
                print ("found 0")
    return phiList

def setWindow(inputList,size):
    windowList=[]
    windowSize=size
    cnt=0

    while True:
        cnt+=1
        if cnt+windowSize>len(inputList):
            break
        selectedList=inputList[cnt:cnt+windowSize]
        sum=0
        for i in selectedList:
            sum+=float(i)
        average=sum/len(selectedList)
        windowList.append(average)
    return windowList

from scipy.stats.mstats import gmean


def cal_mle_windows(sequence,windowSize=10):
    codonList=loadSequence(sequence)
    phiList=method4(codonList)# df['Phi_list'] = np.array(my_windows)
    windowList=setWindow(phiList,windowSize)
    return windowList


def calPhiForGene(sequence):
    codonList=loadSequence(sequence)
    phiList=method4(codonList)
    avg=gmean(phiList)
    return avg
#----------------------------UNMODIFIED----------------------------#

import matplotlib.pyplot as plt
import scipy.stats as stats


def make_histogram_data(data, xlabel, plt_title, bins=None, num_bins=10):
    if bins == None:
        bins = np.linspace(math.floor(min(data)), math.ceil(max(data)), num_bins)

    plt.hist(data, bins)
    plt.title(plt_title)
    plt.xlabel(xlabel)
    plt.ylabel('Count')
    plt.show()

def generate_mlphi_plot(mlphi_list, title, avgRRT=None):
    indexes = list(range(1, len(mlphi_list)+1))
    geo_mean = 1

    for elem in mlphi_list:
        geo_mean *= elem

    geo_mean = geo_mean**(1/float(len(mlphi_list)))

    #print(minmax_list)
    mlphi_list = [None if i==0.0 else i for i in mlphi_list]

    plt.plot(indexes, mlphi_list, label = "ML-Phi", color = '#0000CC')
    if avgRRT != None:
        plt.plot(indexes, avgRRT, label = "Average RRT", color = '#808080')
    plt.plot([0,len(mlphi_list)],[geo_mean, geo_mean],'r--')
    plt.legend(loc = 4)
    plt.xlim(0,len(mlphi_list))
    plt.ylim(0,1.0)
    #plt.fill_between(num, stdevbel0,stdevabv0, color = '#808080', alpha = .1)
    plt.title(title)
    plt.xlabel("Center Of Codon Window")
    plt.ylabel("ML-Phi Windows")

    plt.show()


def add_col_to_df(Phi_list, df):
    Phi_min_list = []
    Phi_max_list = []
    Phi_avg_list = []
    Phi_median_list = []

    for l in Phi_list:
        if not l:
            Phi_min_list.append(None)
            Phi_avg_list.append(None)
            Phi_max_list.append(None)
            Phi_median_list.append(None)
        else:
            Phi_min_list.append(min(l))
            Phi_avg_list.append(sum(l) / len(l))
            Phi_max_list.append(max(l))
            Phi_median_list.append(median(l))

    #print(Phi_min_list)
    df['Phi_min'] = Phi_min_list
    df['Phi_avg'] = Phi_avg_list
    df['Phi_max'] = Phi_max_list
    df['Phi_median'] = Phi_median_list

    return df

#using the above to calculate windows and estimate for each gene
df = pd.read_csv('../Data/Yeast/YeastSeqs.csv')

my_windows = []
my_estimates = []
for seq in df['seq']:
    my_w = cal_mle_windows(seq, windowSize=9)
    my_e = calPhiForGene(seq)
    my_windows.append(my_w)
    my_estimates.append(my_e)

df = add_col_to_df(my_windows, df)
df['phi_estimate'] = my_estimates

df.to_csv('../Files/Phi/Phi_estimates_yeast.csv')

# abundance_list = df['abundance..molecules.cell.'].values.tolist()
# print('Spearman Correlation for Median Phi vs Abundance = ', stats.spearmanr(np.log10(df['Phi_median'].values.tolist()), abundance_list, nan_policy='omit'))
# print('Spearman Correlation for Average Phi vs Abundance = ', stats.spearmanr(np.log10(df['Phi_avg'].values.tolist()), abundance_list, nan_policy='omit'))
# print('Spearman Correlation for Min Phi vs Abundance = ', stats.spearmanr(np.log10(df['Phi_min'].values.tolist()), abundance_list, nan_policy='omit'))

# shuf1 = 'ATGCATCACCATCACCATCACCATAACTATACAAAATTTGATGTAAAAAATTGGGTTCGCCGTGAGCATTTTGAGTTTTATCGGCATCGTTTACCATGTGGTTTTAGCTTAACAAGCAAAATTGATATCACGACGTTAAAAAAGTCATTGGATGATTCAGCGTATAAGTTTTACCCCGTGATGATATACTTAATTGCCCAAGCAGTTAACCAATTTGACGAGCTCCGAATGGCGATTAAAGATGATGAACTCATTGTGTGGGACTCCGTCGATCCGCAATTTACTGTATTCCATCAGGAGACTGAAACTTTTAGTGCGCTATCCTGCCCCTATAGCTCGGACATCGATCAATTTATGGTGAATTACCTGTCGGTCATGGAAAGGTATAAGTCCGACACAAAATTATTTCCCCAGGGCGTGACGCCAGAAAACCATCTAAACATATCCGCGCTGCCATGGGTTAACTTTGACTCGTTCAACTTAAATGTAGCTAACTTTACTGACTATTTTGCGCCGATCATAACCATGGCCAAATACCAACAGGAAGGGGATAGGCTCTTGTTACCCTTGAGCGTCCAGGTGCATCACGCCGTGTGTGACGGCTTTCATGTTGCACGCTTCATTAACAGATTACAGGAGTTATGTAATTCTAAATTGAAGTAA'
#
# windows = cal_mle_windows(shuf1, windowSize=9)
# generate_mlphi_plot(windows, '')
# flatten = [elem for list in my_windows for elem in list]
# #my_avgs = [sum(list)/len(list) for list in my_windows if len(list) != 0]
# #print(flatten)
# make_histogram_data(flatten, 'ML-Phi Window Values', 'Distribution of ML-Phi Windows for All Sequences', num_bins=20)
# make_histogram_data(my_estimates, 'ML-Phi Estimates (one for each Gene)', 'Distribution of ML-Phi Estimates for All Sequences', num_bins=20)

# df['Phi_list'] = np.array(my_windows)
# df['Phi_est'] = my_estimates
#
# df_top_50 = df.sort_values(by='disorder_ratio', ascending=False).head(round(len(df)*.1))
# df_bottom_50 = df.sort_values(by='disorder_ratio', ascending=True).head(round(len(df)*.1))
#
# df_top_50 = df[df['disorder_ratio'] > 0.5]
# df_bottom_50 = df[df['disorder_ratio'] <= 0.5]

# print(df_top_50)
# print(df_bottom_50)

# counts_top_50 = []
# count = 0
#
# for l, est in zip(df_top_50['Phi_list'], df_top_50['Phi_est']):
#     for elem in l:
#         if elem < est:
#             count += 1
#
#     counts_top_50.append(count)
#     count = 0
#
# counts_bottom_50 = []
# count = 0
#
# for l, est in zip(df_bottom_50['Phi_list'], df_bottom_50['Phi_est']):
#     for elem in l:
#         if elem < est:
#             count += 1
#
#     counts_bottom_50.append(count)
#     count = 0

# print(set(counts_top_50))
# print()
# print(set(counts_bottom_50))

# make_histogram_data(counts_top_50, 'Count of Windows below ML-Phi Estimate', 'Atypical ML-Phi Window Count for top 10% Disorder', num_bins=10)
# make_histogram_data(counts_bottom_50, 'Count of Windows below ML-Phi Estimate', 'Atypical ML-Phi Window Count for bottom 10% Disorder', num_bins=10)
#
# print('Spearman Correlation for ML-PHI min vs Disorder = ', stats.spearmanr(df['Phi_est'].values.tolist(), df['disorder_ratio'].values.tolist(), nan_policy='omit'))
# print('Spearman Correlation for ML-PHI min vs Disorder = ', stats.spearmanr(df['Phi_est'].values.tolist(), df['disorder_size'].values.tolist(), nan_policy='omit'))
#
# row = df.loc[df['Gene'] == 'rpsT']
# generate_mlphi_plot(list(row['Phi_list'])[0], 'rpsT')
#
# row = df.loc[df['Gene'] == 'rplS']
# generate_mlphi_plot(list(row['Phi_list'])[0], 'rplS')
#
# row = df.loc[df['Gene'] == 'rplD']
# generate_mlphi_plot(list(row['Phi_list'])[0], 'rplD')
#
# row = df.loc[df['Gene'] == 'rplL']
# generate_mlphi_plot(list(row['Phi_list'])[0], 'rplL')

# df.to_csv('../Files/Phi_estimates_full.csv')
