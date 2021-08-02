import os
import sys
import pandas as pd

def calc_avg_per_cluster(lines, df, colname):

    d = {}

    total_avg = sum(df['MM_avg']) / len(df)

    for line in lines:
        if line.__contains__('Cluster'):
            text = line.split(' ')
            key = text[0]+text[1]
            continue

        list = line.split(', ')
        list.remove('\n')
        #print(list)
        s = 0

        for elem in list:
            #print(elem, df[df[colname] == elem]['MM_avg'].values[0])
            s += df[df[colname] == elem]['MM_avg'].values[0]

        cluster_avg = s / len(list)
        #print(cluster_avg)

        d[key] = cluster_avg

    d['Total'] = total_avg
    return d


#------------------------------------------------------------------------------------------------------#

# DF_ecoli = pd.read_csv('../Files/MM/minMax_full.csv')
# #print(DF_ecoli[DF_ecoli['Gene'] == 'thrA']['MM_avg'].values[0])
# file = open('../Files/Clusters/Hierarchical/Centroid/Full/hclust_HD_ecoli_t0.28.txt', 'r')
# dct = calc_avg_per_cluster(file.readlines(), DF_ecoli, 'Gene')
# print(dct)

FOLDER = 'Full'
#FOLDER = 'Omit10'
#FOLDER = 'Omit25'

LINKAGE = 'Centroid/'
#LINKAGE = 'Single/'

DIRECTORY = '../Files/Clusters/Hierarchical/'+LINKAGE+FOLDER+'/'
#SIG = '0.05'

DF_ecoli = pd.read_csv('../Files/MM/minMax_full.csv')
DF_yeast = pd.read_csv('../Files/MM/minMax_full_yeast.csv')

with os.scandir(DIRECTORY) as d:
    for entry in d:
        if entry.name.endswith('.txt') and entry.is_file():
            input_path = os.path.join(DIRECTORY, entry.name)
            file = open(input_path, 'r')
            lines = file.readlines()

            if entry.name.__contains__('ecoli'):
                df = DF_ecoli
                colname = 'Gene'
            else:
                df = DF_yeast
                colname = 'locus_tag'

        dict = calc_avg_per_cluster(lines, df, colname)
        print(entry.name)
        print(dict)
        print()
