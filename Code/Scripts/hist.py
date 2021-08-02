import pandas as pd
import matplotlib.pyplot as plt
import statistics
import numpy as np

BINS = np.arange(0, 1, 0.01)


def make_histogram_data(X, xlabel, plt_title, bins):
    plt.hist(X, bins)
    plt.title(plt_title)
    plt.xlabel(xlabel)
    plt.ylabel('Count')
    plt.show()


file = open('../Files/Clusters/Ecoli/Kmeans/MM-median_Phi_k5_clusters_ecoli.txt')
lines = file.readlines()

clusters = []
for line in lines:
    if line.__contains__('Cluster'):
        continue

    clust = line.split(',')
    clust = [elem.split('\'')[1] for elem in clust]

    clusters.append(clust)

#print(clusters)
G = pd.read_csv('../Data/Ecoli/full.csv')
genes = G['Gene'].values.tolist()
#[print(elem) for elem in genes if genes.count(elem) > 1]


df = pd.read_csv('../Files/MM/KS_ecoli.csv')
df.drop(df.columns[df.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df.columns = genes
del df['dnaX']
del df['copA']
#print(df)

#for clust in clusters:
avg_pvals = [statistics.mean(df[elem]) for elem in clusters[0]]
make_histogram_data(avg_pvals, 'Average KS Statistic P-value Per Protein', 'Distribution of P-Values By Cluster E. Coli (1)', BINS)
# avg_pvals = [statistics.mean(df[elem]) for elem in clusters[1] if elem in df.columns]
# make_histogram_data(avg_pvals, 'Average KS Statistic P-value Per Protein', 'Distribution of P-Values By Cluster E. Coli (2)', BINS)
# avg_pvals = [statistics.mean(df[elem]) for elem in clusters[2]]
# make_histogram_data(avg_pvals, 'Average KS Statistic P-value Per Protein', 'Distribution of P-Values By Cluster E. Coli (3)', BINS)
# avg_pvals = [statistics.mean(df[elem]) for elem in clusters[3] if elem in df.columns]
# make_histogram_data(avg_pvals, 'Average KS Statistic P-value Per Protein', 'Distribution of P-Values By Cluster E. Coli (4)', BINS)
# avg_pvals = [statistics.mean(df[elem]) for elem in clusters[4]]
# make_histogram_data(avg_pvals, 'Average KS Statistic P-value Per Protein', 'Distribution of P-Values By Cluster E. Coli (5)', BINS)
