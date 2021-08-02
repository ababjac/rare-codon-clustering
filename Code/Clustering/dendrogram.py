import pandas as pd
import numpy as np
import statistics
import math
from sklearn import linear_model
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy import stats
import random
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as col
import seaborn as sns
import plotly.figure_factory as ff
from scipy.cluster.hierarchy import dendrogram,linkage,fcluster,set_link_color_palette, cut_tree
from scipy.spatial.distance import squareform
from scipy_cut_tree_balanced import cut_tree_balanced
from sklearn import metrics


FILENAME = sys.argv[1]
COLNAME = sys.argv[2]
OUTFILE = sys.argv[3]
ERRORPLOT = sys.argv[4]
DENDROGRAMPLOT = sys.argv[5]
LINKAGE = sys.argv[6]
THRESHOLD = float(sys.argv[7])

print('Reading', FILENAME)
vals = pd.read_csv(FILENAME)
genes = vals[COLNAME].values.tolist()

vals.drop(labels=COLNAME, axis=1, inplace=True) #get rid of gene column


matrix = vals.values.tolist()
n = len(matrix)

condensed_matrix = [0]*round(n*(n-1) / 2)
for j in range(0, n):
    for i in range(0, j):
        condensed_matrix[n * i + j - ((i + 2) * (i + 1)) // 2] = matrix[j][i]

print('Computing Linkages...')
Z = linkage(condensed_matrix, method=LINKAGE.lower())

# print('Picking threshold with CH_index')
# thresholds = np.arange(0.05, 0.15, 0.01)
# errors = []
# for t in thresholds:
#     clusters = fcluster(Z, t, criterion='distance')
#     #print(len(set(clusters)))
#     errors.append(metrics.calinski_harabasz_score(matrix, list(clusters)))
#
# THRESHOLD = thresholds[errors.index(max(errors))]
#
# print('Plotting Elbow...')
# #Plot the elbow
# plt.plot(thresholds, errors, 'bx-')
# plt.xlabel('Threshold Value')
# plt.ylabel('CH Index')
# plt.title('Optimal Threshold Elbow Plot')
# plt.savefig(ERRORPLOT)
# plt.close()

print('Plotting Dendrogram...')
colors = cm.gist_ncar(np.arange(0, 1, 0.05))

colorlst=[]# empty list where you will put your colors
for i in range(len(colors)): #get for your color hex instead of rgb
    colorlst.append(col.to_hex(colors[i]))

set_link_color_palette(colorlst)

D = dendrogram(Z,labels=genes, color_threshold=THRESHOLD, above_threshold_color='gray')
plt.axhline(y=THRESHOLD, c='gray', lw=1, linestyle='dashed')

if COLNAME == 'Gene':
    plt.title('E. Coli')
else:
    plt.title('Yeast')

plt.savefig(DENDROGRAMPLOT)
plt.close()

print('Creating clusters using fcluster')
clusters = fcluster(Z, THRESHOLD, criterion='distance')

data = {'Gene' : genes, 'cluster_id' : list(clusters)}
df = pd.DataFrame(data)
length = len(set(clusters))

OUT = OUTFILE+'_t'+str(round(THRESHOLD, 3))+'.txt'

clusters = []
singleton_cluster = []
for i in range(1, length+1):
    l = df['Gene'][df['cluster_id'] == i].values.tolist()

    if len(l) == 1:
        singleton_cluster.append(l[0])
    else:
        clusters.append(l)

print('Writing to', OUT)
file = open(OUT, 'w')
for i in range(1, len(clusters)+1):
    file.write('Cluster '+str(i)+' - \n')
    l = clusters[i-1]

    for item in l:
        file.write(item+', ')

    file.write('\n\n')

file.write('Singletons - \n')
for item in singleton_cluster:
    file.write(item+', ')

file.write('\n\n')


print('\n\n')
