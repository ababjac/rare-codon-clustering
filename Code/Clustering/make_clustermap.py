import pandas as pd
import seaborn as sns
import scipy.spatial as sp
import scipy.cluster.hierarchy as hc
from sklearn.datasets import load_iris
import matplotlib.pyplot as plt
import sys


#EXAMPLE#
# sns.set(font="monospace")
#
# iris = load_iris()
# X, y = iris.data, iris.target
# DF = pd.DataFrame(X, index = ["iris_%d" % (i) for i in range(X.shape[0])], columns = iris.feature_names)
# DF_corr = DF.T.corr()
# DF_dism = 1 - DF_corr
# print(DF_dism)
#
# # distance matrix
# linkage = hc.linkage(sp.distance.squareform(DF_dism), method='average')
# sns.clustermap(DF_dism, row_linkage=linkage, col_linkage=linkage)
# plt.show()

FILENAME = sys.argv[1]
COLNAME = sys.argv[2]
OUTFILE = sys.argv[3]
LINKAGE = sys.argv[4]

print('Reading', FILENAME)
vals = pd.read_csv(FILENAME)
#genes = vals[COLNAME].values.tolist()

vals.set_index(COLNAME, drop=True, inplace=True)
#print(vals)


matrix = vals.values.tolist()
n = len(matrix)

condensed_matrix = [0]*round(n*(n-1) / 2)
for j in range(0, n):
    for i in range(0, j):
        condensed_matrix[n * i + j - ((i + 2) * (i + 1)) // 2] = matrix[j][i]

print('Computing Linkages...')
Z = hc.linkage(condensed_matrix, method=LINKAGE.lower())

print('Creating Clustermap')
sns.clustermap(vals, row_linkage=Z, col_linkage=Z,  cmap='YlGnBu')
plt.savefig(OUTFILE)
plt.close()
#plt.show()

print()
