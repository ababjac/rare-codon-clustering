import pandas as pd
import numpy as np
import statistics
import math
from sklearn import linear_model
from sklearn.cluster import KMeans
import statsmodels.api as sm
from scipy import stats
import random
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as col
import seaborn as sns
import plotly.figure_factory as ff
from scipy.cluster.hierarchy import dendrogram,linkage,fcluster,set_link_color_palette
from sklearn import metrics
from scipy.spatial.distance import cdist
from mpl_toolkits import mplot3d

#-----------------------------------------------------------------------------#
#Average MM versus Phi Est
# f1 = pd.read_csv('../Files/MM/minMax_full_adjusted.csv')
# f2 = pd.read_csv('../Files/Phi/Phi_estimates_full.csv')
#
# sub = pd.DataFrame()
# sub['MM'] = f1['MM_avg']
# sub['Phi'] = f2['Phi_est']
# sub = sub.dropna()
#
#
# #GET ELBOW PLOT
# distortions = []
# K = range(1,20)
# for k in K:
#     kmeanModel = KMeans(n_clusters=k).fit(sub)
#     kmeanModel.fit(sub)
#     distortions.append(sum(np.min(cdist(sub, kmeanModel.cluster_centers_, 'euclidean'), axis=1)) / sub.shape[0])
#
# # Plot the elbow
# plt.plot(K, distortions, 'bx-')
# plt.xlabel('k')
# plt.ylabel('Distortion')
# plt.title('Optimal k Elbow Plot')
# plt.show()

# kmeans = KMeans(n_clusters=5).fit(sub)
# centroids = kmeans.cluster_centers_
# print(centroids)
#
# plt.scatter(sub['MM'], sub['Phi'], c=kmeans.labels_.astype(float), s=50, alpha=0.5)
# plt.scatter(centroids[:, 0], centroids[:, 1], c='red', s=50)
# plt.xlabel('Average %MM')
# plt.ylabel('ML-Phi')
# plt.title('Kmeans Clusters for k=5')
# plt.show()

#-----------------------------------------------------------------------------#
#Median MM versus Phi Est
f1 = pd.read_csv('../Files/MM/minMax_full_adjusted.csv')
f2 = pd.read_csv('../Files/Phi/Phi_estimates_full.csv')

sub = pd.DataFrame()
sub['MM'] = f1['MM_median']
sub['Phi'] = f2['Phi_est']
sub = sub.dropna()


#GET ELBOW PLOT
# distortions = []
# K = range(1,20)
# for k in K:
#     kmeanModel = KMeans(n_clusters=k).fit(sub)
#     kmeanModel.fit(sub)
#     distortions.append(sum(np.min(cdist(sub, kmeanModel.cluster_centers_, 'euclidean'), axis=1)) / sub.shape[0])
#
# # Plot the elbow
# plt.plot(K, distortions, 'bx-')
# plt.xlabel('k')
# plt.ylabel('Distortion')
# plt.title('Optimal k Elbow Plot')
# plt.show()

kmeans = KMeans(n_clusters=5).fit(sub)
sub['labels'] = kmeans.labels_
sub['Gene'] = f1['Gene']
#print(set(kmeans.labels_))

# for i in range(len(set(kmeans.labels_))):
#     clust = sub['Gene'][sub['labels'] == i]
#     print('Cluster', i+1, ':')
#     print(clust.values.tolist())

centroids = kmeans.cluster_centers_
# print(centroids)
#
plt.scatter(sub['MM'], sub['Phi'], c=kmeans.labels_.astype(float), s=50, alpha=0.5)
plt.scatter(centroids[:, 0], centroids[:, 1], c='red', s=50)
plt.xlabel('Median %MM')
plt.ylabel('ML-Phi')
plt.title('Kmeans Clusters for k=5')
plt.legend()
plt.show()


#-----------------------------------------------------------------------------#
#Average MM versus Abundance
# f1 = pd.read_csv('../Files/MM/minMax_full_adjusted.csv')
# f2 = pd.read_csv('../Data/Ecoli/merge.csv')
#
# sub = f1.merge(f2.set_index('Gene'), on='Gene', how='inner')
# sub['log_abundance'] = round(np.log10(sub['abundance..molecules.cell.']), 10)
# sub = sub[['MM_avg', 'log_abundance']]
# #sub.dropna()
# sub=sub[~sub.isin([np.nan, np.inf, -np.inf]).any(1)]
#print(sub)

#GET ELBOW PLOT
# distortions = []
# K = range(1,20)
# for k in K:
#     kmeanModel = KMeans(n_clusters=k).fit(sub)
#     kmeanModel.fit(sub)
#     distortions.append(sum(np.min(cdist(sub, kmeanModel.cluster_centers_, 'euclidean'), axis=1)) / sub.shape[0])
#
# # Plot the elbow
# plt.plot(K, distortions, 'bx-')
# plt.xlabel('k')
# plt.ylabel('Distortion')
# plt.title('Optimal k Elbow Plot')
# plt.show()
#
# kmeans = KMeans(n_clusters=5).fit(sub)
# centroids = kmeans.cluster_centers_
# print(centroids)
#
# plt.scatter(sub['MM_avg'], sub['log_abundance'], c=kmeans.labels_.astype(float), s=50, alpha=0.5)
# plt.scatter(centroids[:, 0], centroids[:, 1], c='red', s=50)
# plt.xlabel('Average %MM')
# plt.ylabel('Log10 Abundance')
# plt.title('Kmeans Clusters for k=5')
# plt.show()

#-----------------------------------------------------------------------------#
#Median MM versus Abundance
# f1 = pd.read_csv('../Files/MM/minMax_full_adjusted.csv')
# f2 = pd.read_csv('../Data/Ecoli/merge.csv')
#
# sub = f1.merge(f2.set_index('Gene'), on='Gene', how='inner')
# sub['log_abundance'] = round(np.log10(sub['abundance..molecules.cell.']), 10)
# sub = sub[['MM_median', 'log_abundance']]
# #sub.dropna()
# sub=sub[~sub.isin([np.nan, np.inf, -np.inf]).any(1)]
#print(sub)

# #GET ELBOW PLOT
# distortions = []
# K = range(1,20)
# for k in K:
#     kmeanModel = KMeans(n_clusters=k).fit(sub)
#     kmeanModel.fit(sub)
#     distortions.append(sum(np.min(cdist(sub, kmeanModel.cluster_centers_, 'euclidean'), axis=1)) / sub.shape[0])
#
# # Plot the elbow
# plt.plot(K, distortions, 'bx-')
# plt.xlabel('k')
# plt.ylabel('Distortion')
# plt.title('Optimal k Elbow Plot')
# plt.show()
#
# kmeans = KMeans(n_clusters=5).fit(sub)
# centroids = kmeans.cluster_centers_
# print(centroids)
#
# plt.scatter(sub['MM_median'], sub['log_abundance'], c=kmeans.labels_.astype(float), s=50, alpha=0.5)
# plt.scatter(centroids[:, 0], centroids[:, 1], c='red', s=50)
# plt.xlabel('Median %MM')
# plt.ylabel('Log10 Abundance')
# plt.title('Kmeans Clusters for k=5')
# plt.show()

#-----------------------------------------------------------------------------#
#Phi versus Abundance
# f1 = pd.read_csv('../Files/MM/minMax_full_adjusted.csv')
# f2 = pd.read_csv('../Data/Ecoli/merge.csv')
#
# sub = f1.merge(f2.set_index('Gene'), on='Gene', how='inner')
# sub['log_abundance'] = round(np.log10(sub['abundance..molecules.cell.']), 10)
# sub = sub[['log10.PHI', 'log_abundance']]
# #sub.dropna()
# sub=sub[~sub.isin([np.nan, np.inf, -np.inf]).any(1)]
# print(sub)

#GET ELBOW PLOT
# distortions = []
# K = range(1,20)
# for k in K:
#     kmeanModel = KMeans(n_clusters=k).fit(sub)
#     kmeanModel.fit(sub)
#     distortions.append(sum(np.min(cdist(sub, kmeanModel.cluster_centers_, 'euclidean'), axis=1)) / sub.shape[0])
#
# # Plot the elbow
# plt.plot(K, distortions, 'bx-')
# plt.xlabel('k')
# plt.ylabel('Distortion')
# plt.title('Optimal k Elbow Plot')
# plt.show()
#
# kmeans = KMeans(n_clusters=5).fit(sub)
# centroids = kmeans.cluster_centers_
# print(centroids)
#
# plt.scatter(sub['log10.PHI'], sub['log_abundance'], c=kmeans.labels_.astype(float), s=50, alpha=0.5)
# plt.scatter(centroids[:, 0], centroids[:, 1], c='red', s=50)
# plt.xlabel('Log10 Phi')
# plt.ylabel('Log10 Abundance')
# plt.title('Kmeans Clusters for k=5')
# plt.show()

#-----------------------------------------------------------------------------#
#Average MM versus Abundance versus Phi
# f1 = pd.read_csv('../Files/MM/minMax_full_adjusted.csv')
# f2 = pd.read_csv('../Data/Ecoli/merge.csv')
#
# sub = f1.merge(f2.set_index('Gene'), on='Gene', how='inner')
# sub['log_abundance'] = round(np.log10(sub['abundance..molecules.cell.']), 10)
# sub = sub[['MM_avg', 'log_abundance', 'log10.PHI']]
# sub=sub[~sub.isin([np.nan, np.inf, -np.inf]).any(1)]

#GET ELBOW PLOT
# distortions = []
# K = range(1,20)
# for k in K:
#     kmeanModel = KMeans(n_clusters=k).fit(sub)
#     kmeanModel.fit(sub)
#     distortions.append(sum(np.min(cdist(sub, kmeanModel.cluster_centers_, 'euclidean'), axis=1)) / sub.shape[0])
#
# # Plot the elbow
# plt.plot(K, distortions, 'bx-')
# plt.xlabel('k')
# plt.ylabel('Distortion')
# plt.title('Optimal k Elbow Plot')
# plt.show()

# kmeans = KMeans(n_clusters=5).fit(sub)
# centroids = kmeans.cluster_centers_
# print(centroids)
#
# #print(kmeans.labels_)
#
# fig = plt.figure()
# ax = plt.axes(projection='3d')
# ax.scatter3D(sub['MM_avg'], sub['log_abundance'], sub['log10.PHI'], c=kmeans.labels_.astype(float), s=50, alpha=0.5)
# #plt.scatter(centroids[:, 0], centroids[:, 1], c='red', s=50)
# ax.set_xlabel('Average %MM')
# ax.set_ylabel('Log10 Abundance')
# ax.set_zlabel('Log10 Phi')
# plt.title('Kmeans Clusters for k=5')
# plt.show()
