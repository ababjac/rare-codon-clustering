import pandas as pd

##########ECOLI############
# df1 = pd.read_csv('../Data/full2.csv')
#
# file = open('../Data/stat.txt', 'r')
# lines = file.readlines()
#
# l = [line.split() for line in lines[1:] if line.split()]
# df2 = pd.DataFrame(l, columns=['NC_ID', 'disorder_size', 'gene_size', 'disorder_ratio'])
# df2['NC_ID'] = df2['NC_ID'].str.replace('>', '')
#
# df = pd.merge(df1, df2, on='NC_ID')
# df.drop(df.columns[df.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
#
# df.to_csv('../Data/disorder_full.csv')

##########YEAST##########
df_1 = pd.read_csv('../Files/Phi/Phi_estimates_yeast.csv')
df_2 = pd.read_csv('../Files/MM/minMax_full_yeast.csv')

df1 = pd.merge(df_1, df_2, on='locus_tag')
df1.drop(df1.columns[df1.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
df1.drop(df1.columns[df1.columns.str.contains('y',case = False)],axis = 1, inplace = True)
df1.columns = ['gene_ID', 'x', 'description', 'locus_tag', 'transcript_id',
       'seq', 'Phi_min', 'Phi_avg', 'Phi_max', 'Phi_median', 'phi_estimate',
       'mRNA_adjusted', 'MM_list', 'MM_min', 'MM_avg', 'MM_max', 'MM_median']


file = open('../Data/Yeast/updated.yeast.txt')
lines = file.readlines()

l = [line.split() for line in lines[1:] if line.split()]
df2 = pd.DataFrame(l, columns=['locus_tag', 'disorder_size', 'gene_size', 'disorder_ratio'])

df = pd.merge(df1, df2, on='locus_tag')
#df.drop(df.columns.str.contains('unnamed', case=False), axis=1, inplace=True)

df.to_csv('../Data/Disorder/disorder_full_yeast.csv')
