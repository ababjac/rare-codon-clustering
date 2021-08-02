import pandas as pd


# df = pd.read_csv('../Data/disorder_full.csv')
#
# file = open('../Data/ecoli_stat.txt', 'r')
# lines = file.readlines()
#
# CAI_list = []
# CBI_list = []
# Fop_list = []
# Nc_list = []
# for line in lines[1:]:
#     _, CAI, CBI, Fop, Nc = line.split(',')
#     CAI_list.append(CAI)
#     CBI_list.append(CBI)
#     Fop_list.append(Fop)
#     Nc_list.append(Nc)
#
# Nc_list = [float('Nan') if elem.__contains__('*') else float(elem) for elem in Nc_list]
#
# df['CAI'] = CAI_list
# df['CBI'] = CBI_list
# df['Fop'] = Fop_list
# df['Nc'] = Nc_list
#
# df.drop(df.columns[df.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
#
# df.to_csv('../Data/disorder_new_full.csv')


# df1 = pd.read_csv('../Data/disorder_new_full.csv')
# df2 = pd.read_csv('../Data/merge.csv')
#
# df = pd.merge(df1, df2, how='inner', on='Gene')
# df.drop(df.columns[df.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)
#
# df.to_csv('../Data/Abundance_new.csv')

df1 = pd.read_csv('../Files/MM/minMax_full_adjusted.csv')
df2 = pd.read_csv('../Data/Ecoli/Abundance_new.csv')

df = pd.merge(df1, df2, how='inner', on='Gene')
df.drop(df.columns[df.columns.str.contains('unnamed',case = False)],axis = 1, inplace = True)

df.to_csv('../Data/Ecoli/ecoli_complete.csv')
