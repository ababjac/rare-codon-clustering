from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt


# file1 = open('Goatools/KS/output_yeast_heir_KS_t50_pvals.txt')
# file2 = open('Goatools/KS/output_yeast_heir_KS_t100_pvals.txt')
# file3 = open('Goatools/KS/output_yeast_heir_KS_t150_pvals.txt')
#
# l_50_BP = []
# l_50_MF = []
# l_50_CC = []
# for line in file1.readlines():
#     if line.__contains__('Cluster'):
#         continue
#     elif line == '\n':
#         continue
#
#     term = line.partition('\t')[0]
#
#     if line.__contains__('BP'):
#         l_50_BP.append(term)
#     elif line.__contains__('MF'):
#         l_50_MF.append(term)
#     elif line.__contains__('CC'):
#         l_50_CC.append(term)
#
# l_100_BP = []
# l_100_MF = []
# l_100_CC = []
# for line in file2.readlines():
#     if line.__contains__('Cluster'):
#         continue
#     elif line == '\n':
#         continue
#
#     term = line.partition('\t')[0]
#
#     if line.__contains__('BP'):
#         l_100_BP.append(term)
#     elif line.__contains__('MF'):
#         l_100_MF.append(term)
#     elif line.__contains__('CC'):
#         l_100_CC.append(term)
#
# l_150_BP = []
# l_150_MF = []
# l_150_CC = []
# for line in file3.readlines():
#     if line.__contains__('Cluster'):
#         continue
#     elif line == '\n':
#         continue
#
#     term = line.partition('\t')[0]
#
#     if line.__contains__('BP'):
#         l_150_BP.append(term)
#     elif line.__contains__('MF'):
#         l_150_MF.append(term)
#     elif line.__contains__('CC'):
#         l_150_CC.append(term)
#
# #BP
# A = len(l_50_BP)
# B = len(l_100_BP)
# C = len(l_150_BP)
# A_B = len(set([elem for elem in l_50_BP if elem in l_100_BP]))
# A_C = len(set([elem for elem in l_50_BP if elem in l_150_BP]))
# B_C = len(set([elem for elem in l_100_BP if elem in l_150_BP]))
# A_B_C = len(set([elem for elem in l_150_BP if elem in l_100_BP and elem in l_50_BP]))
#
# venn3(subsets=(A, B, A_B, C, A_C, B_C, A_B_C), set_labels=('BP-t50', 'BP-t100', 'BP-t150'))
# plt.show()
#
# #MF
# A = len(l_50_MF)
# B = len(l_100_MF)
# C = len(l_150_MF)
# A_B = len(set([elem for elem in l_50_MF if elem in l_100_MF]))
# A_C = len(set([elem for elem in l_50_MF if elem in l_150_MF]))
# B_C = len(set([elem for elem in l_100_MF if elem in l_150_MF]))
# A_B_C = len(set([elem for elem in l_150_MF if elem in l_100_MF and elem in l_50_MF]))
#
# venn3(subsets=(A, B, A_B, C, A_C, B_C, A_B_C), set_labels=('MF-t50', 'MF-t100', 'MF-t150'))
# plt.show()
#
# #CC
# A = len(l_50_CC)
# B = len(l_100_CC)
# C = len(l_150_CC)
# A_B = len(set([elem for elem in l_50_CC if elem in l_100_CC]))
# A_C = len(set([elem for elem in l_50_CC if elem in l_150_CC]))
# B_C = len(set([elem for elem in l_100_CC if elem in l_150_CC]))
# A_B_C = len(set([elem for elem in l_150_CC if elem in l_100_CC and elem in l_50_CC]))
#
# venn3(subsets=(A, B, A_B, C, A_C, B_C, A_B_C), set_labels=('CC-t50', 'CC-t100', 'CC-t150'))
# plt.show()

##############comparing pvals and stats clusters##########################
# file_HD = open('Goatools/HD_revised/output_ecoli_hier_HD_t150.txt')
# file_stats = open('Goatools/KS/output_ecoli_heir_KS_t150_stats.txt')

file_HD = open('Goatools/KS/output_yeast_heir_KS_t50_stats.txt')
file_stats = open('Goatools/Shuf/output_yeast_heir_KS_t50_shuf.txt')

lines_HD = file_HD.readlines()
lines_stats = file_stats.readlines()

goterms_HD = []
for line in lines_HD:
    if line.__contains__('Cluster'):
        continue
    elif line == '\n':
        continue

    term = line.partition('\t')[0]
    goterms_HD.append(term)

goterms_stats = []
for line in lines_stats:
    if line.__contains__('Cluster'):
        continue
    elif line == '\n':
        continue

    term = line.partition('\t')[0]
    goterms_stats.append(term)

sim = [elem for elem in goterms_HD if elem in goterms_stats]
print(sim)
print(len(sim), 'in common out of', len(goterms_HD)+len(goterms_stats) - len(sim))
