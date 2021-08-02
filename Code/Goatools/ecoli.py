from goatools import obo_parser
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from Bio.UniProt.GOA import gafiterator
import gzip
import pandas as pd
from goatools.godag_plot import plot_gos, plot_results, plot_goid2goobj
import sys
import statistics
from scipy.stats import iqr

FILE_IN = sys.argv[1]
FILE_OUT = sys.argv[2]
PVAL_CUT = float(sys.argv[3])

#initialze U_ID -> Gene
df = pd.read_csv('../../Data/Ecoli/uniprot_mapping.csv')
umap = df.groupby('Gene')['Uniprot_ID'].apply(list).to_dict()

go = obo_parser.GODag('../../Data/GO/go-basic.obo')

filename = '../../Data/GO/ecocyc.gaf.gz'
with gzip.open(filename, 'rt') as fp:
    d = {}

    for annotation in gafiterator(fp):
        uniprot_ID = annotation.pop('DB_Object_ID')
        d[uniprot_ID] = annotation

population = d.keys()

association = {}
for elem in d:
    if elem not in association:
        association[elem] = set()
    association[elem].add(str(d[elem]['GO_ID']))

methods = ['fdr']

g = GOEnrichmentStudy(population, association, go, methods=methods)


results = []
results_sig = []
#run on each cluster
file = open(FILE_IN)

for line in file.readlines():
    if line.__contains__('Singleton'):
        break

    if line.__contains__('Cluster'):
        continue

    genes = line.split(', ')
    #genes.remove('\n')

    #get g_res for each cluster
    map = [umap[gene][0] for gene in genes if gene in umap.keys()]
    #print(map)
    h = dict.fromkeys(map)
    study = h.keys()

    g_res = g.run_study(study)
    g_res_sig = [r for r in g_res if r.p_fdr < PVAL_CUT]

    results.append(g_res)
    results_sig.append(g_res_sig)

# file = open('pvals.txt', 'a')
# pvals = [g.p_fdr for r in results_sig for g in r]
# file.write(FILE_IN+': '+str(statistics.median(pvals))+' '+str(iqr(pvals)))

file = open(FILE_OUT, 'w')
for result, i in zip(results_sig, range(1, len(results_sig)+1)):
    file.write('Cluster '+str(i)+' - \n')
    for obj in result:
        file.write(str(obj))
        file.write('\n')
    file.write('\n')

# counter = 0
# for result in results_sig:
#     counter+=1
#     plot_results("cluster"+str(counter)+"_0.15_{NS}.png", result)
