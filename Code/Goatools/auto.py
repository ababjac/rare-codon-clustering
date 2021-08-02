import os
import sys

#-------------------------------------------------------------------------------------------------------------------------------------------------------------#

# ####### ECOLI SHUFFLED ##########
# os.system("python ecoli.py ../../Files/Clusters/Ecoli/Shuf/hier_clust_ecoli_centroid_t50_KS.txt Sig_0.05/Shuf/output_ecoli_hier_KS_t50_shuf.txt 0.05")
# os.system("python ecoli.py ../../Files/Clusters/Ecoli/Shuf/hier_clust_ecoli_centroid_t100_KS.txt Sig_0.05/Shuf/output_ecoli_hier_KS_t100_shuf.txt 0.05")
# os.system("python ecoli.py ../../Files/Clusters/Ecoli/Shuf/hier_clust_ecoli_centroid_t150_KS.txt Sig_0.05/Shuf/output_ecoli_hier_KS_t150_shuf.txt 0.05")
# os.system("python ecoli.py ../../Files/Clusters/Ecoli/Shuf/hier_clust_ecoli_centroid_t50_HD.txt Sig_0.05/Shuf/output_ecoli_hier_HD_t50_shuf.txt 0.05")
# os.system("python ecoli.py ../../Files/Clusters/Ecoli/Shuf/hier_clust_ecoli_centroid_t100_HD.txt Sig_0.05/Shuf/output_ecoli_hier_HD_t100_shuf.txt 0.05")
# os.system("python ecoli.py ../../Files/Clusters/Ecoli/Shuf/hier_clust_ecoli_centroid_t150_HD.txt Sig_0.05/Shuf/output_ecoli_hier_HD_t150_shuf.txt 0.05")
#
# ####### ECOLI RANDOM ##########
# os.system("python ecoli.py ../../Files/Clusters/Ecoli/Rand/hier_clust_ecoli_centroid_t50_KS.txt Sig_0.05/Rand/output_ecoli_hier_KS_t50_rand.txt 0.05")
# os.system("python ecoli.py ../../Files/Clusters/Ecoli/Rand/hier_clust_ecoli_centroid_t100_KS.txt Sig_0.05/Rand/output_ecoli_hier_KS_t100_rand.txt 0.05")
# os.system("python ecoli.py ../../Files/Clusters/Ecoli/Rand/hier_clust_ecoli_centroid_t150_KS.txt Sig_0.05/Rand/output_ecoli_hier_KS_t150_rand.txt 0.05")
# os.system("python ecoli.py ../../Files/Clusters/Ecoli/Rand/hier_clust_ecoli_centroid_t50_HD.txt Sig_0.05/Rand/output_ecoli_hier_HD_t50_rand.txt 0.05")
# os.system("python ecoli.py ../../Files/Clusters/Ecoli/Rand/hier_clust_ecoli_centroid_t100_HD.txt Sig_0.05/Rand/output_ecoli_hier_HD_t100_rand.txt 0.05")
# os.system("python ecoli.py ../../Files/Clusters/Ecoli/Rand/hier_clust_ecoli_centroid_t150_HD.txt Sig_0.05/Rand/output_ecoli_hier_HD_t150_rand.txt 0.05")
#
# ####### ECOLI ACTUAL ##########
# os.system("python ecoli.py ../../Files/Clusters/Ecoli/KS/hier_clust_ecoli_centroid_t50_stats.txt Sig_0.05/KS/output_ecoli_hier_KS_t50_stats.txt 0.05")
# os.system("python ecoli.py ../../Files/Clusters/Ecoli/KS/hier_clust_ecoli_centroid_t100_stats.txt Sig_0.05/KS/output_ecoli_hier_KS_t100_stats.txt 0.05")
# os.system("python ecoli.py ../../Files/Clusters/Ecoli/KS/hier_clust_ecoli_centroid_t150_stats.txt Sig_0.05/KS/output_ecoli_hier_KS_t150_stats.txt 0.05")
# os.system("python ecoli.py ../../Files/Clusters/Ecoli/HD_revised/hier_clust_ecoli_centroid_t50.txt Sig_0.05/HD_revised/output_ecoli_hier_HD_t50.txt 0.05")
# os.system("python ecoli.py ../../Files/Clusters/Ecoli/HD_revised/hier_clust_ecoli_centroid_t100.txt Sig_0.05/HD_revised/output_ecoli_hier_HD_t100.txt 0.05")
# os.system("python ecoli.py ../../Files/Clusters/Ecoli/HD_revised/hier_clust_ecoli_centroid_t150.txt Sig_0.05/HD_revised/output_ecoli_hier_HD_t150.txt 0.05")
#
# ####### YEAST SHUFFLED ##########
# os.system("python yeast.py ../../Files/Clusters/Yeast/Shuf/hier_clust_yeast_centroid_t50_KS.txt Sig_0.05/Shuf/output_yeast_hier_KS_t50_shuf.txt 0.05")
# os.system("python yeast.py ../../Files/Clusters/Yeast/Shuf/hier_clust_yeast_centroid_t100_KS.txt Sig_0.05/Shuf/output_yeast_hier_KS_t100_shuf.txt 0.05")
# os.system("python yeast.py ../../Files/Clusters/Yeast/Shuf/hier_clust_yeast_centroid_t150_KS.txt Sig_0.05/Shuf/output_yeast_hier_KS_t150_shuf.txt 0.05")
# os.system("python yeast.py ../../Files/Clusters/Yeast/Shuf/hier_clust_yeast_centroid_t50_HD.txt Sig_0.05/Shuf/output_yeast_hier_HD_t50_shuf.txt 0.05")
# os.system("python yeast.py ../../Files/Clusters/Yeast/Shuf/hier_clust_yeast_centroid_t100_HD.txt Sig_0.05/Shuf/output_yeast_hier_HD_t100_shuf.txt 0.05")
# os.system("python yeast.py ../../Files/Clusters/Yeast/Shuf/hier_clust_yeast_centroid_t150_HD.txt Sig_0.05/Shuf/output_yeast_hier_HD_t150_shuf.txt 0.05")
#
# ####### YEAST RANDOM ##########
# os.system("python yeast.py ../../Files/Clusters/Yeast/Rand/hier_clust_yeast_centroid_t50_KS.txt Sig_0.05/Rand/output_yeast_hier_KS_t50_rand.txt 0.05")
# os.system("python yeast.py ../../Files/Clusters/Yeast/Rand/hier_clust_yeast_centroid_t100_KS.txt Sig_0.05/Rand/output_yeast_hier_KS_t100_rand.txt 0.05")
# os.system("python yeast.py ../../Files/Clusters/Yeast/Rand/hier_clust_yeast_centroid_t150_KS.txt Sig_0.05/Rand/output_yeast_hier_KS_t150_rand.txt 0.05")
# os.system("python yeast.py ../../Files/Clusters/Yeast/Rand/hier_clust_yeast_centroid_t50_HD.txt Sig_0.05/Rand/output_yeast_hier_HD_t50_rand.txt 0.05")
# os.system("python yeast.py ../../Files/Clusters/Yeast/Rand/hier_clust_yeast_centroid_t100_HD.txt Sig_0.05/Rand/output_yeast_hier_HD_t100_rand.txt 0.05")
# os.system("python yeast.py ../../Files/Clusters/Yeast/Rand/hier_clust_yeast_centroid_t150_HD.txt Sig_0.05/Rand/output_yeast_hier_HD_t150_rand.txt 0.05")
#
# ####### YEAST ACTUAL ##########
# os.system("python yeast.py ../../Files/Clusters/Yeast/KS/hier_clust_yeast_centroid_t50_stats.txt Sig_0.05/KS/output_yeast_hier_KS_t50_stats.txt 0.05")
# os.system("python yeast.py ../../Files/Clusters/Yeast/KS/hier_clust_yeast_centroid_t100_stats.txt Sig_0.05/KS/output_yeast_hier_KS_t100_stats.txt 0.05")
# os.system("python yeast.py ../../Files/Clusters/Yeast/KS/hier_clust_yeast_centroid_t150_stats.txt Sig_0.05/KS/output_yeast_hier_KS_t150_stats.txt 0.05")
# os.system("python yeast.py ../../Files/Clusters/Yeast/HD_revised/hier_clust_yeast_centroid_t50.txt Sig_0.05/HD_revised/output_yeast_hier_HD_t50.txt 0.05")
# os.system("python yeast.py ../../Files/Clusters/Yeast/HD_revised/hier_clust_yeast_centroid_t100.txt Sig_0.05/HD_revised/output_yeast_hier_HD_t100.txt 0.05")
# os.system("python yeast.py ../../Files/Clusters/Yeast/HD_revised/hier_clust_yeast_centroid_t150.txt Sig_0.05/HD_revised/output_yeast_hier_HD_t150.txt 0.05")

#-------------------------------------------------------------------------------------------------------------------------------------------------------------#
FOLDER = 'Full'
#FOLDER = 'Omit10'
#FOLDER = 'Omit25'

LINKAGE = 'Centroid/'
#LINKAGE = 'Single/'

DIRECTORY = '../../Files/Clusters/Hierarchical/'+LINKAGE+FOLDER+'/Cut_0.05/'
SIG = '0.001'


with os.scandir(DIRECTORY) as d:
    for entry in d:
        if entry.name.endswith('.txt') and entry.is_file():
            input_path = os.path.join(DIRECTORY, entry.name)
            output_path = FOLDER+'/Cut_0.05/Sig_'+SIG+'/GO_'+entry.name

            if entry.name.__contains__('ecoli'):
                command = 'python ecoli.py '+input_path+' '+output_path+' '+SIG
            else:
                command = 'python yeast.py '+input_path+' '+output_path+' '+SIG

            os.system(command)
