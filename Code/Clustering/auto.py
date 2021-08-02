import os
#-------------------------------------------------------------------------------------------------------------------------------------------------------------#

#TEST#
#os.system('python dendrogram.py ../Files/MM/Full/Indexed/HD_ecoli.csv Gene ../Files/Clusters/Hierarchical/Full/test.txt ../Images/Clusters/Dendrograms/Full/CH_index/test.png ../Images/Clusters/Dendrograms/Full/Plots/test.png')

#-------------------------------------------------------------------------------------------------------------------------------------------------------------#
FOLDER = 'Full'
OPT = '/Indexed/'
#FOLDER = 'Omit10'
#FOLDER = 'Omit25'
#OPT = '/'

LINKAGE = 'Centroid'
#LINKAGE = 'Single'

IN_DIRECTORY = '../Files/MM/'+FOLDER+OPT
IMG_DIRECTORY = '../Images/Clusters/Dendrograms/'+LINKAGE+'/'+FOLDER+'/Cut_0.05/'
#IMG_DIRECTORY2 = '../Images/Clusters/Clustermaps/'+LINKAGE+'/'+FOLDER+'/'
OUT_DIRECTORY = '../Files/Clusters/Hierarchical/'+LINKAGE+'/'+FOLDER+'/Cut_0.05/'

#based on third pass
# THRESHOLDS = {
#     'Full' : {'HD_ecoli' : 0.28, 'HD_yeast' : 0.25, 'KS_ecoli' : 0.16, 'KS_yeast' : 0.15},
#     'Omit10' : {'HD_ecoli' : 0.26, 'HD_yeast' : 0.26, 'KS_ecoli' : 0.16, 'KS_yeast' : 0.19},
#     'Omit25' : {'HD_ecoli' : 0.27, 'HD_yeast' : 0.27, 'KS_ecoli' : 0.27, 'KS_yeast' : 0.27}
# }

#based on fourth pass
# THRESHOLDS = {
#     'Full' : {'HD_ecoli' : 0.05, 'HD_yeast' : 0.05, 'KS_ecoli' : 0.14, 'KS_yeast' : 0.12},
#     'Omit10' : {'HD_ecoli' : 0.05, 'HD_yeast' : 0.05, 'KS_ecoli' : 0.14, 'KS_yeast' : 0.12},
#     'Omit25' : {'HD_ecoli' : 0.05, 'HD_yeast' : 0.05, 'KS_ecoli' : 0.13, 'KS_yeast' : 0.14}
# }


with os.scandir(IN_DIRECTORY) as d:

    for entry in d:
        if entry.name.endswith('.csv') and entry.is_file():
            input_path = os.path.join(IN_DIRECTORY, entry.name)

            filename = entry.name.split('.')[0]
            output_path = OUT_DIRECTORY+'hclust_'+filename

            if entry.name.__contains__('ecoli'):
                colname = 'Gene'

                # if entry.name.__contains__('KS'):
                #     threshold = THRESHOLDS[FOLDER]['KS_ecoli']
                # else:
                #     threshold = THRESHOLDS[FOLDER]['HD_ecoli']
            else:
                colname = 'locus_tag'

                # if entry.name.__contains__('KS'):
                #     threshold = THRESHOLDS[FOLDER]['KS_yeast']
                # else:
                #     threshold = THRESHOLDS[FOLDER]['HD_yeast']

            errimg_path = IMG_DIRECTORY+'CH_index/ch_'+filename+'.png'
            dendimg_path = IMG_DIRECTORY+'Plots/dend_'+filename+'.png'
            #clustimg_path = IMG_DIRECTORY2+'clustermap_'+filename+'.png'


            command = 'python dendrogram.py '+input_path+' '+colname+' '+output_path+' '+errimg_path+' '+dendimg_path+' '+LINKAGE+' '+str(0.05)
            #command = 'python make_clustermap.py '+input_path+' '+colname+' '+clustimg_path+' '+LINKAGE
            #print(command)
            os.system(command)
