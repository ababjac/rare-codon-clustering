import os
import sys

DIR = sys.argv[1]

def count(l):
    total = 0
    BP = 0
    MF = 0
    CC = 0
    BP_flag = False
    MF_flag = False
    CC_flag = False

    for line in l:
        if line.__contains__('Cluster'):
            total += 1
            BP_flag = False
            MF_flag = False
            CC_flag = False
            continue

        if line.__contains__('BP') and not BP_flag:
            BP += 1
            BP_flag = True
            continue

        if line.__contains__('MF') and not MF_flag:
            MF += 1
            MF_flag = True
            continue

        if line.__contains__('CC') and not CC_flag:
            CC += 1
            CC_flag = True
            continue

    return BP, MF, CC, total

with os.scandir(DIR) as d:
    for entry in d:
        if entry.name.endswith('.txt') and entry.is_file():
            file = open(DIR+entry.name, 'r')
            l = file.readlines()

            BP, MF, CC, total = count(l)

            print(entry.name, 'BP =', BP, 'MF =', MF, 'CC =', CC, 'total =', total)
