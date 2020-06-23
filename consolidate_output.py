#!/usr/bin/env python
# coding: utf-8

import os
import sys

### Consolidate output
protein = sys.argv[1]
pathway = sys.argv[2]

base_path = os.path.join(protein, pathway)
output_path = os.path.join(base_path, "output")
output_files = os.listdir(output_path)
filelist = []
for file in output_files:
    front, back = file.split("_")
    idx = front[len(protein):]
    front = front[:len(protein)]
    new_filename = front + "_" + idx.zfill(3) + "_" + back
    filelist.append(new_filename)
    os.rename(os.path.join(output_path, file), os.path.join(output_path, new_filename))

import glob
filelist = glob.glob(os.path.join(output_path, protein+"_???_output.out"))
filelist.sort()

print(filelist)

consol_output_file = os.path.join(output_path, "{0}_{1}_output.txt".format(protein, pathway))
with open(consol_output_file, 'w') as output:
    for file in filelist:
        print("Writing {0}".format(file.split("_")[1]))
        output.writelines([line for line in open(file).readlines()])