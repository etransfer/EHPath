#!/usr/bin/env python
# coding: utf-8

import os
import sys

### Create pdb_list
protein = sys.argv[1]
pathway = sys.argv[2]

base_path = os.path.join(protein, pathway)
input_file_path = os.path.join(base_path, "input")
input_files = os.listdir(input_file_path)

input_paths = []
for file in input_files:
    fragment = file.split(".")[0].split("_")[-1]
    input_paths.append(fragment)
unique_input_paths = list(set(input_paths))
unique_input_paths.sort()

print(unique_input_paths)

pdb_list_path = os.path.join(base_path, 'pdb_list')
with open(pdb_list_path, 'w') as file:
    for path in unique_input_paths:
        if len(path) > 0:
            file.write(path + '\n')