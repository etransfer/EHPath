#!/usr/bin/env python

import sys
import os
import errno
import subprocess
import multiprocessing

from datetime import datetime

print("System arguments are: {}\n".format(" ".join(sys.argv)))

cutoff_num=int(sys.argv[2])
total_paths=int(sys.argv[3])
type_da=sys.argv[4]
alpha_reorg=int(sys.argv[5])
dir_path=sys.argv[6]
pdb_list_file=os.path.join(dir_path,sys.argv[1])

pdbs = []
with open(pdb_list_file) as file:
	for line in file:
	   pdbs.append(line.rstrip())

print(' ')
print('Start calculations........')

import pandas as pd
import networkx as nx

input_folder = 'input'
output_folder = 'output'

now = datetime.now()
timestamp = now.strftime("%Y%m%d_%H%M%S")

output_logfile = "{0}_{1}.log".format(timestamp,type_da)
output_logpath = os.path.join(dir_path, "log", output_logfile)

if not os.path.exists(os.path.dirname(output_logpath)):
    try:
        os.makedirs(os.path.dirname(output_logpath))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise

with open(output_logpath, 'w') as f:
	f.write("System arguments are: {}\n".format(" ".join(sys.argv)))
	f.write("\n")

no_bridge=[]
no_donor=[]

def do_single_run( pdb, cutoff_num, total_paths, type_da, alpha_reorg, dir_path, num ):
	if type_da not in ('hole', 'electron'):
		print("Invalid type_da for pdb: {0} at iter: {1}".format( pdb, num))
	else:
		com = 'python EHPath_single_{3}.py {0} {1} {2} {3} {4} {5} {6}'.format( pdb, cutoff_num, total_paths, type_da, alpha_reorg, dir_path, num )
		log = subprocess.getoutput(com)
		print(log)
		with open(output_logpath, 'a') as f:
			f.writelines(log)
			f.write("\n")

if __name__ == "__main__": 
    cpu_count = os.cpu_count() 
    pool = multiprocessing.Pool( cpu_count )

    print("Running on {} processes\n".format(cpu_count))

    for k,pdb in enumerate(pdbs):
        pool.apply_async( do_single_run, args=( pdb, cutoff_num, total_paths, type_da, alpha_reorg, dir_path, k+1 ) )
    pool.close()
    pool.join()

print(' ')
print('*****************************************')
print('This version of EHPath was modified first by Hector L. Torres-Vera, and later Xiaochen Du, for large EHPath calculation runs.')
print('If you have any questions or suggestions about EHPath, please contact Dr. Ruijie Teo (rt131@duke.edu) or Dr. Ruobing Wang (ruobing.wang@operasolutions.com).')
print('Citation: Teo, R. D.; Wang, R.; Smithwick, E.; Migliore, A.; Therien, M. J.; Beratan, D. N. Proc. Natl. Acad. Sci. U.S.A., 2019, 116, 15811-15816.')
print('*****************************************')
print(' ')
print('Copyright (C) 2019  Teo, R. D.; Wang, R.; Smithwick, E.; Migliore, A.; Therien, M. J.; Beratan, D. N.')
print(' ')
print('This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.')
print(' ')
print('This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.')
print(' ')
print('You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/gpl.html.')
print(' ')
print('The following PDBs had no bridges and EHPath was not conducted on them:')
print(no_bridge)
print('The following PDBs had no donors and EHPath was not conducted on them:')
print(no_donor)
