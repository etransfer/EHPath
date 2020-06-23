#!/bin/bash
#SBATCH -p scavenger,common,dbchem
#SBATCH --error=ehpath.err
#SBATCH --mem=32G
#SBATCH --output=ehpath.out
#SBATCH --job-name=ehpath_surf

pdb_list=$1
cutoff_num=$2
total_paths=$3
type_da=$4
alpha_reorg=$5
dir_path=$6
module load Anaconda3/3-2019
python EHPath_multirun.py $pdb_list $cutoff_num $total_paths $type_da $alpha_reorg $dir_path
