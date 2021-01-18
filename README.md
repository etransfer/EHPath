# EHPath
EHPath ('Electron-Hole Pathways') is a python module that maps and ranks electron-hole pathways in proteins and nucleic acids according to the mean residence time. Any questions or suggestions about EHPath.py can be directed to Dr. Ruijie Teo (ruijiedarius.teo@ucalgary.ca) or Dr. Ruobing Wang (ruobing.wang@operasolutions.com). 

Citation: Teo, R. D.; Wang, R.; Smithwick, E.; Migliore, A.; Therien, M. J.; Beratan, D. N. Proc. Natl. Acad. Sci. U.S.A. 2019, 116, 15811-15816. Link: https://doi.org/10.1073/pnas.1906394116

## EHPath Multirun
We also created the EHPath Multirun package (EHPath_multirun.py), which is an extension of the EHPath program which allows a user to analyze such pathways with a large number of input files.

A run of EHPath_multirun.py requires a command like the following:

    python path/to/EHPath_multirun.py {pdb_list} {cutoff_num} {total_paths} {hole or electron} {alpha_reorg} {path/to/main_dir}

An example run would be:

    python EHPath_multirun.py pdb_list 4 1 hole 1 1rrq3/HT

The pdb_list is located in the main_dir and each row contains a protein pathway to be analyzed. For each pdb in the pdb_list, EHPath_multirun.py will call either EHPath_single_electron.py or EHPath_single_hole.py depending on user input to run the original EHPath.

The input files should be located in an input folder within the main_dir. The donor node csv file should be named ‘Donor_Nodes_’ + pdb + ‘.csv’, the bridge nodes file should be labelled ‘Bridging_Nodes_’ + pdb + ‘.csv’ and the acceptor node file ‘Acceptor_Nodes_’ + pdb + ‘.csv’. After a successful run, the output from each pdb will be individually saved in the output folder within the main_dir. 

To make it easier for users, a get_input_files.py script can be used to create the pdb_list. In the folder containing the EHPath_multirun.py script, use the bash command:

    python get_input_files.py {protein} {pathway (HT or ET)}

An example is:

    python get_input_files.py 1rrq3 ET

To consolidate all the output files, we can use the consolidate_output.py script:

    python consolidate_output.py {protein} {pathway (HT or ET)}

An example is:

    python consolidate_output.py 1rrq3 HT

EHPath Multirun is also parallelized to enable speed-up when run on an HPC cluster. A Conda environment file (environment.yml) is also supplied for compatibility. Please see (https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file) for installation details.

Citation: Teo, R. D.; Du, X.; Torres Vera, H. L.; Migliore, A.; Beratan, D. N. J. Phys. Chem. B 2021, 125, 1, 17–23. Link: https://doi/10.1021/acs.jpcb.0c08598
