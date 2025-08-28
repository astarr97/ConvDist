Very important note!!!
We need to do prediction on both the input sequence and the reverse complement and then take the average as they do in the typical ChromBPNet pipeline.

The goal of this pipeline is to create sequence inputs to ChromBPNet to be used in the nearest site approach to detecting convergent evolution.

The general strategy is to:
For each convergent/divergent site, predict the accessibility profile and counts of the focal and related allele in the genome context of the LCA of both pairs of focal/related species.
For each of the nearest sites, predict the accessibility profile and counts of the focal and related allele in the genome context of the pair of focal/related species from which that nearest site is drawn.
For example, when comparing orca and otter, you would predict the effect of the Orcinus_orca and Bos_taurus bases in the genome context of the LCA of Orcinus_orca and Bos_taurus.
You would predict the effect of the Enhydra_lutris and Mustela_putorius alleles in the sequence context of their LCA.

Once you have all the predictions in hand, you can compare the average prediction of the convergent/divergent base to the average of the two nearest sites (if both nearest sites are usable)
If only one nearest site is usable, then you can compare the prediction to the average prediction of the convergent/divergent, or only the prediction of the species for which the nearest site is usable.

To set this pipeline up, you need to have a configuration file e.g. Config_GetDL_Aquatic.csv.  It is a comma-delimited file.
There should be a number of rows equal to the number of pairwise comparisons between branches.
For example, if there are 4 branches then there should be 4C2 = 6 rows.  
Each row of the file should consist of two sets (one for each branch in the pairwise comparison) of 4 node/species names. 
The first node in the list should be the last common ancestor of the focal and related species.
The next three should be the focal, related, and outgroups species.
For example, for orca compared to otter, it would be fullTreeAnc191,Orcinus_orca,Bos_taurus,Sus_scrofa,fullTreeAnc211,Enhydra_lutris,Mustela_putorius,Mellivora_capensis

There is a helper script to create the config file from a file that contains each of the comma-delimited sets of 4 node/species names called make_config.py

With this config in hand, you can run, consecutively:
python <config_file.csv> <reference_species> 0
python <config_file.csv> <reference_species> 1

The first argument is the config file, the second is the species to which all the ClosestVar.bed.gz files are referenced, and the third is whether to generate scripts for focal or rel.
This will create three scripts numbered 1, 2, 3 and a single script, shared between rel and foc, numbered 0.
You need to run script 0 first to get the fasta files and sizes of the contigs for later steps

You can submit the first script as a batch job, it extracts the relevant information from the ClosestVar.bed.gz files and then lifts over to the appropriate LCA.
The second script contains code to create jobs to get new TSV files of SNPs similar to the pipeline for PhyloP, but now including the reference species.  
You have to run those in the terminal to submit a bunch of small jobs that should finish in a few hours.
You can run python check_done.py when they all finish to check that there were no errors and all needed files are present for the third script.

The third script uses the TSV files to switch to the bases in those files, expand the single base to a 2114 base region, and pull the sequence from the relevant LCA fasta file.
Finally, it uses the LCA base from the TSV file to determine if the base needs to be reverse complimented.

Lastly, you can use make_rerun.py by running python check_error.py | sort -u, copying all the lines with errors into a file called "to_rerun.txt", and then running make_rerun.py
This will create a shell script that you can run that will resubmit all the errored phase 2 jobs.
