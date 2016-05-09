Phylogenetically aware variant-calling
By Rayaan Khatau and Fiona Rowan 

Prerequisites: In order for this package to run, you must have all the dependencies, as described in dependencies.txt.

Once you have installed dependencies, in the TreeToReads subdirectory, generate sample mutated genome, where each mutated genome is on a leaf of a genetic phylogeny. Do this by running the following shell command:

python treetoreads.py seqsim.cfg

You may change the specifications of these mutated genomes by edited the seqsim.cfg file, and adding your own genetic phylogeny and reference genome in the directory TreeToReads/example. 

The output includes the mutated genome fasta files (in the TreeToReads/example_out/fasta_files directory), as well as some log files, including var_site_matrix. This file lists all the mutations made by taxa (sample genome), including the position of the gene as well as to what it was mutated into. 

For a more detailed explanation on how the TreeToReads package works, you can read the README of the original package in the TreeToReads directory (we have only modified the original version). 

The user may also call variants on raw READ data in FastQ file format using our script raw_reads_preprocessing. The files and directories in that script correspond to example files we have provided, and you only need to change those directories to correspond to the paths of your input reference genome, reads, and supplementary files of known indels and known sites in order to process your own data, and then run the following shell command: 

./raw_reads_preprocessing

Once you have the simulated reads, you can run unique_seqs.py to review the mutated genes by taxa and position, by running the following shell command: 

python unique_seqs.py

Although this file is currently not functional, once it is fully implemented and has no support for .vcf files that are outputted by the GATK pipeline, you can then run the following shell command to call variants based on the given phylogeny:

python call_variants

If the application can successfully call variants of the generated mutated genome, then it will print and log its success. Modifications and further implementation is welcomed by the user. Check back frequently for updates. 
