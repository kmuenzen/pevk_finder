# This code loops through all species fasta files containing
# PEVK exons and performs pairwise alignments for each species
# using megacc

#!/bin/bash

# Loop through all nucleotide files in the nt files folder
for f in *_trimmed.fasta
do
  # Python script that creates fasta file with only pairs of
  #sequences This will add MANY new files to the folder, so it
  # will be important to delete the extraneous files ASAP
  python pairwise_group.py $f

  # Perform pairwise alignments for every combined file in directory
  for q in *_combined.fasta
  do
    # Get name of combination
    pair_name=${q%_combined.fasta}
    # Run distance estimation
    megacc -a clustal_align_coding.mao -d $q -f Fasta -o $pair_name"_alignment" -s -n
  done

  # Remove all combined files
  find . -name "*_combined.fasta" -delete

  #find . -name "*_alignment.fasta" -delete

  # Loop through all Phylip alignment files and estimate divergence

  for z in *_alignment.fasta
  do
    # Get name of combination
    pair_name=${z%_alignment.fasta}
    # Run distance estimation
    megacc -a distance_estimation_overall_mean_coding.mao -d $z -o $pair_name"_divergence"
  done

  # Delete all alignments and divergence summaries
  find . -name "*_alignment.fasta" -delete
  find . -name "*_divergence_summary.txt" -delete

done
