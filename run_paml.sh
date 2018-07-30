#!/bin/bash

# For every exon folder in paml_outputs:
for exon_folder in frame*/
do
  echo $exon_folder
  # Get exon name
  exon_name=${exon_folder%/}

  # Specify path to alignment file
  alignment_file=$exon_folder$exon_name"_seqfile.txt"
  echo $alignment_file

  # Specify path to original tree file
  tree_file=$exon_folder$exon_name".tre"
  echo $tree_file

  # Specify path to M0 outfile
  m0_out=$exon_folder$exon_name"_m0.mlc"
  echo $m0_out

  # Specify path to estimated tree
  estimated_tree_name=$exon_folder$exon_name"_estimated_tree.tre"
  echo $estimated_tree_name

  # Specify path to output of sites models
  final_out=$exon_folder"site/"$exon_name"_site.mlc"
  echo $final_out

  # Run M0 and all sites models
  python paml_codeml_M0_and_sites.py $alignment_file $tree_file $m0_out $estimated_tree_name $final_out

done
