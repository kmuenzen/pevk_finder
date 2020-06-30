#!/bin/bash

species_name=$1 #species name separated by an underscore (i.e. Homo_sapiens)
nucleotide=$2 #the nucleotide sequence
frame_length=$3 #the sliding frame length
percent_pevk=$4 #the desired fraction of pevk (eg: 0.70)
min_length=$5
exon_db=$6

# Name all forward reading frames (species name in second an third positions in header of fast file)

f1="../data/frames/"$species_name"_frame_f1.fasta"
f2="../data/frames/"$species_name"_frame_f2.fasta"
f3="../data/frames/"$species_name"_frame_f3.fasta"


# Run the PEVK finder
python3 pevk_finder_ntsplice.py $f1 $f2 $f3 $nucleotide $frame_length $percent_pevk $min_length

# Define the outputs of the PEVK_finder
protf1=$species_name"_frame_f1_PEVK_protein_final.fasta"
protf2=$species_name"_frame_f2_PEVK_protein_final.fasta"
protf3=$species_name"_frame_f3_PEVK_protein_final.fasta"
nucf1=$species_name"_frame_f1_PEVK_nucleotide_final.fasta"
nucf2=$species_name"_frame_f2_PEVK_nucleotide_final.fasta"
nucf3=$species_name"_frame_f3_PEVK_nucleotide_final.fasta"

# Combine results:
cat $protf1 $protf2 $protf3 > $species_name"_combined_AA_forward.fasta"
cat $nucf1 $nucf2 $nucf3 > $species_name"_combined_NT_forward.fasta"

# Remove unnecessary files:
rm -r *nucleotide_final.fasta
rm -r *protein_final.fasta

# Name combined results:
forward_AA=$species_name"_combined_AA_forward.fasta"
forward_NT=$species_name"_combined_NT_forward.fasta"

# Run local nucleotide blasts:
blastn -query $forward_NT -db $6 -out $1"_fblast_"$3"_"$4"_"$5 -outfmt 6

# Name the blast output
blast_all=$1"_fblast_"$3"_"$4"_"$5

# Remove unnecessary files:
rm -r *forward.fasta

# Pull out 100% matches
grep "100.000" $blast_all > $1"_fblast_"$3"_"$4"_"$5"_perfect"

# Name the exact matches file
blast_exact=$1"_fblast_"$3"_"$4"_"$5"_perfect"

# Convert to csv
tr '\t' ',' < $blast_exact > $species_name"_"$3"_"$4"_"$5"_.csv"

# Name the csv files
blast_csv=$species_name"_"$3"_"$4"_"$5"_.csv"
names_csv="test_names_and_lengths_"$3"_"$4"_"$5".csv"

# Move the csv files to a separate folder
mv $blast_csv $names_csv "../data/mouse_csv/"

# Remove the remaining unnecesary files
rm -r *fblast*

