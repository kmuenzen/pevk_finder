# This code runs PEVK finder on all species found in the pevk_mammals folder

#!/bin/bash

for f in *_data/
do
  # Get species name
  species_name=${f%_data/}

  # Get names of inputs files
  frame_f1=$species_name"_frame_f1.fasta"
  frame_f2=$species_name"_frame_f2.fasta"
  frame_f3=$species_name"_frame_f3.fasta"
  ttn_seq=$species_name"_ttn.fasta"

  # Run PEVK Finder
  python pevk_finder_ntsplice_mammals.py $f/$frame_f1 $f/$frame_f2 $f/$frame_f3 $f/$ttn_seq 10 0.54 12

  mv $species_name*"_final.fasta" $f
done
