# This script performs a BLAST search of a known triplicate sequence in human titin
# Human reference exons 166-174

# Kathleen Muenzen 11/9/18

for f in /Users/kmoney/Desktop/paper_submission_materials/data/exon_libraries/manual_crop_exon_libraries/*.fasta
do
    # Extract query species name
    species_name=${f%_nt_sorted_trimmed.fasta}

    blastn -evalue 1e-6 -query $f -db $f"_db" -out $species_name"_and_"$s"_blast" -outfmt 6

done