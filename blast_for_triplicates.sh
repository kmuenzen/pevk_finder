# This script performs a BLAST search for all possible duplicate and triplicate exon sequences in all 41 mammal species
# Human reference exons 166-174

# Kathleen Muenzen 11/9/18

for f in /Users/kmoney/Desktop/titin_project/pevk_mammals/nt_sorted/*.fasta
do
    # Extract query species name
    species_name=${f%_nt_sorted.fasta}
    species_name=${species_name#/Users/kmoney/Desktop/titin_project/pevk_mammals/nt_sorted/}
    echo $species_name

    blastn -query $f -db '/usr/local/ncbi/blast/db/'$species_name"_db" -out '/Users/kmoney/Desktop/triplicate_analysis/species_triplicate_blast/'$species_name"_and_"$species_name"_blast" -outfmt 6

done

for f in /Users/kmoney/Desktop/triplicate_analysis/species_triplicate_blast/*blast
do
    # Extract blast identifier
    blast_identifier=${f%_blast}
    blast_identifier=${f#/Users/kmoney/Desktop/triplicate_analysis/species_triplicate_blast/}

    # Convert file in csv
    tr '\t' ',' < $f > '/Users/kmoney/Desktop/triplicate_analysis/species_triplicate_csv/'$blast_identifier".csv"
done


#-perc_identity 100