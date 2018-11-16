# This script performs a BLAST search of a known triplicate sequence in human titin
# Human reference exons 166-174

# Kathleen Muenzen 11/9/18

for f in /Users/kmoney/Desktop/titin_project/pevk_mammals/nt_sorted/*.fasta
do
    # Extract query species name
    species_name=${f%_nt_sorted.fasta}
    species_name=${species_name#/Users/kmoney/Desktop/titin_project/pevk_mammals/nt_sorted/}
    echo $species_name

    blastn -evalue 1e-6 -query ~/Desktop/triplicate_analysis/duplication1.fasta  -db '/usr/local/ncbi/blast/db/'$species_name"_db" -out '/Users/kmoney/Desktop/triplicate_analysis/known_triplicate_blast/'$species_name"_blast_triplicate" -outfmt 6

done

for f in /Users/kmoney/Desktop/triplicate_analysis/*triplicate
do
    # Extract blast identifier
    blast_identifier=${f#/Users/kmoney/Desktop/triplicate_analysis/}
    echo $blast_identifier

    # Convert file in csv
    tr '\t' ',' < $f > '/Users/kmoney/Desktop/triplicate_analysis/known_triplicate_csv/'$blast_identifier".csv"
done


#-perc_identity 100

echo ">Homo_sapiens_concatenated" >test.fasta
grep -v ">" Homo_sapiens_nt_sorted.fasta|tr '\n' ' ' | sed -e 's/ //g' >> test.fasta