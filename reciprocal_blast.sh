#!/bin/bash

# Loop through all fasta files in nt_sorted folder
for f in /Users/kmoney/Desktop/paper_submission_materials/data/exon_libraries/manual_crop_exon_libraries/*.fasta
do
    # Extract query species name
    species_name=${f%_nt_sorted_trimmed.fasta}

    # Loop through all species and blast 
    for s in Phascolarctos_cinereus Homo_sapiens Ictidomys_tridecemlineatus Loxodonta_africana Manis_javanica Mus_musculus Mus_pahari Myotis_brandtii Ochotona_princeps Orcinus_orca Oryctolagus_cuniculus Pteropus_vampyrus Rattus_norvegicus Rhinolophus_sinicus Sorex_araneus Sus_scrofa Trichechus_manatus Tupaia_chinensis Tursiops_truncatus Vicugna_pacos Acinonyx_jubatus Bos_taurus Ceratotherium_simum Condylura_cristata Dasypus_novemcinctus Eptesicus_fuscus Equus_caballus Erinaceus_europaeus Felis_catus Galeopterus_variegatus Hipposideros_armiger
    do
        blastn -evalue 1e-6 -query $f -db $s"_db" -out $species_name"_and_"$s"_blast" -outfmt 6
    done
done

# Translate all BLAST inputs to CSV files so they can be read by R
for f in *blast
do
    # Extract blast identifier
    blast_identifier=${f%_blast}

    # Convert file in csv
    tr '\t' ',' < $f > $blast_identifier".csv"
done

# Move blast files into separate folder
mv *blast nt_matches

# Move csv files into separate folder
mv *csv nt_matches_csv





