#!/bin/bash

for a in 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30
do
    for b in 0.45 0.46 0.47 0.48 0.49 0.50 0.51 0.52 0.53 0.54 0.55 0.56 0.57 0.58 0.59 0.60 0.61 0.62 0.63 0.64 0.65 0.66 0.67 0.68 0.69 0.70 0.71 0.72 0.73 0.74 0.75 0.76 0.77 0.78 0.79 0.80 0.81 0.82 0.83
    do
        for c in 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30
        do
            ./pevk_new_mouse.sh Mus_musculus ../data/ttn_seqs/Mus_musculus_ttn.fasta $a $b $c known_mouse_exons_db
        done
    done
done