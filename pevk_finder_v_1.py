
"""
The pevk_finder.py algorithm is designed to annotate the exon-intron structure of 
the PEVK segment of titin, given any given mammalian titin nucleotide sequence.

Inputs:
    1. The complete titin nucleotide sequence for the current species (Fasta file)
    2. The length of the sliding window (int)
    3. The minimum PEVK ratio (float)
    4. The minimum PEVK exon length (int)

Outputs:
    1. [species_name]_PEVK_exons_unbounded_AA.fasta (Fasta file): ALL predicted PEVK exons as
    amino acid sequences, sorted by titin location. Coordinates in sequence descriptions
    are relative to the corresponding reading frame.
    2. [species_name]_PEVK_exons_unbounded_NT.fasta (Fasta file): ALL predicted PEVK exons as
    nucleotide sequences, sorted by titin location. Coordinates in sequence descriptions
    are relative to full titin DNA sequence.
    3. [species_name]_PEVK_exons_bounded_AA.fasta (Fasta file): IQR +- (1.5 x IQR) predicted PEVK exons as
    amino acid sequences, sorted by titin location. Coordinates in sequence descriptions
    are relative to the corresponding reading frame.
    4. [species_name]_PEVK_exons_bounded_NT.fasta (Fasta file): IQR +- (1.5 x IQR) predicted PEVK exons as
    nucleotide sequences, sorted by titin location. Coordinates in sequence descriptions
    are relative to full titin DNA sequence.

Command line instructions (with suggested parameter settings):

1. Navigate to the directory (using cd) where this script and the titin DNA sequence fasta file are saved
2a. To run PEVK_finder on a single DNA sequence, use the following command:
    
    python -W ignore pevk_finder_v_1.py -i [species_name]_ttn.fasta -w window_length -r minimum_pevk_ratio -l minimum_exon_length

    ex: python -W ignore pevk_finder_v_1.py -i ./ttn_seqs/Homo_sapiens_ttn.fasta -w 10 -r 0.54 -l 12

2b. To run PEVK_finder on multiple DNA sequences at once, use the following command:
    
    python -W ignore pevk_finder_v_1.py -i [directory_with_tnn_seqs] -w window_length -r minimum_pevk_ratio -l minimum_exon_length

    ex: python -W ignore pevk_finder_v_1.py -i ./ttn_seqs/ -w 10 -r 0.54 -l 12

Exon libraries will be deposited in the ./data/ directory by default, and separated into bounded/unbounded and translated/untranslated directories.

Fasta files can be viewed with any text editor.

"""


# Import libraries
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_protein
import statistics
import csv
import sys
import argparse
import distutils.util
import numpy as np
import os


def pevk_finder(nt_seq, length, in_ratio, min_length, outpath, optional_outputs):

    # Create nucleotide sequence objects
    ref_seq = SeqIO.read(nt_seq,"fasta") # read the titin nucleotide sequence, will be stored as SeqRecord
    nucleotide_seq = Seq(str(ref_seq.seq)+"*", IUPAC.unambiguous_dna) # create a Seq object for the titin nucleotide sequence


    # Save the description and species name as strings to use later
    description = ref_seq.description
    split = ref_seq.description.split()
    species_name = split[1] + "_" + split[2] + "_"


    # Create Seq objects for each of three forward reading frames
    frame_f1 = Seq(str(ref_seq.seq), generic_dna)
    frame_f2 = Seq(str(ref_seq.seq[1:]), generic_dna)
    frame_f3 = Seq(str(ref_seq.seq[2:]), generic_dna)


    # Translate each reading frame and create SeqRecord objects
    prot_f1 = frame_f1.translate()
    prot_f1_record = SeqRecord(prot_f1, description = description, id="frame_f1")

    prot_f2 = frame_f2.translate()
    prot_f2_record = SeqRecord(prot_f2, description = description, id="frame_f2")

    prot_f3 = frame_f3.translate()
    prot_f3_record = SeqRecord(prot_f3, description = description, id="frame_f3")

    # Compile frames into a list
    frames = [prot_f1_record, prot_f2_record, prot_f3_record]


    # Initialize the lists for the final csv file:
    names = []
    lengths = []
    percentages = []
    starts = []
    ends = []
    frame_names = []
    prot1_regions = []
    prot2_regions = []
    prot3_regions = []
    nuc1_regions = []
    nuc2_regions = []
    nuc3_regions = []


    ################## Loop through each reading frame ###################
    for f in frames:
        
        # Create Seq object for current frame
        prot_seq = Seq(str(f.seq)+"*", IUPAC.protein) #make a seq object


        # Get information on current frame
        identification = f.id # which frame? ex: "frame_f3"
        description = f.description # full object description. ex: 'NC_000002.12:c178807423-178525989 Homo sapiens chromosome 2, GRCh38.p7 Primary Assembly'
        split = f.description.split(" ") # split description to extract species name
        species_name = split[1] + "_" + split[2] + "_" # species name. ex: 'Homo_sapiens_'

        
        # Define paramater inputs
        input_frame_length = int(length) # read in sliding window length
        PEVK_ratio = float(in_ratio) # read in minimum PEVK ratio


        ##### Define Helper Funtions #####
        """
        Helper function pevk_counter accepts a string of amino acids (frame) and returns the int number of
        P, E, V and K (and A) in the string.
        """

        def pevk_counter(frame):
            counter = 0
            for a in frame:
                if a == "P" or a == "E" or a == "V" or a == "K" or a == "A":
                    counter += 1
            return counter


        """
        Helper function frame_define accepts a string of amino acids (seq) and counts the number of steps 
        it takes to reach the first "*" in teh string
        """

        def frame_define(seq):
            counter = 0
            for i in seq:
                if i != "*":
                    counter += 1
                if i == "*":
                    return counter
                    break


        """
        Helper function find_pevk accepts an amino acid sequence without "*"s (seq), the start location of the seq
        within the frame (seq_start) and the end location of the seq within the frame (seq_end)
        """

        def find_pevk(seq, seq_start, seq_end):
            # Find the start and end coordinates of the corresponding nucleotide sequence.
            if identification == "frame_f1":
                nuc_ref_start = seq_start*3
                nuc_ref_end = seq_end*3
            if identification == "frame_f2":
                nuc_ref_start = seq_start*3 + 1
                nuc_ref_end = seq_end*3 + 1
            if identification == "frame_f3":
                nuc_ref_start = seq_start*3 + 2
                nuc_ref_end = seq_end*3 + 2
            

            # Define the nucleotide sequence that correponds with the coordinates found above
            nuc_ref_seq = str(nucleotide_seq[nuc_ref_start:nuc_ref_end])


            # Initialize the list of lists that will hold the results in the following format: [spliced_seq, absolute_seq_start, absolute_seq_end]
            result = []

            ### CASE 1: The length of the sliding frame length is either less than or equal to the seq length
            if input_frame_length > len(seq) or input_frame_length == len(seq):
                result = "x" # result becomes the string "x", which will alow the code to ignore this sequence and move on to the next one
                return result # return the result


            ### CASE 2: The length of the input seq is greater than the sliding frame length
            if len(seq) > input_frame_length:
                primary_positions = [] # list that will hold the start positions of each each PEVK frame in the forward direction
                
                # Fill the primary_positions list the start locations of each frame with PEVK ratio greater than the ratio parameter
                for i in range(0, len(seq) - input_frame_length):
                    frame = seq[i:i+input_frame_length]
                    counts = pevk_counter(frame)
                    ratio = float(counts)/float(input_frame_length)
                    if ratio >= PEVK_ratio:
                        primary_positions.append(i)

                ### CASE 2a: There are no frames with high PEVK content in the sequence
                if primary_positions == []:
                    result = "x" # The result is nullified, and the code will move on the next sequence
                    return result

                ### CASE 2b: There is at least one frame with a high PEVK ratio
                else:
                    # Create a list that contains distinct groupings of high PEVK ratio groupings
                    group_positions = [] # the start points of each grouping
                    groupings = [] # the list that will contain the actual groups
                    temp_min = primary_positions[0]
                    group_positions.append(temp_min)
                    for i in range(1, len(primary_positions)-1):
                        if primary_positions[i] - temp_min <= input_frame_length:
                            temp_min = temp_min
                        else:
                            temp_min = primary_positions[i]
                            group_positions.append(temp_min)
                    
                    # If there is only one starting position:
                    if len(group_positions) == 1:
                        groupings = [primary_positions]
                    
                    # If there are multiple starting positions:
                    else:
                        # Define the groups:
                        for i in range(0,len(group_positions)):
                            if i < len(group_positions) - 1:
                                current_position = group_positions[i]
                                next_position = group_positions[i+1]
                                current_index = primary_positions.index(current_position)
                                next_index = primary_positions.index(next_position)
                                groupings.append(primary_positions[current_index:next_index])
                            
                            # Once we have arrived at the last item in group_positions:
                            else:
                                current_position = group_positions[i]
                                current_index = primary_positions.index(current_position)
                                groupings.append(primary_positions[current_index:])
                    

                    # Go through each sublist in groupings and find donor/acceptor splice sites:
                    for j in range(0, len(groupings)):
                        # Initialize list that will hold candidate donor sites (end of exon)
                        donor_sites = []
                        
                        # Initialize list that will hold candidate acceptor sites (start of exon)
                        acceptor_sites = []
                        
                        # Initialize the list that will hold difference in location between possible donor start sites and the PEVK region start site
                        start_differences = []
                        
                        # Start with the list of primary positions. Find the median of that list.
                        AA_starting_point = int(statistics.median(groupings[j]))
                        
                        # Where is this in the local nucleotide sequence?
                        NT_starting_point = AA_starting_point*3 # This is the first letter of the codon encoding the median AA in the high PEVK region
                        
                        # First go backwards to look for the ACCEPTOR splice site
                        for i in range(NT_starting_point,0,-1):
                            # The rightmost letter (first letter of a codon) should be a G
                            if i%3 == 0 and nuc_ref_seq[i] == "G":
                                
                                # The leftmost letter (last letter of a codon) should be an A
                                if nuc_ref_seq[i-1] == "A":
                                    
                                    # If both of these requirements are met, add the location of G to acceptor_sites
                                    acceptor_sites.append(i)
                        
                        # Is the acceptor site list empty? If so, this is not a real PEVK sequence!
                        if acceptor_sites == []:
                            result.append("x")
                            continue
                        
                        # If there is an acceptor site:
                        else:
                            # We now want to know which acceptor site is closest to the first AA of primary positions. First define the starting position:
                            PEVK_start = NT_starting_point
                            
                            # Then find the closest acceptor site to the start of the PEVK region:
                            for x in acceptor_sites:
                                diff = abs(PEVK_start - x)
                                start_differences.append(diff)
                            
                            # The best guess G of the start splice site is the min of the differences
                            acceptor_g_location = acceptor_sites[start_differences.index(min(start_differences))]
                            
                            # Now go forward to look for the DONOR splice site
                            for i in range(NT_starting_point, len(nuc_ref_seq)-3):
                                # The leftmost letter (third letter of a codon) should be an A
                                if i%3 == 2 and nuc_ref_seq[i] == "A":
                                    
                                    # The second letter (first letter of a codon) should be a G
                                    if (i+1)%3 == 0 and nuc_ref_seq[i+1] == "G":
                                        
                                        # The third letter (second letter of a codon) should be a G
                                        if (i+2)%3 == 1 and nuc_ref_seq[i+2] == "G":
                                            
                                            # The fourth letter (third letter of a codon) should be a T
                                            if (i+3)%3 == 2 and nuc_ref_seq[i+3] == "T":
                                                
                                                # If all of these requirements are met, add the location of the first G to donor sites:
                                                donor_sites.append(i+1)
                            

                            # If there are no donor sites:
                            if donor_sites == []:
                                result.append("x")
                                continue
                            

                            # If there is at least one donor site
                            if donor_sites != []:

                                # If there is only one possible donor site:
                                if len(donor_sites) == 1:
                                    
                                    # This becomes the acceptor location
                                    donor_g1_location = donor_sites[0]
                                
                                # If there is more than one possible acceptor site:
                                if len(donor_sites) > 1:
                                    end_differences = []
                                    PEVK_end = (groupings[j][-1]+input_frame_length)*3
                                    for x in donor_sites:
                                        diff = abs(PEVK_end - x)
                                        end_differences.append(diff)
                                    
                                    # The best guess G of the start splice site is the min of differences
                                    donor_g1_location = donor_sites[end_differences.index(min(end_differences))]
                                
                                # Now extract the region between the donor and acceptor sites
                                final_seq =  str(Seq(nuc_ref_seq[acceptor_g_location+3:donor_g1_location]).translate())

                                # Define the relative start and end positions (within the local amino acid sequence) of the final sequence
                                relative_seq_start = (acceptor_g_location/3)+1
                                relative_seq_end = donor_g1_location/3
                                
                                # Define the relative start and end positions of the sequence within the whole amino acid sequence
                                absolute_seq_start = seq_start + relative_seq_start
                                absolute_seq_end = seq_start + relative_seq_end

                                # Define the full sequence with no splice sites
                                final_seq = seq[relative_seq_start:relative_seq_end]
                                # Put the sequence and the start and end locations in the result list
                                result.append((final_seq, absolute_seq_start, absolute_seq_end))

                # Woohoo!
                return result

        ################################## Main Code ##############################
        # Loop through the entire protein sequence, extracting each region between stop codons and running it through find_pevk
        
        # Create empty lists for storing final PEVK exons
        #prot_regions = []
        #nuc_regions = []
        
        # Loop through entire frame
        for i in range(len(prot_seq)-1):
            if prot_seq[i] == "*": # when we get to a stop codon
                
                seq_start = i+1 # get the position of the first AA in the exon
                # If the next amino acid is a "*", keep moving through the sequence
                if prot_seq[seq_start] == "*":
                    continue
                
                # if the next amino acid is an actual letter:
                else:
                    seq_length = frame_define(str(prot_seq)[seq_start:]) # find the length of the seq between 2 stop codons using frame_define
                    seq_end = seq_start + seq_length # get the position of the last AA in the exon (don't include stop codon)
                    seq = str(prot_seq[seq_start:seq_end]) # get the actual AA sequence of the frame


                    ### Run the sequence through the find_pevk helper function to find PEVK exons!
                    pevk_seqs = find_pevk(seq, seq_start, seq_end)



                    # At this point, we have a result in the form of [(final_seq, absolute_seq_start, absolute_seq_end)], or a result of "x". These are compiled into the result list
                    # First check to see if the result is null:
                    if pevk_seqs != "x" and pevk_seqs != ["x"]:
                        # Initialize result_final, which will hold the contiguous and high % PEVK results
                        result_final = []
                        
                        ### Case 1. If there is only one item in result:
                        if len(pevk_seqs) == 1:
                            current_seq = pevk_seqs[0][0]
                            if current_seq != "x":
                                # Use the pevk_counter helper function to determine % PEVK
                                final_counts = pevk_counter(current_seq)
                                
                                # Shorter exons will be subject to the input minimum % PEVK
                                # Longer exons will only be required to pass a 50% PEVK threshold
                                if len(current_seq) > 0:
                                    final_ratio = float(final_counts)/float(len(current_seq))
                                    if len(current_seq) <= 30:
                                        if final_ratio >= PEVK_ratio:
                                            result_final.append(pevk_seqs[0])
                                    if len(current_seq) > 30:
                                        if final_ratio >= 0.50:
                                            result_final.append(pevk_seqs[0])
                        

                        ### Case 2. If there are multiple items in result:
                        if len(pevk_seqs) > 1:
                            # Create another list of results that will hold contiguous and non-contiguous sequences
                            result_with_contigs = []
                            
                            # Create a list that will hold the initial contiguous ranges, that will be used to create result_with_contigs
                            range_result = []
                            
                            # To look for contiguous ranges, first create a new list that holds the distinct ranges:
                            ranges = []
                            
                            # First, only add "real" results to ranges. Forget about PEVK seqs after this point
                            for i in range(0, len(pevk_seqs)):
                                if pevk_seqs[i] != "x":
                                    ranges.append((pevk_seqs[i][1],pevk_seqs[i][2]))
                            

                            # Now create a list with overlapping ranges combined
                            # Code source: https://codereview.stackexchange.com/questions/21307/consolidate-list-of-ranges-that-overlap
                            # Winston Ewert, stack exchange
                            if ranges != []:
                                current_start = -1
                                current_stop = -1
                                for start, stop in sorted(ranges):
                                    if start > current_stop:
                                        range_result.append((start,stop))
                                        current_start, current_stop = start, stop
                                    else:
                                        range_result[-1] = (current_start, stop)
                                        current_stop = max(current_stop, stop)
                            #### End of borrowed code #####
                            

                            # Create result_with_contigs with these ranges:
                            if range_result != []:
                                for i in range(0, len(range_result)):
                                    current_start = range_result[i][0]
                                    current_stop = range_result[i][1]
                                    current_seq = str(prot_seq[current_start:current_stop])
                                    result_with_contigs.append((current_seq, current_start, current_stop))
                            

                            # Now loop through result_with_contigs and check to see if each sequence contains at least 50% PEVK (for larger sequences)
                            if result_with_contigs != []:
                                for i in range(0, len(result_with_contigs)):
                                    current_seq = result_with_contigs[i][0]
                                    
                                    # Check to see if the current sequence is >50 PEVK:
                                    final_counts = pevk_counter(current_seq)
                                    if len(current_seq) > 0:
                                        final_ratio = float(final_counts)/float(len(current_seq))
                                        if len(current_seq) <= 30:
                                            if final_ratio >= PEVK_ratio:
                                                result_final.append(result_with_contigs[i])
                                        if len(current_seq) > 30:
                                            if final_ratio >= 0.50:
                                                result_final.append(result_with_contigs[i])
                            

                            # Now we have our result in result_final. If this list isn't empty:
                            if result_final != []:

                                for i in range(0, len(result_final)):
                                    # Check to see if the length of the results is greater than the min length parameter
                                    if len(result_final[i][0]) >= int(min_length):
                                        prot_start = result_final[i][1]
                                        prot_end = result_final[i][2]
                                        final_seq = prot_seq[prot_start:prot_end]
                                        prot_seq_record = SeqRecord(final_seq, id=identification+ "_" + str(prot_start) + ":" + str(prot_end), description=species_name)
                                        
                                        #prot_regions.append(prot_seq_record)

                                        # Find PEVK ratio, add to percentages
                                        final_counts = pevk_counter(final_seq)
                                        final_ratio = float(final_counts)/float(len(final_seq))
                                        percentages.append(final_ratio)

                                        # Finish up based on what frame we're in
                                        if identification == "frame_f1":
                                            nuc_start = prot_start*3
                                            nuc_end = prot_end*3
                                            nt_seq = str(nucleotide_seq[nuc_start:nuc_end]) #extract the full nucleotide sequence defined by that region
                                            nt_seq_record = SeqRecord(Seq(nt_seq, IUPAC.unambiguous_dna), id=identification+ "_" + str(nuc_start) + ":" + str(nuc_end-1), description='')
                                            #nuc_regions.append(nt_seq_record)
                                            names.append(identification+ "_" + str(nuc_start) + ":" + str(nuc_end-1))
                                            lengths.append(nuc_end - nuc_start)
                                            starts.append(nuc_start)
                                            ends.append(nuc_end-1)
                                            frame_names.append(identification)
                                            prot1_regions.append(prot_seq_record)
                                            nuc1_regions.append(nt_seq_record)
                                        
                                        if identification == "frame_f2":
                                            nuc_start = prot_start*3 + 1
                                            nuc_end = prot_end*3 + 1
                                            nt_seq = str(nucleotide_seq[nuc_start:nuc_end]) #extract the full nucleotide sequence defined by that region
                                            nt_seq_record = SeqRecord(Seq(nt_seq, IUPAC.unambiguous_dna), id=identification+ "_" + str(nuc_start) + ":" + str(nuc_end-1), description='')
                                            #nuc_regions.append(nt_seq_record)
                                            names.append(identification+ "_" + str(nuc_start) + ":" + str(nuc_end-1))
                                            lengths.append(nuc_end - nuc_start)
                                            starts.append(nuc_start)
                                            ends.append(nuc_end-1)
                                            frame_names.append(identification)
                                            prot2_regions.append(prot_seq_record)
                                            nuc2_regions.append(nt_seq_record)
                                        
                                        if identification == "frame_f3":
                                            nuc_start = prot_start*3 + 2
                                            nuc_end = prot_end*3 + 2
                                            nt_seq = str(nucleotide_seq[nuc_start:nuc_end]) #extract the full nucleotide sequence defined by that region
                                            nt_seq_record = SeqRecord(Seq(nt_seq, IUPAC.unambiguous_dna), id=identification+ "_" + str(nuc_start) + ":" + str(nuc_end-1), description='')
                                            #nuc_regions.append(nt_seq_record)
                                            names.append(identification+ "_" + str(nuc_start) + ":" + str(nuc_end-1))
                                            lengths.append(nuc_end - nuc_start)
                                            starts.append(nuc_start)
                                            ends.append(nuc_end-1)
                                            frame_names.append(identification)
                                            prot3_regions.append(prot_seq_record)
                                            nuc3_regions.append(nt_seq_record)


    # Combine results from frames, and sort by titin coordinates:

    # Amino acid:
    combined_prot = prot1_regions + prot2_regions + prot3_regions
    sorted_prot = sorted(combined_prot, key=lambda x: int(x.id.split("_")[2].split(":")[0]))

    # Nucleotide:
    combined_nuc = nuc1_regions + nuc2_regions + nuc3_regions
    sorted_nuc = sorted(combined_nuc, key=lambda x: int(x.id.split("_")[2].split(":")[0]))   
    #print(sorted_nuc)

    nuc_starts = [int(x.id.split("_")[2].split(":")[0]) for x in sorted_nuc]
    nuc_ends = [int(x.id.split("_")[2].split(":")[1]) for x in sorted_nuc]

    # Find lower and upper bounds of PEVK region
    q3, q1 = np.percentile(nuc_starts, [75,25])
    iqr = q3 - q1

    lower = q1 - (iqr*1.5)
    upper = q3 + (iqr*1.5)

    lower_bound_index = 0
    upper_bound_index = 0

    for i in range(0,len(nuc_starts)-1):
        if nuc_starts[i] < lower and nuc_starts[i+1] > lower:
            lower_bound_value = nuc_starts[i+1]
            lower_bound_index = i+1
        if nuc_starts[i] < upper and nuc_starts[i+1] > upper:
            upper_bound_value = nuc_starts[i]
            upper_bound_index = i
    
    # If lower bound is less than the first sequence start
    if lower < nuc_starts[0]:
        lower_bound_value = nuc_starts[0]
        lower_bound_index = 0
    # If upper bound is greater than the last sequence start
    if upper > nuc_starts[-1]:
        upper_bound_value = nuc_starts[-1]
        upper_bound_index = -2 # because 1 will be added later

    # 'Best guess' exon boundaries output:
    # Create an output file with the amino acid sequences for all found PEVK exons (location relative to full titin DNA sequence)
    bounded_aa_outdir = outpath + "bounded/translated/"
    if not os.path.exists(bounded_aa_outdir):
        os.makedirs(bounded_aa_outdir)
    bounded_nt_outdir = outpath + "bounded/untranslated/"
    if not os.path.exists(bounded_nt_outdir):
        os.makedirs(bounded_nt_outdir)
    unbounded_aa_outdir = outpath + "unbounded/translated/"
    if not os.path.exists(unbounded_aa_outdir):
        os.makedirs(unbounded_aa_outdir)
    unbounded_nt_outdir = outpath + "unbounded/untranslated/"
    if not os.path.exists(unbounded_nt_outdir):
        os.makedirs(unbounded_nt_outdir)


    SeqIO.write(sorted_prot[lower_bound_index:upper_bound_index+1], bounded_aa_outdir + species_name +"PEVK_exons_AA_bounded_"+str(length)+"_"+str(in_ratio)+"_"+str(min_length)+".fasta", "fasta")
    # Create an output file with the nucleotide sequences for all found PEVK exons (location relative to full titin DNA sequence)
    SeqIO.write(sorted_nuc[lower_bound_index:upper_bound_index+1], bounded_nt_outdir + species_name +"PEVK_exons_NT_bounded_"+str(length)+"_"+str(in_ratio)+"_"+str(min_length)+".fasta", "fasta")


    # All exons (no boundaries) output:
    # Create an output file with the amino acid sequences for all found PEVK exons (location relative to full titin DNA sequence)
    SeqIO.write(sorted_prot, unbounded_aa_outdir + species_name +"PEVK_exons_AA_unbounded_"+str(length)+"_"+str(in_ratio)+"_"+str(min_length)+".fasta", "fasta")

    # Create an output file with the nucleotide sequences for all found PEVK exons (location relative to full titin DNA sequence)
    SeqIO.write(sorted_nuc, unbounded_nt_outdir + species_name +"PEVK_exons_NT_unbounded_"+str(length)+"_"+str(in_ratio)+"_"+str(min_length)+".fasta", "fasta")



    if optional_outputs == True:
        # Optional output 1: Create a csv file that contains the exon names and their corresponding lengths and PEVK ratios
        lengths_and_ratios_outdir = outpath + "exon_lengths_and_ratios/"
        if not os.path.exists(lengths_and_ratios_outdir):
            os.makedirs(lengths_and_ratios_outdir)
        with open(lengths_and_ratios_outdir + species_name + 'exon_lengths_and_ratios'+"_"+str(length)+"_"+str(in_ratio)+"_"+str(min_length)+".csv", 'w') as csv_file:
            wr = csv.writer(csv_file, delimiter=',')
            wr.writerow(names)
            wr.writerow(lengths)
            wr.writerow(percentages)


        # Optional output 2: Create a csv file with frames, starts, ends and lengths
        coordinates_outdir = outpath + "exon_coordinates/"
        if not os.path.exists(coordinates_outdir):
            os.makedirs(coordinates_outdir)
        with open(coordinates_outdir + species_name + 'exon_coordinates'+"_"+str(length)+"_"+str(in_ratio)+"_"+str(min_length)+".csv", 'w') as csv_file:
            wr = csv.writer(csv_file, delimiter=',')
            wr.writerow(frame_names)
            wr.writerow(starts)
            wr.writerow(ends)
            wr.writerow(lengths)


def main():

    # get input/output/filename
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", default="./ttn_seqs/Homo_sapiens_ttn.fasta",
                   help="Path to the file(s) that contain the full titin nucleotide sequence, the default is ./ttn_seqs/Homo_sapiens_ttn.fasta",
                   type=str)
    ap.add_argument("-w", "--window_length", default=10,
                   help="The integer value of the sliding window length that PEVK Finder will use to find PEVK exons, the default is 10",
                   type=int)
    ap.add_argument("-r", "--pevk_ratio", default=0.54,
                   help="The float value of the minimum PEVK ratio that PEVK Finder will use to find PEVK exons, the default is 0.54",
                   type=float)
    ap.add_argument("-l", "--exon_length", default=12,
                   help="The integer value of the minimum PEVK exon length that PEVK Finder will use to find PEVK exons, the default is 12",
                   type=int)
    ap.add_argument("-o", "--outpath", default="./data/",
                   help="Path to the directory where the output files will be stored, the default is ./data/",
                   type=str)
    ap.add_argument("-v", "--verbose", default=True,
                   help="When verbose is True, will emit messages about script progress, the default is True",
                   type=lambda x:bool(distutils.util.strtobool(x)))
    ap.add_argument("-p", "--optional_outputs", default=False,
                   help="When optional_outputs is True, will create the two optional output files containing exon length, PEVK ratio and coordinate information, the default is False",
                   type=lambda x:bool(distutils.util.strtobool(x)))

    args = ap.parse_args()

    if args.verbose:
        print("\n")
        print "Titin Sequence(s): " + str(args.input)
        print "Sliding Window Length: " + str(args.window_length)
        print "Minimum PEVK Ratio: " + str(args.pevk_ratio)
        print "Minimum Exon Length: " + str(args.exon_length)
        print("\n")
    
    # Run PEVK Finder
    if os.path.isdir(args.input):
        for filename in os.listdir(args.input):
            if filename.endswith('.fasta'):
                print "RUNNING: " + filename
                filepath = args.input + filename
                pevk_finder(filepath, args.window_length, args.pevk_ratio, args.exon_length, args.outpath, args.optional_outputs)
    if not os.path.isdir(args.input):
        print "RUNNING: " + args.input
        pevk_finder(args.input, args.window_length, args.pevk_ratio, args.exon_length, args.outpath, args.optional_outputs)

    # Write output files

if __name__ == "__main__":
   main()

