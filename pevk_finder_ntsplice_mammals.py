# Import libraries

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
import statistics
import csv
import sys

# Define reading frame objects

frame_f1_in = sys.argv[1] #frame f1 sequence
frame_f2_in = sys.argv[2] #frame f2 sequence
frame_f3_in = sys.argv[3] #frame f3 sequence

# Define complete nucleotide sequence and input parameters

nt_seq = sys.argv[4] #the original nucleotide sequence
length = sys.argv[5] #frame length used as sliding frame
in_ratio = sys.argv[6] #desired PEVK ratio
min_length = sys.argv[7]

# Create a Seq object for the forward and reverse nucleotide sequences

ref_seq = SeqIO.read(nt_seq,"fasta") #read the nucleotide sequence into biopython
nucleotide_seq = Seq(str(ref_seq.seq)+"*", IUPAC.unambiguous_dna) #make a seq object
rev_nucleotide_seq = Seq(str(ref_seq.seq)[::-1], IUPAC.unambiguous_dna) #make a seq object for the reversed sequence

# Compile each frame into a list

frames = [frame_f1_in, frame_f2_in, frame_f3_in]

# Initialize the lists for the final csv file:

names = []
lengths = []
percentages = []
starts = []
ends = []
frame_names = []

##### Begin the for loop #####
for f in frames:
    ##### Get frame sequence #####

    handle = SeqIO.read(f,"fasta") #read the protein sequence into biopython
    prot_seq = Seq(str(handle.seq)+"*", IUPAC.protein) #make a seq object

    ##### Get Seq Info #####
    identification = handle.id #which frame
    description = handle.description #store the object name as a string for later identification
    split = handle.description.split()
    species_name = split[2] + "_" + split[3] + "_"

    ##### Set Search Parameters #####
    input_frame_length = int(length) #how big is the sliding window
    PEVK_ratio = float(in_ratio) #what percentage PEVK are we looking for?

    ##### Helper Funtions #####
    """
    helper function pevk_counter accepts a string of
    amino acids (frame) and determines whether the frame
    contains a sufficient percentage of P, E, V and K (and A)
    """

    def pevk_counter(frame):
        counter = 0
        for a in frame:
            if a == "P" or a == "E" or a == "V" or a == "K" or a == "A":
                counter += 1
        return counter

    # Helper function frame_define fefines the between stop codon frame:
    def frame_define(seq):
        counter = 0
        for i in seq:
            if i != "*":
                counter += 1
            if i == "*":
                return counter
                break


    ##### Main Function: pevk_finder #####
    """
    This function takes in an amino acid sequence between stop codons
    """
    ########## START OF NEW CODE 6/12/17 #########
    ## First loop forward to find all possible PEVK regions
    def find_pevk(seq, seq_start, seq_end):
        # Define the start and end of the corresponding nucleotide sequence
        # This depends on which reading frame we're working with:
        if identification == "frame_f1":
            nuc_ref_start = seq_start*3
            nuc_ref_end = seq_end*3
        if identification == "frame_f2":
            nuc_ref_start = seq_start*3 + 1
            nuc_ref_end = seq_end*3 + 1
        if identification == "frame_f3":
            nuc_ref_start = seq_start*3 + 2
            nuc_ref_end = seq_end*3 + 2
        # Define the actual sequence (as a string) of the corresponding nucleotide sequence
        nuc_ref_seq = str(nucleotide_seq[nuc_ref_start:nuc_ref_end])
        # Initialize a list that will hold the results in format (spliced_seq, absolute_seq_start, absolute_seq_end)
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
                # If there is only one starting position
                if len(group_positions) == 1:
                    groupings = [primary_positions]
                # If there are multiple starting positions
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
                #print groupings
                for j in range(0, len(groupings)):
                    # Initialize list that will hold candidate donor sites
                    donor_sites = []
                    # Initialize list that will hold candidate acceptor sites
                    acceptor_sites = []
                    # Initialize the list that will hold difference in location between possible donor start sites and the PEVK region start site
                    start_differences = []
                    # Start with the list of primary positions. Find the median of that list.
                    AA_starting_point = int(statistics.median(groupings[j]))
                    # Where is this in the local nucleotide sequence?
                    NT_starting_point = AA_starting_point*3 # This is the first letter of the codon encoding the median AA in the high PEVK region
                    # First go backwards to look for the DONOR splice site
                    for i in range(NT_starting_point,0,-1):
                        # The rightmost letter (first letter of a codon) should be a G
                        if i%3 == 0 and nuc_ref_seq[i] == "G":
                            # The leftmost letter (last letter of a codon) should be an A
                            if nuc_ref_seq[i-1] == "A":
                                # If both of these requirements are met, add the location of G to donor_sites
                                donor_sites.append(i)
                    # Is the donor site list empty? If so, this is not a real PEVK sequence!
                    if donor_sites == []:
                        result.append("x")
                        continue
                    # If there is a donor site:
                    else:
                        # We now want to know which donor site is closest to the first AA of primary positions
                        # First define the starting position of the high PEVK region
                        #PEVK_start = groupings[j][0]*3
                        PEVK_start = NT_starting_point
                        # Then find the closest donor site to the start of the PEVK region
                        for x in donor_sites:
                            diff = abs(PEVK_start - x)
                            start_differences.append(diff)
                        # The best guess G of the start splice site is the min of differences
                        donor_g_location = donor_sites[start_differences.index(min(start_differences))]
                        # Now go forwards to look for the ACCEPTOR splice site
                        for i in range(NT_starting_point, len(nuc_ref_seq)-3):
                            # The leftmost letter (third letter of a codon) should be an A
                            if i%3 == 2 and nuc_ref_seq[i] == "A":
                                # The second letter (first letter of a codon) should be a G
                                if (i+1)%3 == 0 and nuc_ref_seq[i+1] == "G":
                                    # The third letter (second letter of a codon) should be a G
                                    if (i+2)%3 == 1 and nuc_ref_seq[i+2] == "G":
                                        # The fourth letter (third letter of a codon) should be a T
                                        if (i+3)%3 == 2 and nuc_ref_seq[i+3] == "T":
                                            # If all of these requirements are met, add the location of the first G to
                                            acceptor_sites.append(i+1)
                        # If there are no acceptor sites:
                        if acceptor_sites == []:
                            result.append("x")
                            continue
                        # If there is at least one acceptor site
                        if acceptor_sites != []:

                            # If there is only one possible acceptor site:
                            if len(acceptor_sites) == 1:
                                # This becomes the acceptor location
                                acceptor_g1_location = acceptor_sites[0]
                            # If there is more than one possible acceptor site:
                            if len(acceptor_sites) > 1:
                                end_differences = []
                                PEVK_end = (groupings[j][-1]+input_frame_length)*3
                                for x in acceptor_sites:
                                    diff = abs(PEVK_end - x)
                                    end_differences.append(diff)
                                # The best guess G of the start splice site is the min of differences
                                acceptor_g1_location = acceptor_sites[end_differences.index(min(end_differences))]
                            # Now extract the region between the donor and acceptor sites
                            final_seq =  str(Seq(nuc_ref_seq[donor_g_location+3:acceptor_g1_location]).translate())
                            # Define the relative start and end positions (within the local amino acid sequence) of the final sequence
                            relative_seq_start = (donor_g_location/3)+1
                            relative_seq_end = acceptor_g1_location/3
                            # Define the relative start and end positions of the sequence within the whole amino acid sequence
                            absolute_seq_start = seq_start + relative_seq_start
                            absolute_seq_end = seq_start + relative_seq_end
                            ###############
                            #print identification
                            #print absolute_seq_start
                            #print absolute_seq_end
                            #print final_seq
                            ##############
                            # Define the full sequence with no splice sites
                            final_seq = seq[relative_seq_start:relative_seq_end]
                            # Put the sequence and the start and end locations in the result list
                            result.append((final_seq, absolute_seq_start, absolute_seq_end))
            #print result
            return result

############ END OF NEW CODE 6/12/17 ###########

    # Loop through the entire protein sequence, extracting each region between stop codons and running it through find_pevk
    prot_regions = []
    nuc_regions = []
    for i in range(len(prot_seq)-1):
        if prot_seq[i] == "*": #when we get to a stop codon
            seq_start = i+1 #define the position of the first AA in the exon
            if prot_seq[seq_start] == "*":
                continue
            else:
                seq_length = frame_define(str(prot_seq)[seq_start:]) #find the length of the seq between 2 stop codons
                seq_end = seq_start + seq_length #define the position of the last AA in the exon (don't include stop codon)
                seq = str(prot_seq[seq_start:seq_end]) # get the actual AA sequence of the frame
                #print seq
                #print seq_start
                #print seq_end
                #####################
                pevk_seqs = find_pevk(seq, seq_start, seq_end)

                #print pevk_seqs
                ####################
                # At this point, we have a result in the form of [(final_seq, absolute_seq_start, absolute_seq_end)], or a result of "x". These are compiled into the result list
                # First check to see if the result is null:
                if pevk_seqs != "x" and pevk_seqs != ["x"]:
                    # Initialize result_final, which will hold the contiguous and high PEVK results
                    result_final = []
                    # Loop through the initial result. If there is only one item in result, check to see if the sequence is greater than 50% PEVK, then add it to result_final
                    if len(pevk_seqs) == 1:
                        current_seq = pevk_seqs[0][0]
                        if current_seq != "x":
                            # Check to see if the current sequence is >50 PEVK:
                            final_counts = pevk_counter(current_seq)
                            if len(current_seq) > 0:
                                final_ratio = float(final_counts)/float(len(current_seq))
                                if len(current_seq) <= 30:
                                    if final_ratio >= PEVK_ratio:
                                        result_final.append(pevk_seqs[0])
                                if len(current_seq) > 30:
                                    if final_ratio >= 0.50:
                                        result_final.append(pevk_seqs[0])
                    # If there are multiple items in result:
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
                                    prot_regions.append(prot_seq_record)


                                    # Find PEVK ratio, add to percentages
                                    final_counts = pevk_counter(final_seq)
                                    final_ratio = float(final_counts)/float(len(final_seq))
                                    percentages.append(final_ratio)
                                    #names.append(identification+ "_" + str(prot_start) + ":" + str(prot_end))
                                    #lengths.append(prot_end - prot_start + 1
                                    if identification == "frame_f1":
                                        nuc_start = prot_start*3
                                        nuc_end = prot_end*3
                                        nt_seq = str(nucleotide_seq[nuc_start:nuc_end]) #extract the full nucleotide sequence defined by that region
                                        nt_seq_record = SeqRecord(Seq(nt_seq, IUPAC.unambiguous_dna), id=identification+ "_" + str(nuc_start) + ":" + str(nuc_end-1), description=species_name)
                                        nuc_regions.append(nt_seq_record)
                                        names.append(identification+ "_" + str(nuc_start) + ":" + str(nuc_end-1))
                                        lengths.append(nuc_end - nuc_start)
                                        starts.append(nuc_start)
                                        ends.append(nuc_end-1)
                                        frame_names.append(identification)
                                    if identification == "frame_f2":
                                        nuc_start = prot_start*3 + 1
                                        nuc_end = prot_end*3 + 1
                                        nt_seq = str(nucleotide_seq[nuc_start:nuc_end]) #extract the full nucleotide sequence defined by that region
                                        nt_seq_record = SeqRecord(Seq(nt_seq, IUPAC.unambiguous_dna), id=identification+ "_" + str(nuc_start) + ":" + str(nuc_end-1), description=species_name)
                                        nuc_regions.append(nt_seq_record)
                                        names.append(identification+ "_" + str(nuc_start) + ":" + str(nuc_end-1))
                                        lengths.append(nuc_end - nuc_start)
                                        starts.append(nuc_start)
                                        ends.append(nuc_end-1)
                                        frame_names.append(identification)
                                    if identification == "frame_f3":
                                        nuc_start = prot_start*3 + 2
                                        nuc_end = prot_end*3 + 2
                                        nt_seq = str(nucleotide_seq[nuc_start:nuc_end]) #extract the full nucleotide sequence defined by that region
                                        nt_seq_record = SeqRecord(Seq(nt_seq, IUPAC.unambiguous_dna), id=identification+ "_" + str(nuc_start) + ":" + str(nuc_end-1), description=species_name)
                                        nuc_regions.append(nt_seq_record)
                                        names.append(identification+ "_" + str(nuc_start) + ":" + str(nuc_end-1))
                                        lengths.append(nuc_end - nuc_start)
                                        starts.append(nuc_start)
                                        ends.append(nuc_end-1)
                                        frame_names.append(identification)
    # Create an output file for all found PEVK sequences
    SeqIO.write(prot_regions, species_name+identification + "_" +"PEVK_protein_final.fasta", "fasta")

    # Make another fasta file with the corresponding nucleotide sequence
    SeqIO.write(nuc_regions, species_name+identification + "_" +"PEVK_nucleotide_final.fasta", "fasta")

    # Create a csv file that contains the exon names and their corresponding lengths and PEVK ratios
    with open(species_name + 'test_lengths_and_ratios'+"_"+length+"_"+in_ratio+"_"+min_length+".csv", 'w') as csv_file:
        wr = csv.writer(csv_file, delimiter=',')
        wr.writerow(names)
        wr.writerow(lengths)
        wr.writerow(percentages)
    # File with frames, starts, ends and lengths
    with open(species_name + 'test_locations'+"_"+length+"_"+in_ratio+"_"+min_length+".csv", 'w') as csv_file:
        wr = csv.writer(csv_file, delimiter=',')
        wr.writerow(frame_names)
        wr.writerow(starts)
        wr.writerow(ends)
        wr.writerow(lengths)
