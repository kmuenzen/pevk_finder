# This file contains code for figures representing the divergence and dnds results of the 
# discrete PEVK exons, relative to the human exon-intron organization 

# Figure 1: Extended exon-intron plots with significance asterisks over codons under
# positive selection, according to models 7 vs. 8 in PAML, with a sig value under 0.05
# Approach: split segment into 8000 NT increments (total 64000 NT long)
# For each line of exons, have another line on top that has the locations of where a
# "*" should be to signify positive selection
# Have non-orthologous exons in grey, orthologous exons in color
# May have to make a separate figure for each line, and then combine them all afterwards

# Import R color brewer
library(RColorBrewer)
library(randomcoloR)
library(stringr)

# 1. Read CSV of all human exons into R
all_exon_data <- read.csv("Homo_sapiens_test_locations_10_0.54_12.csv", header=FALSE)
all_exons <- t(all_exon_data)
all_exons <- data.frame(all_exons)
colnames(all_exons) <- c("frame","start","end","length")
all_exons$start <- as.numeric(as.character(all_exons$start))
all_exons$end <- as.numeric(as.character(all_exons$end))
all_exons$length <- as.numeric(as.character(all_exons$length))
all_exons <- all_exons[order(all_exons$start),]
rownames(all_exons) <- 1:nrow(all_exons)
all_exons <- as.data.frame(all_exons)

# Get starts and stops for all exons
all_starts <- as.numeric(as.character(all_exons$start))
all_ends <- as.numeric(as.character(all_exons$end))



# 2. Read CSV of the 37 orthologous exons into R
#ortholog_data <- read.csv("ortholog_data_7_8.csv")
ortholog_data <- read.csv("ortholog_dnds_data_cleaned.csv")
ortholog_data <- as.data.frame(ortholog_data)

# Initialize list to hold orthologous exon identifiers
ortholog_ids <- c()

# Initialize list to hold ortholog starts
ortholog_starts <- c()

# Initialize list to hold ortholog ends
ortholog_ends <- c()

# Loop through csv file and get unique exon identifiers
for (i in 1:nrow(ortholog_data)) {
  current_exon <- as.character(ortholog_data$ortholog[i])
  
  # Is the exon already in the exon id list?
  if (is.element(current_exon, ortholog_ids) == FALSE && current_exon != "") {
    # Add the exon to the id list
    ortholog_ids <- c(ortholog_ids, current_exon)
    
    # Get exon starts and ends
    ortholog_start <- as.numeric(strsplit(current_exon, "_")[[1]][3])
    ortholog_end <- as.numeric(strsplit(current_exon, "_")[[1]][4])
    
    # Add exon starts and end to their respective lists
    ortholog_starts <- c(ortholog_starts, ortholog_start)
    ortholog_ends <- c(ortholog_ends, ortholog_end)
  }
}

# 3. Make a list of "grey" exons that are not orthologs
# Initialize list of non_ortholog_starts
non_ortholog_starts <- c()
non_ortholog_ends <- c()

for (i in 1:length(all_starts)) {
  exon_start <- all_starts[i]
  exon_end <- all_ends[i]
  
  # If the current exon IS NOT an ortholog (i.e. not found in ortholog starts)
  if (is.element(exon_start, ortholog_starts) == FALSE) {
    # Add this exon's starts and ends to the list of non-ortholog starts and ends
    non_ortholog_starts <- c(exon_start, non_ortholog_starts)
    non_ortholog_ends <- c(exon_end, non_ortholog_ends)
  }
}

# 3a. Make a list of "empty" exons that were missed by PEVK Finder
missed_starts <- c(112749, 113417, 137978, 141965, 165108)
missed_ends <- c(112823, 113500, 138049, 142045, 165185)

# 4. Figure out where to put asterisks
# Will plot all sites under selection, and color code for significance...red for high signficance, purple for mid, and blue for low
# To calculate location of the codon under selection...extract number part of the codon under selection, multiply by 3, and add to the start location of the exon
# Initialize the vectors that will hold the horizontal asterisk coordinates for each significance level
sig1_coordinates <- c()
sig2_coordinates <- c()
sig3_coordinates <- c()

for (i in 1:nrow(ortholog_data)) {
  # Get exon name, start and stop
  exon_name <- as.character(ortholog_data[i,1])
  exon_start <- as.numeric(strsplit(exon_name, "_")[[1]][3])
  exon_end <- as.numeric(strsplit(exon_name, "_")[[1]][4])
  
  sig1_codons <- as.character(ortholog_data$"sites...0.95"[i])
  sig2_codons <- as.character(ortholog_data$"sites...0.75"[i])
  sig3_codons <- as.character(ortholog_data$"sites...0.5"[i])
  
  if (sig1_codons != "" | sig2_codons != "" | sig3_codons != "") {
    sig1_vector <- strsplit(sig1_codons, ", ")
    sig2_vector <- strsplit(sig2_codons, ", ")
    sig3_vector <- strsplit(sig3_codons, ", ")
    
    sig1_numvec <- c()
    sig2_numvec <- c()
    sig3_numvec <- c()
    
    # Loop through each of these vectors and create vectors with only numeric locations
    for (j in 1:length(sig1_vector[[1]])) {
      current_codon <- sig1_vector[[1]][j]
      numeric_value <- as.numeric(str_extract(current_codon, "[0-9]+"))
      sig1_numvec <- c(sig1_numvec, numeric_value)
    }
    
    for (j in 1:length(sig2_vector[[1]])) {
      current_codon <- sig2_vector[[1]][j]
      numeric_value <- as.numeric(str_extract(current_codon, "[0-9]+"))
      sig2_numvec <- c(sig2_numvec, numeric_value)
    }
    
    for (j in 1:length(sig3_vector[[1]])) {
      current_codon <- sig3_vector[[1]][j]
      numeric_value <- as.numeric(str_extract(current_codon, "[0-9]+"))
      sig3_numvec <- c(sig3_numvec, numeric_value)
    }
    
    # Now that the relative locations are specified, we need to specify the location in reference to the nucleotide location
    # and add to selection coordinates vectors, only if's not an empty string
    
    for (k in 1:length(sig1_numvec)) {
      relative_location <- sig1_numvec[k]*3 + exon_start - 2
      if (is.na(relative_location) == FALSE) {
        sig1_coordinates <- c(sig1_coordinates, relative_location)
      }
    }
    
    for (k in 1:length(sig2_numvec)) {
      relative_location <- sig2_numvec[k]*3 + exon_start - 2
      if (is.na(relative_location) == FALSE) {
        sig2_coordinates <- c(sig2_coordinates, relative_location)
      }
    }
    
    for (k in 1:length(sig3_numvec)) {
      relative_location <- sig3_numvec[k]*3 + exon_start - 2
      if (is.na(relative_location) == FALSE) {
        sig3_coordinates <- c(sig3_coordinates, relative_location)
      }
    }
  }
}

# AA 45, NT 133
# AA 49 NT 145
sig1_coordinates
sig2_coordinates
sig3_coordinates

# 5. Plot each row of exons separately

# SEGMENT 1: 106000 - 114000
pdf(file="segment1", width=50, height=7, pointsize=22)
op <- par(bg = "white", bty="n", mar=c(0.5,0.5,0.5,10))
plot(c(106000,114000), c(0,0.75), type="n", xlab="", ylab="", yaxt="n", xaxt="n")

rect(ortholog_starts, 0, ortholog_ends, 0.5, col="gray35", border=NA)
rect(non_ortholog_starts, 0, non_ortholog_ends, 0.5, col="lightgrey", border=NA)
rect(missed_starts, 0, missed_ends, 0.5, col=NA, border= 'black')

abline(h=0.25)

# Plot selection coordinates above 
points(sig1_coordinates, y=rep(0.55, times=length(sig1_coordinates)), pch=8, col="red", cex=3)
points(sig2_coordinates, y=rep(0.65, times=length(sig2_coordinates)), pch=8, col="green", cex=3)
points(sig3_coordinates, y=rep(0.75, times=length(sig3_coordinates)), pch=8, col="blue", cex=3)

mtext("114000bp", side=4, cex=2, las=1, padj=2.7)

dev.off()

# SEGMENT 2: 114001 - 122000
pdf(file="segment2",width=50, height=7, pointsize=22)
op <- par(bg = "white", bty="n", mar=c(0.5,0.5,0.5,10))
plot(c(114500,122000), c(0,0.75), type="n", xlab="", ylab="", yaxt="n", xaxt="n")

rect(ortholog_starts, 0, ortholog_ends, 0.5, col="gray35", border=NA)
rect(non_ortholog_starts, 0, non_ortholog_ends, 0.5, col="lightgrey", border=NA)
rect(missed_starts, 0, missed_ends, 0.5, col=NA, border= 'black')

abline(h=0.25)

# Plot selection coordinates above 
points(sig1_coordinates, y=rep(0.55, times=length(sig1_coordinates)), pch=8, col="red", cex=3)
points(sig2_coordinates, y=rep(0.67, times=length(sig2_coordinates)), pch=8, col="green", cex=3)
points(sig3_coordinates, y=rep(0.75, times=length(sig3_coordinates)), pch=8, col="blue", cex=3)

mtext("122000bp", side=4, cex=2, las=1, padj=2.7)

dev.off()

# SEGMENT 3: 122001 - 130000
pdf(file="segment3", width=50, height=7, pointsize=22)
op <- par(bg = "white", bty="n", mar=c(0.5,0.5,0.5,10))
plot(c(122500,130000), c(0,0.75), type="n", xlab="", ylab="", yaxt="n", xaxt="n")

rect(ortholog_starts, 0, ortholog_ends, 0.5, col="gray35", border=NA)
rect(non_ortholog_starts, 0, non_ortholog_ends, 0.5, col="lightgrey", border=NA)
rect(missed_starts, 0, missed_ends, 0.5, col=NA, border= 'black')

abline(h=0.25)

# Plot selection coordinates above 
points(sig1_coordinates, y=rep(0.55, times=length(sig1_coordinates)), pch=8, col="red", cex=3)
points(sig2_coordinates, y=rep(0.65, times=length(sig2_coordinates)), pch=8, col="green", cex=3)
points(sig3_coordinates, y=rep(0.75, times=length(sig3_coordinates)), pch=8, col="blue", cex=3)

mtext("130000bp", side=4, cex=2, las=1, padj=2.7)

dev.off()

# SEGMENT 4: 130001 - 138000
pdf(file="segment4", width=50, height=7, pointsize=22)
op <- par(bg = "white", bty="n", mar=c(0.5,0.5,0.5,10))
plot(c(130550,138000), c(0,0.75), type="n", xlab="", ylab="", yaxt="n", xaxt="n")

rect(ortholog_starts, 0, ortholog_ends, 0.5, col="gray35", border=NA)
rect(non_ortholog_starts, 0, non_ortholog_ends, 0.5, col="lightgrey", border=NA)
rect(missed_starts, 0, missed_ends, 0.5, col=NA, border= 'black')

abline(h=0.25)

# Plot selection coordinates above 
points(sig1_coordinates, y=rep(0.55, times=length(sig1_coordinates)), pch=8, col="red", cex=3)
points(sig2_coordinates, y=rep(0.65, times=length(sig2_coordinates)), pch=8, col="green", cex=3)
points(sig3_coordinates, y=rep(0.75, times=length(sig3_coordinates)), pch=8, col="blue", cex=3)

mtext("138000bp", side=4, cex=2, las=1, padj=2.7)

dev.off()

# SEGMENT 5: 138001 - 146000
pdf(file="segment5", width=50, height=7, pointsize=22)
op <- par(bg = "white", bty="n", mar=c(0.5,0.5,0.5,10))
plot(c(138500,146000), c(0,0.75), type="n", xlab="", ylab="", yaxt="n", xaxt="n")

rect(ortholog_starts, 0, ortholog_ends, 0.5, col="gray35", border=NA)
rect(non_ortholog_starts, 0, non_ortholog_ends, 0.5, col="lightgrey", border=NA)
rect(missed_starts, 0, missed_ends, 0.5, col=NA, border= 'black')

abline(h=0.25)

# Plot selection coordinates above 
points(sig1_coordinates, y=rep(0.55, times=length(sig1_coordinates)), pch=8, col="red", cex=3)
points(sig2_coordinates, y=rep(0.65, times=length(sig2_coordinates)), pch=8, col="green", cex=3)
points(sig3_coordinates, y=rep(0.75, times=length(sig3_coordinates)), pch=8, col="blue", cex=3)

mtext("146000bp", side=4, cex=2, las=1, padj=2.7)

dev.off()

# SEGMENT 6: 146001 - 154000
pdf(file="segment6", width=50, height=7, pointsize=22)
op <- par(bg = "white", bty="n", mar=c(0.5,0.5,0.5,10))
plot(c(146500,154000), c(0,0.75), type="n", xlab="", ylab="", yaxt="n", xaxt="n")

rect(ortholog_starts, 0, ortholog_ends, 0.5, col="gray35", border=NA)
rect(non_ortholog_starts, 0, non_ortholog_ends, 0.5, col="lightgrey", border=NA)
rect(missed_starts, 0, missed_ends, 0.5, col=NA, border= 'black')

abline(h=0.25)

# Plot selection coordinates above 
points(sig1_coordinates, y=rep(0.55, times=length(sig1_coordinates)), pch=8, col="red", cex=3)
points(sig2_coordinates, y=rep(0.65, times=length(sig2_coordinates)), pch=8, col="green", cex=3)
points(sig3_coordinates, y=rep(0.75, times=length(sig3_coordinates)), pch=8, col="blue", cex=3)

mtext("154000bp", side=4, cex=2, las=1, padj=2.7)

dev.off()

# SEGMENT 7: 154001 - 162000
pdf(file="segment7", width=50, height=7, pointsize=22)
op <- par(bg = "white", bty="n", mar=c(0.5,0.5,0.5,10))
plot(c(154500,162000), c(0,0.75), type="n", xlab="", ylab="", yaxt="n", xaxt="n")

rect(ortholog_starts, 0, ortholog_ends, 0.5, col="gray35", border=NA)
rect(non_ortholog_starts, 0, non_ortholog_ends, 0.5, col="lightgrey", border=NA)
rect(missed_starts, 0, missed_ends, 0.5, col=NA, border= 'black')

abline(h=0.25)

# Plot selection coordinates above 
points(sig1_coordinates, y=rep(0.55, times=length(sig1_coordinates)), pch=8, col="red", cex=3)
points(sig2_coordinates, y=rep(0.65, times=length(sig2_coordinates)), pch=8, col="green", cex=3)
points(sig3_coordinates, y=rep(0.75, times=length(sig3_coordinates)), pch=8, col="blue", cex=3)

mtext("162000bp", side=4, cex=2, las=1, padj=2.7)

dev.off()

# SEGMENT 8: 162001 - 170000
pdf(file="segment8", width=50, height=7, pointsize=22)
op <- par(bg = "white", bty="n", mar=c(0.5,0.5,0.5,10))
plot(c(162500,170000), c(0,0.75), type="n", xlab="", ylab="", yaxt="n", xaxt="n")

rect(ortholog_starts, 0, ortholog_ends, 0.5, col="gray35", border=NA)
rect(non_ortholog_starts, 0, non_ortholog_ends, 0.5, col="lightgrey", border=NA)
rect(missed_starts, 0, missed_ends, 0.5, col=NA, border= 'black')

abline(h=0.25)

# Plot selection coordinates above 
points(sig1_coordinates, y=rep(0.55, times=length(sig1_coordinates)), pch=8, col="red", cex=3)
points(sig2_coordinates, y=rep(0.65, times=length(sig2_coordinates)), pch=8, col="green", cex=3)
points(sig3_coordinates, y=rep(0.75, times=length(sig3_coordinates)), pch=8, col="blue", cex=3)

mtext("170000bp", side=4, cex=2, las=1, padj=2.7)

dev.off()

