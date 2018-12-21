##### Handling Exon Libraries from Unannotated Species ######
# Goals: 
#	a. Organize exons visually, compare with controls and each other
#	b. Make a data frame that collects/displays data on different PEVK characteristics:
#		1. Overall length of region
#		2. Number of exons
#		3. Number of longer exons (>=30 AAs)
#		4. Number of shorter exons (< 30 AAs)
#		5. % PEVK (P vs E vs V vs K)
#		6. Muscle characteristics/phenotype at the species level
#	c. Incorporate trait values into a phylogeny (or multiple phylogenies), do ancestral state reconstructions for different PEVK characteristics

##### Set working directory
setwd("~/Desktop/titin_project/pevk_mammals")
# a. Visual exon organization. This code will be similar to the code used for exon plots during control testing

library(plotrix)
library(SDMTools)
library(gplots)
library(RColorBrewer)
library(MASS)

# Initilize data frame that will hold start/stop index info
starts_and_stops <- data.frame(matrix(ncol=4, nrow=0))
colnames(starts_and_stops) <- c("taxon","start_index","stop_index","hv_start") 

################### Pig ####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Sus_scrofa_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Sus_scrofa_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Sus_scrofa_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Sus_scrofa_start <- 3
Sus_scrofa_end <- 92
Sus_scrofa_hv_start <- 54

starts_and_stops[1,1] <- "Sus_scrofa"
starts_and_stops[1,2] <- Sus_scrofa_start
starts_and_stops[1,3] <- Sus_scrofa_end
starts_and_stops[1,4] <- Sus_scrofa_hv_start

####################### Horse ####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Equus_caballus_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Equus_caballus_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Equus_caballus_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(100000,170000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()
Equus_caballus_start <- 6
Equus_caballus_end <- 105
Equus_caballus_hv_start <- 58

starts_and_stops[2,1] <- "Equus_caballus"
starts_and_stops[2,2] <- Equus_caballus_start
starts_and_stops[2,3] <- Equus_caballus_end
starts_and_stops[2,4] <- Equus_caballus_hv_start
######################## Cat ####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Felis_catus_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Felis_catus_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Felis_catus_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Felis_catus_start <- 4
Felis_catus_end <- 91
Felis_catus_hv_start <- 52

starts_and_stops[3,1] <- "Felis_catus"
starts_and_stops[3,2] <- Felis_catus_start
starts_and_stops[3,3] <- Felis_catus_end
starts_and_stops[3,4] <- Felis_catus_hv_start


#################### Cheetah ##########


# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Acinonyx_jubatus_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Acinonyx_jubatus_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Acinonyx_jubatus_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Acinonyx_jubatus_start <- 4
Acinonyx_jubatus_end <- 90
Acinonyx_jubatus_hv_start <- 53

starts_and_stops[4,1] <- "Acinonyx_jubatus"
starts_and_stops[4,2] <- Acinonyx_jubatus_start
starts_and_stops[4,3] <- Acinonyx_jubatus_end
starts_and_stops[4,4] <- Acinonyx_jubatus_hv_start


##################### Cattle ###############

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Bos_taurus_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Bos_taurus_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Bos_taurus_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()
Bos_taurus_start <- 4
Bos_taurus_end <- 98
Bos_taurus_hv_start <- 56

starts_and_stops[5,1] <- "Bos_taurus"
starts_and_stops[5,2] <- Bos_taurus_start
starts_and_stops[5,3] <- Bos_taurus_end
starts_and_stops[5,4] <- Bos_taurus_hv_start

################ Rabbit ################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Oryctolagus_cuniculus_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Oryctolagus_cuniculus_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Oryctolagus_cuniculus_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Oryctolagus_cuniculus_start <- 6
Oryctolagus_cuniculus_end <- 98
Oryctolagus_cuniculus_hv_start <- 56

starts_and_stops[6,1] <- "Oryctolagus_cuniculus"
starts_and_stops[6,2] <- Oryctolagus_cuniculus_start
starts_and_stops[6,3] <- Oryctolagus_cuniculus_end
starts_and_stops[6,4] <- Oryctolagus_cuniculus_hv_start


################### Rat ##############

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Rattus_norvegicus_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Rattus_norvegicus_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Rattus_norvegicus_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Rattus_norvegicus_start <- 3
Rattus_norvegicus_end <- 91
Rattus_norvegicus_hv_start <- 50

starts_and_stops[7,1] <- "Rattus_norvegicus"
starts_and_stops[7,2] <- Rattus_norvegicus_start
starts_and_stops[7,3] <- Rattus_norvegicus_end
starts_and_stops[7,4] <- Rattus_norvegicus_hv_start


##################### Pangolin #####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Manis_javanica_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Manis_javanica_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Manis_javanica_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(100000,170000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Manis_javanica_start <- 3 
Manis_javanica_end <- 94
Manis_javanica_hv_start <- 57

starts_and_stops[8,1] <- "Manis_javanica"
starts_and_stops[8,2] <- Manis_javanica_start
starts_and_stops[8,3] <- Manis_javanica_end
starts_and_stops[8,4] <- Manis_javanica_hv_start


##################### Koala #####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Phascolarctos_cinereus_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Phascolarctos_cinereus_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Phascolarctos_cinereus_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(135000,205000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()
Phascolarctos_cinereus_start <- 4
Phascolarctos_cinereus_end <- 104
Phascolarctos_cinereus_hv_start <- 52

starts_and_stops[9,1] <- "Phascolarctos_cinereus"
starts_and_stops[9,2] <- Phascolarctos_cinereus_start
starts_and_stops[9,3] <- Phascolarctos_cinereus_end
starts_and_stops[9,4] <- Phascolarctos_cinereus_hv_start


##################### Manatee #####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Trichechus_manatus_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Trichechus_manatus_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Trichechus_manatus_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(100000,170000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()
Trichechus_manatus_start <- 6
Trichechus_manatus_end <- 91
Trichechus_manatus_hv_start <- 59

starts_and_stops[10,1] <- "Trichechus_manatus"
starts_and_stops[10,2] <- Trichechus_manatus_start
starts_and_stops[10,3] <- Trichechus_manatus_end
starts_and_stops[10,4] <- Trichechus_manatus_hv_start


##################### Shrew #####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Mus_pahari_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Mus_pahari_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Mus_pahari_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))
dev.off()
Mus_pahari_start <- 5
Mus_pahari_end <- 79
Mus_pahari_hv_start <- 53

starts_and_stops[11,1] <- "Mus_pahari"
starts_and_stops[11,2] <- Mus_pahari_start
starts_and_stops[11,3] <- Mus_pahari_end
starts_and_stops[11,4] <- Mus_pahari_hv_start

##################### Mole #####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Condylura_cristata_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Condylura_cristata_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Condylura_cristata_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Condylura_cristata_start <- 7
Condylura_cristata_end <- 80
Condylura_cristata_hv_start <- 52

starts_and_stops[12,1] <- "Condylura_cristata"
starts_and_stops[12,2] <- Condylura_cristata_start
starts_and_stops[12,3] <- Condylura_cristata_end
starts_and_stops[12,4] <- Condylura_cristata_hv_start

##################### Chinese Shrew #####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Tupaia_chinensis_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Tupaia_chinensis_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Tupaia_chinensis_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(115000,185000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Tupaia_chinensis_start <- 6
Tupaia_chinensis_end <- 88
Tupaia_chinensis_hv_start <- 54

starts_and_stops[13,1] <- "Tupaia_chinensis"
starts_and_stops[13,2] <- Tupaia_chinensis_start
starts_and_stops[13,3] <- Tupaia_chinensis_end
starts_and_stops[13,4] <- Tupaia_chinensis_hv_start


##################### Armadillo #####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Dasypus_novemcinctus_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Dasypus_novemcinctus_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Dasypus_novemcinctus_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Dasypus_novemcinctus_start <- 5
Dasypus_novemcinctus_end <- 101
Dasypus_novemcinctus_hv_start <- 52

starts_and_stops[14,1] <- "Dasypus_novemcinctus"
starts_and_stops[14,2] <- Dasypus_novemcinctus_start
starts_and_stops[14,3] <- Dasypus_novemcinctus_end
starts_and_stops[14,4] <- Dasypus_novemcinctus_hv_start


##################### Lemur #####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Galeopterus_variegatus_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Galeopterus_variegatus_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Galeopterus_variegatus_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Galeopterus_variegatus_start <- 5
Galeopterus_variegatus_end <- 108
Galeopterus_variegatus_hv_start <- 56

starts_and_stops[15,1] <- "Galeopterus_variegatus"
starts_and_stops[15,2] <- Galeopterus_variegatus_start
starts_and_stops[15,3] <- Galeopterus_variegatus_end
starts_and_stops[15,4] <- Galeopterus_variegatus_hv_start

##################### European Shrew #####################

# Read in csv containing info on length and % pevk
#pevk_lengths_and_ratios <- read.csv("Sorex_araneus_test_lengths_and_ratios_10_0.54_12.csv", header=F)
#pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
#pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
#rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
#colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
#new_starts = c()
#new_ends = c()
#for (i in 1:nrow(pevk_lengths_and_ratios))	{
#	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
#	split <- strsplit(exon_name, "_")
#	location <- split[[1]][3]
#	location <- strsplit(location, ":")
#	match_start <- as.numeric(location[[1]][1])
#	match_end <- as.numeric(location[[1]][2])
#	new_starts <- append(new_starts, match_start)
#	new_ends <- append(new_ends, match_end)
#}

# Add a new column called test_location to the data frame
#pevk_lengths_and_ratios$start <- new_starts
#pevk_lengths_and_ratios$end <- new_ends
#pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
#rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
#test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
#pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
#test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

#my_exons <- read.csv("Sorex_araneus_test_locations_10_0.54_12.csv", header=F)
#my_exons <- t(my_exons)
#my_exons <- data.frame(my_exons)
#colnames(my_exons) <- c("frame","start","end","length")
#my_exons$start <- as.numeric(as.character(my_exons$start))
#my_exons$end <- as.numeric(as.character(my_exons$end))
#my_exons$length <- as.numeric(as.character(my_exons$length))
#my_exons <- my_exons[order(my_exons$start),]
#rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

#my_starts <- as.numeric(as.character(my_exons$start))
#my_ends <- as.numeric(as.character(my_exons$end))

#pdf(file="european_shrew_exon_plot.pdf", width=25, height=10, pointsize=22)

#op <- par(bg = "white")
#plot(c(107000,170000), c(1.7, 3.3), type="n", xlab="Titin Coordinates", ylab="", main="PEVK Exon Distribution (European Shrew)", yaxt="n")

#rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

#dev.off()

#Sorex_araneus_start <- 5
#Sorex_araneus_end <- 109

#starts_and_stops[16,1] <- "Sorex_araneus"
#starts_and_stops[16,2] <- Sorex_araneus_start
#starts_and_stops[16,3] <- Sorex_araneus_end

##################### Ground Squirrel #####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Ictidomys_tridecemlineatus_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Ictidomys_tridecemlineatus_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Ictidomys_tridecemlineatus_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(100000,170000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Ictidomys_tridecemlineatus_start <- 5
Ictidomys_tridecemlineatus_end <- 82
Ictidomys_tridecemlineatus_hv_start <- 53

starts_and_stops[17,1] <- "Ictidomys_tridecemlineatus"
starts_and_stops[17,2] <- Ictidomys_tridecemlineatus_start
starts_and_stops[17,3] <- Ictidomys_tridecemlineatus_end
starts_and_stops[17,4] <- Ictidomys_tridecemlineatus_hv_start


##################### Pika #####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Ochotona_princeps_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Ochotona_princeps_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Ochotona_princeps_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(100000,170000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Ochotona_princeps_start <- 7
Ochotona_princeps_end <- 99
Ochotona_princeps_hv_start <- 58

starts_and_stops[18,1] <- "Ochotona_princeps"
starts_and_stops[18,2] <- Ochotona_princeps_start
starts_and_stops[18,3] <- Ochotona_princeps_end
starts_and_stops[18,4] <- Ochotona_princeps_hv_start

##################### Hedgehog #####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Erinaceus_europaeus_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Erinaceus_europaeus_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Erinaceus_europaeus_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(126000,201000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Erinaceus_europaeus_start <- 6
Erinaceus_europaeus_end <- 108
Erinaceus_europaeus_hv_start <- 53

starts_and_stops[19,1] <- "Erinaceus_europaeus"
starts_and_stops[19,2] <- Erinaceus_europaeus_start
starts_and_stops[19,3] <- Erinaceus_europaeus_end
starts_and_stops[19,4] <- Erinaceus_europaeus_hv_start

##################### Killer Whale #####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Orcinus_orca_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Orcinus_orca_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Orcinus_orca_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(110000,180000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Orcinus_orca_start <- 3
Orcinus_orca_end <- 92
Orcinus_orca_hv_start <- 54

starts_and_stops[20,1] <- "Orcinus_orca"
starts_and_stops[20,2] <- Orcinus_orca_start
starts_and_stops[20,3] <- Orcinus_orca_end
starts_and_stops[20,4] <- Orcinus_orca_hv_start


##################### Flying Fox #####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Pteropus_vampyrus_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Pteropus_vampyrus_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Pteropus_vampyrus_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Pteropus_vampyrus_start <- 5
Pteropus_vampyrus_end <- 87
Pteropus_vampyrus_hv_start <- 54

starts_and_stops[21,1] <- "Pteropus_vampyrus"
starts_and_stops[21,2] <- Pteropus_vampyrus_start
starts_and_stops[21,3] <- Pteropus_vampyrus_end
starts_and_stops[21,4] <- Pteropus_vampyrus_hv_start


##################### Rhinoceros #####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Ceratotherium_simum_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Ceratotherium_simum_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Ceratotherium_simum_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Ceratotherium_simum_start <- 3
Ceratotherium_simum_end <- 97
Ceratotherium_simum_hv_start <- 55

starts_and_stops[22,1] <- "Ceratotherium_simum"
starts_and_stops[22,2] <- Ceratotherium_simum_start
starts_and_stops[22,3] <- Ceratotherium_simum_end
starts_and_stops[22,4] <- Ceratotherium_simum_hv_start

##################### Elephant #####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Loxodonta_africana_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Loxodonta_africana_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Loxodonta_africana_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Loxodonta_africana_start <- 6
Loxodonta_africana_end <- 104
Loxodonta_africana_hv_start <- 58

starts_and_stops[23,1] <- "Loxodonta_africana"
starts_and_stops[23,2] <- Loxodonta_africana_start
starts_and_stops[23,3] <- Loxodonta_africana_end
starts_and_stops[23,4] <- Loxodonta_africana_hv_start


##################### Horseshoe Bat #####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Rhinolophus_sinicus_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Rhinolophus_sinicus_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Rhinolophus_sinicus_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Rhinolophus_sinicus_start <- 5
Rhinolophus_sinicus_end <- 91
Rhinolophus_sinicus_hv_start <- 53

starts_and_stops[24,1] <- "Rhinolophus_sinicus"
starts_and_stops[24,2] <- Rhinolophus_sinicus_start
starts_and_stops[24,3] <- Rhinolophus_sinicus_end
starts_and_stops[24,4] <- Rhinolophus_sinicus_hv_start


##################### Dolphin #####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Tursiops_truncatus_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Tursiops_truncatus_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Tursiops_truncatus_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Tursiops_truncatus_start <- 3
Tursiops_truncatus_end <- 76
Tursiops_truncatus_hv_start <- 54

starts_and_stops[25,1] <- "Tursiops_truncatus"
starts_and_stops[25,2] <- Tursiops_truncatus_start
starts_and_stops[25,3] <- Tursiops_truncatus_end
starts_and_stops[25,4] <- Tursiops_truncatus_hv_start

##################### Alpaca #####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Vicugna_pacos_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Vicugna_pacos_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Vicugna_pacos_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Vicugna_pacos_start <- 4
Vicugna_pacos_end <- 91
Vicugna_pacos_hv_start <- 56

starts_and_stops[26,1] <- "Vicugna_pacos"
starts_and_stops[26,2] <- Vicugna_pacos_start
starts_and_stops[26,3] <- Vicugna_pacos_end
starts_and_stops[26,4] <- Vicugna_pacos_hv_start


##################### Brown Bat #####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Eptesicus_fuscus_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Eptesicus_fuscus_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Eptesicus_fuscus_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Eptesicus_fuscus_start <- 2
Eptesicus_fuscus_end <- 83
Eptesicus_fuscus_hv_start <- 52


starts_and_stops[27,1] <- "Eptesicus_fuscus"
starts_and_stops[27,2] <- Eptesicus_fuscus_start
starts_and_stops[27,3] <- Eptesicus_fuscus_end
starts_and_stops[27,4] <- Eptesicus_fuscus_hv_start

##################### Brandts Bat #####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Myotis_brandtii_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Myotis_brandtii_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Myotis_brandtii_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(110000,180000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Myotis_brandtii_start <- 4
Myotis_brandtii_end <- 89
Myotis_brandtii_hv_start <- 53

starts_and_stops[28,1] <- "Myotis_brandtii"
starts_and_stops[28,2] <- Myotis_brandtii_start
starts_and_stops[28,3] <- Myotis_brandtii_end
starts_and_stops[28,4] <- Myotis_brandtii_hv_start

##################### Roundleaf Bat #####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Hipposideros_armiger_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Hipposideros_armiger_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Hipposideros_armiger_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Hipposideros_armiger_start <- 4
Hipposideros_armiger_end <- 89
Hipposideros_armiger_hv_start <- 53

starts_and_stops[29,1] <- "Hipposideros_armiger"
starts_and_stops[29,2] <- Hipposideros_armiger_start
starts_and_stops[29,3] <- Hipposideros_armiger_end
starts_and_stops[29,4] <- Hipposideros_armiger_hv_start


##################### Human #####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Homo_sapiens_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Homo_sapiens_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Homo_sapiens_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Homo_sapiens_start <- 6
Homo_sapiens_end <- 114
Homo_sapiens_hv_start <- 56

starts_and_stops[30,1] <- "Homo_sapiens"
starts_and_stops[30,2] <- Homo_sapiens_start
starts_and_stops[30,3] <- Homo_sapiens_end
starts_and_stops[30,4] <- Homo_sapiens_hv_start

##################### Mouse #####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Mus_musculus_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
	exon_name <- as.character(pevk_lengths_and_ratios[i,1])
	split <- strsplit(exon_name, "_")
	location <- split[[1]][3]
	location <- strsplit(location, ":")
	match_start <- as.numeric(location[[1]][1])
	match_end <- as.numeric(location[[1]][2])
	new_starts <- append(new_starts, match_start)
	new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Mus_musculus_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Mus_musculus_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Mus_musculus_start <- 4
Mus_musculus_end <- 96
Mus_musculus_hv_start <- 53

starts_and_stops[31,1] <- "Mus_musculus"
starts_and_stops[31,2] <- Mus_musculus_start
starts_and_stops[31,3] <- Mus_musculus_end
starts_and_stops[31,4] <- Mus_musculus_hv_start

###################### Beaver ###################### 

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Castor_canadensis_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
  exon_name <- as.character(pevk_lengths_and_ratios[i,1])
  split <- strsplit(exon_name, "_")
  location <- split[[1]][3]
  location <- strsplit(location, ":")
  match_start <- as.numeric(location[[1]][1])
  match_end <- as.numeric(location[[1]][2])
  new_starts <- append(new_starts, match_start)
  new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Castor_canadensis_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Castor_canadensis_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Castor_canadensis_start <- 4
Castor_canadensis_end <- 93
Castor_canadensis_hv_start <- 51

starts_and_stops[32,1] <- "Castor_canadensis"
starts_and_stops[32,2] <- Castor_canadensis_start
starts_and_stops[32,3] <- Castor_canadensis_end
starts_and_stops[32,4] <- Castor_canadensis_hv_start


################### Orangutan ####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Pongo_abelii_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
  exon_name <- as.character(pevk_lengths_and_ratios[i,1])
  split <- strsplit(exon_name, "_")
  location <- split[[1]][3]
  location <- strsplit(location, ":")
  match_start <- as.numeric(location[[1]][1])
  match_end <- as.numeric(location[[1]][2])
  new_starts <- append(new_starts, match_start)
  new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Pongo_abelii_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Pongo_abelii_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Pongo_abelii_start <- 6
Pongo_abelii_end <- 105
Pongo_abelii_hv_start <- 54

starts_and_stops[33,1] <- "Pongo_abelii"
starts_and_stops[33,2] <- Pongo_abelii_start
starts_and_stops[33,3] <- Pongo_abelii_end
starts_and_stops[33,4] <- Pongo_abelii_hv_start

####################### Little brown bat #################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Myotis_lucifugus_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
  exon_name <- as.character(pevk_lengths_and_ratios[i,1])
  split <- strsplit(exon_name, "_")
  location <- split[[1]][3]
  location <- strsplit(location, ":")
  match_start <- as.numeric(location[[1]][1])
  match_end <- as.numeric(location[[1]][2])
  new_starts <- append(new_starts, match_start)
  new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Myotis_lucifugus_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Myotis_lucifugus_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Myotis_lucifugus_start <- 3
Myotis_lucifugus_end <- 85
Myotis_lucifugus_hv_start <- 45

starts_and_stops[34,1] <- "Myotis_lucifugus"
starts_and_stops[34,2] <- Myotis_lucifugus_start
starts_and_stops[34,3] <- Myotis_lucifugus_end
starts_and_stops[34,4] <- Myotis_lucifugus_hv_start


####################### Monk Seal #######################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Neomonachus_schauinslandi_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
  exon_name <- as.character(pevk_lengths_and_ratios[i,1])
  split <- strsplit(exon_name, "_")
  location <- split[[1]][3]
  location <- strsplit(location, ":")
  match_start <- as.numeric(location[[1]][1])
  match_end <- as.numeric(location[[1]][2])
  new_starts <- append(new_starts, match_start)
  new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Neomonachus_schauinslandi_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Neomonachus_schauinslandi_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()


Neomonachus_schauinslandi_start <- 3
Neomonachus_schauinslandi_end <- 80
Neomonachus_schauinslandi_hv_start <- 51

starts_and_stops[35,1] <- "Neomonachus_schauinslandi"
starts_and_stops[35,2] <- Neomonachus_schauinslandi_start
starts_and_stops[35,3] <- Neomonachus_schauinslandi_end
starts_and_stops[35,4] <- Neomonachus_schauinslandi_hv_start



############################### Walrus #####################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Odobenus_rosmarus_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
  exon_name <- as.character(pevk_lengths_and_ratios[i,1])
  split <- strsplit(exon_name, "_")
  location <- split[[1]][3]
  location <- strsplit(location, ":")
  match_start <- as.numeric(location[[1]][1])
  match_end <- as.numeric(location[[1]][2])
  new_starts <- append(new_starts, match_start)
  new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Odobenus_rosmarus_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Odobenus_rosmarus_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Odobenus_rosmarus_start <- 3
Odobenus_rosmarus_end <- 80
Odobenus_rosmarus_hv_start <- 53

starts_and_stops[36,1] <- "Odobenus_rosmarus"
starts_and_stops[36,2] <- Odobenus_rosmarus_start
starts_and_stops[36,3] <- Odobenus_rosmarus_end
starts_and_stops[36,4] <- Odobenus_rosmarus_hv_start

########################## Platypus ###################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Ornithorhynchus_anatinus_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
  exon_name <- as.character(pevk_lengths_and_ratios[i,1])
  split <- strsplit(exon_name, "_")
  location <- split[[1]][3]
  location <- strsplit(location, ":")
  match_start <- as.numeric(location[[1]][1])
  match_end <- as.numeric(location[[1]][2])
  new_starts <- append(new_starts, match_start)
  new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Ornithorhynchus_anatinus_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Ornithorhynchus_anatinus_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(120000,190000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Ornithorhynchus_anatinus_start <- 3
Ornithorhynchus_anatinus_end <- 93
Ornithorhynchus_anatinus_hv_start <- 52

starts_and_stops[37,1] <- "Ornithorhynchus_anatinus"
starts_and_stops[37,2] <- Ornithorhynchus_anatinus_start
starts_and_stops[37,3] <- Ornithorhynchus_anatinus_end
starts_and_stops[37,4] <- Ornithorhynchus_anatinus_hv_start


####################### Chimpanzee #################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Pan_troglodytes_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
  exon_name <- as.character(pevk_lengths_and_ratios[i,1])
  split <- strsplit(exon_name, "_")
  location <- split[[1]][3]
  location <- strsplit(location, ":")
  match_start <- as.numeric(location[[1]][1])
  match_end <- as.numeric(location[[1]][2])
  new_starts <- append(new_starts, match_start)
  new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Pan_troglodytes_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Pan_troglodytes_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Pan_troglodytes_start <- 6
Pan_troglodytes_end <- 92
Pan_troglodytes_hv_start <- 51

starts_and_stops[38,1] <- "Pan_troglodytes"
starts_and_stops[38,2] <- Pan_troglodytes_start
starts_and_stops[38,3] <- Pan_troglodytes_end
starts_and_stops[38,4] <- Pan_troglodytes_hv_start

####################### Panda ###################
# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Ailuropoda_melanoleuca_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
  exon_name <- as.character(pevk_lengths_and_ratios[i,1])
  split <- strsplit(exon_name, "_")
  location <- split[[1]][3]
  location <- strsplit(location, ":")
  match_start <- as.numeric(location[[1]][1])
  match_end <- as.numeric(location[[1]][2])
  new_starts <- append(new_starts, match_start)
  new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Ailuropoda_melanoleuca_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Ailuropoda_melanoleuca_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,170000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()

Ailuropoda_melanoleuca_start <- 5
Ailuropoda_melanoleuca_end <- 100
Ailuropoda_melanoleuca_hv_start <- 54

starts_and_stops[39,1] <- "Ailuropoda_melanoleuca"
starts_and_stops[39,2] <- Ailuropoda_melanoleuca_start
starts_and_stops[39,3] <- Ailuropoda_melanoleuca_end
starts_and_stops[39,4] <- Ailuropoda_melanoleuca_hv_start


####################### Davids Bat ###################
# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Myotis_davidii_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
  exon_name <- as.character(pevk_lengths_and_ratios[i,1])
  split <- strsplit(exon_name, "_")
  location <- split[[1]][3]
  location <- strsplit(location, ":")
  match_start <- as.numeric(location[[1]][1])
  match_end <- as.numeric(location[[1]][2])
  new_starts <- append(new_starts, match_start)
  new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Myotis_davidii_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Myotis_davidii_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()


Myotis_davidii_start <- 2
Myotis_davidii_end <- 86
Myotis_davidii_hv_start <- 51

starts_and_stops[40,1] <- "Myotis_davidii"
starts_and_stops[40,2] <- Myotis_davidii_start
starts_and_stops[40,3] <- Myotis_davidii_end
starts_and_stops[40,4] <- Myotis_davidii_hv_start

######################### Tasmanian Devil ################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Sarcophilus_harrisii_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
  exon_name <- as.character(pevk_lengths_and_ratios[i,1])
  split <- strsplit(exon_name, "_")
  location <- split[[1]][3]
  location <- strsplit(location, ":")
  match_start <- as.numeric(location[[1]][1])
  match_end <- as.numeric(location[[1]][2])
  new_starts <- append(new_starts, match_start)
  new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Sarcophilus_harrisii_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Sarcophilus_harrisii_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(128000,198000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()


Sarcophilus_harrisii_start <- 5
Sarcophilus_harrisii_end <- 92
Sarcophilus_harrisii_hv_start <- 48

starts_and_stops[41,1] <- "Sarcophilus_harrisii"
starts_and_stops[41,2] <- Sarcophilus_harrisii_start
starts_and_stops[41,3] <- Sarcophilus_harrisii_end
starts_and_stops[41,4] <- Sarcophilus_harrisii_hv_start

########################## Polar Bear ##################

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("Ursus_maritimus_test_lengths_and_ratios_10_0.54_12.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
  exon_name <- as.character(pevk_lengths_and_ratios[i,1])
  split <- strsplit(exon_name, "_")
  location <- split[[1]][3]
  location <- strsplit(location, ":")
  match_start <- as.numeric(location[[1]][1])
  match_end <- as.numeric(location[[1]][2])
  new_starts <- append(new_starts, match_start)
  new_ends <- append(new_ends, match_end)
}

# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))

# Assign colors to each PEVK ratio
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

# Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("Ursus_maritimus_test_locations_10_0.54_12.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Create rectangle objects for each exon coordinate

my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

pdf(file="Ursus_maritimus_exon_plot.pdf", width=50, height=10, pointsize=22)

op <- par(bg = "white", mar=c(.05,.05,.05,.05), bty="n")
plot(c(105000,175000), c(1.7, 3.3), type="n", xlab="", ylab="", main="", yaxt="n", xaxt="n")

rect(my_starts,2.1,my_ends,2.9, col=test_exon_cols, border=NA)
#abline(h=2.5) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(2.5), labels=c("PEVK Finder Exons"))

dev.off()



Ursus_maritimus_start <- 3
Ursus_maritimus_end <- 93
Ursus_maritimus_hv_start <- 51

starts_and_stops[42,1] <- "Ursus_maritimus"
starts_and_stops[42,2] <- Ursus_maritimus_start
starts_and_stops[42,3] <- Ursus_maritimus_end
starts_and_stops[42,4] <- Ursus_maritimus_hv_start



#************************************************************************************************************
################################ PEVK Characterization #############################

library(ape)
library(phytools)
library(geiger)
library(nlme)
library(plotrix)

# Read in mammal tree from timetree.org
#mammal.tree <- read.tree(url("http://schmitzlab.info/mammalia_species.tre"))
mammal.tree <- read.tree("~/Desktop/titin_project/exon_alignments_named/mammalia_species.tre")

# Create list of the names/lengths ad locations files for each species
lengths_and_ratios = list.files(pattern="test_lengths_and_ratios")
locations = list.files(pattern="test_locations")

# Read files into R
for (i in 1:length(lengths_and_ratios)) assign(lengths_and_ratios[i], read.csv(lengths_and_ratios[i]))
for (i in 1:length(locations)) assign(locations[i], read.csv(locations[i]))


# Create data frame that will hold info on PEVK composition in each species 
# Columns: taxon, [region length], # PEVK, # large exons, # small exons,  
col_names <- c("taxon", "region_length", "total_pevk", "pevk_density","small_exons","large_exons","mean_ratio","mean_error","mean_sd","mean_hv","hv_exons","conserved_exons")
pevk_data <- data.frame(matrix(ncol=12, nrow=length(lengths_and_ratios)))
colnames(pevk_data) <- col_names

# Fill first column of table with species names
for (i in 1:length(lengths_and_ratios)) {
	file_name <- as.character(lengths_and_ratios[i])
	split <- strsplit(file_name, "_")
	names <- c(split[[1]][1],split[[1]][2])
	species_name <- paste(names, collapse='_')
	
	# Add taxon name to pevk_data
	#if (species_name != "Sorex_araneus" && species_name != "Physeter_catodon" && species_name != "Leptonychotes_weddellii") {
	pevk_data[i,1] <- species_name
}

# Make the following species names compatible with the tree
pevk_data[24,1] <- "Monachus_schauinslandi"
pevk_data[41,1] <- "Tupaia_belangeri"

#if (species_name != "Sorex_araneus" && species_name != "Physeter_catodon" && species_name != "Leptonychotes_weddellii") {
# remove Sorex_araneus, Physeter_catodon, Leptonychotes_weddellii (little to no data for these species)
pevk_data <- rbind(pevk_data[1:15,], pevk_data[17:31,], pevk_data[33:37,], pevk_data[39:44,])

# Set row names to taxon names
rownames(pevk_data) <- pevk_data$taxon

# Match tree to data
compare <- treedata(mammal.tree, pevk_data, sort=TRUE)
pevk_data <- as.data.frame(compare$data)
mammal.tree <- ladderize(compare$phy)

# Scale down branch lengths
mammal.tree$edge.length <- mammal.tree$edge.length*0.25

# Plot the plain mammal tree
#pdf(file="mammal_tree.pdf", width=15, height=30, pointsize=30)
#plot(mammal.tree)
#dev.off()


# Switch those species back to the original names
pevk_data$taxon <- as.character(pevk_data$taxon)
pevk_data[8,1] <- "Neomonachus_schauinslandi"
pevk_data[37,1] <- "Tupaia_chinensis"
rownames(pevk_data)[8] <- "Neomonachus_schauinslandi"
rownames(pevk_data)[37] <- "Tupaia_chinensis"

# Initialize data frame
pevk_data$region_length <- integer(nrow(pevk_data))
pevk_data$total_pevk <- integer(nrow(pevk_data))
pevk_data$pevk_density <- integer(nrow(pevk_data))
pevk_data$large_exons <- integer(nrow(pevk_data))
pevk_data$small_exons <- integer(nrow(pevk_data))
pevk_data$mean_ratio <- integer(nrow(pevk_data))
pevk_data$mean_error <- integer(nrow(pevk_data))
pevk_data$mean_sd <- integer(nrow(pevk_data))
pevk_data$mean_hv <- integer(nrow(pevk_data))
pevk_data$hv_exons <- integer(nrow(pevk_data))
pevk_data$conserved_exons <- integer(nrow(pevk_data))
#################### Get exon count per species  and hy/conserved count #######
pevk_exon_counts = c()
pevk_conserved_counts = c()
pevk_hypervariable_counts = c()

conserved_exon_lengths = c()
conserved_exon_ratios = c()
hv_exon_lengths = c()
hv_exon_ratios = c()

# Collect data on each species
for (z in 1:nrow(pevk_data)) {
	
	# Extract species name
	taxon <- toString(pevk_data$taxon[z])
	
	# Get the corresponding lengths and ratios data file
	for (j in 1:length(lengths_and_ratios)) {
		if (grepl(taxon,lengths_and_ratios[j]) == TRUE) {
			LAR = read.csv(lengths_and_ratios[j], header=F)
		}
	}
	
	# Get corresponding PEVK ranges
	for (l in 1:nrow(starts_and_stops)) {
		if (grepl(taxon,starts_and_stops$taxon[l]) == TRUE) {
			start_index <- starts_and_stops$start_index[l]
			stop_index <- starts_and_stops$stop_index[l]
			hv_index <- starts_and_stops$hv_start[l]
		}
	}
		
	###################### Transform LAR #######################
	# Purpose: be able to access data on lengths and ratios of PEVK exons
	LAR <- t(LAR)
	LAR <- data.frame(LAR)
	rownames(LAR) <- 1:nrow(LAR)
	colnames(LAR) <- c("name","length","percent_pevk")

	# Extract locations and append to LAR
	new_starts = c()
	new_ends = c()
	for (i in 1:nrow(LAR))	{
		exon_name <- as.character(LAR[i,1])
		split <- strsplit(exon_name, "_")
		location <- split[[1]][3]
		location <- strsplit(location, ":")
		match_start <- as.numeric(location[[1]][1])
		match_end <- as.numeric(location[[1]][2])
		new_starts <- append(new_starts, match_start)
		new_ends <- append(new_ends, match_end)
	}
	
	# Add starts and ends to the data frame
	LAR$start <- new_starts
	LAR$end <- new_ends
	LAR <- LAR[order(LAR$start),]
	rownames(LAR) <- 1:nrow(LAR)
###### Changed
	# Extract best-guess PEVK region from data frame
	LAR_range <- LAR[1:nrow(LAR),]
  
	pevk_exon_counts = c(pevk_exon_counts, nrow(LAR))
	pevk_conserved_counts = c(pevk_conserved_counts, hv_index)
	pevk_hypervariable_counts = c(pevk_hypervariable_counts, nrow(LAR) - hv_index)
	
	################## Add data to pevk_data ############ 
	
	# 1. Total span of the PEVK region
	region_start <- LAR$start[start_index]
	region_end <- LAR$end[stop_index]
	region_length <- region_end - region_start
	
	pevk_data$region_length[z] <- region_length
	
	# 2. Total number of exons
	total_pevk <- stop_index - start_index + 1
	
	pevk_data$total_pevk[z] <- total_pevk

	# 3. Number of exons: region length
	pevk_density <- total_pevk/region_length	
	
	pevk_data$pevk_density[z] <- pevk_density
	
	# 4. Mean PEVK/SE
	pevk_ratios <- as.vector(as.numeric(as.character(LAR$percent_pevk)))
	mean_ratio <- mean(pevk_ratios)
	se_ratio <- std.error(pevk_ratios)
	sd_ratio <- sd(pevk_ratios)
	
	pevk_data$mean_ratio[z] <- mean_ratio
	pevk_data$mean_error[z] <- se_ratio
	pevk_data$mean_sd[z] <- sd_ratio
	
	# 5. Mean PEVK for HV region
	hv_ratios <- as.vector(as.numeric(as.character(LAR$percent_pevk[hv_index:stop_index])))
	mean_ratio <- mean(hv_ratios)
	pevk_data$mean_hv[z] <- mean_ratio
	
	# 6b. Number of exons in conserved region
	conserved_exon_number <- hv_index - start_index + 1
	pevk_data$conserved_exons[z] <- conserved_exon_number

	# 6a. Number of exons in hv region
	hv_exon_number <- (stop_index - start_index + 1) - conserved_exon_number
	pevk_data$hv_exons[z] <- hv_exon_number
	
	# 7. Length of exons in conserved
	conserved_exon_lengths = c(conserved_exon_lengths, as.numeric(as.character(LAR$length[start_index:hv_index-1])))
	
	# 8. Length of exons in hv
	conserved_exon_ratios = c(conserved_exon_ratios, as.numeric(as.character(LAR$percent_pevk[start_index:hv_index-1])))
	
	# 9. Ratio of exons in conserved
	hv_exon_lengths = c(hv_exon_lengths, as.numeric(as.character(LAR$length[hv_index:stop_index])))
	
	# 10. Ratio of exons in hv
	hv_exon_ratios = c(hv_exon_ratios, as.numeric(as.character(LAR$percent_pevk[hv_index:stop_index])))
	
}

######## Compare phylogeny and data
#is_tip <- mammal.tree$edge[,2] <= length(mammal.tree$tip.label)
#ordered_tips <- mammal.tree$edge[is_tip,2]
#mammal.tree$tip.label <- mammal.tree$tip.label[rev(ordered_tips)]

# Match tree to data
#compare <- treedata(mammal.tree, pevk_data, sort=TRUE)
#pevk_data <- as.data.frame(compare$data)
#mammal.tree <- compare$phy


# Visualizing data
density.data <- (pevk_data$pevk_density)*1500
# make narrow margins
par(mar=c(2,0,0,1))
# plot
#plot(mammal.tree)
#tiplabels(pch=20, col="black", cex=density.data)

length.data <- (pevk_data$region_length)/30000
# make narrow margins
par(mar=c(2,0,0,1))
# plot
#plot(mammal.tree)
#tiplabels(pch=20, col="black", cex=length.data)

number.data <- (pevk_data$total_pevk)/50
# make narrow margins
par(mar=c(2,0,0,1))
# plot
#plot(mammal.tree)
#tiplabels(pch=20, col="black", cex=number.data)

#### Plot mean PEVk ratios #######
pdf(file="mammal_tree_with_mean.pdf", width=18, height=30, pointsize=30)
# Mean and SE
mean <- pevk_data$mean_ratio
# split the plotting panel
layout(matrix(c(1,2),1,2), c(0.8, 0.2))
#make narrow margins
par(mar=c(2,0,0,1))
# plot tree
plot(mammal.tree, show.tip.label=FALSE)
#plot barplot (with trait ordered by tip label)
names(mean) <- mammal.tree$tip.label
par(mar=c(2,0,0,1))
barplot(mean[mammal.tree$tip.label], xlim=c(0.7,0.85), horiz=TRUE, space=0, ylim=c(1,length(mammal.tree$tip.label))-0.5, names="", col="steel blue")
#barploterrbar()
dev.off()

##### Plot region length (number of exons)
pdf(file="mammal_tree_with_lengths.pdf", width=18, height=30, pointsize=30)
# Mean and SE
lengths <- as.vector(as.numeric(as.character(pevk_data$total_pevk)))
# split the plotting panel
layout(matrix(c(1,2),1,2), c(0.8, 0.2))
#make narrow margins
par(mar=c(2,0,0,1))
# plot tree
plot(mammal.tree)
#plot barplot (with trait ordered by tip label)
names(lengths) <- mammal.tree$tip.label
par(mar=c(2,0,0,1))
barplot(lengths, xlim=c(70,110), horiz=TRUE, space=0, ylim=c(1,length(mammal.tree$tip.label))-0.5, names="", col="pink")
#barploterrbar()
dev.off()



#################### Start figure code ###################

# 1. Overall PEVK length
pevk_lengths <- as.vector(as.numeric(as.character(pevk_data$total_pevk)))
pevk_lengths_mean <- as.integer(mean(pevk_lengths))
pevk_lengths_mean_differences <- vector(length=29)

for (i in 1:length(pevk_lengths)) {
	diff <- pevk_lengths[i] - pevk_lengths_mean
	pevk_lengths_mean_differences[i] <- diff
}
names(pevk_lengths_mean_differences) <- mammal.tree$tip.label

# Reorder the vector so it matches tip labels
is_tip <- mammal.tree$edge[,2] <= length(mammal.tree$tip.label)
ordered_tips <- mammal.tree$edge[is_tip,2]
tip.order <- mammal.tree$tip.label[mammal.tree$edge[mammal.tree$edge[,2] <= Ntip(mammal.tree),2]]

lengths_final <- vector(length=29)
for (i in 1:length(ordered_tips)) {
	lengths_final[i] <- pevk_lengths_mean_differences[ordered_tips[i]]
}
names(lengths_final) <- tip.order

##### Plot region length (number of exons)
pdf(file="mammal_tree_with_lengths.pdf", width=18, height=30, pointsize=30)
layout(matrix(c(1,2),1,2), c(0.8, 0.2))
#make narrow margins
par(mar=c(2,0,0,1))
# plot tree
plot(mammal.tree)
#plot barplot (with trait ordered by tip label)
par(mar=c(2,0,0,1))
barplot(lengths_final, xlim=c(-20,20), horiz=TRUE, space=0.2, ylim=c(1,length(mammal.tree$tip.label)+6)-0.5, names="", col="darkgrey", border=FALSE)
#barploterrbar()
dev.off()





############################# HV region mean lengths ###############################

# 1. Overall PEVK length
hv_lengths <- as.vector(as.numeric(as.character(pevk_data$hv_exons)))
hv_lengths_mean <- as.integer(mean(hv_lengths))
hv_lengths_mean_differences <- vector(length=29)

for (i in 1:length(hv_lengths)) {
	diff <- hv_lengths[i] - hv_lengths_mean
	hv_lengths_mean_differences[i] <- diff
}






############################# PEVK-C/N Exon Number/Ratios (7/13/18) ##############
pevk_lengths_updated <- as.numeric(as.character(pevk_data$total_pevk))
hv_lengths_updated <- as.numeric(as.character(pevk_data$hv_exons))
conserved_lengths_updated <- as.numeric(as.character(pevk_data$conserved_exons))
#conserved_lengths_updated <- c()
#for (i in 1:length(pevk_lengths_updated)){
#  conserved_length = pevk_lengths_updated[i] - hv_lengths_updated[i]
#  conserved_lengths_updated <- c(conserved_lengths_updated, conserved_length)
#}
# Region lengths:
conserved_lengths_mean <- mean(conserved_lengths_updated)
conserved_lengths_se <- std.error(conserved_lengths_updated)

hv_lengths_mean <- mean(hv_lengths_updated)
hv_lengths_se <- std.error(hv_lengths_updated)


# Per-Exon lengths:
conserved_exon_lengths_mean = mean(conserved_exon_lengths)
conserved_exon_lengths_se = std.error(conserved_exon_lengths)

hv_exon_lengths_mean = mean(hv_exon_lengths)
hv_exon_lengths_se = std.error(hv_exon_lengths)

# Per-Exon ratios:

conserved_exon_ratios_mean = mean(conserved_exon_ratios)
conserved_exon_ratios_se = std.error(conserved_exon_ratios)

hv_exon_ratios_mean = mean(hv_exon_ratios)
hv_exon_ratios_se = std.error(hv_exon_ratios)


# Stat tests: overall length comparison
conserved_hv_exon_lengths_comparison = t.test(conserved_exon_lengths, hv_exon_lengths, paired = FALSE)
conserved_hv_exon_ratios_comparison = t.test(hv_exon_ratios, conserved_exon_ratios, paired = FALSE)

# Variance of lengths
conserved_hv_exon_lengths_var_test = var.test(conserved_exon_lengths, hv_exon_lengths)


#################################################################






names(hv_lengths_mean_differences) <- mammal.tree$tip.label

# Reorder the vector so it matches tip labels
is_tip <- mammal.tree$edge[,2] <= length(mammal.tree$tip.label)
ordered_tips <- mammal.tree$edge[is_tip,2]
tip.order <- mammal.tree$tip.label[mammal.tree$edge[mammal.tree$edge[,2] <= Ntip(mammal.tree),2]]

hv_lengths_final <- vector(length=29)
for (i in 1:length(ordered_tips)) {
	hv_lengths_final[i] <- hv_lengths_mean_differences[ordered_tips[i]]
}
names(hv_lengths_final) <- tip.order

##### Plot region length (number of exons)
pdf(file="mammal_tree_with_hv_lengths.pdf", width=18, height=30, pointsize=30)
layout(matrix(c(1,2),1,2), c(0.8, 0.2))
#make narrow margins
par(mar=c(2,0,0,1))
# plot tree
plot(mammal.tree)
#plot barplot (with trait ordered by tip label)
par(mar=c(2,0,0,1))
barplot(hv_lengths_final, xlim=c(-20,20), horiz=TRUE, space=0.2, ylim=c(1,length(mammal.tree$tip.label)+6)-0.5, names="", col="black", border=FALSE)
#barploterrbar()
dev.off()

############################# HV region ratios ##########################################

# 1. Overall PEVK length
hv_ratios <- as.vector(as.numeric(as.character(pevk_data$mean_hv)))
hv_ratios_mean <- mean(hv_ratios)
hv_ratios_mean_differences <- vector(length=29)

for (i in 1:length(hv_ratios)) {
	diff <- hv_ratios[i] - hv_ratios_mean
	hv_ratios_mean_differences[i] <- diff
}
names(hv_ratios_mean_differences) <- mammal.tree$tip.label

# Reorder the vector so it matches tip labels
is_tip <- mammal.tree$edge[,2] <= length(mammal.tree$tip.label)
ordered_tips <- mammal.tree$edge[is_tip,2]
tip.order <- mammal.tree$tip.label[mammal.tree$edge[mammal.tree$edge[,2] <= Ntip(mammal.tree),2]]

hv_ratios_final <- vector(length=29)
for (i in 1:length(ordered_tips)) {
	hv_ratios_final[i] <- hv_ratios_mean_differences[ordered_tips[i]]
}
names(hv_ratios_final) <- tip.order

##### Plot region length (number of exons)
pdf(file="mammal_tree_with_hv_ratios.pdf", width=18, height=30, pointsize=30)
layout(matrix(c(1,2),1,2), c(0.8, 0.2))
#make narrow margins
par(mar=c(2,0,0,1))
# plot tree
plot(mammal.tree)
#plot barplot (with trait ordered by tip label)
par(mar=c(2,0,0,1))
barplot(hv_ratios_final, xlim=c(-0.06,0.06), horiz=TRUE, space=0.2, ylim=c(1,length(mammal.tree$tip.label)+6)-0.5, names="", col="grey", border=FALSE)
#barploterrbar()
dev.off()





