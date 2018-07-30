library(ape)
library(phytools)
library(geiger)
library(nlme)
library(plotrix)

# Read in mammal tree from timetree.org
#mammal.tree <- read.tree(url("http://schmitzlab.info/mammalia_species.tre"))
mammal.tree <- read.tree("mammalia_species.tre")

# Create list of the names/lengths ad locations files for each species
lengths_and_ratios = list.files(pattern="test_lengths_and_ratios")
locations = list.files(pattern="test_locations")

# Read files into R
for (i in 1:length(lengths_and_ratios)) assign(lengths_and_ratios[i], read.csv(lengths_and_ratios[i]))
for (i in 1:length(locations)) assign(locations[i], read.csv(locations[i]))

# Create data frame that will hold info on PEVK composition in each species 
# Columns: taxon, [region length], # PEVK, # large exons, # small exons,  
col_names <- c("taxon", "region_length", "total_pevk", "pevk_density","small_exons","large_exons","mean_ratio","mean_error","mean_sd","mean_hv","hv_exons")
pevk_data <- data.frame(matrix(ncol=11, nrow=length(lengths_and_ratios)))
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
pevk_data[24,1] <- "Monachus_schauinslandi"
pevk_data[41,1] <- "Tupaia_belangeri"

#if (species_name != "Sorex_araneus" && species_name != "Physeter_catodon" && species_name != "Leptonychotes_weddellii") {
pevk_data <- rbind(pevk_data[1:15,], pevk_data[17:21,], pevk_data[24:31,], pevk_data[33:37,], pevk_data[39:44,])

# Set row names to taxon names
rownames(pevk_data) <- pevk_data$taxon

# Match tree to data
compare <- treedata(mammal.tree, pevk_data, sort=TRUE)
pevk_data <- as.data.frame(compare$data)
mammal.tree <- ladderize(compare$phy)

# Scale down branch lengths
mammal.tree$edge.length <- mammal.tree$edge.length*0.25

# Plot the plain mammal tree
pdf(file="mammal_tree.pdf", width=15, height=30, pointsize=30)
plot(mammal.tree)
dev.off()