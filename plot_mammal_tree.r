library(ape)
library(phytools)
library(geiger)
library(nlme)
library(plotrix)

# Read in mammal tree from timetree.org
#mammal.tree <- read.tree(url("http://schmitzlab.info/mammalia_species.tre"))
mammal.tree <- read.tree("/Users/kathleenmuenzen/Desktop/titin_project/exon_alignments_named/mammalia_species.tre")

# Create list of the names/lengths ad locations files for each species
setwd('/Users/kathleenmuenzen/Desktop/titin_project/pevk_mammals')
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
#pevk_data <- rbind(pevk_data[1:15,], pevk_data[17:21,], pevk_data[24:31,], pevk_data[33:37,], pevk_data[39:44,])
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
#pdf(file="mammal_tree.pdf", width=11, height=200, pointsize=30)
#plot(mammal.tree)
#dev.off()

####### Plot repeats with tree ########
# This script creates a phylogenetic tree with information about tandem repeats and LINEs plotted alongside the tree
# Kathleen Muenzen, 12/21/18, kmuenzen@uw.edu

# Read in repeat data
repeat_data <- read.csv("/Users/kathleenmuenzen/Desktop/final_ttn_submission/dup_trip_str_classifications.csv",stringsAsFactors = FALSE)
repeat_data[23,1] <- "Monachus_schauinslandi"
repeat_data[38,1] <- "Tupaia_belangeri"
rownames(repeat_data) <- repeat_data$species

repeat_data[match(mammal.tree$tip.label, repeat_data$species),]

# Create clean vector for plotting
species_order = c("Ornithorhynchus_anatinus","Phascolarctos_cinereus","Sarcophilus_harrisii","Dasypus_novemcinctus",
                  "Loxodonta_africana","Trichechus_manatus","Tupaia_belangeri","Galeopterus_variegatus",
                  "Pongo_abelii","Homo_sapiens","Pan_troglodytes","Ochotona_princeps","Oryctolagus_cuniculus",
                  "Ictidomys_tridecemlineatus","Castor_canadensis","Rattus_norvegicus","Mus_pahari","Mus_musculus",
                  "Condylura_cristata","Erinaceus_europaeus","Pteropus_vampyrus","Rhinolophus_sinicus",
                  "Hipposideros_armiger","Eptesicus_fuscus","Myotis_davidii","Myotis_lucifugus","Myotis_brandtii",
                  "Vicugna_pacos","Sus_scrofa","Bos_taurus","Orcinus_orca","Tursiops_truncatus","Ceratotherium_simum",
                  "Equus_caballus","Manis_javanica","Acinonyx_jubatus","Felis_catus","Ursus_maritimus",
                  "Ailuropoda_melanoleuca","Monachus_schauinslandi","Odobenus_rosmarus")

# Match data frame order with species order
repeat_data = repeat_data[match(species_order, repeat_data$species),]

# Create clean word vector for plotting repeat data
repeat_text = c()
for (p in 1:nrow(repeat_data)){
  repeat_type = repeat_data$repeat_type[p]
  repeat_type_subbed = gsub(";",", ",repeat_type)
  if (is.na(repeat_type_subbed)){
    repeat_text = c(repeat_text,' ')
  }
  else{
    repeat_text = c(repeat_text,repeat_type_subbed)
  }
}
repeat_text = rev(repeat_text)

# Create text vector for LINE analysis
human_lines <- c(rep(' ',1),'Identical Seqs x 2',rep(' ',6),'LINE-1 x 1','LINE-1 x 2','LINE-1 x 1',rep(' ',30))
human_lines = rev(human_lines)

pdf(file="~/Desktop/mammal_tree.pdf", width=15, height=25, pointsize=20)
plot(mammal.tree,y.lim=50,x.lim=150)
text(repeat_text,x=rep.int(100,41),y=seq(from=1,to=42,by=1))
text(human_lines,x=rep.int(130,41),y=seq(from=1,to=42,by=1))
dev.off()






