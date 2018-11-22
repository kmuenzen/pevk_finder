# This code reads the csv files generated by MEGA divergence calculations and creates divergence heatmaps
# and calculates overall divergence for each quadrant

library(RColorBrewer)
library(seqinr)
library(ggplot2)
library(reshape2)
library(plyr)
library(scales)

#library(geom_tile)
# Species pevk bundaries
species_pevk_boundaries <- as.data.frame(read.csv("~/Desktop/new_materials/species_pevk_boundaries.csv",sep=','))
mean(species_pevk_boundaries$pevkcb_length)
sd(species_pevk_boundaries$pevkcb_length)

setwd("~/Desktop/titin_project/megacc_tests")
# Loop through all folders in the directory
all_folders <- list.dirs(recursive=FALSE)

# Create data frame that will hold info on quadrant stats
quadrant_stats <- data.frame(matrix(nrow=length(all_folders)-1, ncol=4))
colnames(quadrant_stats) <- c("quad1_mean", "quad2_mean", "quad3_mean", "quad4_mean")

# Create data frame that will hold info for all 9 segments
quad_9_stats <- data.frame(matrix(nrow=length(all_folders)-1, ncol=9))
colnames(quad_9_stats) <- c("quadi_mean", "quadii_mean", "quadiii_mean", "quadiv_mean", "quadv_mean", "quadvi_mean", "quadvii_mean", "quadviii_mean", "quadix_mean")

# Loop through all divergence folders
for (i in 1:(length(all_folders)-1)){
  # Get name of current folder
  current_folder = all_folders[i]
  
  # Set working directroy
  setwd(current_folder)
  
  # Get name of current species
  species_name <- sub("./","",current_folder)
  species_name <- sub("_divergence","", species_name)
  
  # Read in CSV files
  divergence_files <- list.files(pattern=".csv")
  
  # Create a data frame for file names to organize in ascending order
  names <- data.frame(matrix(ncol=3, nrow=1))
  colnames(names) <- c("frame", "start", "end")
  
  counter = 0
  start_vec = c()
  for (j in (1:length(divergence_files))){
    split=strsplit(divergence_files[j], "_")
    start_end=strsplit(split[[1]][3], ":")
    frame=split[[1]][2]
    start=as.numeric(start_end[[1]][1])
    end=as.numeric(start_end[[1]][2])
    
    # If the current exon isn't already in the matrix...
    if (is.element(start,start_vec) == FALSE) {
      counter = counter + 1
      start_vec <- c(start_vec, start)
      names[counter, "frame"] = frame
      names[counter, "start"] = start
      names[counter, "end"] = end
    }
  }
  
  names_ordered <- names[order(names$start),]
  rownames(names_ordered) <- c(1:nrow(names_ordered))
  
  # Create vector to hold ordered names
  names_vec <- c()
  for (k in 1:nrow(names_ordered)) {
    current_frame <- as.character(names_ordered[k,"frame"])
    current_start <- as.character(names_ordered[k,"start"])
    current_end <- as.character(names_ordered[k,"end"])
    current_exon <- paste(current_frame, current_start, current_end, sep="_")
    names_vec <- c(names_vec, current_exon)
  }
  
  # Initialize matrix that will hold divergence data
  # option 1
  #heatmap_matrix <- data.frame(matrix(ncol=length(names_vec)+1,nrow=length(names_vec)))
  #colnames(heatmap_matrix) <- c("exon_id",names_vec)
  #rownames(heatmap_matrix) <- c(1:length(names_vec))
  #heatmap_matrix$exon_id <- names_vec
  
  # option 2
  heatmap_matrix <- data.frame(matrix(ncol=length(names_vec),nrow=length(names_vec)))
  colnames(heatmap_matrix) <- names_vec
  rownames(heatmap_matrix) <- names_vec
  
  # Loop through csv files and fill in heatmap matrix
  for (m in 1:length(divergence_files)){
    current_file <- as.data.frame(read.csv(divergence_files[m]))
    div_and_se <- as.character(current_file[[1]])
    div_and_se <- strsplit(div_and_se, " ")
    if (div_and_se[[1]][3] != ""){
      div <- as.numeric(div_and_se[[1]][3])
    }
    if (div_and_se[[1]][3] == ""){
      div <- as.numeric(div_and_se[[1]][4])
    }
    if (div_and_se[[1]][3] == "" && div_and_se[[1]][4] == ""){
      div <- as.numeric(div_and_se[[1]][5])
    }
    
    # Get exon info
    name_split=strsplit(divergence_files[m], "_")
    
    exon1_start_end=strsplit(name_split[[1]][3], ":")
    exon2_start_end=strsplit(name_split[[1]][6], ":")
    
    exon1_frame=name_split[[1]][2]
    exon2_frame=name_split[[1]][5]
    
    exon1_start=exon1_start_end[[1]][1]
    exon2_start=exon2_start_end[[1]][1]
    
    exon1_end=exon1_start_end[[1]][2]
    exon2_end=exon2_start_end[[1]][2]
    
    # Combine info
    exon1 <- paste(exon1_frame, exon1_start, exon1_end, sep="_")
    exon2 <- paste(exon2_frame, exon2_start, exon2_end, sep="_")
    
    # Get location of div value in data frame
    row_location <- match(exon1, names_vec)
    col_location <- match(exon2, names_vec)
    
    # Add div value to data frame
    heatmap_matrix[row_location, col_location] <- div
  }
  
  # Prepare data frame for plotting
  #heatmap_matrix_melted <- melt(heatmap_matrix, measure.vars = names_vec)
  #heatmap_matrix_rescaled <- ddply(heatmap_matrix_melted, .(exon_id), rescale = rescale(value))
  
  # Convert to matrix
  heatmap_matrix <- data.matrix(heatmap_matrix)
  
  ##################### Get quadrant stats ##############
  # Add species name to data frame
  rownames(quadrant_stats)[i] <- species_name
  # Get species pevk boundaries
  for (row in 1:nrow(species_pevk_boundaries)){
    current_species_name = species_pevk_boundaries$taxon[row]
    if (species_name == current_species_name){
      pevkn_length = species_pevk_boundaries$pevkn_length[row]
      pevkcb_length = species_pevk_boundaries$pevkcb_length[row]
    }
  }
  
  pevk_length <- length(names_vec)
  quad_1 <- heatmap_matrix[1:pevkn_length,][,(pevkn_length+1):pevk_length]
  quad_2 <- heatmap_matrix[(pevkn_length+1):pevk_length,][,(pevkn_length+1):pevk_length]
  quad_3 <- heatmap_matrix[1:pevkn_length,][,1:pevkn_length]
  quad_4 <- heatmap_matrix[(pevkn_length+1):pevk_length,][,1:pevkn_length]
  
  # Separate substitutions matrix into 9 seqments
  pevk_n_end <- pevkn_length
  pevk_ca_start <- pevkn_length+1
  pevk_ca_end <- nrow(heatmap_matrix)-pevkcb_length
  pevk_cb_start <- pevk_ca_end + 1
  pevk_cb_end <- nrow(heatmap_matrix)
  
  quad_ix <- heatmap_matrix[1:pevk_n_end,][,pevk_cb_start:pevk_cb_end]
  quad_vi <- heatmap_matrix[pevk_ca_start:pevk_ca_end,][,pevk_cb_start:pevk_cb_end]
  quad_iii <- heatmap_matrix[pevk_cb_start:pevk_cb_end,][,pevk_cb_start:pevk_cb_end]
  quad_viii <- heatmap_matrix[1:pevk_n_end,][,pevk_ca_start:pevk_ca_end]
  quad_v <- heatmap_matrix[pevk_ca_start:pevk_ca_end,][,pevk_ca_start:pevk_ca_end]
  quad_ii <- heatmap_matrix[pevk_cb_start:pevk_cb_end,][,pevk_ca_start:pevk_ca_end]
  quad_vii <- heatmap_matrix[1:pevk_n_end,][,1:pevk_n_end]
  quad_iv <- heatmap_matrix[pevk_ca_start:pevk_ca_end,][,1:pevk_n_end]
  quad_i <- heatmap_matrix[pevk_cb_start:pevk_cb_end,][,1:pevk_n_end]
  
  # Add means to data frame
  quadrant_stats$quad1_mean[i] <- mean(quad_1)
  quadrant_stats$quad2_mean[i] <- mean(quad_2)
  quadrant_stats$quad3_mean[i] <- mean(quad_3)
  quadrant_stats$quad4_mean[i] <- mean(quad_4)
  
  
  # Add means to larger data frame
  quad_9_stats$quadi_mean[i] <- mean(quad_i)
  quad_9_stats$quadii_mean[i] <- mean(quad_ii)
  quad_9_stats$quadiii_mean[i] <- mean(quad_iii)
  quad_9_stats$quadiv_mean[i] <- mean(quad_iv)
  quad_9_stats$quadv_mean[i] <- mean(quad_v)
  quad_9_stats$quadvi_mean[i] <- mean(quad_vi)
  quad_9_stats$quadvii_mean[i] <- mean(quad_vii)
  quad_9_stats$quadviii_mean[i] <- mean(quad_viii)
  quad_9_stats$quadix_mean[i] <- mean(quad_ix)
  
  # Scale matrix
  matrix_test <- apply(heatmap_matrix, MARGIN=2, FUN=function(X) (X - min(X))/diff(range(X)))
  # Color palette
  pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
  
  # Make heatmap
  pdf(paste(species_name,"_divergence.pdf",sep=''))
  test_heatmap <- heatmap(matrix_test, Rowv=NA, Colv=NA, col=pal, scale="column", labRow=c(1:nrow(names_ordered)), labCol=c(1:nrow(names_ordered)), cexRow=0.6, cexCol=0.6, xlab=paste(species_name," PEVK Finder Exons", sep=''), ylab=paste(species_name," PEVK Finder Exons", sep=''), margins=c(5,5))
  dev.off()
  #(ggplot_heatmap <- ggplot(heatmap_matrix_rescaled, aes(rev(variable), exon_id)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient(low = "steelblue", high = "white"))
  setwd("~/Desktop/titin_project/megacc_tests")
}

# Process quadrant stats and make figure

# Make separate matrix for mean and stdev values:
overall_quad_matrix <- data.frame(matrix(nrow=2, ncol=4))
colnames(overall_quad_matrix) <- c("I","II","III","IV")
rownames(overall_quad_matrix) <- c("mean","stdev")

# one way ANOVA
subs <- c(quadrant_stats$quad1_mean, quadrant_stats$quad2_mean, quadrant_stats$quad3_mean, quadrant_stats$quad4_mean)
box <- c(rep("quad1_mean",41), rep("quad2_mean",41), rep("quad3_mean",41), rep("quad4_mean",41))
box_means <- data.frame(subs,box)
plot(subs ~ box, data=box_means)
results <- aov(subs ~ box, data=box_means)
summary(results)
tukey.test <- TukeyHSD(results)
tukey.test

means <- tapply(subs,box,mean)
se <- tapply(subs,box,std.error)
means_and_ses <- data.frame(means,se)
means_and_ses

pdf(file="quadrant_comparison_updated.pdf", width=7, height=8, pointsize=22)
par(par(mar = c(5, 6, 4, 5) + 0.1))
barCenters <- barplot(height=means_and_ses$means,ylim=c(0,1.5), space = 0.5,col="darkgrey", beside=T, names.arg=c('I','II','III','IV'), xlab="Quadrant", ylab="Mean Substitutions Per Exon")
segments(barCenters, means_and_ses$means - means_and_ses$se, barCenters, means_and_ses$means + means_and_ses$se, lwd=1.5)
arrows(barCenters, means_and_ses$means - means_and_ses$se, barCenters, means_and_ses$means + means_and_ses$se, lwd=1.5, angle=90, code=3, length=0.05)
dev.off()


# T-tests
I_v_II = t.test(quadrant_stats$quad1_mean, quadrant_stats$quad2_mean, paired = TRUE)
I_v_III = t.test(quadrant_stats$quad1_mean, quadrant_stats$quad3_mean, paired = TRUE)
I_v_IV = t.test(quadrant_stats$quad1_mean, quadrant_stats$quad4_mean, paired = TRUE)
II_v_III = t.test(quadrant_stats$quad2_mean, quadrant_stats$quad3_mean, paired = TRUE)
II_v_IV = t.test(quadrant_stats$quad2_mean, quadrant_stats$quad4_mean, paired = TRUE)
III_v_IV = t.test(quadrant_stats$quad3_mean, quadrant_stats$quad4_mean, paired = TRUE)



# Get values
mean_I <- mean(quadrant_stats$quad1_mean)
mean_II <- mean(quadrant_stats$quad2_mean)
mean_III <- mean(quadrant_stats$quad3_mean)
mean_IV <- mean(quadrant_stats$quad4_mean)

stdev_I <- sd(quadrant_stats$quad1_mean)
stdev_II <- sd(quadrant_stats$quad2_mean)
stdev_III <- sd(quadrant_stats$quad3_mean)
stdev_IV <- sd(quadrant_stats$quad4_mean)

overall_quad_matrix[1,1] <- mean_I
overall_quad_matrix[1,2] <- mean_II
overall_quad_matrix[1,3] <- mean_III
overall_quad_matrix[1,4] <- mean_IV

overall_quad_matrix[2,1] <- stdev_I
overall_quad_matrix[2,2] <- stdev_II
overall_quad_matrix[2,3] <- stdev_III
overall_quad_matrix[2,4] <- stdev_IV

#
# x axis: I, II, III, IV (Divergence Heatmap Quadrant)
# y axis: Divergence
# bars with stdev plotted
means <- as.numeric(as.character(overall_quad_matrix[1,]))
names <- colnames(overall_quad_matrix)
errors <- as.numeric(as.character(overall_quad_matrix[2,]))

pdf(file="quadrant_comparison.pdf", width=10, height=10, pointsize=22)
par(par(mar = c(5, 6, 4, 5) + 0.1))
barCenters <- barplot(height=means, names.arg=names, beside=TRUE, ylim=c(0.8,1.4), xpd=FALSE, xlab="Quadrant", ylab="Mean Divergence")
segments(barCenters, means - errors, barCenters, means + errors, lwd=1.5)
arrows(barCenters, means - errors, barCenters, means + errors, lwd=1.5, angle=90, code=3, length=0.05)
dev.off()





################################## Variance tests ##############################

# Check for normal distribution with shapiro-wilk test:
shapiro.test(quadrant_stats$quad1_mean) #--> p>0.05
shapiro.test(quadrant_stats$quad2_mean) #--> p<0.05
shapiro.test(quadrant_stats$quad3_mean) #--> p>0.05
shapiro.test(quadrant_stats$quad4_mean) #--> p>0.05


# Because most of our data is not normally distributed, we should use either 
# Levene's test of the Fligner-Killeen test (http://www.sthda.com/english/wiki/f-test-compare-two-variances-in-r)

# Get individual variances
I_var = var(quadrant_stats$quad1_mean)
II_var = var(quadrant_stats$quad2_mean)
III_var = var(quadrant_stats$quad3_mean)
IV_var = var(quadrant_stats$quad4_mean)

# F-tests
I_v_II_var = var.test(quadrant_stats$quad1_mean, quadrant_stats$quad2_mean)
I_v_III_var = var.test(quadrant_stats$quad1_mean, quadrant_stats$quad3_mean)
I_v_IV_var = var.test(quadrant_stats$quad1_mean, quadrant_stats$quad4_mean)
II_v_III_var = var.test(quadrant_stats$quad2_mean, quadrant_stats$quad3_mean)
II_v_IV_var = var.test(quadrant_stats$quad2_mean, quadrant_stats$quad4_mean)
III_v_IV_var = var.test(quadrant_stats$quad3_mean, quadrant_stats$quad4_mean)



# x axis: I, II, III, IV (Divergence Heatmap Quadrant)
# y axis: Variance

means <- c(I_var,II_var,III_var,IV_var)
names <- colnames(overall_quad_matrix)

pdf(file="quadrant_comparison_variance.pdf", width=10, height=10, pointsize=22)
par(par(mar = c(5, 6, 4, 5) + 0.1))
barCenters <- barplot(height=means, names.arg=names, beside=TRUE, ylim=c(0.0,0.03), xpd=FALSE, xlab="Quadrant", ylab="Mean Variance")
dev.off()


############################ Repeat for 9 quad matrix #############
# Make separate matrix for mean and stdev values:
overall_quad_9_matrix <- data.frame(matrix(nrow=2, ncol=9))
colnames(overall_quad_9_matrix) <- c("i","ii","iii","iv","v","vi","vii","viii","ix")
rownames(overall_quad_9_matrix) <- c("mean","stdev")
############################## Used for paper ###############################################

# One-way ANOVA
subs <- c(quad_9_stats$quadi_mean, quad_9_stats$quadii_mean, quad_9_stats$quadiii_mean, quad_9_stats$quadiv_mean, quad_9_stats$quadv_mean, quad_9_stats$quadvi_mean, quad_9_stats$quadvii_mean, quad_9_stats$quadviii_mean, quad_9_stats$quadix_mean)
box <- c(rep("quadi_mean",41), rep("quadii_mean",41), rep("quadiii_mean",41), rep("quadiv_mean",41), rep("quadv_mean",41), rep("quadvi_mean",41), rep("quadvii_mean",41), rep("quadviii_mean",41), rep("quadix_mean",41))
box_means <- data.frame(subs,box)
plot(subs ~ box, data=box_means)
results <- aov(subs ~ box, data=box_means)
summary(results)
tukey.test <- TukeyHSD(results)

tukey.test <- TukeyHSD(results)
tukey.test

means <- tapply(subs,box,mean)
se <- tapply(subs,box,std.error)
means_and_ses <- data.frame(means,se)
means_and_ses


# Plot substitutions
pdf(file="quadrant_comparison_9.pdf", width=10, height=8, pointsize=22)
par(par(mar = c(5, 6, 4, 5) + 0.1))
barCenters <- barplot(height=means_and_ses$means,ylim=c(0,2), col="darkgrey", beside=T, names.arg=c('i','ii','iii','iv','v','vi','vii','viii','ix'), xlab="Box", ylab="Mean substitutions per exon")
segments(barCenters, means_and_ses$means - means_and_ses$se, barCenters, means_and_ses$means + means_and_ses$se, lwd=1.5)
arrows(barCenters, means_and_ses$means - means_and_ses$se, barCenters, means_and_ses$means + means_and_ses$se, lwd=1.5, angle=90, code=3, length=0.05)
dev.off()

#############################################################################

# Means and SDs
overall_quadi_mean_sd <- c(mean(quad_9_stats$quadi_mean),sd(quad_9_stats$quadi_mean))
overall_quadii_mean_sd <- c(mean(quad_9_stats$quadii_mean),sd(quad_9_stats$quadii_mean))
overall_quadiii_mean_sd <- c(mean(quad_9_stats$quadiii_mean),sd(quad_9_stats$quadiii_mean))
overall_quadiv_mean_sd <- c(mean(quad_9_stats$quadiv_mean),sd(quad_9_stats$quadiv_mean))
overall_quadv_mean_sd <- c(mean(quad_9_stats$quadv_mean),sd(quad_9_stats$quadv_mean))
overall_quadvi_mean_sd <- c(mean(quad_9_stats$quadvi_mean),sd(quad_9_stats$quadvi_mean))
overall_quadvii_mean_sd <- c(mean(quad_9_stats$quadvii_mean),sd(quad_9_stats$quadvii_mean))
overall_quadviii_mean_sd <- c(mean(quad_9_stats$quadviii_mean),sd(quad_9_stats$quadviii_mean))
overall_quadix_mean_sd <- c(mean(quad_9_stats$quadix_mean),sd(quad_9_stats$quadix_mean))




# T-tests
#i_v_ii = t.test(quad_9_stats$quadi_mean, quad_9_stats$quadii_mean, paired=TRUE)
#i_v_iii = t.test(quad_9_stats$quadi_mean, quad_9_stats$quadiii_mean, paired = TRUE)
#i_v_iv = t.test(quad_9_stats$quadi_mean, quad_9_stats$quadiv_mean, paired = TRUE)
#i_v_v = t.test(quad_9_stats$quadi_mean, quad_9_stats$quadv_mean, paired = TRUE)
#i_v_vi = t.test(quad_9_stats$quadi_mean, quad_9_stats$quadvi_mean, paired = TRUE)
#i_v_vii = t.test(quad_9_stats$quadi_mean, quad_9_stats$quadvii_mean, paired = TRUE)
#i_v_viii = t.test(quad_9_stats$quadi_mean, quad_9_stats$quadviii_mean, paired = TRUE)
#i_v_ix = t.test(quad_9_stats$quadi_mean, quad_9_stats$quadix_mean, paired = TRUE)
#ii_v_iii = t.test(quad_9_stats$quadii_mean, quad_9_stats$quadiii_mean, paired = TRUE)
#ii_v_iv = t.test(quad_9_stats$quadii_mean, quad_9_stats$quadiv_mean, paired = TRUE)
#ii_v_v = t.test(quad_9_stats$quadii_mean, quad_9_stats$quadv_mean, paired = TRUE)
#ii_v_vi = t.test(quad_9_stats$quadii_mean, quad_9_stats$quadvi_mean, paired = TRUE)
#ii_v_vii = t.test(quad_9_stats$quadii_mean, quad_9_stats$quadvii_mean, paired = TRUE)
#ii_v_viii = t.test(quad_9_stats$quadii_mean, quad_9_stats$quadviii_mean, paired = TRUE)
#ii_v_ix = t.test(quad_9_stats$quadii_mean, quad_9_stats$quadix_mean, paired = TRUE)
#iii_v_iv = t.test(quad_9_stats$quadiii_mean, quad_9_stats$quadiv_mean, paired = TRUE)
#iii_v_v = t.test(quad_9_stats$quadiii_mean, quad_9_stats$quadv_mean, paired = TRUE)
#iii_v_vi = t.test(quad_9_stats$quadiii_mean, quad_9_stats$quadvi_mean, paired = TRUE)
#iii_v_vii = t.test(quad_9_stats$quadiii_mean, quad_9_stats$quadvii_mean, paired = TRUE)
#iii_v_viii = t.test(quad_9_stats$quadiii_mean, quad_9_stats$quadviii_mean, paired = TRUE)
#iii_v_ix = t.test(quad_9_stats$quadiii_mean, quad_9_stats$quadix_mean, paired = TRUE)
#iv_v_v = t.test(quad_9_stats$quadiv_mean, quad_9_stats$quadv_mean, paired = TRUE)
#iv_v_vi = t.test(quad_9_stats$quadiv_mean, quad_9_stats$quadvi_mean, paired = TRUE)
#iv_v_vii = t.test(quad_9_stats$quadiv_mean, quad_9_stats$quadvii_mean, paired = TRUE)
#iv_v_viii = t.test(quad_9_stats$quadiv_mean, quad_9_stats$quadviii_mean, paired = TRUE)
#iv_v_ix = t.test(quad_9_stats$quadiv_mean, quad_9_stats$quadix_mean, paired = TRUE)
#v_v_vi = t.test(quad_9_stats$quadv_mean, quad_9_stats$quadvi_mean, paired = TRUE)
#v_v_vii = t.test(quad_9_stats$quadv_mean, quad_9_stats$quadvii_mean, paired = TRUE)
#v_v_viii = t.test(quad_9_stats$quadv_mean, quad_9_stats$quadviii_mean, paired = TRUE)
#v_v_ix = t.test(quad_9_stats$quadv_mean, quad_9_stats$quadix_mean, paired = TRUE)
#vi_v_vii = t.test(quad_9_stats$quadvi_mean, quad_9_stats$quadvii_mean, paired = TRUE)
#vi_v_viii = t.test(quad_9_stats$quadvi_mean, quad_9_stats$quadviii_mean, paired = TRUE)
#vi_v_ix = t.test(quad_9_stats$quadvi_mean, quad_9_stats$quadix_mean, paired = TRUE)
#vii_v_viii = t.test(quad_9_stats$quadvii_mean, quad_9_stats$quadviii_mean, paired = TRUE)
#vii_v_ix = t.test(quad_9_stats$quadvii_mean, quad_9_stats$quadix_mean, paired = TRUE)
#viii_v_ix = t.test(quad_9_stats$quadviii_mean, quad_9_stats$quadix_mean, paired = TRUE)


# Get values
mean_i <- mean(quad_9_stats$quadi_mean)
mean_ii <- mean(quad_9_stats$quadii_mean)
mean_iii <- mean(quad_9_stats$quadiii_mean)
mean_iv <- mean(quad_9_stats$quadiv_mean)
mean_v <- mean(quad_9_stats$quadv_mean)
mean_vi <- mean(quad_9_stats$quadvi_mean)
mean_vii <- mean(quad_9_stats$quadvii_mean)
mean_viii <- mean(quad_9_stats$quadviii_mean)
mean_ix <- mean(quad_9_stats$quadix_mean)

stdev_i <- sd(quad_9_stats$quadi_mean)
stdev_ii <- sd(quad_9_stats$quadii_mean)
stdev_iii <- sd(quad_9_stats$quadiii_mean)
stdev_iv <- sd(quad_9_stats$quadiv_mean)
stdev_v <- sd(quad_9_stats$quadv_mean)
stdev_vi <- sd(quad_9_stats$quadvi_mean)
stdev_vii <- sd(quad_9_stats$quadvii_mean)
stdev_viii <- sd(quad_9_stats$quadviii_mean)
stdev_ix <- sd(quad_9_stats$quadix_mean)

overall_quad_9_matrix[1,1] <- mean_i
overall_quad_9_matrix[1,2] <- mean_ii
overall_quad_9_matrix[1,3] <- mean_iii
overall_quad_9_matrix[1,4] <- mean_iv
overall_quad_9_matrix[1,5] <- mean_v
overall_quad_9_matrix[1,6] <- mean_vi
overall_quad_9_matrix[1,7] <- mean_vii
overall_quad_9_matrix[1,8] <- mean_viii
overall_quad_9_matrix[1,9] <- mean_ix

overall_quad_9_matrix[2,1] <- stdev_i
overall_quad_9_matrix[2,2] <- stdev_ii
overall_quad_9_matrix[2,3] <- stdev_iii
overall_quad_9_matrix[2,4] <- stdev_iv
overall_quad_9_matrix[2,5] <- stdev_v
overall_quad_9_matrix[2,6] <- stdev_vi
overall_quad_9_matrix[2,7] <- stdev_vii
overall_quad_9_matrix[2,8] <- stdev_viii
overall_quad_9_matrix[2,9] <- stdev_ix

#
# x axis: I, II, III, IV (Divergence Heatmap Quadrant)
# y axis: Divergence
# bars with stdev plotted
means <- as.numeric(as.character(overall_quad_9_matrix[1,]))
names <- colnames(overall_quad_9_matrix)
errors <- as.numeric(as.character(overall_quad_9_matrix[2,]))

#pdf(file="quadrant_comparison.pdf", width=10, height=10, pointsize=22)
par(par(mar = c(5, 6, 4, 5) + 0.1))
barCenters <- barplot(height=means, names.arg=names, beside=TRUE, ylim=c(0,2), xpd=FALSE, xlab="Quadrant", ylab="Mean Divergence")
segments(barCenters, means - errors, barCenters, means + errors, lwd=1.5)
arrows(barCenters, means - errors, barCenters, means + errors, lwd=1.5, angle=90, code=3, length=0.05)
#dev.off()





################################## Variance tests ##############################

# Check for normal distribution with shapiro-wilk test:
shapiro.test(quad_9_stats$quadi_mean) #--> p<0.05
shapiro.test(quad_9_stats$quadii_mean) #--> p>0.05
shapiro.test(quad_9_stats$quadiii_mean) #--> p>0.05
shapiro.test(quad_9_stats$quadiv_mean) #--> p>0.05
shapiro.test(quad_9_stats$quadv_mean) #--> p>0.05
shapiro.test(quad_9_stats$quadvi_mean) #--> p>0.05
shapiro.test(quad_9_stats$quadvii_mean) #--> p>0.05
shapiro.test(quad_9_stats$quadviii_mean) #--> p>0.05
shapiro.test(quad_9_stats$quadix_mean) #--> p>0.05


# Because most of our data is not normally distributed, we should use either 
# Levene's test of the Fligner-Killeen test (http://www.sthda.com/english/wiki/f-test-compare-two-variances-in-r)

# Get individual variances
i_var = var(quad_9_stats$quadi_mean)
ii_var = var(quad_9_stats$quadii_mean)
iii_var = var(quad_9_stats$quadiii_mean)
iv_var = var(quad_9_stats$quadiv_mean)
v_var = var(quad_9_stats$quadv_mean)
vi_var = var(quad_9_stats$quadvi_mean)
vii_var = var(quad_9_stats$quadvii_mean)
viii_var = var(quad_9_stats$quadviii_mean)
ix_var = var(quad_9_stats$quadix_mean)

# F-tests
i_v_ii = var.test(quad_9_stats$quadi_mean, quad_9_stats$quadii_mean, paired=TRUE)
i_v_iii = var.test(quad_9_stats$quadi_mean, quad_9_stats$quadiii_mean, paired = TRUE)
i_v_iv = var.test(quad_9_stats$quadi_mean, quad_9_stats$quadiv_mean, paired = TRUE)
i_v_v = var.test(quad_9_stats$quadi_mean, quad_9_stats$quadv_mean, paired = TRUE)
i_v_vi = var.test(quad_9_stats$quadi_mean, quad_9_stats$quadvi_mean, paired = TRUE)
i_v_vii = var.test(quad_9_stats$quadi_mean, quad_9_stats$quadvii_mean, paired = TRUE)
i_v_viii = var.test(quad_9_stats$quadi_mean, quad_9_stats$quadviii_mean, paired = TRUE)
i_v_ix = var.test(quad_9_stats$quadi_mean, quad_9_stats$quadix_mean, paired = TRUE)
ii_v_iii = var.test(quad_9_stats$quadii_mean, quad_9_stats$quadiii_mean, paired = TRUE)
ii_v_iv = var.test(quad_9_stats$quadii_mean, quad_9_stats$quadiv_mean, paired = TRUE)
ii_v_v = var.test(quad_9_stats$quadii_mean, quad_9_stats$quadv_mean, paired = TRUE)
ii_v_vi = var.test(quad_9_stats$quadii_mean, quad_9_stats$quadvi_mean, paired = TRUE)
ii_v_vii = var.test(quad_9_stats$quadii_mean, quad_9_stats$quadvii_mean, paired = TRUE)
ii_v_viii = var.test(quad_9_stats$quadii_mean, quad_9_stats$quadviii_mean, paired = TRUE)
ii_v_ix = var.test(quad_9_stats$quadii_mean, quad_9_stats$quadix_mean, paired = TRUE)
iii_v_iv = var.test(quad_9_stats$quadiii_mean, quad_9_stats$quadiv_mean, paired = TRUE)
iii_v_v = var.test(quad_9_stats$quadiii_mean, quad_9_stats$quadv_mean, paired = TRUE)
iii_v_vi = var.test(quad_9_stats$quadiii_mean, quad_9_stats$quadvi_mean, paired = TRUE)
iii_v_vii = var.test(quad_9_stats$quadiii_mean, quad_9_stats$quadvii_mean, paired = TRUE)
iii_v_viii = var.test(quad_9_stats$quadiii_mean, quad_9_stats$quadviii_mean, paired = TRUE)
iii_v_ix = var.test(quad_9_stats$quadiii_mean, quad_9_stats$quadix_mean, paired = TRUE)
iv_v_v = var.test(quad_9_stats$quadiv_mean, quad_9_stats$quadv_mean, paired = TRUE)
iv_v_vi = var.test(quad_9_stats$quadiv_mean, quad_9_stats$quadvi_mean, paired = TRUE)
iv_v_vii = var.test(quad_9_stats$quadiv_mean, quad_9_stats$quadvii_mean, paired = TRUE)
iv_v_viii = var.test(quad_9_stats$quadiv_mean, quad_9_stats$quadviii_mean, paired = TRUE)
iv_v_ix = var.test(quad_9_stats$quadiv_mean, quad_9_stats$quadix_mean, paired = TRUE)
v_v_vi = var.test(quad_9_stats$quadv_mean, quad_9_stats$quadvi_mean, paired = TRUE)
v_v_vii = var.test(quad_9_stats$quadv_mean, quad_9_stats$quadvii_mean, paired = TRUE)
v_v_viii = var.test(quad_9_stats$quadv_mean, quad_9_stats$quadviii_mean, paired = TRUE)
v_v_ix = var.test(quad_9_stats$quadv_mean, quad_9_stats$quadix_mean, paired = TRUE)
vi_v_vii = var.test(quad_9_stats$quadvi_mean, quad_9_stats$quadvii_mean, paired = TRUE)
vi_v_viii = var.test(quad_9_stats$quadvi_mean, quad_9_stats$quadviii_mean, paired = TRUE)
vi_v_ix = var.test(quad_9_stats$quadvi_mean, quad_9_stats$quadix_mean, paired = TRUE)
vii_v_viii = var.test(quad_9_stats$quadvii_mean, quad_9_stats$quadviii_mean, paired = TRUE)
vii_v_ix = var.test(quad_9_stats$quadvii_mean, quad_9_stats$quadix_mean, paired = TRUE)
viii_v_ix = var.test(quad_9_stats$quadviii_mean, quad_9_stats$quadix_mean, paired = TRUE)




# x axis: I, II, III, IV (Divergence Heatmap Quadrant)
# y axis: Variance

means <- c(i_var,ii_var,iii_var,iv_var,v_var,vi_var,vii_var,viii_var,ix_var)
names <- colnames(overall_quad_9_matrix)

#pdf(file="quadrant_comparison_variance.pdf", width=10, height=10, pointsize=22)
par(par(mar = c(5, 6, 4, 5) + 0.1))
barCenters <- barplot(height=means, names.arg=names, beside=TRUE, ylim=c(0.0,0.06), xpd=FALSE, xlab="Quadrant", ylab="Mean Variance")
#dev.off()


