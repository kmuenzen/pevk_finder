# This script calculates the min, mean and max exon length and PEVK ratio values for human and mouse PEVK exons


library(gplots)
library(RColorBrewer)
library(MASS)
library(plot3D)


# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios_known <- read.csv("human_known_stats.csv", header=F)
pevk_lengths_and_ratios_known <- t(pevk_lengths_and_ratios_known)
pevk_lengths_and_ratios_known <- data.frame(pevk_lengths_and_ratios_known)
rownames(pevk_lengths_and_ratios_known) <- 1:nrow(pevk_lengths_and_ratios_known)
colnames(pevk_lengths_and_ratios_known) <- c("name","length","percent_pevk")


min_length = min(pevk_lengths_and_ratios_known$length)
min_ratio = min(pevk_lengths_and_ratios_known$percent_pevk)

max_length = max(pevk_lengths_and_ratios_known$length)
max_ratio = max(pevk_lengths_and_ratios_known$percent_pevk)

mean_length = mean(pevk_lengths_and_ratios_known$length)
mean_ratio = mean(pevk_lengths_and_ratios_known$percent_pevk)



print(paste("Minimum Exon Length:", min_length))
print(paste("Maximum Exon Length:", max_length))
print(paste("Mean Exon Length:", mean_length))

print(paste("Minimum PEVK Ratio:", min_ratio))
print(paste("Maximum PEVK Ratio:", max_ratio))
print(paste("Mean PEVK Ratio:", mean_ratio))



# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios_known <- read.csv("mouse_known_stats.csv", header=F)
pevk_lengths_and_ratios_known <- t(pevk_lengths_and_ratios_known)
pevk_lengths_and_ratios_known <- data.frame(pevk_lengths_and_ratios_known)
rownames(pevk_lengths_and_ratios_known) <- 1:nrow(pevk_lengths_and_ratios_known)
colnames(pevk_lengths_and_ratios_known) <- c("name","length","percent_pevk")


min_length = min(as.numeric(as.vector(pevk_lengths_and_ratios_known$length)))
min_ratio = min(as.numeric(as.vector(pevk_lengths_and_ratios_known$percent_pevk)))

max_length = max(as.numeric(as.vector(pevk_lengths_and_ratios_known$length)))
max_ratio = max(as.numeric(as.vector(pevk_lengths_and_ratios_known$percent_pevk)))

mean_length = mean(as.numeric(as.vector(pevk_lengths_and_ratios_known$length)))
mean_ratio = mean(as.numeric(as.vector(pevk_lengths_and_ratios_known$percent_pevk)))



print(paste("Minimum Exon Length:", min_length))
print(paste("Maximum Exon Length:", max_length))
print(paste("Mean Exon Length:", mean_length))

print(paste("Minimum PEVK Ratio:", min_ratio))
print(paste("Maximum PEVK Ratio:", max_ratio))
print(paste("Mean PEVK Ratio:", mean_ratio))