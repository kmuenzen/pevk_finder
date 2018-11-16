# This code creates plots for visualizing evolutaionry analyses of the orthologous exons, based on the human exon-intron organiation determined by tool

# Set workng directory:
setwd("~/Desktop/titin_project/ortholog_data")

# Read CSV file into R and convert to CSV
exon_data <- read.csv("exon_evolutionary_analysis.csv")
exon_data <- as.data.frame(exon_data)

# Extract data on evolutionary distance, dn.ds values, and hit count
nt_dist <- exon_data$mean_distance_nucleotide
#dn_ds <- exon_data$dn_ds_stat # star values that are not statistically significant...can do manually
hits <- exon_data$hit_count

# Assign colors to different parameters
pal1 <- brewer.pal(11, "PiYG")[c(1:4,8:11)]
pal2 <-rev( brewer.pal(11, "PuOr")[c(1:4,8:11)])
pal3 <- brewer.pal(11, "RdGy")[c(1:4,8:11)]
pal4 <- rev(brewer.pal(11, "RdYlBu")[c(1:4,8:11)])


###### Create color scale ######
library(neurobase)

pdf(file="~/Desktop/titin_project/ortholog_data/ortholog_divergence_colorbar.pdf", width=25, height=13, pointsize=22)
colorbar(c(0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1.0), pal4, maxleft=0.5)
#axis(at = c(1,2,3,4,5,6,7,8),side=4, labels=c(0.125,0.25,0.375,0.5,0.625,0.75,0.875,1.0))
dev.off()

10-28
################################


nt_dist_cols <- pal1[as.numeric(cut(nt_dist, breaks=8))]
#dn_ds_cols <- pal1[as.numeric(cut(dn_ds, breaks=8))]
hit_cols <- pal4[as.numeric(cut(hits, breaks=8))]

# Extract starts and ends
exon_starts <- as.numeric(as.character(exon_data$start))
exon_ends <- as.numeric(as.character(exon_data$end))

# Create rectangle plot
pdf(file="human_hit_plot.pdf", width=30, height=10, pointsize=22)


op <- par(bg = "white")
plot(c(100000,167500), c(0,2), type="n", xlab="Titin Coordinates", ylab="", yaxt="n", bty='n')

# Plot hit count on top:
rect(exon_starts, 0.5, exon_ends, 1.5, col=hit_cols, border=NA)
abline(h=1.0)

# Plot divergence underneath
#rect(exon_starts, 2.5, exon_ends, 3.5, col=nt_dist_cols, border=NA)
#abline(h=3.0)

# Plot dn/ds on bottom
#rect(exon_starts, 0.5, exon_ends, 1.5, col=dn_ds_cols, border=NA)
#abline(h=1.0)

#axis(side=2, at=c(1.0, 3.0, 5.0), labels=c("dN/dS", "Divergence","Hits"))

axis(side=3, at=c(1.0, 3.0, 5.0), labels=c("dN/dS", "Divergence","Hits"))


dev.off()


