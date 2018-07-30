# This script performs a statistical comparison of PEVK Finder and cDNA/gnomon titin annotations
setwd("~/Desktop/titin_project/ttn_annotations")

# cdna only anotations
cdna_only_annotations = read.csv("cdna_gnomon_annotations_updated.csv")
cdna_only_annotations = data.frame(cdna_only_annotations)
rownames(cdna_only_annotations) <- 1:nrow(cdna_only_annotations)

# gnomon only anotations
gnomon_only_annotations = read.csv("gnomon_annotations_updated.csv")
gnomon_only_annotations = data.frame(gnomon_only_annotations)
rownames(gnomon_only_annotations) <- 1:nrow(gnomon_only_annotations)



############## Updated code 07/27/18 ##############
cdna_cols = c(1:7,10:11)
cDNA_annotations_subset = cdna_only_annotations[,cdna_cols]
gnomon_annotations_subset = gnomon_only_annotations[,1:9]

all_annotations <- rbind(cDNA_annotations_subset, gnomon_annotations_subset)

finder_mean_overall <- mean(all_annotations$all_pevk_finder_exons)
cdna_gnomon_mean_overall <- mean(all_annotations$gnomon_exons)

finder_mean_se <- std.error(all_annotations$all_pevk_finder_exons)
cdna_gnomon_mean_se <- std.error(all_annotations$gnomon_exons)

annotation_mean_ttest <- t.test(all_annotations$all_pevk_finder_exons,all_annotations$gnomon_exons)


##### Proportion of union
finder_union_mean_cdna <- mean(cDNA_annotations_subset$X._pevk_of_union)
finder_union_mean_gnomon <- mean(gnomon_annotations_subset$X._pevk_of_unions)

finder_union_mean_cdna_se <- std.error(cDNA_annotations_subset$X._pevk_of_union)
finder_union_mean_gnomon_se <- std.error(gnomon_annotations_subset$X._pevk_of_unions)

finder_union_ttest <- t.test(gnomon_annotations_subset$X._pevk_of_unions,cDNA_annotations_subset$X._pevk_of_union)

#####
cdna_union_mean_cdna <- mean(cDNA_annotations_subset$X._gnomon_of_union)
gnomon_union_mean_gnomon <- mean(gnomon_annotations_subset$X._gnomon_of_union)

cdna_union_mean_cdna_se <- std.error(cDNA_annotations_subset$X._gnomon_of_union)
gnomon_union_mean_gnomon_se <- std.error(gnomon_annotations_subset$X._gnomon_of_union)

ttest_finder_vs_cdna <- t.test(cDNA_annotations_subset$X._pevk_of_union,cDNA_annotations_subset$X._gnomon_of_union)
ttest_finder_vs_gnomon <- t.test(gnomon_annotations_subset$X._pevk_of_unions,gnomon_annotations_subset$X._gnomon_of_union)

###################################################



##### stat tests
# 1. PEVK-N exons
mean_pevkn_pevk_finder = mean(all_annotations$pevk_n_exons)
mean_pevkn_pevk_finder_se = std.error(all_annotations$pevk_n_exons)

# 1. PEVK-C exons
mean_pevkc_pevk_finder = mean(all_annotations$pevk_c_exons)
mean_pevkc_pevk_finder_se = std.error(all_annotations$pevk_c_exons)

# Var test
var_test_regions_exons = var.test(all_annotations$pevk_n_exon, all_annotations$pevk_c_exons)
