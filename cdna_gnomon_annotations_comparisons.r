# This script performs a statistical comparison of PEVK Finder and cDNA/gnomon titin annotations


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





# PEVK Finder means and sd
finder_mean <- mean(gnomon_and_cdna_annotations$pevk_finder)
finder_se <- sd(gnomon_and_cdna_annotations$pevk_finder)/sqrt(length(gnomon_and_cdna_annotations$pevk_finder))

# Gnomon + cdna mean
gnomon_cdna_mean <- mean(gnomon_and_cdna_annotations$annotated)
gnomon_cdna_se <- sd(gnomon_and_cdna_annotations$annotated)/sqrt(length(gnomon_and_cdna_annotations$annotated))


# overall % identification mean
id_mean <- mean(gnomon_and_cdna_annotations$ident)
id_se <- sd(gnomon_and_cdna_annotations$ident)/sqrt(length(gnomon_and_cdna_annotations$ident))


# cdna only id mean
cdna_id_mean <- mean(cdna_only_annotations$ident)
cdna_id_se <- sd(cdna_only_annotations$ident)/sqrt(length(cdna_only_annotations$ident))


# gnomon only id mean
gnomon_id_mean <- mean(gnomon_only_annotations$ident)
gnomon_id_se <- sd(gnomon_only_annotations$ident)/sqrt(length(gnomon_only_annotations$ident))


# missing exons
missing_mean <- mean(gnomon_and_cdna_annotations$missing)
missing_se <- sd(gnomon_and_cdna_annotations$missing)/sqrt(length(gnomon_and_cdna_annotations$missing))


# novel exons
novel_mean <- mean(gnomon_and_cdna_annotations$novel)
novel_se <- sd(gnomon_and_cdna_annotations$novel)/sqrt(length(gnomon_and_cdna_annotations$novel))


# cdna novel exons
cdna_novel_mean <- mean(cdna_only_annotations$novel)
cdna_novel_se <- sd(cdna_only_annotations$novel)/sqrt(length(cdna_only_annotations$novel))

# gnomo nnovel exons
gnomon_novel_mean <- mean(gnomon_only_annotations$novel)
gnomon_novel_se <- sd(gnomon_only_annotations$novel)/sqrt(length(gnomon_only_annotations$novel))



##### stat tests
# 1. total exons identified
comp_1 = t.test(gnomon_and_cdna_annotations$pevk_finder, gnomon_and_cdna_annotations$annotated, paired=TRUE)

# 2. novel exons
comp_2 = t.test(gnomon_and_cdna_annotations$pevk_finder, gnomon_and_cdna_annotations$annotated, paired=TRUE)


