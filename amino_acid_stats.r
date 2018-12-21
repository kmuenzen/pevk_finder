###### Performs statistical analyses and creates chart for overall amino acid counts #######

##### Set working directory
setwd("~/Desktop/titin_project")

# Read in csv file
aa_stats <- as.data.frame(read.csv("residue_overall_stats.csv"))

# Define column names
column_names <- c("species",	"exon_lengths1",	"A_count",	"I_count",	"L_count",	"M_count",	"V_count",	"F_count",	"W_count",	"Y_count",	"N_count",	"C_count",	"Q_count",	"S_count",	"T_count",	"D_count",	"E_count",	"R_count",	"H_count",	"K_count",	"G_count",	"P_count")

# Create vector of conserved means
as.character(aa_stats[43,])

# Create vector of variable means

# Create matrix for conserved value
conserved_stats <- data.frame(aa_stats[3:44,])
colnames(conserved_stats) <- column_names

# Create matrix for hypervariable stats
variable_stats <- data.frame(aa_stats[48:89,])
colnames(variable_stats) <- column_names

# Make data frame to hold all test statistics and p-values
t_tests <- data.frame(matrix(ncol=21, nrow=11))
colnames(t_tests) <- c("exon_length", "A",	"I",	"L",	"M",	"V",	"F",	"W",	"Y",	"N",	"C",	"Q",	"S",	"T",	"D",	"E",	"R",	"H",	"K",	"G",	"P")
rownames(t_tests) <- c("test_stat","p-value","Conserved","Conserved_SE", "Hypervariable","Hypervariable_SE", "Overall","Conserved_Var","Hypervariable_Var","f_test_stat","f_test_p_value")

# Do t-tests for each category
for (i in 2:length(column_names)) {
  conserved_col <- as.numeric(as.character(conserved_stats[1:41,][,i]))
  variable_col <- as.numeric(as.character(variable_stats[1:41,][,i]))
  test <- t.test(conserved_col,variable_col,paired=TRUE)
  conserved_mean <- mean(as.numeric(as.character(conserved_stats[1:41,][,i])))
  conserved_se <- std.error(as.numeric(as.character(conserved_stats[1:41,][,i])))
  variable_mean <- mean(as.numeric(as.character(variable_stats[1:41,][,i])))
  variable_se <- std.error(as.numeric(as.character(variable_stats[1:41,][,i])))
  overall_mean <- mean(conserved_mean, variable_mean)
  conserved_var <- var(as.numeric(as.character(conserved_stats[1:41,][,i])))
  t_tests[1,i-1] <- test$statistic
  t_tests[2,i-1] <- test$p.value
  t_tests[3,i-1] <- conserved_mean
  t_tests[4,i-1] <- conserved_se
  t_tests[5,i-1] <- variable_mean
  t_tests[6,i-1] <- variable_se
  t_tests[7,i-1] <- overall_mean
}


# Calculate variance
# E:





# Order by overall mean
mean_sorted <- t_tests[order(-t_tests[3,])]

get_rows = c(3,5)
new_rownames = c('PEVK-N','PEVK-C')
conserved_variable_means <- mean_sorted[get_rows,][,2:5]
all_conserved_variable_means <- mean_sorted[get_rows,][,2:21]
rownames(all_conserved_variable_means) = new_rownames
rownames(conserved_variable_means) = new_rownames

#### Find errors ####
get_error_rows = c(4,6)
errors <- mean_sorted[get_error_rows,][,2:21]
rownames(errors) <- new_rownames

pdf(file="amino_acid_stats.pdf", width=30, height=10, pointsize=22)
par(par(mar = c(5, 6, 4, 5) + 0.1))
barCenters <- barplot(as.matrix(all_conserved_variable_means),ylim=c(0,7), col=c("black","darkgrey"), beside=T, names.arg=colnames(all_conserved_variable_means), xlab="Amino Acid", ylab="Mean Count Per Exon", legend = rownames(all_conserved_variable_means), args.legend = list(x ='topright', bty='n', inset=c(0.035,0.001)))
segments(barCenters, data.matrix(all_conserved_variable_means) - data.matrix(errors), barCenters, data.matrix(all_conserved_variable_means) + data.matrix(errors), lwd=1.5)
arrows(barCenters, data.matrix(all_conserved_variable_means) - data.matrix(errors), barCenters, data.matrix(all_conserved_variable_means) + data.matrix(errors), lwd=1.5, angle=90, code=3, length=0.05)

#####
dev.off()

# Order into P, E, V, K
ordered_cols <- c(1,3,2,4)
conserved_variable_means_ordered = conserved_variable_means[1:2,][,ordered_cols]
conserved_variable_errors <- mean_sorted[get_error_rows,][,2:5]
conserved_variable_errors <- conserved_variable_errors[,ordered_cols]

pdf(file="amino_acid_stats_updated.pdf", width=10, height=8, pointsize=22)

par(mar=c(5,5,4,8))
barCenters = barplot(as.matrix(conserved_variable_means_ordered),ylim=c(0,7.5), col=c("darkgrey","gray25"), beside=T, names.arg=colnames(conserved_variable_means_ordered), xlab="Amino Acid", ylab="Mean number of residues per exon")
segments(barCenters, data.matrix(conserved_variable_means_ordered) - data.matrix(conserved_variable_errors), barCenters, data.matrix(conserved_variable_means_ordered) + data.matrix(conserved_variable_errors), lwd=1.5)
arrows(barCenters, data.matrix(conserved_variable_means_ordered) - data.matrix(conserved_variable_errors), barCenters, data.matrix(conserved_variable_means_ordered) + data.matrix(conserved_variable_errors), lwd=1.5, angle=90, code=3, length=0.05)

dev.off()


# Read in Human residue stats
human_stats <- as.data.frame(read.csv("Homo_sapiens_residue_stats.csv"))

exon_lengths = human_stats$exon_lengths
e_count = human_stats$E_count

percent_e = c()

for (i in 1:length(e_count)){
  current_exon_length = exon_lengths[i]
  current_e_count = e_count[i]
  current_percent_e = current_e_count/current_exon_length
  percent_e = c(percent_e, current_percent_e)
}

percent_e

###### Threshold tests #####
sampled_means = c()
counter = 1
while (counter <= 100){
  sampled_mean = mean(sample(percent_e,111,replace=TRUE))
  sampled_means = c(sampled_means,sampled_mean)
  counter  = counter + 1
}
mean(sampled_means) #--> threshold is 0.1722969


for (i in (1:length(percent_e))){
  if (percent_e[i] > 0.1723562){
    print(i)
  }
}


#pdf(file="percent_e_stats.pdf", width=30, height=10, pointsize=22)

par(mar=c(5,5,4,8))
barplot(percent_e,ylim=c(0, 0.7),xlab="Human PEVK Exons", ylab="Percent Glutamate (E) Per Exon")

#dev.off()
