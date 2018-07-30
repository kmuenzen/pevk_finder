# Bar plot for overall ortholog divergence

library(ggplot2)

# Import csv
divergence_data <- as.data.frame(read.csv("mean_alignment_distances.csv"))

values <- as.numeric(as.character(divergence_data$mean_dist_nt))
names <- c(1:nrow(divergence_data))
variances <- as.numeric(as.character(divergence_data$mean_dist_nt_variance))

# Get regional stats
conserved_values = values[1:22]
hypervariable_values = values[23:35]

# Means and SEs
conserved_mean = mean(conserved_values)
conserved_se = std.error(conserved_values)

hv_mean = mean(hypervariable_values)
hv_se = std.error(hypervariable_values)


# T-test
diffs = t.test(hypervariable_values, conserved_values, paired=FALSE)


colors <- c("red", "black","black","black","black","red","black","red","red","black","black","black","black","red","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black","black")

pdf(file="overall_ortholog_divergence.pdf", width=20, height=10, pointsize=22)

par(mar = c(5, 5, 2, 0.05))
plot(x=names, y=values, pch=19, bty="n", col= colors, ylim=c(0,0.25), xlab="Ortholog", ylab="Divergence")
#segments(barCenters, values - variances, barCenters, values + variances, lwd=1.5)
arrows(names, values - variances, names, values + variances, lwd=1.5, angle=90, code=3, length=0.05)
#abline(h=overall_mean, lty=2)
dev.off()

# 1, 6, 8, 9, 14


