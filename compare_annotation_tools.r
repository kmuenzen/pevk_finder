##### Set working directory

library(ggplot2)

# Read in csv
tool_stats <- data.frame(read.csv("tool_comparison_stats.csv"))
tool_stats_compressed <- cbind(tool_stats[,1:2], tool_stats[,5])
colnames(tool_stats_compressed) <- c("Species","tool","ident")

# Mouse and human plot
c <- ggplot(tool_stats_compressed, aes(x=tool, y=ident, fill=Species))
bars <- geom_bar(stat = "identity", position="dodge", width=0.5)
x_labels <- xlab("Annotation Tool")
y_labels <- ylab("% Exons Identified")
colors <- scale_fill_manual("Species", values = c("human" = "black", "mouse" = "darkgrey"))

pdf(file="tool_comparison.pdf", width=7, height=5, pointsize=22)
c + bars + x_labels + y_labels + colors + theme(axis.line=element_blank(),
                                                panel.background=element_blank(),
                                                axis.ticks=element_blank(),
                                                panel.border=element_blank(),
                                                panel.grid.major=element_blank(),
                                                panel.grid.minor=element_blank(),
                                                plot.background=element_blank())
dev.off()

# Human only plot

human_only_labels = c(levels(droplevels(tool_stats$program[9])), levels(droplevels(tool_stats$program[1])), levels(droplevels(tool_stats$program[7])), levels(droplevels(tool_stats$program[3])), levels(droplevels(tool_stats$program[5])))
human_only_values = c(tool_stats$X._identified[9], tool_stats$X._identified[1], tool_stats$X._identified[8], tool_stats$X._identified[3], tool_stats$X._identified[6])

pdf(file="tool_comparison_human.pdf", width=10, height=10, pointsize=19)
barplot(height=human_only_values, names.arg=human_only_labels, xlab="Annotation Tool", ylab="% Exons Identified", ylim=c(0,1), width=0.3,space=0.5, border=NA, beside=TRUE)
dev.off()

# Mouse only plot

mouse_only_labels = c(levels(droplevels(tool_stats$program[4])), levels(droplevels(tool_stats$program[10])), levels(droplevels(tool_stats$program[2])), levels(droplevels(tool_stats$program[7])), levels(droplevels(tool_stats$program[5])))
mouse_only_labels = c(mouse_only_labels[2],mouse_only_labels[3],mouse_only_labels[4],mouse_only_labels[1],mouse_only_labels[5])
mouse_only_values = c(tool_stats$X._identified[4], tool_stats$X._identified[10], tool_stats$X._identified[2], tool_stats$X._identified[7], tool_stats$X._identified[5])
mouse_only_values = c(mouse_only_values[2],mouse_only_values[3],mouse_only_values[4],mouse_only_values[1],mouse_only_values[5])

pdf(file="tool_comparison_mouse.pdf", width=10, height=10, pointsize=19)
barplot(height=mouse_only_values, names.arg=mouse_only_labels, xlab="Annotation Tool", ylab="% Exons Identified", ylim=c(0,1), width=0.3,space=0.5, border=NA, beside=TRUE)
dev.off()


# pevk finder, genscan, augustus, fgenesh, geneid