

##### Install new packages
library(gplots)
library(RColorBrewer)
library(MASS)
library(plot3D)
#library(plotly)


##################################### DATA HANDLING AND HEATMAPS ####################################

##### Read the csv files into R
human_results <- read.csv("human_final_results")
#human_results_ordered <- human_results[order(human_results$min_length),]



##### Separate results into min exon length
list <- split(human_results,human_results$min_length)

min_10_results <- data.frame(list[1])
rownames(min_10_results) <- 1:nrow(min_10_results)
colnames(min_10_results) <- colnames(human_results)
min_11_results <- data.frame(list[2])
rownames(min_11_results) <- 1:nrow(min_11_results)
colnames(min_11_results) <- colnames(human_results)
min_12_results <- data.frame(list[3])
rownames(min_12_results) <- 1:nrow(min_12_results)
colnames(min_12_results) <- colnames(human_results)
min_13_results <- data.frame(list[4])
rownames(min_13_results) <- 1:nrow(min_13_results)
colnames(min_13_results) <- colnames(human_results)
min_14_results <- data.frame(list[5])
rownames(min_14_results) <- 1:nrow(min_14_results)
colnames(min_14_results) <- colnames(human_results)
min_15_results <- data.frame(list[6])
rownames(min_15_results) <- 1:nrow(min_15_results)
colnames(min_15_results) <- colnames(human_results)
min_16_results <- data.frame(list[7])
rownames(min_16_results) <- 1:nrow(min_16_results)
colnames(min_16_results) <- colnames(human_results)
min_17_results <- data.frame(list[8])
rownames(min_17_results) <- 1:nrow(min_17_results)
colnames(min_17_results) <- colnames(human_results)
min_18_results <- data.frame(list[9])
rownames(min_18_results) <- 1:nrow(min_18_results)
colnames(min_18_results) <- colnames(human_results)
min_19_results <- data.frame(list[10])
rownames(min_19_results) <- 1:nrow(min_19_results)
colnames(min_19_results) <- colnames(human_results)
min_20_results <- data.frame(list[11])
rownames(min_20_results) <- 1:nrow(min_20_results)
colnames(min_20_results) <- colnames(human_results)
min_21_results <- data.frame(list[12])
rownames(min_21_results) <- 1:nrow(min_21_results)
colnames(min_21_results) <- colnames(human_results)
min_22_results <- data.frame(list[13])
rownames(min_22_results) <- 1:nrow(min_22_results)
colnames(min_22_results) <- colnames(human_results)
min_23_results <- data.frame(list[14])
rownames(min_23_results) <- 1:nrow(min_23_results)
colnames(min_23_results) <- colnames(human_results)
min_24_results <- data.frame(list[15])
rownames(min_24_results) <- 1:nrow(min_24_results)
colnames(min_24_results) <- colnames(human_results)
min_25_results <- data.frame(list[16])
rownames(min_25_results) <- 1:nrow(min_25_results)
colnames(min_25_results) <- colnames(human_results)
min_26_results <- data.frame(list[17])
rownames(min_26_results) <- 1:nrow(min_26_results)
colnames(min_26_results) <- colnames(human_results)
min_27_results <- data.frame(list[18])
rownames(min_27_results) <- 1:nrow(min_27_results)
colnames(min_27_results) <- colnames(human_results)
min_28_results <- data.frame(list[19])
rownames(min_28_results) <- 1:nrow(min_28_results)
colnames(min_28_results) <- colnames(human_results)
min_29_results <- data.frame(list[20])
rownames(min_29_results) <- 1:nrow(min_29_results)
colnames(min_29_results) <- colnames(human_results)
min_30_results <- data.frame(list[21])
rownames(min_30_results) <- 1:nrow(min_30_results)
colnames(min_30_results) <- colnames(human_results)


###### get_score: function that scores each match based on exon matches, extraneous exons and perfect matches

get_score <- function(df) {
  
  for (i in 1:nrow(df)) {
    
    # Score total exon matches
    total_exon_matches <- as.numeric(df[i,5])
    total_exon_matches_score <- (total_exon_matches/113)*0.7
    
    # Score extraneous exons:
    extraneous_exons <- as.numeric(df[i,6])
    # If the number of extraneous exons is greater than or equal to 100:
    if (extraneous_exons >= 100) {
      extraneous_exons_score <- 0
    }
    # If the number of extraneous exons is less than 100:
    if (extraneous_exons < 100) {
      extraneous_exons_score <- (100 - extraneous_exons)*0.002
    }
    
    # Score perfect_matches
    perfect_matches <- as.numeric(df[i,7])
    perfect_matches_score <- (perfect_matches/total_exon_matches)*0.1
    
    score <- total_exon_matches_score + extraneous_exons_score + perfect_matches_score
    df$score[i] <- score
  }
  df
}



##### Apply get score function to separated results 

min_10_results <- get_score(min_10_results)
min_11_results <- get_score(min_11_results)
min_12_results <- get_score(min_12_results)
min_13_results <- get_score(min_13_results)
min_14_results <- get_score(min_14_results)
min_15_results <- get_score(min_15_results)
min_16_results <- get_score(min_16_results)
min_17_results <- get_score(min_17_results)
min_18_results <- get_score(min_18_results)
min_19_results <- get_score(min_19_results)
min_20_results <- get_score(min_20_results)
min_21_results <- get_score(min_21_results)
min_22_results <- get_score(min_22_results)
min_23_results <- get_score(min_23_results)
min_24_results <- get_score(min_24_results)
min_25_results <- get_score(min_25_results)
min_26_results <- get_score(min_26_results)
min_27_results <- get_score(min_27_results)
min_28_results <- get_score(min_28_results)
min_29_results <- get_score(min_29_results)
min_30_results <- get_score(min_30_results)



##### Create list of scored results

frame_list <- list(min_10_results,min_11_results,min_12_results,min_13_results,min_14_results,min_15_results,min_16_results,min_17_results,min_18_results,min_19_results,min_20_results,min_21_results,min_22_results,min_23_results,min_24_results,min_25_results,min_26_results,min_27_results,min_28_results,min_29_results,min_30_results)



##### Order scores from greatest to least

frame_list_scored <- lapply(frame_list, function(df) {
  df <- df[order(-df$score),]
  rownames(df) <- 1:nrow(df)
  df
})



##### Find mins and maxes of all score frames

min_10_results_scored_scores <- frame_list_scored[[1]]$score
min_10 <- min_10_results_scored_scores[819]
max_10 <- min_10_results_scored_scores[1]

min_11_results_scored_scores <- frame_list_scored[[2]]$score
min_11 <- min_11_results_scored_scores[819]
max_11 <- min_11_results_scored_scores[1]

min_12_results_scored_scores <- frame_list_scored[[3]]$score
min_12 <- min_12_results_scored_scores[819]
max_12 <- min_12_results_scored_scores[1]

min_13_results_scored_scores <- frame_list_scored[[4]]$score
min_13 <- min_13_results_scored_scores[819]
max_13 <- min_13_results_scored_scores[1]

min_14_results_scored_scores <- frame_list_scored[[5]]$score
min_14 <- min_14_results_scored_scores[819]
max_14 <- min_14_results_scored_scores[1]

min_15_results_scored_scores <- frame_list_scored[[6]]$score
min_15 <- min_15_results_scored_scores[819]
max_15 <- min_15_results_scored_scores[1]

min_16_results_scored_scores <- frame_list_scored[[7]]$score
min_16 <- min_16_results_scored_scores[819]
max_16 <- min_16_results_scored_scores[1]

min_17_results_scored_scores <- frame_list_scored[[8]]$score
min_17 <- min_17_results_scored_scores[819]
max_17 <- min_17_results_scored_scores[1]

min_18_results_scored_scores <- frame_list_scored[[9]]$score
min_18 <- min_18_results_scored_scores[819]
max_18 <- min_18_results_scored_scores[1]

min_19_results_scored_scores <- frame_list_scored[[10]]$score
min_19 <- min_19_results_scored_scores[819]
max_19 <- min_19_results_scored_scores[1]

min_20_results_scored_scores <- frame_list_scored[[11]]$score
min_20 <- min_20_results_scored_scores[819]
max_20 <- min_20_results_scored_scores[1]

min_21_results_scored_scores <- frame_list_scored[[12]]$score
min_21 <- min_21_results_scored_scores[819]
max_21 <- min_21_results_scored_scores[1]

min_22_results_scored_scores <- frame_list_scored[[13]]$score
min_22 <- min_22_results_scored_scores[819]
max_22 <- min_22_results_scored_scores[1]

min_23_results_scored_scores <- frame_list_scored[[14]]$score
min_23 <- min_23_results_scored_scores[819]
max_23 <- min_23_results_scored_scores[1]

min_24_results_scored_scores <- frame_list_scored[[15]]$score
min_24 <- min_24_results_scored_scores[819]
max_24 <- min_24_results_scored_scores[1]

min_25_results_scored_scores <- frame_list_scored[[16]]$score
min_25 <- min_25_results_scored_scores[819]
max_25 <- min_25_results_scored_scores[1]

min_26_results_scored_scores <- frame_list_scored[[17]]$score
min_26 <- min_26_results_scored_scores[819]
max_26 <- min_26_results_scored_scores[1]

min_27_results_scored_scores <- frame_list_scored[[18]]$score
min_27 <- min_27_results_scored_scores[819]
max_27 <- min_27_results_scored_scores[1]

min_28_results_scored_scores <- frame_list_scored[[19]]$score
min_28 <- min_28_results_scored_scores[819]
max_28 <- min_28_results_scored_scores[1]

min_29_results_scored_scores <- frame_list_scored[[20]]$score
min_29 <- min_29_results_scored_scores[819]
max_29 <- min_29_results_scored_scores[1]

min_30_results_scored_scores <- frame_list_scored[[21]]$score
min_30 <- min_30_results_scored_scores[819]
max_30 <- min_30_results_scored_scores[1]



##### List of heatmap matrices

ratios <- c(0.45,0.46,0.47,0.48,0.49,0.50,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.60,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.70,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.80,0.81,0.82,0.83)



##### Create empty data frames for heatmap matrices

heatmap_matrix_10 <- data.frame(matrix(ncol=21,nrow=39))
colnames(heatmap_matrix_10) <- 10:30
rownames(heatmap_matrix_10) <- ratios

heatmap_matrix_11 <- data.frame(matrix(ncol=21,nrow=39))
colnames(heatmap_matrix_11) <- 10:30
rownames(heatmap_matrix_11) <- ratios

heatmap_matrix_12 <- data.frame(matrix(ncol=21,nrow=39))
colnames(heatmap_matrix_12) <- 10:30
rownames(heatmap_matrix_12) <- ratios

heatmap_matrix_13 <- data.frame(matrix(ncol=21,nrow=39))
colnames(heatmap_matrix_13) <- 10:30
rownames(heatmap_matrix_13) <- ratios

heatmap_matrix_14 <- data.frame(matrix(ncol=21,nrow=39))
colnames(heatmap_matrix_14) <- 10:30
rownames(heatmap_matrix_14) <- ratios

heatmap_matrix_15 <- data.frame(matrix(ncol=21,nrow=39))
colnames(heatmap_matrix_15) <- 10:30
rownames(heatmap_matrix_15) <- ratios

heatmap_matrix_16 <- data.frame(matrix(ncol=21,nrow=39))
colnames(heatmap_matrix_16) <- 10:30
rownames(heatmap_matrix_16) <- ratios

heatmap_matrix_17 <- data.frame(matrix(ncol=21,nrow=39))
colnames(heatmap_matrix_17) <- 10:30
rownames(heatmap_matrix_17) <- ratios

heatmap_matrix_18 <- data.frame(matrix(ncol=21,nrow=39))
colnames(heatmap_matrix_18) <- 10:30
rownames(heatmap_matrix_18) <- ratios

heatmap_matrix_19 <- data.frame(matrix(ncol=21,nrow=39))
colnames(heatmap_matrix_19) <- 10:30
rownames(heatmap_matrix_19) <- ratios

heatmap_matrix_20 <- data.frame(matrix(ncol=21,nrow=39))
colnames(heatmap_matrix_20) <- 10:30
rownames(heatmap_matrix_20) <- ratios

heatmap_matrix_21 <- data.frame(matrix(ncol=21,nrow=39))
colnames(heatmap_matrix_21) <- 10:30
rownames(heatmap_matrix_21) <- ratios

heatmap_matrix_22 <- data.frame(matrix(ncol=21,nrow=39))
colnames(heatmap_matrix_22) <- 10:30
rownames(heatmap_matrix_22) <- ratios

heatmap_matrix_23 <- data.frame(matrix(ncol=21,nrow=39))
colnames(heatmap_matrix_23) <- 10:30
rownames(heatmap_matrix_23) <- ratios

heatmap_matrix_24 <- data.frame(matrix(ncol=21,nrow=39))
colnames(heatmap_matrix_24) <- 10:30
rownames(heatmap_matrix_24) <- ratios

heatmap_matrix_25 <- data.frame(matrix(ncol=21,nrow=39))
colnames(heatmap_matrix_25) <- 10:30
rownames(heatmap_matrix_25) <- ratios

heatmap_matrix_26 <- data.frame(matrix(ncol=21,nrow=39))
colnames(heatmap_matrix_26) <- 10:30
rownames(heatmap_matrix_26) <- ratios

heatmap_matrix_27 <- data.frame(matrix(ncol=21,nrow=39))
colnames(heatmap_matrix_27) <- 10:30
rownames(heatmap_matrix_27) <- ratios

heatmap_matrix_28 <- data.frame(matrix(ncol=21,nrow=39))
colnames(heatmap_matrix_28) <- 10:30
rownames(heatmap_matrix_28) <- ratios

heatmap_matrix_29 <- data.frame(matrix(ncol=21,nrow=39))
colnames(heatmap_matrix_29) <- 10:30
rownames(heatmap_matrix_29) <- ratios

heatmap_matrix_30 <- data.frame(matrix(ncol=21,nrow=39))
colnames(heatmap_matrix_30) <- 10:30
rownames(heatmap_matrix_30) <- ratios



##### Put empty matrices into a lsit 

matrix_list <- list(heatmap_matrix_10,heatmap_matrix_11,heatmap_matrix_12,heatmap_matrix_13,heatmap_matrix_14,heatmap_matrix_15,heatmap_matrix_16,heatmap_matrix_17,heatmap_matrix_18,heatmap_matrix_19,heatmap_matrix_20,heatmap_matrix_21,heatmap_matrix_22,heatmap_matrix_23,heatmap_matrix_24,heatmap_matrix_25,heatmap_matrix_26,heatmap_matrix_27,heatmap_matrix_28,heatmap_matrix_29,heatmap_matrix_30)



##### Fill: function that puts the score data into matrix format for a specified min frame length

fill <- function(df,z) {
  min_frame <- frame_list_scored[[z]]
  
  rownames <- row.names(df)
  colnames <- colnames(df)
  for (i in 1:nrow(df)) {
    for (j in 1:ncol(df)) {
      
      ratio <- as.numeric(rownames[i])
      length <- as.numeric(colnames[j])
      
      for (n in 1:nrow(min_frame)) {
        frame_length <- as.numeric(min_frame[n,2])
        pevk_ratio <- as.numeric(min_frame[n,3])
        score <- as.numeric(min_frame[n,10])
        if (length == frame_length && ratio == pevk_ratio) {
          found_score <- score
        }
      }
      df[i,j] <- found_score
    }
  }
  df
}



##### Create matrices for heatmaps:
heatmap_10 <- fill(matrix_list[[1]],1)
heatmap_10 <- data.matrix(heatmap_10)

heatmap_11 <- fill(matrix_list[[2]],2)
heatmap_11 <- data.matrix(heatmap_11)

heatmap_12 <- fill(matrix_list[[3]],3)
heatmap_12 <- data.matrix(heatmap_12)

heatmap_13 <- fill(matrix_list[[4]],4)
heatmap_13 <- data.matrix(heatmap_13)

heatmap_14 <- fill(matrix_list[[5]],5)
heatmap_14 <- data.matrix(heatmap_14)

heatmap_15 <- fill(matrix_list[[6]],6)
heatmap_15 <- data.matrix(heatmap_15)

heatmap_16 <- fill(matrix_list[[7]],7)
heatmap_16 <- data.matrix(heatmap_16)

heatmap_17 <- fill(matrix_list[[8]],8)
heatmap_17 <- data.matrix(heatmap_17)

heatmap_18 <- fill(matrix_list[[9]],9)
heatmap_18 <- data.matrix(heatmap_18)

heatmap_19 <- fill(matrix_list[[10]],10)
heatmap_19 <- data.matrix(heatmap_19)

heatmap_20 <- fill(matrix_list[[11]],11)
heatmap_20 <- data.matrix(heatmap_20)

heatmap_21 <- fill(matrix_list[[12]],12)
heatmap_21 <- data.matrix(heatmap_21)

heatmap_22 <- fill(matrix_list[[13]],13)
heatmap_22 <- data.matrix(heatmap_22)

heatmap_23 <- fill(matrix_list[[14]],14)
heatmap_23 <- data.matrix(heatmap_23)

heatmap_24 <- fill(matrix_list[[15]],15)
heatmap_24 <- data.matrix(heatmap_24)

heatmap_25 <- fill(matrix_list[[16]],16)
heatmap_25 <- data.matrix(heatmap_25)

heatmap_26 <- fill(matrix_list[[17]],17)
heatmap_26 <- data.matrix(heatmap_26)

heatmap_27 <- fill(matrix_list[[18]],18)
heatmap_27 <- data.matrix(heatmap_27)

heatmap_28 <- fill(matrix_list[[19]],19)
heatmap_28 <- data.matrix(heatmap_28)

heatmap_29 <- fill(matrix_list[[20]],21)
heatmap_29 <- data.matrix(heatmap_29)

heatmap_30 <- fill(matrix_list[[21]],21)
heatmap_30 <- data.matrix(heatmap_30)



##### Make contour maps #####
setwd("./PEVK_DATA/I/4/A")
library("RColorBrewer")

min_10 # 0.5306785
max_10 # 0.8347061
min_30 # 0.2772325
max_30 # 0.6137815
#Optimal contour map:

z_range <- c(0.2,0.3,0.4,0.5,0.6,0.70,0.75,0.80,0.85,0.87,0.89,0.91)
mypalette <- brewer.pal(11,"PiYG")
test_palette <- colorRampPalette(colors=c('#8e0152','#c51b7d','#de77ae','#f1b6da','#fde0ef','#f7f7f7','#e6f5d0','#b8e186','#7fbc41','#4d9221','#276419'), space="rgb", bias=0.3)

opt_x = as.numeric(as.character(colnames(heatmap_10)))
opt_y = as.numeric(as.character(rownames(heatmap_10)))
opt_z = heatmap_10
#png(filename="opt_contour_human.png", width=700, height=700, res=80, pointsize=17) #write to png
#pdf(file="opt_contour_human.pdf",width=16,height=16,pointsize=32) #write to pdf

#filled.contour(opt_y,opt_x,opt_z, col=mypalette, levels=z_range, xlab="PEVK Ratio", ylab="Sliding Frame Length") #main="Optimal Contour Plot (Human)")
#mtext("Minimum Frame Length 13",adj=0.2, padj=-0.5)
#dev.off() #close png

# Get key
#pdf(file="3d_contour_human_key.pdf",width=20,height=16,pointsize=20) #write to pdf
#contour3D(x=as.numeric(as.character(rownames(heatmap_10))), y=as.numeric(as.character(colnames(heatmap_10))), z=10, colvar=heatmap_10, col=test_palette(50), zlim=c(0,40), addbox=FALSE, alpha=0.5, theta=40, phi=35, nticks=10, ticktype="detailed")
#contour3D(x=as.numeric(as.character(rownames(heatmap_11))), y=as.numeric(as.character(colnames(heatmap_11))), z=11, colvar=heatmap_11, col=test_palette(50), add=TRUE, zlim=c(0,40), addbox=FALSE, colorkey = FALSE)
#dev.off()

# 3d
#pdf(file="3d_contour_human.pdf",width=20,height=16,pointsize=20) #write to pdf

op <- par(mar = c(5,7,4,2) + 0.1, mgp=c(2,1,0))
contour3D(x=as.numeric(as.character(rownames(heatmap_10))), y=as.numeric(as.character(colnames(heatmap_10))), z=10, colvar=heatmap_10, col=test_palette(50), zlim=c(0,40), addbox=FALSE, theta=40, phi=35, nticks=10, ticktype="detailed", xlab="PEVK Ratio",ylab="Window Length",zlab="Minimum Exon Length", colkey=FALSE)
contour3D(x=as.numeric(as.character(rownames(heatmap_11))), y=as.numeric(as.character(colnames(heatmap_11))), z=11, colvar=heatmap_11, col=test_palette(50), add=TRUE, zlim=c(0,40), addbox=FALSE, alpha=0.6, colkey = TRUE)
contour3D(x=as.numeric(as.character(rownames(heatmap_12))), y=as.numeric(as.character(colnames(heatmap_12))), z=12, colvar=heatmap_12, col=test_palette(50), add=TRUE, zlim=c(0,40), addbox=FALSE, alpha=0.6, colkey = FALSE)
contour3D(x=as.numeric(as.character(rownames(heatmap_13))), y=as.numeric(as.character(colnames(heatmap_13))), z=13, colvar=heatmap_13, col=test_palette(50), add=TRUE, zlim=c(0,40), addbox=FALSE, alpha=0.6, colkey = FALSE)
contour3D(x=as.numeric(as.character(rownames(heatmap_14))), y=as.numeric(as.character(colnames(heatmap_14))), z=14, colvar=heatmap_14, col=test_palette(50), add=TRUE, zlim=c(0,40), addbox=FALSE, alpha=0.6, colkey = FALSE)
contour3D(x=as.numeric(as.character(rownames(heatmap_15))), y=as.numeric(as.character(colnames(heatmap_15))), z=15, colvar=heatmap_15, col=test_palette(50), add=TRUE, zlim=c(0,40), addbox=FALSE, alpha=0.6, colkey = FALSE)
contour3D(x=as.numeric(as.character(rownames(heatmap_16))), y=as.numeric(as.character(colnames(heatmap_16))), z=16, colvar=heatmap_16, col=test_palette(50), add=TRUE, zlim=c(0,40), addbox=FALSE, alpha=0.6, colkey = FALSE)
contour3D(x=as.numeric(as.character(rownames(heatmap_17))), y=as.numeric(as.character(colnames(heatmap_17))), z=17, colvar=heatmap_17, col=test_palette(50), add=TRUE, zlim=c(0,40), addbox=FALSE, alpha=0.6, colkey = FALSE)
contour3D(x=as.numeric(as.character(rownames(heatmap_18))), y=as.numeric(as.character(colnames(heatmap_18))), z=18, colvar=heatmap_18, col=test_palette(50), add=TRUE, zlim=c(0,40), addbox=FALSE, alpha=0.6, colkey = FALSE)
contour3D(x=as.numeric(as.character(rownames(heatmap_19))), y=as.numeric(as.character(colnames(heatmap_19))), z=19, colvar=heatmap_19, col=test_palette(50), add=TRUE, zlim=c(0,40), addbox=FALSE, alpha=0.6, colkey = FALSE)
contour3D(x=as.numeric(as.character(rownames(heatmap_20))), y=as.numeric(as.character(colnames(heatmap_20))), z=20, colvar=heatmap_20, col=test_palette(50), add=TRUE, zlim=c(0,40), addbox=FALSE, alpha=0.6, colkey = FALSE)
contour3D(x=as.numeric(as.character(rownames(heatmap_21))), y=as.numeric(as.character(colnames(heatmap_21))), z=21, colvar=heatmap_21, col=test_palette(50), add=TRUE, zlim=c(0,40), addbox=FALSE, alpha=0.6, colkey = FALSE)
contour3D(x=as.numeric(as.character(rownames(heatmap_22))), y=as.numeric(as.character(colnames(heatmap_22))), z=22, colvar=heatmap_22, col=test_palette(50), add=TRUE, zlim=c(0,40), addbox=FALSE, alpha=0.6, colkey = FALSE)
contour3D(x=as.numeric(as.character(rownames(heatmap_23))), y=as.numeric(as.character(colnames(heatmap_23))), z=23, colvar=heatmap_23, col=test_palette(50), add=TRUE, zlim=c(0,40), addbox=FALSE, alpha=0.6, colkey = FALSE)
contour3D(x=as.numeric(as.character(rownames(heatmap_24))), y=as.numeric(as.character(colnames(heatmap_24))), z=24, colvar=heatmap_24, col=test_palette(50), add=TRUE, zlim=c(0,40), addbox=FALSE, alpha=0.6, colkey = FALSE)
contour3D(x=as.numeric(as.character(rownames(heatmap_25))), y=as.numeric(as.character(colnames(heatmap_25))), z=25, colvar=heatmap_25, col=test_palette(50), add=TRUE, zlim=c(0,40), addbox=FALSE, alpha=0.6, colkey = FALSE)
contour3D(x=as.numeric(as.character(rownames(heatmap_26))), y=as.numeric(as.character(colnames(heatmap_26))), z=26, colvar=heatmap_26, col=test_palette(50), add=TRUE, zlim=c(0,40), addbox=FALSE, alpha=0.6, colkey = FALSE)
contour3D(x=as.numeric(as.character(rownames(heatmap_27))), y=as.numeric(as.character(colnames(heatmap_27))), z=27, colvar=heatmap_27, col=test_palette(50), add=TRUE, zlim=c(0,40), addbox=FALSE, alpha=0.6, colkey = FALSE)
contour3D(x=as.numeric(as.character(rownames(heatmap_28))), y=as.numeric(as.character(colnames(heatmap_28))), z=28, colvar=heatmap_28, col=test_palette(50), add=TRUE, zlim=c(0,40), addbox=FALSE, alpha=0.6, colkey = FALSE)
contour3D(x=as.numeric(as.character(rownames(heatmap_29))), y=as.numeric(as.character(colnames(heatmap_29))), z=29, colvar=heatmap_29, col=test_palette(50), add=TRUE, zlim=c(0,40), addbox=FALSE, alpha=0.6, colkey = FALSE)
contour3D(x=as.numeric(as.character(rownames(heatmap_30))), y=as.numeric(as.character(colnames(heatmap_30))), z=30, colvar=heatmap_30, col=test_palette(50), add=TRUE, zlim=c(0,40), addbox=FALSE, alpha=0.6, colkey=FALSE)



#dev.off()


# Bad contour map:
bad_x = as.numeric(as.character(colnames(heatmap_30)))
bad_y = as.numeric(as.character(rownames(heatmap_30)))
bad_z = heatmap_30
#png(filename="bad_contour_human.png", width=700, height=700, res=80, pointsize=17) #write to png
#pdf(file="bad_contour_human.pdf",width=16,height=16, pointsize=32) #write to pdf
#filled.contour(bad_y,bad_x,bad_z, col=mypalette, levels=z_range, xlab="PEVK Ratio", ylab="Sliding Frame Length")#main="Sub-optimal Contour Plot (Human)")
#mtext("Minimum Frame Length 30",adj=0.2, padj=-0.5)
#dev.off()


######################################### EXON PLOTS ############################################

library(plotrix)
library(SDMTools)

###### Get test PEVK ratios

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios <- read.csv("human_stats.csv", header=F)
pevk_lengths_and_ratios <- t(pevk_lengths_and_ratios)
pevk_lengths_and_ratios <- data.frame(pevk_lengths_and_ratios)
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)
colnames(pevk_lengths_and_ratios) <- c("name","length","percent_pevk")

# Extract locations and append to pevk_lengths_and_ratios
new_starts = c()
new_ends = c()
for (i in 1:nrow(pevk_lengths_and_ratios))	{
  exon_name <- as.character(pevk_lengths_and_ratios[i,1])
  split <- strsplit(exon_name, "_")
  location <- split[[1]][3]
  location <- strsplit(location, ":")
  match_start <- as.numeric(location[[1]][1])
  match_end <- as.numeric(location[[1]][2])
  new_starts <- append(new_starts, match_start)
  new_ends <- append(new_ends, match_end)
}
# Add a new column called test_location to the data frame
pevk_lengths_and_ratios$start <- new_starts
pevk_lengths_and_ratios$end <- new_ends
pevk_lengths_and_ratios <- pevk_lengths_and_ratios[order(pevk_lengths_and_ratios$start),]
rownames(pevk_lengths_and_ratios) <- 1:nrow(pevk_lengths_and_ratios)

#lengths_and_ratios <- data.frame(matrix(ncol=4,nrow=0))
#colnames(lengths_and_ratios) <- c("known_exon", "test_exon", "pevk_ratio", "aa_length")
#lengths_and_ratios$known_exon

# Extract ratios from the data frame
test_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios$percent_pevk))
#ratio_mean <- mean(test_pevk_ratios)

## Assign colors to each PEVK ratio
#rbpal <- colorRampPalette(c("blue","red"))
#test_exon_cols <- rbpal(200)[as.numeric(cut(test_pevk_ratios, breaks=200))]
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
test_exon_cols <- pal[as.numeric(cut(test_pevk_ratios, breaks=8))]

###### Get known PEVK ratios

# Read in csv containing info on length and % pevk
pevk_lengths_and_ratios_known <- read.csv("human_known_stats.csv", header=F)
pevk_lengths_and_ratios_known <- t(pevk_lengths_and_ratios_known)
pevk_lengths_and_ratios_known <- data.frame(pevk_lengths_and_ratios_known)
rownames(pevk_lengths_and_ratios_known) <- 1:nrow(pevk_lengths_and_ratios_known)
colnames(pevk_lengths_and_ratios_known) <- c("name","length","percent_pevk")


# Extract ratios from the data frame
known_pevk_ratios <- as.numeric(as.character(pevk_lengths_and_ratios_known$percent_pevk))
#ratio_mean <- mean(test_pevk_ratios)

## Assign colors to each PEVK ratio
#rbpal <- colorRampPalette(c("blue","red"))
#known_exon_cols <- rbpal(200)[as.numeric(cut(known_pevk_ratios, breaks=200))]
pal <- brewer.pal(11,"PiYG")[c(1:4,8:11)]
known_exon_cols <- pal[as.numeric(cut(known_pevk_ratios, breaks=8))]

#known_pevk_ratios,c(0,1,1),c(1,1,0),0
########### Read in the file containing the locations of the found exons with optimized parameters

my_exons <- read.csv("test_locations_10_0.54_10.csv", header=F)
my_exons <- t(my_exons)
my_exons <- data.frame(my_exons)
colnames(my_exons) <- c("frame","start","end","length")
my_exons$start <- as.numeric(as.character(my_exons$start))
my_exons$end <- as.numeric(as.character(my_exons$end))
my_exons$length <- as.numeric(as.character(my_exons$length))
my_exons <- my_exons[order(my_exons$start),]
rownames(my_exons) <- 1:nrow(my_exons)

# Read in the file containing the locations of the known exons, relative to the gene start
known_exons <- read.csv("known_locations_human.csv")
known_exons <- data.frame(known_exons)

############ Reading in the files containing the locations and PEVK ratios of the GENSCAN Exons

# Starts and ends
genscan_exons <- read.csv("genscan_locations_human.csv")
genscan_exons <- data.frame(genscan_exons)

# PEVK ratios
genscan_ratios <- read.csv("genscan_human_ratios.csv", header=F)
genscan_ratios <- t(genscan_ratios)
colnames(genscan_ratios) <- c("ratio")
genscan_cols <- pal[as.numeric(cut(genscan_ratios, breaks=8))]

####### Create rectangle objects for each exon coordinate #######


# Known exons
known_starts <- as.character(known_exons$rel_start)
known_starts <-as.numeric(gsub(",","",known_starts))
known_ends <- as.character(known_exons$rel_end)
known_ends <-as.numeric(gsub(",","",known_ends))

# My exons
my_starts <- as.numeric(as.character(my_exons$start))
my_ends <- as.numeric(as.character(my_exons$end))

# GENSCAN exons
genscan_starts <- as.numeric(as.character(genscan_exons$start))
genscan_ends <- as.numeric(as.character(genscan_exons$end))

#png(filename="human_exon_plot.png", width=1000, height=700)
#pdf(file="~/Desktop/titin_project/pevk_ntsplice_human/human_analysis/human_exon_plot2.pdf", width=25, height=13, pointsize=22)

# Set up plot region
op <- par(bg = "white", bty="n")
plot(c(106000,171000), c(0, 4.5), type="n", xlab="Titin Coordinates", ylab="", yaxt="n", axes=FALSE)

# Add known exons to the diagram (above)
rect(known_starts,3.3,known_ends,4.2, col=known_exon_cols, border=NA)
abline(h=3.75) # tranverse line

# Add my exons to the diagram (below)
rect(my_starts,1.8,my_ends,2.7, col=test_exon_cols, border=NA)
abline(h=2.25) # transverse line

# Add genscan exons to the diagram
rect(genscan_starts,0.3,genscan_ends,1.2, col=genscan_cols, border=NA)
abline(h=0.75) # transverse line

# Add labels to the x axis:
axis(side=1, at=c(105000,110000,115000,120000,125000,130000,135000,140000,145000,150000,155000,160000,165000,170000))

# Add labels to the y axis:
axis(side=2, lty=0, at=c(0.75, 2.25, 3.75), labels=c("GENSCAN Exons", "PEVK Finder Exons","Annotated Exons"))

# Add color legend:
#legend(100000,100000, legend=c("Low PEVK", "High PEVK"), col=c("blue","red"))

#dev.off()


####### With GENSCAN Exons ######

#op <- par(bg = "white")
#plot(c(108000,170500), c(0, 2.5), type="n", xlab="Titin Coordinates", ylab="", main="PEVK Exon Distribution (Human)", yaxt="n")

# Add known exons to the diagram (above)
#rect(known_starts,1.5,known_ends,2, col=known_exon_cols, border=NA)
#abline(h=1.75) # tranverse line

# Add my exons to the diagram (below)
#rect(genscan_starts,0.5,genscan_ends,1, col="black", border=NA)
#abline(h=0.75) # transverse line

# Add labels to the y axis:
#axis(side=2, at=c(0.75, 1.75), labels=c("GENSCAN Exons", "Annotated Exons"))

############################# Plotting differences between subject and query #######################

# Read in the csv file containing blast data
blast.data <- read.csv("final_blast.csv", header=F)

# Set up the csv file
data_cols <- c("test_exon", "known_exon", "%identity", "alignment_length", "gaps", "misalignments", "test_start", "test_end", "known_start", "known_end", "bit_score", "certainty")
colnames(blast.data) <- data_cols

# Split each tets exon name to get start and end info
test_starts = c()
test_ends = c()
for (i in 1:nrow(blast.data))	{
  exon_name <- as.character(blast.data[i,1])
  split <- strsplit(exon_name, "_")
  location <- split[[1]][3]
  location <- strsplit(location, ":")
  match_start <- as.numeric(location[[1]][1])
  match_end <- as.numeric(location[[1]][2])
  test_starts <- append(test_starts, match_start)
  test_ends <- append(test_ends, match_end)
}
# Add a new column called test_location to the data frame
blast.data$start <- test_starts
blast.data$end <- test_ends
blast.data <- blast.data[order(blast.data$start, blast.data$bit_score),]
rownames(blast.data) <- 1:nrow(blast.data)

# Make an empty data frame to hold the exon matches
more_col_names <- c("known_exon", "test_exon", "alignment_length", "test_start", "test_end", "known_start", "known_end", "bit_score", "test_start","test_end") # specify the column names of exon_matches
exon_matches <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(exon_matches) <- more_col_names

###### Filling the exon_matches data frame ######

# Initialize exon_matches:
exon_matches[1,1] <- blast.data[1,2]
exon_matches[1,2] <- toString(blast.data[1,1])
exon_matches[1,3] <- blast.data[1,4]
exon_matches[1,4] <- blast.data[1,7]
exon_matches[1,5] <- blast.data[1,8]
exon_matches[1,6] <- blast.data[1,9]
exon_matches[1,7] <- blast.data[1,10]
exon_matches[1,8] <- blast.data[1,11]
exon_matches[1,9] <- blast.data[1,13]
exon_matches[1,10] <- blast.data[1,14]


# Get the length of the csv file:
file_length <- nrow(blast.data)

# Initialize the row number for exon_matches entries
current_row <- 1

# Loop through each row of blast.data:
for (i in 2:(file_length-1)) {
  
  # Get the exon number:
  exon_number <- blast.data[i,2]
  
  # Get the start location:
  test_start <- blast.data[i,13]
  
  # Get the end location:
  test_end <- blast.data[i,14]
  
  # If the start location is greater than the previously entered start location:
  if (test_start > blast.data[i-1,13] && test_start > blast.data[i-1,14] && test_start != blast.data[i+1,13]) {
    if (exon_number > exon_matches[current_row,1]) {
      
      # Increment the row number by 1
      current_row = current_row + 1
      
      # Add all relevant data
      exon_matches[current_row,1] <- exon_number
      exon_matches[current_row,2] <- toString(blast.data[i,1])
      exon_matches[current_row,3] <- blast.data[i,4]
      exon_matches[current_row,4] <- blast.data[i,7]
      exon_matches[current_row,5] <- blast.data[i,8]
      exon_matches[current_row,6] <- blast.data[i,9]
      exon_matches[current_row,7] <- blast.data[i,10]
      exon_matches[current_row,8] <- blast.data[i,11]
      exon_matches[current_row,9] <- blast.data[i,13]
      exon_matches[current_row,10] <- blast.data[i,14]
    }
    if (exon_number == exon_matches[current_row,1] && test_start != blast.data[i-1,13]) {
      bit1 <- blast.data[i-1,11]
      bit2 <- blast.data[i,11]
      if (bit2 > bit1) {
        exon_matches[current_row,1] <- exon_number
        exon_matches[current_row,2] <- toString(blast.data[i,1])
        exon_matches[current_row,3] <- blast.data[i,4]
        exon_matches[current_row,4] <- blast.data[i,7]
        exon_matches[current_row,5] <- blast.data[i,8]
        exon_matches[current_row,6] <- blast.data[i,9]
        exon_matches[current_row,7] <- blast.data[i,10]
        exon_matches[current_row,8] <- blast.data[i,11]
        exon_matches[current_row,9] <- blast.data[i,13]
        exon_matches[current_row,10] <- blast.data[i,14]
      }
      else {
        next
      }
    }
  }
  
  # If the start location is equal to the previously entered start location;
  if (test_start > blast.data[i-1,13] && test_start == blast.data[i+1,13]) {
    
    # make a list of all the exons that matched that location
    row_start <- i
    subset <- c()
    for (q in row_start:file_length) {
      if (blast.data[q,13] == test_start) {
        subset <- c(subset, blast.data[q,2])
      }
    }
    # find differences between the possible exons and the previously entered exon
    differences <- c()
    for (r in subset) {
      diff <- r - exon_matches[current_row,1]
      if (diff > 0) {
        differences <- c(differences, diff)
      }
      if (diff < 0) {
        differences <- c(differences, 200)
      }
    }
    # find the location of the closest exon
    min_location <- which.min(differences)
    min_row <- row_start + min_location-1
    
    # Increment the row number by 1
    current_row = current_row + 1
    
    # add the info for that exon into the current row of exon_matches
    exon_matches[current_row,1] <- blast.data[min_row,2]
    exon_matches[current_row,2] <- toString(blast.data[min_row,1])
    exon_matches[current_row,3] <- blast.data[min_row,4]
    exon_matches[current_row,4] <- blast.data[min_row,7]
    exon_matches[current_row,5] <- blast.data[min_row,8]
    exon_matches[current_row,6] <- blast.data[min_row,9]
    exon_matches[current_row,7] <- blast.data[min_row,10]
    exon_matches[current_row,8] <- blast.data[min_row,11]
    exon_matches[current_row,9] <- blast.data[min_row,13]
    exon_matches[current_row,10] <- blast.data[min_row,14]
    
  }
  if (test_start == blast.data[i-1,13] && test_start == blast.data[i+1,13]) {
    next
  }
  if (test_start == blast.data[i-1,13] && test_start < blast.data[i+1,13]) {
    next
  }
}
if (blast.data[file_length,13] > blast.data[file_length-1,13]) {
  exon_matches[current_row+1,1] <- blast.data[file_length,2]
  exon_matches[current_row+1,2] <- toString(blast.data[file_length,1])
  exon_matches[current_row+1,3] <- blast.data[file_length,4]
  exon_matches[current_row+1,4] <- blast.data[file_length,7]
  exon_matches[current_row+1,5] <- blast.data[file_length,8]
  exon_matches[current_row+1,6] <- blast.data[file_length,9]
  exon_matches[current_row+1,7] <- blast.data[file_length,10]
  exon_matches[current_row+1,8] <- blast.data[file_length,11]
  exon_matches[current_row+1,9] <- blast.data[file_length,13]
  exon_matches[current_row+1,10] <- blast.data[file_length,14]
}

######## Lengths:

# Extract test lengths from the data frame:
#test_aa_lengths <- as.numeric(as.character(pevk_lengths_and_ratios$length))

# Extract known lengths from the known data frame:
#known_aa_lengths <- as.numeric(as.character(pevk_lengths_and_ratios_known$length))

############## % Diff in Subj - Query: Overestimates, Underestimates and Perfect Matches ###########
matched_exons <- c(exon_matches$known_exon+111)
differences <- c()
known_lengths <- c()
test_lengths <- c()
for (i in 1:nrow(known_exons)) {
  exon <- as.numeric(known_exons[i,1])
  known_length <- as.numeric(known_exons[i,6])
  
  if (is.element(exon, matched_exons) == TRUE) {
    for (j in 1:nrow(exon_matches)) {
      test_exon <- as.numeric(exon_matches[j,1]) + 111
      test_length <- as.numeric(exon_matches[j,10]) - as.numeric(exon_matches[j,9]) + 1
      if (test_exon == exon) {
        differences[i] <- (test_length-known_length)/((test_length+known_length)/2)
        known_lengths[i] <- known_length
        test_lengths[i] <- test_length
      }
    }	
  }
  if (is.element(exon, matched_exons) == FALSE) {
    differences[i] <- 0
    known_lengths[i] <- known_length
    test_lengths[i] <- 0
  }
}

# 1) Plot differences
# Specify rectangle starts and ends:
rec_starts <- seq(from=111.5, to=223.5, by=1)
rec_ends <- seq(from=112.5, to=224.5, by=1)

# Make plotting area
#png(filename="lengths_diffs_human.png", width=1000, height=700, res=80)
#pdf(file="lengths_diffs_human.pdf", width=25, height=11, pointsize=22)

par(mar=c(5,5,4,4) + 0.1)
op <- par(bg = "white")
plot(c(112,224), c(-1.5, 1.5), type="n", xlab="Exon Number", ylab="% Difference from Subject", main="Subject-Query Lengths and Percent Differences (Human)")

rect(rec_starts,differences, rec_ends, 0, col=rgb(197,27,125, maxColorValue=200))

# Plot known lengths:
par(new=TRUE)
plot(as.numeric(c(112:224)),known_lengths, pch=16, cex=1, col=rgb(77,146,33,maxColorValue=200), axes=FALSE, bty="n", xlab="",ylab="", ylim=c(-400,400))
axis(side=4, c(0,200,400))
mtext("Exon Length (nt)", side=4,line=2)
# Plot test lengths:
points(as.numeric(c(112:224)),test_lengths, pch=1, cex=1, col="black")
#dev.off()