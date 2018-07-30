# Reconstructs ancestral states of amino acids in selected codons

# Load packages
library(ape)
library(geiger)
library(phytools)

# Write functiont to assign charcteristics to each amino acid
assignCharacteristics <- function(codon_data){
  new_codon_data <- factor(vector(mode="character", length=length(codon_data)))
  levels(new_codon_data)=c('hydrophobic','polar','neg_charged','pos_charged')
  names(new_codon_data) <- names(codon_data)
  for (i in 1:length(codon_data)){
    if (as.character(codon_data[i]) %in% c('A','I','L','M','F','V','P','G','W')){
      new_codon_data[i] <- 'hydrophobic'
    }
    if (as.character(codon_data[i]) %in% c('Q','N','S','T','Y','C')){
      new_codon_data[i] <- 'polar'
    }
    if (as.character(codon_data[i]) %in% c('D','E')){
      new_codon_data[i] <- 'neg_charged'
    }
    if (as.character(codon_data[i]) %in% c('H','K','R')){
      new_codon_data[i] <- 'pos_charged'
    }
  }
  return(new_codon_data)
}

###################################### 114 #######################################

# Read tree
exon_114_tree = read.tree("frame_f1_111156_111557.tre")

# Get codon data
exon_114_data = read.csv("frame_f1_111156_111557_E114_codons.csv", header = TRUE)
rownames(exon_114_data) = exon_114_data$species

# Compare tree and data
compare <- treedata(exon_114_tree, exon_114_data, sort=TRUE) 

# Revise data/tree
exon_114_data <- as.data.frame(compare$data)
exon_114_data[,2:3] <- lapply(exon_114_data[,2:3], as.character)
exon_114_tree <- compare$phy

# save new tree and load it into r
write.tree(exon_114_tree, "exon_114.tre")
exon_114_tree <- read.tree("exon_114.tre")

############# Exon 114, Codon 1
z <- as.factor(exon_114_data$codon1)
names(z) <- exon_114_tree$tip.label
assigned_z <- assignCharacteristics(z)

fitER <- rerootingMethod(exon_114_tree, assigned_z, model = "ER")

pdf(file="E114_C1_selection.pdf", width=20, height=20, pointsize=30)
par(mar=c(1,0,0,6), xpd = T)
plot(exon_114_tree, label.offset = 1.5)
nodelabels(node = as.numeric(rownames(fitER$marginal.anc)), pie = fitER$marginal.anc, 
           piecol = c("blue", "yellow", "red","orange"), cex = 0.6)
tiplabels(pie = to.matrix(assigned_z, sort(unique(assigned_z))), piecol = c("blue", "yellow", "red","orange"), 
          cex = 0.3)
xlim = par("usr")[2]
ylim= par("usr")[4]
add.simmap.legend(colors=setNames(c("blue", "yellow", "red"),c("Hydrophobic","Polar","Pos. Charged")),
                  prompt=FALSE,shape="circle",x=xlim-5,y=ylim-1.1)
dev.off()

############# Exon 114, Codon 2 ---------> all hydrophobic
z <- as.factor(exon_114_data$codon2)
names(z) <- exon_114_tree$tip.label
assigned_z <- assignCharacteristics(z)

fitER <- rerootingMethod(exon_114_tree, assigned_z, model = "ER")
par(mar=c(1,0,0,6), xpd = T)
plot(exon_114_tree, label.offset = 1.5)
nodelabels(node = as.numeric(rownames(fitER$marginal.anc)), pie = fitER$marginal.anc, 
           piecol = c("blue", "red", "yellow","green"), cex = 0.6)
tiplabels(pie = to.matrix(assigned_z, sort(unique(assigned_z))), piecol = c("blue", "red", "yellow","green"), 
          cex = 0.3)
xlim = par("usr")[2]
ylim= par("usr")[4]
add.simmap.legend(colors=setNames(c("blue","red"),c("hydrophobic","polar")),
                  prompt=FALSE,shape="circle",x=xlim-5,y=ylim-0.8)










###################################### 135 #######################################

# Read tree
exon_135_tree = read.tree("frame_f2_124495_124725.tre")

# Get codon data
exon_135_data = read.csv("frame_f2_124495_124725_E135_codons.csv", header = TRUE)
rownames(exon_135_data) = exon_135_data$species

# Compare tree and data
compare <- treedata(exon_135_tree, exon_135_data, sort=TRUE) 

# Revise data/tree
exon_135_data <- as.data.frame(compare$data)
exon_135_data[,2:3] <- lapply(exon_135_data[,2:3], as.character)
exon_135_tree <- compare$phy

# save new tree and load it into r
write.tree(exon_135_tree, "exon_135.tre")
exon_135_tree <- read.tree("exon_135.tre")



############# Exon 135, Codon 1
z <- as.factor(exon_135_data$codon1)
names(z) <- exon_135_tree$tip.label
assigned_z <- assignCharacteristics(z)

fitER <- rerootingMethod(exon_135_tree, assigned_z, model = "ER")

pdf(file="E135_C1_selection.pdf", width=20, height=20, pointsize=30)
par(mar=c(1,0,0,6), xpd = T)
plot(exon_135_tree, label.offset = 1.5)
nodelabels(node = as.numeric(rownames(fitER$marginal.anc)), pie = fitER$marginal.anc, 
           piecol = c("blue", "red", "yellow","green"), cex = 0.6)
tiplabels(pie = to.matrix(assigned_z, sort(unique(assigned_z))), piecol = c("blue", "red", "yellow","green"), 
          cex = 0.3)
# Get xlim and y lim for legend justification
xlim = par("usr")[2]
ylim= par("usr")[4]
add.simmap.legend(colors=setNames(c("blue","red"),c("Hydrophobic","Polar")),
                  prompt=FALSE,shape="circle",x=xlim-5,y=ylim-0.8)
dev.off()

############# Exon 135, Codon 2
z <- as.factor(exon_135_data$codon2)
names(z) <- exon_135_tree$tip.label
assigned_z <- assignCharacteristics(z)

fitER <- rerootingMethod(exon_135_tree, assigned_z, model = "ER")

pdf(file="E135_C2_selection.pdf", width=20, height=20, pointsize=30)
par(mar=c(1,0,0,6), xpd = T)
plot(exon_135_tree, label.offset = 1.5)
nodelabels(node = as.numeric(rownames(fitER$marginal.anc)), pie = fitER$marginal.anc, 
           piecol = c("blue", "yellow", "red","green"), cex = 0.6)
tiplabels(pie = to.matrix(assigned_z, sort(unique(assigned_z))), piecol = c("blue", "yellow", "red","green"), 
          cex = 0.3)
# Get xlim and y lim for legend justification
xlim = par("usr")[2]
ylim= par("usr")[4]
add.simmap.legend(colors=setNames(c("blue", "yellow", "red"),c("Hydrophobic","Polar","Pos. Charged")),
                  prompt=FALSE,shape="circle",x=xlim-5,y=ylim-0.8)
dev.off()

############# Exon 135, Codon 3
z <- as.factor(exon_135_data$codon3)
names(z) <- exon_135_tree$tip.label
assigned_z <- assignCharacteristics(z)

fitER <- rerootingMethod(exon_135_tree, assigned_z, model = "ER")

pdf(file="E135_C3_selection.pdf", width=20, height=20, pointsize=30)
par(mar=c(1,0,0,6), xpd = T)
plot(exon_135_tree, label.offset = 1.5)
nodelabels(node = as.numeric(rownames(fitER$marginal.anc)), pie = fitER$marginal.anc, 
           piecol = c("blue", "red", "yellow","green"), cex = 0.6)
tiplabels(pie = to.matrix(assigned_z, sort(unique(assigned_z))), piecol = c("blue", "red", "yellow","green"), 
          cex = 0.3)
# Get xlim and y lim for legend justification
xlim = par("usr")[2]
ylim= par("usr")[4]
add.simmap.legend(colors=setNames(c("blue","red"),c("Hydrophobic","Polar")),
                  prompt=FALSE,shape="circle",x=xlim-5,y=ylim-0.8)
dev.off()

############# Exon 135, Codon 4
z <- as.factor(exon_135_data$codon4)
names(z) <- exon_135_tree$tip.label
assigned_z <- assignCharacteristics(z)

fitER <- rerootingMethod(exon_135_tree, assigned_z, model = "ER")

pdf(file="E135_C4_selection.pdf", width=20, height=20, pointsize=30)
par(mar=c(1,0,0,6), xpd = T)
plot(exon_135_tree, label.offset = 1.5)
nodelabels(node = as.numeric(rownames(fitER$marginal.anc)), pie = fitER$marginal.anc, 
           piecol = c("blue", "red", "yellow","green"), cex = 0.6)
tiplabels(pie = to.matrix(assigned_z, sort(unique(assigned_z))), piecol = c("blue", "red", "yellow","green"), 
          cex = 0.3)
# Get xlim and y lim for legend justification
xlim = par("usr")[2]
ylim= par("usr")[4]
add.simmap.legend(colors=setNames(c("blue","red"),c("Hydrophobic","Polar")),
                  prompt=FALSE,shape="circle",x=xlim-5,y=ylim-0.8)
dev.off()












###################################### 137 #######################################
# Read tree
exon_137_tree = read.tree("frame_f3_125975_126046.tre")

# Get codon data
exon_137_data = read.csv("frame_f3_125975_126046_E137_codons.csv", header = TRUE)
rownames(exon_137_data) = exon_137_data$species

# Compare tree and data
compare <- treedata(exon_137_tree, exon_137_data, sort=TRUE) 

# Revise data/tree
exon_137_data <- as.data.frame(compare$data)
#exon_137_data[,2:3] <- lapply(exon_137_data[,2], as.character)
exon_137_tree <- compare$phy

# save new tree and load it into r
write.tree(exon_137_tree, "exon_137.tre")
exon_137_tree <- read.tree("exon_137.tre")


############ Exon 137, Codon 1 --------------------> all hydrophobic
z <- as.factor(exon_137_data$codon1)
names(z) <- exon_137_tree$tip.label
assigned_z <- assignCharacteristics(z)

fitER <- rerootingMethod(exon_137_tree, assigned_z, model = "ER")
par(mar=c(1,0,0,6), xpd = T)
plot(exon_137_tree, label.offset = 1.5)
nodelabels(node = as.numeric(rownames(fitER$marginal.anc)), pie = fitER$marginal.anc, 
           piecol = c("blue", "red", "yellow","green"), cex = 0.6)
tiplabels(pie = to.matrix(assigned_z, sort(unique(assigned_z))), piecol = c("blue", "red", "yellow","green"), 
          cex = 0.3)
# Get xlim and y lim for legend justification
xlim = par("usr")[2]
ylim= par("usr")[4]
add.simmap.legend(colors=setNames(c("blue","red","yellow","green"),c("hydrophobic","polar","pos_charged","neg_charged")),
                  prompt=FALSE,shape="circle",x=xlim-5,y=ylim-0.8)












###################################### 138 #######################################
# Read tree
exon_138_tree = read.tree("frame_f3_126254_126343.tre")

# Get codon data
exon_138_data = read.csv("frame_f3_126254_126343_E138_codons.csv", header = TRUE)
rownames(exon_138_data) = exon_138_data$species

# Compare tree and data
compare <- treedata(exon_138_tree, exon_138_data, sort=TRUE) 

# Revise data/tree
exon_138_data <- as.data.frame(compare$data)
#exon_137_data[,2:3] <- lapply(exon_137_data[,2], as.character)
exon_138_tree <- compare$phy

# save new tree and load it into r
write.tree(exon_138_tree, "exon_138.tre")
exon_138_tree <- read.tree("exon_138.tre")


############ Exon 138, Codon 1
z <- as.factor(exon_138_data$codon1)
names(z) <- exon_138_tree$tip.label
assigned_z <- assignCharacteristics(z)

fitER <- rerootingMethod(exon_138_tree, assigned_z, model = "ER")

pdf(file="E138_C1_selection.pdf", width=20, height=20, pointsize=30)
par(mar=c(1,0,0,6), xpd = T)
plot(exon_138_tree, label.offset = 1.5)
nodelabels(node = as.numeric(rownames(fitER$marginal.anc)), pie = fitER$marginal.anc, 
           piecol = c("blue", "red", "yellow","green"), cex = 0.6)
tiplabels(pie = to.matrix(assigned_z, sort(unique(assigned_z))), piecol = c("blue", "red", "yellow","green"), 
          cex = 0.3)
# Get xlim and y lim for legend justification
xlim = par("usr")[2]
ylim= par("usr")[4]
add.simmap.legend(colors=setNames(c("blue","red"),c("Hydrophobic","Polar")),
                  prompt=FALSE,shape="circle",x=xlim-5,y=ylim-0.65)
dev.off()


############ Exon 138, Codon 2
z <- as.factor(exon_138_data$codon2)
names(z) <- exon_138_tree$tip.label
assigned_z <- assignCharacteristics(z)

fitER <- rerootingMethod(exon_138_tree, assigned_z, model = "ER")

pdf(file="E138_C2_selection.pdf", width=20, height=20, pointsize=30)
par(mar=c(1,0,0,6), xpd = T)
plot(exon_138_tree, label.offset = 1.5)
nodelabels(node = as.numeric(rownames(fitER$marginal.anc)), pie = fitER$marginal.anc, 
           piecol = c("blue", "red", "yellow","green"), cex = 0.6)
tiplabels(pie = to.matrix(assigned_z, sort(unique(assigned_z))), piecol = c("blue", "red", "yellow","green"), 
          cex = 0.3)
# Get xlim and y lim for legend justification
xlim = par("usr")[2]
ylim= par("usr")[4]
add.simmap.legend(colors=setNames(c("blue","red"),c("Hydrophobic","Polar")),
                  prompt=FALSE,shape="circle",x=xlim-5,y=ylim-0.7)
dev.off()

############ Exon 138, Codon 3 ------------------> all positively charged!
z <- as.factor(exon_138_data$codon3)
names(z) <- exon_138_tree$tip.label
assigned_z <- assignCharacteristics(z)

fitER <- rerootingMethod(exon_138_tree, assigned_z, model = "ER")
par(mar=c(1,0,0,6), xpd = T)
plot(exon_138_tree, label.offset = 1.5)
nodelabels(node = as.numeric(rownames(fitER$marginal.anc)), pie = fitER$marginal.anc, 
           piecol = c("blue", "red", "yellow","green"), cex = 0.6)
tiplabels(pie = to.matrix(assigned_z, sort(unique(assigned_z))), piecol = c("blue", "red", "yellow","green"), 
          cex = 0.3)
# Get xlim and y lim for legend justification
xlim = par("usr")[2]
ylim= par("usr")[4]
add.simmap.legend(colors=setNames(c("blue","red","yellow","green"),c("hydrophobic","polar","pos_charged","neg_charged")),
                  prompt=FALSE,shape="circle",x=xlim-5,y=ylim-0.8)

############ Exon 138, Codon 4
z <- as.factor(exon_138_data$codon4)
names(z) <- exon_138_tree$tip.label
assigned_z <- assignCharacteristics(z)

fitER <- rerootingMethod(exon_138_tree, assigned_z, model = "ER")

pdf(file="E138_C4_selection.pdf", width=20, height=20, pointsize=30)
par(mar=c(1,0,0,6), xpd = T)
plot(exon_138_tree, label.offset = 1.5)
nodelabels(node = as.numeric(rownames(fitER$marginal.anc)), pie = fitER$marginal.anc, 
           piecol = c("blue", "red", "yellow","green"), cex = 0.6)
tiplabels(pie = to.matrix(assigned_z, sort(unique(assigned_z))), piecol = c("blue", "red", "yellow","green"), 
          cex = 0.3)
# Get xlim and y lim for legend justification
xlim = par("usr")[2]
ylim= par("usr")[4]
add.simmap.legend(colors=setNames(c("blue","red"),c("Hydrophobic","Polar")),
                  prompt=FALSE,shape="circle",x=xlim-5,y=ylim-0.7)
dev.off()

############ Exon 138, Codon 5
z <- as.factor(exon_138_data$codon5)
names(z) <- exon_138_tree$tip.label
assigned_z <- assignCharacteristics(z)

fitER <- rerootingMethod(exon_138_tree, assigned_z, model = "ER")

pdf(file="E138_C5_selection.pdf", width=20, height=20, pointsize=30)
par(mar=c(1,0,0,6), xpd = T)
plot(exon_138_tree, label.offset = 1.5)
nodelabels(node = as.numeric(rownames(fitER$marginal.anc)), pie = fitER$marginal.anc, 
           piecol = c("blue", "red", "yellow","green"), cex = 0.6)
tiplabels(pie = to.matrix(assigned_z, sort(unique(assigned_z))), piecol = c("blue", "red", "yellow","green"), 
          cex = 0.3)
# Get xlim and y lim for legend justification
xlim = par("usr")[2]
ylim= par("usr")[4]
add.simmap.legend(colors=setNames(c("blue","red"),c("Hydrophobic","Polar")),
                  prompt=FALSE,shape="circle",x=xlim-5,y=ylim-0.7)
dev.off()












###################################### 145 #######################################
# Read tree
exon_145_tree = read.tree("frame_f2_129217_129297.tre")

# Get codon data
exon_145_data = read.csv("frame_f2_129217_129297_E145_codons.csv", header = TRUE)
rownames(exon_145_data) = exon_145_data$species

# Compare tree and data
compare <- treedata(exon_145_tree, exon_145_data, sort=TRUE) 

# Revise data/tree
exon_145_data <- as.data.frame(compare$data)
#exon_137_data[,2:3] <- lapply(exon_137_data[,2], as.character)
exon_145_tree <- compare$phy

# save new tree and load it into r
write.tree(exon_145_tree, "exon_145.tre")
exon_145_tree <- read.tree("exon_145.tre")

############# Exon 145, Codon 1 --------------------------> all hydrophobic!
z <- as.factor(exon_145_data$codon1)
names(z) <- exon_145_tree$tip.label
assigned_z <- assignCharacteristics(z)

fitER <- rerootingMethod(exon_145_tree, assigned_z, model = "ER")

par(mar=c(1,0,0,6), xpd = T)
plot(exon_145_tree, label.offset = 1.5)
nodelabels(node = as.numeric(rownames(fitER$marginal.anc)), pie = fitER$marginal.anc, 
           piecol = c("blue", "red", "yellow","green"), cex = 0.6)
tiplabels(pie = to.matrix(assigned_z, sort(unique(assigned_z))), piecol = c("blue", "red", "yellow","green"), 
          cex = 0.3)
# Get xlim and y lim for legend justification
xlim = par("usr")[2]
ylim= par("usr")[4]
add.simmap.legend(colors=setNames(c("blue","red","yellow","green"),c("hydrophobic","polar","pos_charged","neg_charged")),
                  prompt=FALSE,shape="circle",x=xlim-5,y=ylim-0.8)



############## Exon 145, Codon 2 ------------------------------> all hydrophobic!
z <- as.factor(exon_145_data$codon2)
names(z) <- exon_145_tree$tip.label
assigned_z <- assignCharacteristics(z)

fitER <- rerootingMethod(exon_145_tree, assigned_z, model = "ER")
par(mar=c(1,0,0,6), xpd = T)
plot(exon_145_tree, label.offset = 1.5)
nodelabels(node = as.numeric(rownames(fitER$marginal.anc)), pie = fitER$marginal.anc, 
           piecol = c("blue", "red", "yellow","green"), cex = 0.6)
tiplabels(pie = to.matrix(assigned_z, sort(unique(assigned_z))), piecol = c("blue", "red", "yellow","green"), 
          cex = 0.3)
# Get xlim and y lim for legend justification
xlim = par("usr")[2]
ylim= par("usr")[4]
add.simmap.legend(colors=setNames(c("blue","red","yellow","green"),c("hydrophobic","polar","pos_charged","neg_charged")),
                  prompt=FALSE,shape="circle",x=xlim-5,y=ylim-0.8)


############# Exon 145, Codon 3
z <- as.factor(exon_145_data$codon3)
names(z) <- exon_145_tree$tip.label
assigned_z <- assignCharacteristics(z)

fitER <- rerootingMethod(exon_145_tree, assigned_z, model = "ER")

pdf(file="E145_C3_selection.pdf", width=20, height=20, pointsize=30)
par(mar=c(1,0,0,6), xpd = T)
plot(exon_145_tree, label.offset = 1.5)
nodelabels(node = as.numeric(rownames(fitER$marginal.anc)), pie = fitER$marginal.anc, 
           piecol = c("blue", "red", "yellow","green"), cex = 0.6)
tiplabels(pie = to.matrix(assigned_z, sort(unique(assigned_z))), piecol = c("blue", "red", "yellow","green"), 
          cex = 0.3)
# Get xlim and y lim for legend justification
xlim = par("usr")[2]
ylim= par("usr")[4]
add.simmap.legend(colors=setNames(c("blue","red"),c("Hydrophobic","Polar")),
                  prompt=FALSE,shape="circle",x=xlim-5,y=ylim-1.0)
dev.off()

################ Exon 145, Codon 4
z <- as.factor(exon_145_data$codon4)
names(z) <- exon_145_tree$tip.label
assigned_z <- assignCharacteristics(z)

fitER <- rerootingMethod(exon_145_tree, assigned_z, model = "ER")

pdf(file="E145_C4_selection.pdf", width=20, height=20, pointsize=30)
par(mar=c(1,0,0,6), xpd = T)
plot(exon_145_tree, label.offset = 1.5)
nodelabels(node = as.numeric(rownames(fitER$marginal.anc)), pie = fitER$marginal.anc, 
           piecol = c("blue", "red", "yellow","green"), cex = 0.6)
tiplabels(pie = to.matrix(assigned_z, sort(unique(assigned_z))), piecol = c("blue", "red", "yellow","green"), 
          cex = 0.3)
# Get xlim and y lim for legend justification
xlim = par("usr")[2]
ylim= par("usr")[4]
add.simmap.legend(colors=setNames(c("blue","red"),c("Hydrophobic","Polar")),
                  prompt=FALSE,shape="circle",x=xlim-5,y=ylim-1.0)
dev.off()







###################################### 199 #######################################
# Read tree
exon_199_tree = read.tree("frame_f3_154307_154381.tre")

# Get codon data
exon_199_data = read.csv("frame_f3_154307_154381_E199_codons.csv", header = TRUE)
rownames(exon_199_data) = exon_199_data$species

# Compare tree and data
compare <- treedata(exon_199_tree, exon_199_data, sort=TRUE) 

# Revise data/tree
exon_199_data <- as.data.frame(compare$data)
#exon_137_data[,2:3] <- lapply(exon_137_data[,2], as.character)
exon_199_tree <- compare$phy

# save new tree and load it into r
write.tree(exon_199_tree, "exon_199.tre")
exon_199_tree <- read.tree("exon_199.tre")

############## Exon 199, Codon 1
z <- as.factor(exon_199_data$codon1)
names(z) <- exon_199_tree$tip.label
assigned_z <- assignCharacteristics(z)

fitER <- rerootingMethod(exon_199_tree, assigned_z, model = "ER")

pdf(file="E199_C1_selection.pdf", width=20, height=20, pointsize=30)
par(mar=c(1,0,0,6), xpd = T)
plot(exon_199_tree, label.offset = 1.5)
nodelabels(node = as.numeric(rownames(fitER$marginal.anc)), pie = fitER$marginal.anc, 
           piecol = c("blue", "red", "yellow","green"), cex = 0.6)
tiplabels(pie = to.matrix(assigned_z, sort(unique(assigned_z))), piecol = c("blue", "red", "yellow","green"), 
          cex = 0.3)
# Get xlim and y lim for legend justification
xlim = par("usr")[2]
ylim= par("usr")[4]
add.simmap.legend(colors=setNames(c("blue","red"),c("Hydrophobic","Polar")),
                  prompt=FALSE,shape="circle",x=xlim-5,y=ylim-0.35)
dev.off()

############## Exon 199, Codon 2
z <- as.factor(exon_199_data$codon2)
names(z) <- exon_199_tree$tip.label
assigned_z <- assignCharacteristics(z)

fitER <- rerootingMethod(exon_199_tree, assigned_z, model = "ER")

pdf(file="E199_C2_selection.pdf", width=20, height=20, pointsize=30)
par(mar=c(1,0,0,6), xpd = T)
plot(exon_199_tree, label.offset = 1.5)
nodelabels(node = as.numeric(rownames(fitER$marginal.anc)), pie = fitER$marginal.anc, 
           piecol = c("blue", "red", "yellow","green"), cex = 0.6)
tiplabels(pie = to.matrix(assigned_z, sort(unique(assigned_z))), piecol = c("blue", "red", "yellow","green"), 
          cex = 0.3)
# Get xlim and y lim for legend justification
xlim = par("usr")[2]
ylim= par("usr")[4]
add.simmap.legend(colors=setNames(c("blue","red"),c("Hydrophobic","Polar")),
                  prompt=FALSE,shape="circle",x=xlim-5,y=ylim-0.35)
dev.off()

############## Exon 199, Codon 3
z <- as.factor(exon_199_data$codon3)
names(z) <- exon_199_tree$tip.label
assigned_z <- assignCharacteristics(z)

fitER <- rerootingMethod(exon_199_tree, assigned_z, model = "ER")

pdf(file="E199_C3_selection.pdf", width=20, height=20, pointsize=30)
par(mar=c(1,0,0,6), xpd = T)
plot(exon_199_tree, label.offset = 1.5)
nodelabels(node = as.numeric(rownames(fitER$marginal.anc)), pie = fitER$marginal.anc, 
           piecol = c("blue", "red", "yellow","green"), cex = 0.6)
tiplabels(pie = to.matrix(assigned_z, sort(unique(assigned_z))), piecol = c("blue", "red", "yellow","green"), 
          cex = 0.3)
# Get xlim and y lim for legend justification
xlim = par("usr")[2]
ylim= par("usr")[4]
add.simmap.legend(colors=setNames(c("blue","red"),c("Hydrophobic","Polar")),
                  prompt=FALSE,shape="circle",x=xlim-5,y=ylim-0.35)

dev.off()

############## Exon 199, Codon 4 ----------------> all hydrophobic!
z <- as.factor(exon_199_data$codon4)
names(z) <- exon_199_tree$tip.label
assigned_z <- assignCharacteristics(z)

fitER <- rerootingMethod(exon_199_tree, assigned_z, model = "ER")
par(mar=c(1,0,0,6), xpd = T)
plot(exon_199_tree, label.offset = 1.5)
nodelabels(node = as.numeric(rownames(fitER$marginal.anc)), pie = fitER$marginal.anc, 
           piecol = c("blue", "red", "yellow","green"), cex = 0.6)
tiplabels(pie = to.matrix(assigned_z, sort(unique(assigned_z))), piecol = c("blue", "red", "yellow","green"), 
          cex = 0.3)
# Get xlim and y lim for legend justification
xlim = par("usr")[2]
ylim= par("usr")[4]
add.simmap.legend(colors=setNames(c("blue","red","yellow","green"),c("hydrophobic","polar","pos_charged","neg_charged")),
                  prompt=FALSE,shape="circle",x=xlim-5,y=ylim-0.35)












