# This script creates dot plots of any full nucleotide PEVK-C sequence
# Kathleen Muenzen, 12/12/18
# kmuenzen@uw.edu

# Import dotplot package from evolvedmicrobe/dotplot
# Note: in order to install dotplot with install_github("evolvedmicrobe/dotplot", build_vignettes = FALSE),
# I needed to run R in the command line with "sudo R" and run the installation commands there

library(devtools)
library(dotplot)
library(Biostrings)
library(ggplot2)

# The default setttings in the package don't work for very dense plots
# Create modified functions from https://github.com/evolvedmicrobe/dotplot/tree/master/R
mkDotPlotDataFrame <- function(seq1, seq2, wsize, wstep, nmatch) {
  .Call('dotplot_mkDotPlotDataFrame', PACKAGE = 'dotplot', seq1, seq2, wsize, wstep, nmatch)
}

dotPlotg <- function (seq1, seq2, wsize = 10, wstep = 1, nmatch = -1)
{
  if (length(seq1[1]) > 1)
    stop("seq1 should be provided as a single string")
  if (length(seq2[1]) > 1)
    stop("seq2 should be provided as a single string")
  if (wsize < 1)
    stop("non allowed value for wsize")
  if (wstep < 1)
    stop("non allowed value for wstep")
  if (nmatch < 1)
    nmatch = wsize
  if (nmatch > wsize)
    stop("nmatch > wsize is not allowed")
  xy <- mkDotPlotDataFrame(seq1, seq2, wsize, wstep, nmatch)
  ggplot2::ggplot(xy, ggplot2::aes(x=x, y=y)) + ggplot2::geom_point(shape=15,size=0.01)
}

# Cheetah dot
cheetah_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Acinonyx_jubatus_ttn.fasta")))
seq1 = substr(cheetah_ttn,142239,162172)
seq2 = substr(cheetah_ttn,142239,162172)
cheetah_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/cheetah_dot.png")
cheetah_dot + labs(x="Cheetah PEVK-C", y = "Cheetah PEVK-C")
dev.off()


# Panda dot
panda_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Ailuropoda_melanoleuca_ttn.fasta")))
seq1 = substr(panda_ttn,136454,157867)
seq2 = substr(panda_ttn,136454,157867)
panda_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/panda_dot.png")
panda_dot + labs(x="Panda PEVK-C", y = "Panda PEVK-C")
dev.off()

# Cow dot
cow_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Bos_taurus_ttn.fasta")))
seq1 = substr(cow_ttn,144080,163728)
seq2 = substr(cow_ttn,144080,163728)
cow_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/cow_dot.png")
cow_dot + labs(x="Cow PEVK-C", y = "Cow PEVK-C")
dev.off()


# Beaver dot
beaver_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Castor_canadensis_ttn.fasta")))
seq1 = substr(beaver_ttn,137911,159291)
seq2 = substr(beaver_ttn,137911,159291)
beaver_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/beaver_dot.png")
beaver_dot + labs(x="Beaver PEVK-C", y = "Beaver PEVK-C")
dev.off()


# Rhino dot
rhino_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Ceratotherium_simum_ttn.fasta")))
seq1 = substr(rhino_ttn,139636,161248)
seq2 = substr(rhino_ttn,139636,161248)
rhino_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/rhino_dot.png")
rhino_dot + labs(x="Rhino PEVK-C", y = "Rhino PEVK-C")
dev.off()


# star_nosed mole dot
snmole_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Condylura_cristata_ttn.fasta")))
seq1 = substr(snmole_ttn,143353,161139)
seq2 = substr(snmole_ttn,143353,161139)
snmole_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/star_nosed_mole_dot.png")
snmole_dot + labs(x="Star-Nosed Mole PEVK-C", y = "Star-Nosed Mole PEVK-C")
dev.off()

# armadillo dot
armadillo_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Dasypus_novemcinctus_ttn.fasta")))
seq1 = substr(armadillo_ttn,139417,160228)
seq2 = substr(armadillo_ttn,139417,160228)
armadillo_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/armadillo_dot.png")
armadillo_dot + labs(x="Armadillo PEVK-C", y = "Armadillo PEVK-C")
dev.off()

# big brown bat dot
brwnbat_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Eptesicus_fuscus_ttn.fasta")))
seq1 = substr(brwnbat_ttn,146203,165152)
seq2 = substr(brwnbat_ttn,146203,165152)
brwnbat_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/big_brown_bat_dot.png")
brwnbat_dot + labs(x="Big Brown Bat PEVK-C", y = "Big Brown Bat PEVK-C")
dev.off()

# horse dot
horse_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Equus_caballus_ttn.fasta")))
seq1 = substr(horse_ttn,138535,159465)
seq2 = substr(horse_ttn,138535,159465)
horse_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/horse_dot.png")
horse_dot + labs(x="Horse PEVK-C", y = "Horse PEVK-C")
dev.off()

# hedgehog dot
hedgehog_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Erinaceus_europaeus_ttn.fasta")))
seq1 = substr(hedgehog_ttn,164558,202489)
seq2 = substr(hedgehog_ttn,164558,202489)
hedgehog_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/hedgehog_dot.png")
hedgehog_dot + labs(x="Hedgehog PEVK-C", y = "Hedgehog PEVK-C")
dev.off()


# cat dot
cat_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Felis_catus_ttn.fasta")))
seq1 = substr(cat_ttn,143776,163830)
seq2 = substr(cat_ttn,143776,163830)
cat_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/cat_dot.png")
cat_dot + labs(x="Cat PEVK-C", y = "Cat PEVK-C")
dev.off()


# lemur dot
lemur_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Galeopterus_variegatus_ttn.fasta")))
seq1 = substr(lemur_ttn,143143,168215)
seq2 = substr(lemur_ttn,143143,168215)
lemur_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/lemur_dot.png")
lemur_dot + labs(x="Lemur PEVK-C", y = "Lemur PEVK-C")
dev.off()

# roundleaf bat dot
roundleaf_bat_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Hipposideros_armiger_ttn.fasta")))
seq1 = substr(roundleaf_bat_ttn,143335,160971)
seq2 = substr(roundleaf_bat_ttn,143335,160971)
roundleaf_bat_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/roundleaf_bat_dot.png")
roundleaf_bat_dot + labs(x="Roundleaf Bat PEVK-C", y = "Roundleaf Bat PEVK-C")
dev.off()

# Human dot
human_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Homo_sapiens_ttn.fasta")))
seq1 = substr(human_ttn,141634,170054)
seq2 = substr(human_ttn,141634,170054)
human_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/human_dot.png")
human_dot + labs(x="Human PEVK-C", y = "Human PEVK-C")
dev.off()

# squirrel dot
squirrel_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Ictidomys_tridecemlineatus_ttn.fasta")))
seq1 = substr(squirrel_ttn,137660,157040)
seq2 = substr(squirrel_ttn,137660,157040)
squirrel_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/squirrel_dot.png")
squirrel_dot + labs(x="Squirrel PEVK-C", y = "Squirrel PEVK-C")
dev.off()

# elephant dot
elephant_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Loxodonta_africana_ttn.fasta")))
seq1 = substr(elephant_ttn,142725,165277)
seq2 = substr(elephant_ttn,142725,165277)
elephant_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/elephant_dot.png")
elephant_dot + labs(x="Elephant PEVK-C", y = "Elephent PEVK-C")
dev.off()

# pangolin dot
pangolin_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Manis_javanica_ttn.fasta")))
seq1 = substr(pangolin_ttn,136601,155080)
seq2 = substr(pangolin_ttn,136601,155080)
pangolin_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/pangolin_dot.png")
pangolin_dot + labs(x="Elephant PEVK-C", y = "Elephent PEVK-C")
dev.off()


# mouse dot
mouse_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Mus_musculus_ttn.fasta")))
seq1 = substr(mouse_ttn,148463,167028)
seq2 = substr(mouse_ttn,148463,167028)
mouse_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/mouse_dot.png")
mouse_dot + labs(x="Mouse PEVK-C", y = "Mouse PEVK-C")
dev.off()

# Shrewmouse dot
shrewmouse_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Mus_pahari_ttn.fasta")))
seq1 = substr(shrewmouse_ttn,143001,161955)
seq2 = substr(shrewmouse_ttn,143001,161955)
shrewmouse_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/shrewmouse_dot.png")
shrewmouse_dot + labs(x="Shrewmouse PEVK-C", y = "Shrewmouse PEVK-C")
dev.off()


# brandts bat dot
brandtsbat_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Myotis_brandtii_ttn.fasta")))
seq1 = substr(brandtsbat_ttn,151210,170828)
seq2 = substr(brandtsbat_ttn,151210,170828)
brandtsbat_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/brandts_bat_dot.png")
brandtsbat_dot + labs(x="Brandt's Bat PEVK-C", y = "Brandt's Bat PEVK-C")
dev.off()

# davids bat dot
davidsbat_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Myotis_davidii_ttn.fasta")))
seq1 = substr(davidsbat_ttn,150196,172225)
seq2 = substr(davidsbat_ttn,150196,172225)
davidsbat_ttn_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/davidsbat_dot.png")
davidsbat_ttn_dot + labs(x="David's Bat PEVK-C", y = "David's Bat PEVK-C")
dev.off()

# little brown bat dot
lilbrwnbat_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Myotis_lucifugus_ttn.fasta")))
seq1 = substr(lilbrwnbat_ttn,142041,166467)
seq2 = substr(lilbrwnbat_ttn,142041,166467)
lilbrwnbat_ttn_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/lil_brown_bat_dot.png")
lilbrwnbat_ttn_dot + labs(x="Little Brown Bat PEVK-C", y = "Little Brown Bat PEVK-C")
dev.off()


# monk seal dot
monkseal_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Neomonachus_schauinslandi_ttn.fasta")))
seq1 = substr(monkseal_ttn,138731,159340)
seq2 = substr(monkseal_ttn,138731,159340)
monkseal_ttn_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/monk_seal_dot.png")
monkseal_ttn_dot + labs(x="Monk Seal PEVK-C", y = "Monk Seal PEVK-C")
dev.off()


# pika dot
pika_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Ochotona_princeps_ttn.fasta")))
seq1 = substr(pika_ttn,133647,150992)
seq2 = substr(pika_ttn,133647,150992)
pika_ttn_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/pika_dot.png")
pika_ttn_dot + labs(x="Pika PEVK-C", y = "Pika PEVK-C")
dev.off()

# walrus dot
walrus_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Odobenus_rosmarus_ttn.fasta")))
seq1 = substr(walrus_ttn,137525,166049)
seq2 = substr(walrus_ttn,137525,166049)
walrus_ttn_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/walrus_dot.png")
walrus_ttn_dot + labs(x="Walrus PEVK-C", y = "Walrus PEVK-C")
dev.off()

# Orca dot
orca_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Orcinus_orca_ttn.fasta")))
seq1 = substr(orca_ttn,148206,168011)
seq2 = substr(orca_ttn,148206,168011)
orca_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/orca_dot.png")
orca_dot + labs(x="Orca PEVK-C", y = "Orca PEVK-C")
dev.off()

# Platypus dot
platypus_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Ornithorhynchus_anatinus_ttn.fasta")))
seq1 = substr(platypus_ttn,160208,185739)
seq2 = substr(platypus_ttn,160208,185739)
platypus_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/platypus_dot.png")
platypus_dot + labs(x="Platypus PEVK-C", y = "Platypus PEVK-C")
dev.off()

# Rabbit dot
rabbit_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Oryctolagus_cuniculus_ttn.fasta")))
seq1 = substr(rabbit_ttn,139116,160407)
seq2 = substr(rabbit_ttn,139116,160407)
rabbit_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/rabbit_dot.png")
rabbit_dot + labs(x="Rabbit PEVK-C", y = "Rabbit PEVK-C")
dev.off()

# Chimp dot
chimp_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Pan_troglodytes_ttn.fasta")))
seq1 = substr(chimp_ttn,137457,175029)
seq2 = substr(chimp_ttn,137457,175029)
chimp_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/chimp_dot.png")
chimp_dot + labs(x="Chimpanzee PEVK-C", y = "Chimpanzee PEVK-C")
dev.off()


# Koala dot
koala_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Phascolarctos_cinereus_ttn.fasta")))
seq1 = substr(koala_ttn,176127,205810)
seq2 = substr(koala_ttn,176127,205810)
koala_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/koala_dot.png")
koala_dot + labs(x="Koala PEVK-C", y = "Koala PEVK-C")
dev.off()

# Orangutan dot
orangutan_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Pongo_abelii_ttn.fasta")))
seq1 = substr(orangutan_ttn,141562,169209)
seq2 = substr(orangutan_ttn,141562,169209)
orangutan_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/orangutan_dot.png")
orangutan_dot + labs(x="Orangutan PEVK-C", y = "Orangutan PEVK-C")
dev.off()

# Flying fox dot
flyingfox_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Pteropus_vampyrus_ttn.fasta")))
seq1 = substr(flyingfox_ttn,139070,158160)
seq2 = substr(flyingfox_ttn,139070,158160)
flyingfox_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/flying_fox_dot.png")
flyingfox_dot + labs(x="Flying Fox PEVK-C", y = "Flying Fox PEVK-C")
dev.off()

# Rat dot
rat_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Rattus_norvegicus_ttn.fasta")))
seq1 = substr(rat_ttn,144135,162127)
seq2 = substr(rat_ttn,144135,162127)
rat_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/rat_dot.png")
rat_dot + labs(x="Rat PEVK-C", y = "Rat PEVK-C")
dev.off()

# Horseshoe bat dot
horseshoebat_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Rhinolophus_sinicus_ttn.fasta")))
seq1 = substr(horseshoebat_ttn,141681,160526)
seq2 = substr(horseshoebat_ttn,141681,160526)
horseshoebat_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/horseshoe_bat_dot.png")
horseshoebat_dot + labs(x="Horseshoe Bat PEVK-C", y = "Horseshoe Bat PEVK-C")
dev.off()

# Tasmanian devil dot
tasmaniandevil_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Sarcophilus_harrisii_ttn.fasta")))
seq1 = substr(tasmaniandevil_ttn,161169,194288)
seq2 = substr(tasmaniandevil_ttn,161169,194288)
tasmaniandevil_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/tasmanian_devil_dot.png")
tasmaniandevil_dot + labs(x="Tasmanian Devil PEVK-C", y = "Tasmanian Devil PEVK-C")
dev.off()

# Pig dot
pig_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Sus_scrofa_ttn.fasta")))
seq1 = substr(pig_ttn,144210,163437)
seq2 = substr(pig_ttn,144210,163437)
pig_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/pig_dot.png")
pig_dot + labs(x="Pig PEVK-C", y = "Pig PEVK-C")
dev.off()

# Manatee dot
manatee_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Trichechus_manatus_ttn.fasta")))
seq1 = substr(manatee_ttn,141453,158040)
seq2 = substr(manatee_ttn,141453,158040)
manatee_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/manatee_dot.png")
manatee_dot + labs(x="Manatee PEVK-C", y = "Manatee PEVK-C")
dev.off()


# Chinese shrew dot
chineseshrew_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Tupaia_chinensis_ttn.fasta")))
seq1 = substr(chineseshrew_ttn,159164,180463)
seq2 = substr(chineseshrew_ttn,159164,180463)
chineseshrew_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/chinese_shrew_dot.png")
chineseshrew_dot + labs(x="Chinese Shrew PEVK-C", y = "Chinese Shrew PEVK-C")
dev.off()

# Dolphin dot
dolphin_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Tursiops_truncatus_ttn.fasta")))
seq1 = substr(dolphin_ttn,145527,163558)
seq2 = substr(dolphin_ttn,145527,163558)
dolphin_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/dolphin_dot.png")
dolphin_dot + labs(x="Dolphin PEVK-C", y = "Dolphin PEVK-C")
dev.off()

# Polar bear dot
polarbear_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Ursus_maritimus_ttn.fasta")))
seq1 = substr(polarbear_ttn,136078,156826)
seq2 = substr(polarbear_ttn,136078,156826)
polarbear_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/polar_bear_dot.png")
polarbear_dot + labs(x="Polar Bear PEVK-C", y = "Polar Bear PEVK-C")
dev.off()

# Alpaca dot
alpaca_ttn = toString(as.character(readDNAStringSet("~/Desktop/paper_submission_materials/data/ttn_seqs/Vicugna_pacos_ttn.fasta")))
seq1 = substr(alpaca_ttn,138814,161356)
seq2 = substr(alpaca_ttn,138814,161356)
alpaca_dot <- dotPlotg(seq1, seq2, wsize=7)
png(filename = "~/Desktop/titin_project/pevkc_dotplots/alpaca_dot.png")
alpaca_dot + labs(x="Alpaca PEVK-C", y = "Alpaca PEVK-C")
dev.off()


####### LINE analysis #########


