library(compositions)

phylum <- read.table("phylum.pcl", header=T, row.names=1, stringsAsFactors=F)
family <- read.table("family.pcl", header=T, row.names=1, stringsAsFactors=F)              
genus <- read.table("genus.pcl", header=T, row.names=1, stringsAsFactors=F)            
species <- read.table("species.pcl", header=T, row.names=1, stringsAsFactors=F) 
meta <- read.csv("MICROBIOME_SAMPLES_02JUL2018.csv", header = T, row.names = 1, stringsAsFactors=F)

ra_phylum <- read.table("ra_phylum.pcl", header=T, row.names=1, stringsAsFactors=F)
ra_family <- read.table("ra_family.pcl", header=T, row.names=1, stringsAsFactors=F)
ra_genus <- read.table("ra_genus.pcl", header=T, row.names=1, stringsAsFactors=F)
ra_species <- read.table("ra_species.pcl", header=T, row.names=1, stringsAsFactors=F)

clr_phylum <- clr(t(phylum+1))
clr_family <- clr(t(family+1))
clr_genus <- clr(t(genus+1))
clr_species <- clr(t(species+1))   

d <- cbind(clr_phylum, clr_family, clr_genus, clr_species)
rownames(d) <- colnames(phylum)
d <- as.data.frame(d)
meta <- meta[as.character(rownames(d)),]
white <- subset(meta, Race==1)
black <- subset(meta, Race==0)

ra <- cbind(t(ra_phylum), t(ra_family), t(ra_genus), t(ra_species))
rownames(ra) <- colnames(phylum)
ra <- as.data.frame(ra)
ra_white <- ra[as.character(rownames(white)),]
ra_black <- ra[as.character(rownames(black)),]

hhincome <- meta$HHIncome
hhincome[hhincome==3] <- 2
hhincome[hhincome==4] <- 3
hhincome[hhincome==5] <- 3
hhincome[hhincome==8888] <- 4
hhincome[hhincome==9999] <- 4 
hhincome[is.na(hhincome)] <- 4 

drink <- meta$drink
drink[(drink != 1 | is.na(drink))] <- 2 
drink <- (drink == 1) ** 2

toothloss <- meta$ToothLossF1
toothloss[toothloss == 3] <- 2
toothloss[toothloss == 9999] <- 6
toothloss[is.na(toothloss)] <- 6


energy <- meta$FFQ_KCal_Tertile
energy[is.na(energy)] <- 4

state <- meta$Enrollment_State

batch <- c()
for (i in 1:nrow(meta)) {
  if (meta$CRC_SetID[i] == "") {
    batch <- c(batch, 1)
  }
  else {
    batch <- c(batch, 0)
  }
}

smoking <- (meta$SmokingStatus == 1) ** 2

race <- ((meta$Race == 0) ** 2)
age <- as.numeric(meta$Enrollment_Age)
sex <- meta$Gender
bmi <- as.numeric(meta$BMI)
disease <- meta$Disease_Status

cova <- cbind(race, bmi, age, sex, smoking, drink, energy, toothloss, hhincome, state, disease, batch)
rownames(cova) <- rownames(meta)
cova <- as.data.frame(cova)

batch1_cova <- subset(cova, batch==1)
batch2_cova <- subset(cova, batch==0) 
batch1_d <- d[as.character(rownames(batch1_cova)),]
batch2_d <- d[as.character(rownames(batch2_cova)),]
batch1_ra <- ra[as.character(rownames(batch1_cova)),]
batch2_ra <- ra[as.character(rownames(batch2_cova)),]

batch1_black <- subset(batch1_cova, race==1)
batch1_white <- subset(batch1_cova, race==0)
batch1_ra_black <- batch1_ra[as.character(rownames(batch1_black)),]
batch1_ra_white <- batch1_ra[as.character(rownames(batch1_white)),]

batch2_black <- subset(batch2_cova, race==1)
batch2_white <- subset(batch2_cova, race==0)
batch2_ra_black <- batch2_ra[as.character(rownames(batch2_black)),]
batch2_ra_white <- batch2_ra[as.character(rownames(batch2_white)),]


for (i in 1:ncol(d)) {

  batch1_median_white <- toString(round(median(batch1_ra_white[,i]), 4))
  batch1_median_black <- toString(round(median(batch1_ra_black[,i]), 4))
  batch1_model <- lm(as.numeric(batch1_d[,i]) ~ batch1_cova$race + as.numeric(batch1_cova$bmi) + as.numeric(batch1_cova$age) + batch1_cova$sex + batch1_cova$smoking + batch1_cova$drink + factor(batch1_cova$energy) + factor(batch1_cova$toothloss) + factor(batch1_cova$hhincome) + factor(batch1_cova$state) + batch1_cova$disease)
  batch1_effect <- toString(coef(summary(batch1_model))[2,1])
  batch1_p <- toString(coef(summary(batch1_model))[2,4])

  batch2_median_white <- toString(round(median(batch2_ra_white[,i]), 4))
  batch2_median_black <- toString(round(median(batch2_ra_black[,i]), 4))
  batch2_model <- lm(as.numeric(batch2_d[,i]) ~ batch2_cova$race + as.numeric(batch2_cova$bmi) + as.numeric(batch2_cova$age) + batch2_cova$sex + batch2_cova$smoking + batch2_cova$drink + factor(batch2_cova$energy) + factor(batch2_cova$toothloss) + factor(batch2_cova$hhincome) + factor(batch2_cova$state) + batch2_cova$disease)
  batch2_effect <- toString(coef(summary(batch2_model))[2,1])
  batch2_p <- toString(coef(summary(batch2_model))[2,4])

  cat(colnames(d)[i], batch1_median_white, batch1_median_black, batch1_effect, batch1_p, batch2_median_white, batch2_median_black, batch2_effect, batch2_p, sep="\t", end="\n")
}

 
