phylum <- read.table("phylum.pcl", header=T, row.names=1, stringsAsFactors=F)
family <- read.table("family.pcl", header=T, row.names=1, stringsAsFactors=F)              
genus <- read.table("genus.pcl", header=T, row.names=1, stringsAsFactors=F)            
species <- read.table("species.pcl", header=T, row.names=1, stringsAsFactors=F) 
meta <- read.csv("MICROBIOME_SAMPLES_02JUL2018.csv", header = T, row.names = 1, stringsAsFactors=F)

d <- cbind(t(phylum), t(family), t(genus), t(species))
rownames(d) <- colnames(phylum)
d <- as.data.frame(d)
meta <- meta[as.character(rownames(d)),]

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
african_fr <- meta$AFR

cova <- cbind(race, bmi, age, sex, smoking, drink, energy, toothloss, hhincome, state, disease, batch, african_fr)
rownames(cova) <- rownames(meta)
cova <- as.data.frame(cova)

batch1_cova <- subset(cova, batch==1)
batch2_cova <- subset(cova, batch==0) 
batch1_ra <- d[as.character(rownames(batch1_cova)),]
batch2_ra <- d[as.character(rownames(batch2_cova)),]

batch1_black <- subset(batch1_cova, race==1)
batch1_white <- subset(batch1_cova, race==0)
batch1_ra_black <- batch1_ra[as.character(rownames(batch1_black)),]
batch1_ra_white <- batch1_ra[as.character(rownames(batch1_white)),]

batch2_black <- subset(batch2_cova, race==1)
batch2_white <- subset(batch2_cova, race==0)
batch2_ra_black <- batch2_ra[as.character(rownames(batch2_black)),]
batch2_ra_white <- batch2_ra[as.character(rownames(batch2_white)),]

for (i in 1:ncol(d)) {
  if (!(0 %in% as.numeric(d[,i]))) {next}
  batch1_pre_white <- subset(batch1_ra_white[,i], batch1_ra_white[,i] != 0) 
  batch1_prevalence_white <- length(batch1_pre_white)/length(batch1_ra_white[,i])
  batch1_prevalence_white <- toString(round(batch1_prevalence_white, 4))

  batch1_pre_black <- subset(batch1_ra_black[,i], batch1_ra_black[,i] != 0)
  batch1_prevalence_black <- length(batch1_pre_black)/length(batch1_ra_black[,i])
  batch1_prevalence_black <- toString(round(batch1_prevalence_black, 4))

  y <- (as.numeric(batch1_ra[,i]) !=0 ) ** 2
  batch1_model <- glm(y ~ batch1_cova$race + as.numeric(batch1_cova$bmi) + as.numeric(batch1_cova$age) + batch1_cova$sex + batch1_cova$smoking + batch1_cova$drink + factor(batch1_cova$energy) + factor(batch1_cova$toothloss) + factor(batch1_cova$hhincome) + factor(batch1_cova$state) + batch1_cova$disease, family = binomial("logit"))
  batch1_effect <- toString(suppressMessages(coef(summary(batch1_model)))[2,1])
  batch1_p <- toString(suppressMessages(coef(summary(batch1_model)))[2,4])

  batch2_pre_white <- subset(batch2_ra_white[,i], batch2_ra_white[,i] != 0)
  batch2_prevalence_white <- length(batch2_pre_white)/length(batch2_ra_white[,i])
  batch2_prevalence_white <- toString(round(batch2_prevalence_white, 4))

  batch2_pre_black <- subset(batch2_ra_black[,i], batch2_ra_black[,i] != 0)
  batch2_prevalence_black <- length(batch2_pre_black)/length(batch2_ra_black[,i])
  batch2_prevalence_black <- toString(round(batch2_prevalence_black, 4))

  y1 <- (as.numeric(batch2_ra[,i]) !=0 ) ** 2
  batch2_model <- glm(y1 ~ batch2_cova$race + as.numeric(batch2_cova$bmi) + as.numeric(batch2_cova$age) + batch2_cova$sex + batch2_cova$smoking + batch2_cova$drink + factor(batch2_cova$energy) + factor(batch2_cova$toothloss) + factor(batch2_cova$hhincome) + factor(batch2_cova$state) + batch2_cova$disease, family = binomial("logit"))
  batch2_effect <- toString(suppressMessages(coef(summary(batch2_model)))[2,1])
  batch2_p <- toString(suppressMessages(coef(summary(batch2_model)))[2,4])

  cat(colnames(d)[i], batch1_prevalence_white, batch1_prevalence_black, batch1_effect, batch1_p, batch2_prevalence_white, batch2_prevalence_black, batch2_effect, batch2_p, sep="\t", end="\n")
  
}

 
