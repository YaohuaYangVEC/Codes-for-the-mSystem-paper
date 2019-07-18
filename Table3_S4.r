phylum <- read.table("phylum.pcl", header=T, row.names=1, stringsAsFactors=F)
family <- read.table("family.pcl", header=T, row.names=1, stringsAsFactors=F)              
genus <- read.table("genus.pcl", header=T, row.names=1, stringsAsFactors=F)            
species <- read.table("species.pcl", header=T, row.names=1, stringsAsFactors=F) 
meta <- read.csv("MICROBIOME_SAMPLES_02JUL2018.csv", header = T, row.names = 1, stringsAsFactors=F)

d <- cbind(t(phylum), t(family), t(genus), t(species))
rownames(d) <- colnames(phylum)
d <- as.data.frame(d)
meta <- meta[as.character(rownames(d)),]
white <- subset(meta, Race==1)
black <- subset(meta, Race==0)
ra_white <- d[as.character(rownames(white)),]
ra_black <- d[as.character(rownames(black)),]

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

afr_sample <- c(rownames(subset(meta, !is.na(AFR))))
for (i in 1:nrow(white)) {
  if (!(rownames(white)[i] %in% afr_sample)) {afr_sample <- c(afr_sample, rownames(white)[i])}
}
afr_d <- d[as.character(afr_sample),]
afr_cova <- cova[as.character(afr_sample),]
afr <- as.character(afr_cova$african_fr)
afr[is.na(afr)] <- 0.0
afr <- as.numeric(afr)

for (i in 1:ncol(d)) {
  if (!(0 %in% as.numeric(d[,i]))) {next}
  pre_white <- subset(ra_white[,i], ra_white[,i] != 0) 
  prevalence_white <- length(pre_white)/length(ra_white[,i])
  prevalence_white <- toString(round(prevalence_white, 4))

  pre_black <- subset(ra_black[,i], ra_black[,i] != 0)
  prevalence_black <- length(pre_black)/length(ra_black[,i])
  prevalence_black <- toString(round(prevalence_black, 4))

  y <- (as.numeric(d[,i]) !=0 ) ** 2
  model_race <- glm(y ~ race + as.numeric(bmi) + as.numeric(age) + sex + smoking + drink + factor(energy) + factor(toothloss) + factor(hhincome) + factor(state) + disease + batch, family = binomial("logit"))
  effect_race <- toString(suppressMessages(coef(summary(model_race)))[2,1])
  effect_bmi <- toString(suppressMessages(coef(summary(model_race)))[3,1])
  effect_age <- toString(suppressMessages(coef(summary(model_race)))[4,1])
  effect_sex <- toString(suppressMessages(coef(summary(model_race)))[5,1])
  effect_smoking <- toString(suppressMessages(coef(summary(model_race)))[6,1])
  effect_drink <- toString(suppressMessages(coef(summary(model_race)))[7,1])
  effect_energy2 <- toString(suppressMessages(coef(summary(model_race)))[8,1])
  effect_energy3 <- toString(suppressMessages(coef(summary(model_race)))[9,1])
  effect_energy4 <- toString(suppressMessages(coef(summary(model_race)))[10,1])
  effect_toothloss2 <- toString(suppressMessages(coef(summary(model_race)))[11,1])
  effect_toothloss3 <- toString(suppressMessages(coef(summary(model_race)))[12,1])
  effect_toothloss4 <- toString(suppressMessages(coef(summary(model_race)))[13,1])
  effect_toothloss5 <- toString(suppressMessages(coef(summary(model_race)))[14,1])
  effect_hh1 <- toString(suppressMessages(coef(summary(model_race)))[15,1])
  effect_hh2 <- toString(suppressMessages(coef(summary(model_race)))[16,1])
  effect_hh3 <- toString(suppressMessages(coef(summary(model_race)))[17,1])
  effect_ar <- toString(suppressMessages(coef(summary(model_race)))[18,1])
  effect_fl <- toString(suppressMessages(coef(summary(model_race)))[19,1])
  effect_ga <- toString(suppressMessages(coef(summary(model_race)))[20,1])
  effect_ky <- toString(suppressMessages(coef(summary(model_race)))[21,1])
  effect_la <- toString(suppressMessages(coef(summary(model_race)))[22,1])
  effect_ms <- toString(suppressMessages(coef(summary(model_race)))[23,1])
  effect_nc <- toString(suppressMessages(coef(summary(model_race)))[24,1])
  effect_sc <- toString(suppressMessages(coef(summary(model_race)))[25,1])
  effect_tn <- toString(suppressMessages(coef(summary(model_race)))[26,1])
  effect_va <- toString(suppressMessages(coef(summary(model_race)))[27,1])
  effect_wv <- toString(suppressMessages(coef(summary(model_race)))[28,1])
  effect_disease <- toString(coef(summary(model_race))[29,1])
  effect_batch <- toString(coef(summary(model_race))[30,1])

  p_race <- toString(suppressMessages(coef(summary(model_race)))[2,4])
  bonfer_p_race <- p_race
  if (p_race == "NaN") {
    bonfer_p_race <- "NaN"
  }
  else {
    bonfer_p_race <- coef(summary(model_race))[2,4] * 71
    if (bonfer_p_race >= 1.0) {
        bonfer_p_race <- "1.0"
    }
    else {bonfer_p_race <- toString(bonfer_p_race)}
  }
  y1 <- (as.numeric(afr_d[,i]) !=0 ) ** 2
  model_afr <- glm(y1 ~ afr + as.numeric(afr_cova$bmi) + as.numeric(afr_cova$age) + afr_cova$sex + afr_cova$smoking + afr_cova$drink + factor(afr_cova$energy) + factor(afr_cova$toothloss) + factor(afr_cova$hhincome) + factor(afr_cova$state) + afr_cova$disease + afr_cova$batch, family = binomial("logit"))

  effect_afr <- toString(suppressMessages(coef(summary(model_afr)))[2,1])
  p_afr <- toString(suppressMessages(coef(summary(model_afr)))[2,4])

  cat(colnames(d)[i], prevalence_white, prevalence_black, effect_race, p_race, bonfer_p_race, effect_afr, p_afr, effect_bmi, effect_age, effect_sex, effect_smoking, effect_drink, effect_energy2, effect_energy3, effect_energy4, effect_toothloss2, effect_toothloss3, effect_toothloss4, effect_toothloss5, effect_hh1, effect_hh2, effect_hh3, effect_ar, effect_fl, effect_ga, effect_ky, effect_la, effect_ms, effect_nc, effect_sc, effect_tn, effect_va, effect_wv, effect_disease, effect_batch, sep="\t", end="\n")
  
}

 
