library(MiRKAT)
library(GUniFrac)

otu <- read.table("merge_from_raw.lcp", header = T, row.names = 1, check.names=F)
meta <- read.csv("MICROBIOME_SAMPLES_02JUL2018.csv", header = T, row.names = 1, stringsAsFactors=F)
tree <- read.tree("rep_set.tre") 

meta <- meta[as.character(rownames(otu)),]

white <- subset(meta, Race==1)
black <- subset(meta, Race==0)

hhincome <- meta$HHIncome
hhincome[hhincome==3] <- 2
hhincome[hhincome==4] <- 3
hhincome[hhincome==5] <- 3
hhincome[hhincome==8888] <- 4
hhincome[hhincome==9999] <- 4 
hhincome[is.na(hhincome)] <- 4  
hh1 <- (hhincome == 1) ** 2
hh2 <- (hhincome == 2) ** 2
hh3 <- (hhincome == 3) ** 2


drink <- meta$drink
drink[(drink != 1 | is.na(drink))] <- 2 
drink <- (drink == 1) ** 2

toothloss <- meta$ToothLossF1
toothloss[toothloss == 3] <- 2
toothloss[toothloss == 9999] <- 6
toothloss[is.na(toothloss)] <- 6
toothloss1 <- (toothloss == 1) ** 2
toothloss2 <- (toothloss == 2) ** 2
toothloss3 <- (toothloss == 4) ** 2
toothloss4 <- (toothloss == 5) ** 2

energy <- meta$FFQ_KCal_Tertile
energy[is.na(energy)] <- 4
energy1 <- (energy == 1) ** 2
energy2 <- (energy == 2) ** 2
energy3 <- (energy == 3) ** 2

state <- meta$Enrollment_State
al <- (state=="AL") ** 2
ar <- (state=="AR") ** 2
fl <- (state=="FL") ** 2
ga <- (state=="GA") ** 2
ky <- (state=="KY") ** 2
la <- (state=="LA") ** 2
ms <- (state=="MS") ** 2
nc <- (state=="NC") ** 2
sc <- (state=="SC") ** 2
tn <- (state=="TN") ** 2
va <- (state=="VA") ** 2


batch <- c()
for (i in 1:nrow(meta)) {
  if (meta$CRC_SetID[i] == "") {
    batch <- c(batch, 1)
  }
  else {
    batch <- c(batch, 0)
  }
}

race <- ((meta$Race == 0) ** 2)
age <- as.numeric(meta$Enrollment_Age)
sex <- meta$Gender
bmi <- as.numeric(meta$BMI)
disease <- meta$Disease_Status
smoking <- ((meta$SmokingStatus == 1) ** 2)

cova <- cbind(bmi, age, sex, smoking, drink, energy1, energy2, energy3, toothloss1, toothloss2, toothloss3, toothloss4, hh1, hh2, hh3, al, ar, fl, ga, ky, la, ms, nc, sc, tn, va, disease, batch)

otu_rff <- Rarefy(otu)$otu.tab.rff

unifracs <- GUniFrac(otu_rff, tree, alpha=c(0, 0.5, 1))$unifracs
D.weighted = unifracs[,,"d_1"]
D.unweighted = unifracs[,,"d_UW"]
D.BC= as.matrix(vegdist(otu_rff, method="bray"))

K.weighted = D2K(D.weighted)
K.unweighted = D2K(D.unweighted)
K.BC = D2K(D.BC)

set.seed(123)
MiRKAT(y = race, Ks = K.weighted, X = cova, out_type = "D", nperm=1000, method = "davies")
MiRKAT(y = race, Ks = K.unweighted, X = cova, out_type = "D", nperm=1000, method = "davies")
MiRKAT(y = race, Ks = K.BC, X = cova, out_type = "D", nperm=1000, method = "davies")
