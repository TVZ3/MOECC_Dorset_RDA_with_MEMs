library(ade4)
library(vegan)
library(packfor)
library(MASS)
library(gclus)
library(ape)
library(ellipse)
library(lmom)
library(dplyr)
library(tidyr)

setwd("C:\\Users\\adminuser\\Desktop\\EVERYTHING THOMAS\\gitprojects\\MOECC_Dorset_RDA_with_MEMs")

dorset <- read.csv("03Dorset.csv")

names(dorset)
dorset <- dorset[, -c(16:19, 32:35, 44:55, 84:87, 100,101), drop=FALSE]
dorset <- na.omit(dorset)

morphology <- cbind(dorset[, c(11:13)]) # Lake morphology variables
chemistry <- cbind(dorset[, c(16:63)])  # Lake chemistry variables
climate <- cbind(dorset[, c(83:121)])   # Climate data from Climatic Research Unit (CRU)
algae <- cbind(dorset[, c(64:75)])      # Algae species biomass and relative percentage data
mem <- cbind(dorset[, c(122:128)])      # Moran Eigenvector cycle variables
year <- cbind(dorset[, c(15)])          # Column for year
ID <- cbind(dorset[, c(1:5)])           # Lake ID by time, and names

IF_chem <- select(chemistry, ends_with("IF"))
IF_clim <- cbind(climate[, c(13, 26, 39)])

rel_al <- select(algae, ends_with("REL"))
rel_spe_hel = decostand(rel_al, "hellinger")      # Algal biomass data Hellinger-transformed
morpho_std = decostand(morphology, "standardize") # Morphology standardized
IF_chem_std = decostand(IF_chem, "standardize")   # Ice free chemistry standardized
IF_clim_std = decostand(IF_clim, "standardize")   # Ice free climate standardized

icefree <- data.frame(morpho_std, IF_chem_std, IF_clim_std, mem)

forward.sel(rel_spe_hel, icefree)

icefree <- data.frame(ID, mem, morpho_std, IF_chem_std, IF_clim_std)

par(mfrow=c(1,1))
mod <- rda(rel_spe_hel ~ zmax+pH_IF+Cl_IF+A0+NH4_IF+TN_IF+NO3_IF+Cond_IF+IF_Cloud+IF_Precip+IF_MeanTemp+PCNM2+PCNM4+PCNM9+PCNM3+PCNM8+PCNM17, data=icefree, scale=FALSE)
plot(mod, type = "n", scaling = 2)
colvec <- c("blue","dodgerblue3","cyan","darkseagreen1","yellow2","goldenrod1","orange","tomato","red2")
with(icefree, levels(Bin_Time))
text(mod, display = "species", cex = 0.8, col = "black")
with(icefree, points(mod, display = "sites", col = adjustcolor(colvec[Bin_Time], alpha=0.8),
                     scaling = 2, pch = 18, cex = 0.7, bg = colvec[Bin_Time]))
text(mod, scaling = 2, display = "bp", col = "blue", cex = 0.8) 
with(icefree, legend("topright", legend = levels(Bin_Time), bty = "n",
                     col = colvec, pch = 18, cex = 1.0, pt.bg = colvec))

(R2 <- RsquareAdj(mod)$r.squared)
(R2adj <- RsquareAdj(mod)$adj.r.squared)
vif.cca(mod)

colnames(icefree)
v.morph <- cbind(icefree[, c(13, 15)])             #X1 - morphology variables
v.chem <- cbind(icefree[, c(17, 18, 21:23, 27)])   #X2 - lake chemistry variables
v.clim <- cbind(icefree[, c(28:30)])               #X3 - climate variables
v.mem <- cbind(icefree[, c(6:8, 10:12)])           #X4 - Eigenvector variables

(x <- varpart(rel_spe_hel,v.morph,v.chem,v.clim,v.mem))
showvarparts(4)
