##########################################################
# Script to analyse the relationships between the community-weighted means
# of 20 plant traits and climate regulation processes
# script by: Dr. Stephan Kambach (stephan.kambach@gmail.com)
# initial date: 05.06.2023
########################################################

# 1. set session -------------------------------------------------------

# clean workspace
rm(list=ls())
gc()
options(na.action = "na.fail")

# load packages
library(tidyverse)
library(tidyselect)
library(vroom)
library(lubridate)
library(FactoMineR)
library(ggforce)
library(corrplot)
library(dismo)
library(viridis)
library(grid)
library(gridExtra)
library(randomForest)
library(scales)
library(randomForestExplainer)
library(ggeffects)
library(ggnewscale)
library(rstatix)
library(ggpubr)
library(ggtext)
library(vegan)
library(ggthemes)
library(edarf)
library(multcomp)
library(multcompView)
library(psych)
library(GGally)
library(ggrepel)

# set working directory that contains the data 
# and the additional functions
setwd("C:/Users/Stephan/OneDrive/current tasks/postdoc feedbacks/scripts/feedbacks_project2/manuscript/11. submission to Nat Clim Change") 

# 2. load and prepare data ------------------------------------------------------------
dat.raw = read_delim("eva_all_merged.csv", delim = "\t")
source("kambach - additional functions for analysis.R")

# log-transform variables
dat.log = dat.raw %>%
  mutate(sp_richness = log(sp.richness),
         et = log(et),
         chelsa_gsp = log(chelsa_gsp),
         
         CWM_LDMC = log(CWM_LDMC),
         CWM_Leaf_delta_15N = CWM_Leaf_delta_15N,
         CWM_LeafArea_leaf_undef = log(CWM_LeafArea_leaf_undef),
         CWM_LeafC_perdrymass = log(CWM_LeafC_perdrymass),
         CWM_LeafCN_ratio = log(CWM_LeafCN_ratio),
         CWM_LeafDryMass_single = log(CWM_LeafDryMass_single),
         CWM_Leaffreshmass = log(CWM_Leaffreshmass),
         CWM_LeafN = log(CWM_LeafN),
         CWM_LeafP = log(CWM_LeafP),
         CWM_LeafThickness = log(CWM_LeafThickness),
         CWM_LeafWaterCont = log(CWM_LeafWaterCont),
         CWM_PlantHeight = log(CWM_PlantHeight),
         CWM_RootingDepth = log(CWM_RootingDepth),
         CWM_SeedMass = log(CWM_SeedMass),
         CWM_SLA = log(CWM_SLA),
         CWM_SpecificRootLength = log(CWM_SpecificRootLength),
         CWM_Stem.cond.dens = log(CWM_Stem.cond.dens),
         CWM_StemConduitDiameter = log(CWM_StemConduitDiameter),
         CWM_StemDens = log(CWM_StemDens),
         CWM_StemDiam = log(CWM_StemDiam),
         
         
         eunis_class1 = factor(eunis_class1),
         eunis_class2 = factor(eunis_class2),
         eunis_class3 = factor(eunis_class3)) %>%
  drop_na(eunis_class2) 

# scale variables within eunis level 1 classes
dat.scaled.eunis2 = dat.log %>%
  group_by(eunis_class2) %>%
  mutate(sp_richness = scale_this(sp.richness),
         cover_total_cwm = scale_this(cover_total_cwm),
         et = scale_this(et),
         npp = scale_this(npp),
         CWM_LDMC = scale_this(CWM_LDMC),
         CWM_Leaf_delta_15N = scale_this(CWM_Leaf_delta_15N),
         CWM_LeafArea_leaf_undef = scale_this(CWM_LeafArea_leaf_undef),
         CWM_LeafC_perdrymass = scale_this(CWM_LeafC_perdrymass),
         CWM_LeafCN_ratio = scale_this(CWM_LeafCN_ratio),
         CWM_LeafDryMass_single = scale_this(CWM_LeafDryMass_single),
         CWM_Leaffreshmass = scale_this(CWM_Leaffreshmass),
         CWM_LeafN = scale_this(CWM_LeafN),
         CWM_LeafP = scale_this(CWM_LeafP),
         CWM_LeafThickness = scale_this(CWM_LeafThickness),
         CWM_LeafWaterCont = scale_this(CWM_LeafWaterCont),
         CWM_PlantHeight = scale_this(CWM_PlantHeight),
         CWM_RootingDepth = scale_this(CWM_RootingDepth),
         CWM_SeedMass = scale_this(CWM_SeedMass),
         CWM_SLA = scale_this(CWM_SLA),
         CWM_SpecificRootLength = scale_this(CWM_SpecificRootLength),
         CWM_Stem.cond.dens = scale_this(CWM_Stem.cond.dens),
         CWM_StemConduitDiameter = scale_this(CWM_StemConduitDiameter),
         CWM_StemDens = scale_this(CWM_StemDens),
         CWM_StemDiam = scale_this(CWM_StemDiam),
         
         chelsa_cmi = scale_this(chelsa_cmi), chelsa_gsl = scale_this(chelsa_gsl),
         chelsa_gsp = scale_this(chelsa_gsp), chelsa_gst = scale_this(chelsa_gst),
         chelsa_npp = scale_this(chelsa_npp), chelsa_pet = scale_this(chelsa_pet),
         chelsa_rsds = scale_this(chelsa_rsds), chelsa_sfcWind = scale_this(chelsa_sfcWind),
         chelsa_swb = scale_this(chelsa_swb),
         
         bio_01 = scale_this(bio_01), bio_02 = scale_this(bio_02),
         bio_03 = scale_this(bio_03), bio_04 = scale_this(bio_04),
         bio_05 = scale_this(bio_05), bio_06 = scale_this(bio_06),
         bio_07 = scale_this(bio_07), bio_08 = scale_this(bio_08),
         bio_09 = scale_this(bio_09), bio_10 = scale_this(bio_10),
         bio_11 = scale_this(bio_11), bio_12 = scale_this(bio_12),
         bio_13 = scale_this(bio_13), bio_14 = scale_this(bio_14),
         bio_15 = scale_this(bio_15), bio_16 = scale_this(bio_16),
         bio_17 = scale_this(bio_17), bio_18 = scale_this(bio_18),
         bio_19 = scale_this(bio_19)) %>%
  ungroup() %>%
  drop_na(eunis_class2)

# set the colouring scheme
color.gradient = c("#a6c6e1", "#2171b5", "#13436c",
                   "#fed98e", "#fe9929", "#993404",
                   "#acdaba", "#45ac65", "#22723a", "#134121")

# create empty list that will be used to save analysis results
empty.list.for.results = list("Forest_coniferous" = list("prop_sis_reflected" = NA, "et" = NA, "npp" = NA),
                              "Forest_deciduous" = list("prop_sis_reflected" = NA, "et" = NA, "npp" = NA),
                              "Forest_evergreen" = list("prop_sis_reflected" = NA, "et" = NA, "npp" = NA),
                              "Shrubland_alpine" = list("prop_sis_reflected" = NA, "et" = NA, "npp" = NA),
                              "Shrubland_heathland" = list("prop_sis_reflected" = NA, "et" = NA, "npp" = NA),
                              "Shrubland_temperate" = list("prop_sis_reflected" = NA, "et" = NA, "npp" = NA),
                              "Grassland_alpine" = list("prop_sis_reflected" = NA, "et" = NA, "npp" = NA),
                              "Grassland_dry" = list("prop_sis_reflected" = NA, "et" = NA, "npp" = NA),
                              "Grassland_mesic" = list("prop_sis_reflected" = NA, "et" = NA, "npp" = NA),
                              "Grassland_wet" = list("prop_sis_reflected" = NA, "et" = NA, "npp" = NA))

# create vector with the trait axes used for analysis
traits.selected = c("CWM_pca1", "CWM_pca2", "CWM_pca3", "CWM_pca4", 
                    "CWM_pca5", "CWM_pca6", "CWM_pca7", "CWM_pca8")

# create function to shorten trait names
recode.trait.names = function(trait.vector.temp){
  trait.vector.temp[trait.vector.temp == "CWM_LDMC"] = "LDMC"
  trait.vector.temp[trait.vector.temp == "CWM_Leaf_delta_15N"] = "L15N"
  trait.vector.temp[trait.vector.temp == "CWM_LeafArea_leaf_undef"] = "LArea"
  trait.vector.temp[trait.vector.temp == "CWM_LeafC_perdrymass"] = "LCcont"
  trait.vector.temp[trait.vector.temp == "CWM_LeafCN_ratio"] = "LCNratio"
  trait.vector.temp[trait.vector.temp == "CWM_LeafDryMass_single"] = "LMdry"
  trait.vector.temp[trait.vector.temp == "CWM_Leaffreshmass"] = "LMfresh"
  trait.vector.temp[trait.vector.temp == "CWM_LeafN"] = "LNcont"
  trait.vector.temp[trait.vector.temp == "CWM_LeafP"] = "LPcont"
  trait.vector.temp[trait.vector.temp == "CWM_LeafThickness"] = "LThic"
  trait.vector.temp[trait.vector.temp == "CWM_LeafWaterCont"] = "LWcont"
  trait.vector.temp[trait.vector.temp == "CWM_PlantHeight"] = "PlantHeight"
  trait.vector.temp[trait.vector.temp == "CWM_RootingDepth"] = "RDepth"
  trait.vector.temp[trait.vector.temp == "CWM_SeedMass"] = "SM"
  trait.vector.temp[trait.vector.temp == "CWM_SLA"] = "SLA"
  trait.vector.temp[trait.vector.temp == "CWM_SpecificRootLength"] = "SRL"
  trait.vector.temp[trait.vector.temp == "CWM_Stem.cond.dens"] = "SCdens"
  trait.vector.temp[trait.vector.temp == "CWM_StemConduitDiameter"] = "SCdiam"
  trait.vector.temp[trait.vector.temp == "CWM_StemDens"] = "Sdens"
  trait.vector.temp[trait.vector.temp == "CWM_StemDiam"] = "Sdiam"
  
  return(trait.vector.temp)
}

# 3. determine the principal components of trait community-weighted means  -----------------------------------------------
# show correlation
trait.cor = cor(dat.scaled.eunis2 %>%
                  dplyr::select(CWM_LDMC:CWM_StemDiam))
ggcorrplot::ggcorrplot(trait.cor, method = "square") + 
  scale_fill_stepsn(colours = c("blue", "white", "white", "red"))

# create empty list to store the results of separate PCAs
all.trait.pcas = list()

# add empty colums to the data to store the position of 
# each plot along the first 8 principal component axes
dat.scaled.eunis2$CWM_pca1 = as.numeric(NA)
dat.scaled.eunis2$CWM_pca2 = as.numeric(NA)
dat.scaled.eunis2$CWM_pca3 = as.numeric(NA)
dat.scaled.eunis2$CWM_pca4 = as.numeric(NA)
dat.scaled.eunis2$CWM_pca5 = as.numeric(NA)
dat.scaled.eunis2$CWM_pca6 = as.numeric(NA)
dat.scaled.eunis2$CWM_pca7 = as.numeric(NA)
dat.scaled.eunis2$CWM_pca8 = as.numeric(NA)

# calculate varimax rotated trait pcas with 8 axes
# for all forest habitats
all.trait.pcas[["Forest"]] = make.a.nice.pca.plot.test(
  data.temp = dat.scaled.eunis2 %>% filter(eunis_class1 == "Forest") %>%
    dplyr::select(CWM_LDMC:CWM_StemDiam),
  nfactors.temp = 8, rotate.temp = "varimax", print.plot = T,
  recode.trait.names = recode.trait.names, switch.axes.temp = "Dim.5")

# for all shrubland habitats
all.trait.pcas[["Shrubland"]] = make.a.nice.pca.plot.test(
  data.temp = dat.scaled.eunis2 %>% filter(eunis_class1 == "Shrubland") %>%
    dplyr::select(CWM_LDMC:CWM_StemDiam),
  nfactors.temp = 8, rotate.temp = "varimax", print.plot = T,
  recode.trait.names = recode.trait.names)

# for all grassland habitats
all.trait.pcas[["Grassland"]] = make.a.nice.pca.plot.test(
  data.temp = dat.scaled.eunis2 %>% filter(eunis_class1 == "Grassland") %>%
    dplyr::select(CWM_LDMC:CWM_StemDiam),
  nfactors.temp = 8, rotate.temp = "varimax", print.plot = T,
  recode.trait.names = recode.trait.names)

# check amount of captured variation
all.trait.pcas[["Forest"]]$pca_rotated$Vaccounted
all.trait.pcas[["Shrubland"]]$pca_rotated$Vaccounted
all.trait.pcas[["Grassland"]]$pca_rotated$Vaccounted

# attach trait axes scores to dat.scaled.eunis2
for(eunis_class1.temp in as.character(unique(dat.scaled.eunis2$eunis_class1))){
  dat.scaled.eunis2[dat.scaled.eunis2$eunis_class1 == eunis_class1.temp, "CWM_pca1"] =
    all.trait.pcas[[eunis_class1.temp]]$data_to_plot$ind_temp$Dim.1
  dat.scaled.eunis2[dat.scaled.eunis2$eunis_class1 == eunis_class1.temp, "CWM_pca2"] =
    all.trait.pcas[[eunis_class1.temp]]$data_to_plot$ind_temp$Dim.2
  dat.scaled.eunis2[dat.scaled.eunis2$eunis_class1 == eunis_class1.temp, "CWM_pca3"] =
    all.trait.pcas[[eunis_class1.temp]]$data_to_plot$ind_temp$Dim.3
  dat.scaled.eunis2[dat.scaled.eunis2$eunis_class1 == eunis_class1.temp, "CWM_pca4"] =
    all.trait.pcas[[eunis_class1.temp]]$data_to_plot$ind_temp$Dim.4
  dat.scaled.eunis2[dat.scaled.eunis2$eunis_class1 == eunis_class1.temp, "CWM_pca5"] =
    all.trait.pcas[[eunis_class1.temp]]$data_to_plot$ind_temp$Dim.5
  dat.scaled.eunis2[dat.scaled.eunis2$eunis_class1 == eunis_class1.temp, "CWM_pca6"] =
    all.trait.pcas[[eunis_class1.temp]]$data_to_plot$ind_temp$Dim.6
  dat.scaled.eunis2[dat.scaled.eunis2$eunis_class1 == eunis_class1.temp, "CWM_pca7"] =
    all.trait.pcas[[eunis_class1.temp]]$data_to_plot$ind_temp$Dim.7
  dat.scaled.eunis2[dat.scaled.eunis2$eunis_class1 == eunis_class1.temp, "CWM_pca8"] =
    all.trait.pcas[[eunis_class1.temp]]$data_to_plot$ind_temp$Dim.8}


# 4. check relationship between bioclim and BIOCLIM+ variables -----------------
pca.climate = PCA(dplyr::select(dat.scaled.eunis2,
                                chelsa_cmi, chelsa_gsl, chelsa_gsp, chelsa_gst,
                                chelsa_npp, chelsa_pet, chelsa_rsds, chelsa_sfcWind, chelsa_swb,
                                bio_01, bio_02, bio_03, bio_04, bio_05, bio_06,
                                bio_07, bio_08, bio_09, bio_10, bio_11, bio_12,
                                bio_13, bio_14, bio_15, bio_16, bio_17, bio_18, bio_19))

# amount of explaine variation in bioclim variables with 6 bioclim+ variables
rda(dat.log %>% dplyr::select(bio_01, bio_02, bio_03, bio_04, bio_05, bio_06,
                              bio_07, bio_08, bio_09, bio_10, bio_11, bio_12,
                              bio_13, bio_14, bio_15, bio_16, bio_17, bio_18, bio_19) ~
      chelsa_cmi + chelsa_gsl + chelsa_gsp + chelsa_gst + chelsa_rsds + chelsa_sfcWind ,
    data = dat.log)
# 94% of variation

# -> use  the following climate variables
# chelsa_cmi - climate moisture index
# chelsa_gsl - growthing season length
# chelsa_gsp - growing season precipitation
# chelsa_gst - growing season temperature
# chelsa_rsds - Surface downwelling shortwave radiation
# chelsa_sfcWind -  near-surface wind speed

# 5. calculate ANOVAs between habitat types ----------------------------------------------------------------
data.for.cfp.dist.eunis2 = data.frame("eunis_class2" = names(empty.list.for.results)) %>%
  left_join(dat.log, by = "eunis_class2", multiple = "all")

list.cfp.dist.eunis2 = list(list("prop_sis_reflected" = NA, "et" = NA, "npp" = NA))
list.cfp.dist.letters.eunis.2 = list.cfp.dist.eunis2

# linear models
list.cfp.dist.eunis2$prop_sis_reflected = 
  lm(prop_sis_reflected ~ eunis_class2, data = data.for.cfp.dist.eunis2)
list.cfp.dist.eunis2$et = 
  lm(et ~ eunis_class2, data = data.for.cfp.dist.eunis2)
list.cfp.dist.eunis2$npp = 
  lm(npp ~ eunis_class2, data = data.for.cfp.dist.eunis2)

# anova
summary(list.cfp.dist.eunis2$prop_sis_reflected)
summary(list.cfp.dist.eunis2$et)
summary(list.cfp.dist.eunis2$npp)

# get multicomp letters, i.e. posthoc comparisons of mean corrected for multiple comparison
letters.prop.sis.reflected = 
  multcompLetters4(aov(prop_sis_reflected ~ eunis_class2, data = data.for.cfp.dist.eunis2),
                   TukeyHSD(aov(prop_sis_reflected ~ eunis_class2, data = data.for.cfp.dist.eunis2)))$eunis_class2$Letters
letters.et = 
  multcompLetters4(aov(et ~ eunis_class2, data = data.for.cfp.dist.eunis2),
                   TukeyHSD(aov(et ~ eunis_class2, data = data.for.cfp.dist.eunis2)))$eunis_class2$Letters
letters.npp = 
  multcompLetters4(aov(npp ~ eunis_class2, data = data.for.cfp.dist.eunis2),
                   TukeyHSD(aov(npp ~ eunis_class2, data = data.for.cfp.dist.eunis2)))$eunis_class2$Letters

# check results
letters.prop.sis.reflected
letters.et
letters.npp

# 6. account for effects of climate on CRPs ---------------------------
# create empty lists to store results
lm.eunis2 = empty.list.for.results
lm.expl.var.eunis2 = empty.list.for.results

# add columns to the data to store the residuals after accounting for climate
dat.scaled.eunis2$resid_prop_sis_reflected = NA
dat.scaled.eunis2$resid_et = NA
dat.scaled.eunis2$resid_npp = NA

# for each habitat type, run linear models for each CRP 
for(i in 1:length(lm.eunis2)){
  
  lines.temp = which(dat.scaled.eunis2$eunis_class2 == names(lm.eunis2)[i])
  
  lm.eunis2[[i]]$prop_sis_reflected = 
    lm(prop_sis_reflected ~ chelsa_gsl + chelsa_gsp + chelsa_gst + chelsa_cmi + chelsa_rsds + chelsa_sfcWind, 
       data = dat.scaled.eunis2[lines.temp,])
  
  lm.eunis2[[i]]$et = 
    lm(et ~ chelsa_gsl + chelsa_gsp + chelsa_gst + chelsa_cmi + chelsa_rsds + chelsa_sfcWind, 
       data = dat.scaled.eunis2[lines.temp,])
  
  lm.eunis2[[i]]$npp = 
    lm(npp ~ chelsa_gsl + chelsa_gsp + chelsa_gst + chelsa_cmi + chelsa_rsds + chelsa_sfcWind, 
       data = dat.scaled.eunis2[lines.temp,])
  
  dat.scaled.eunis2$resid_prop_sis_reflected[lines.temp] = lm.eunis2[[i]]$prop_sis_reflected$residuals
  dat.scaled.eunis2$resid_et[lines.temp] = lm.eunis2[[i]]$et$residuals
  dat.scaled.eunis2$resid_npp[lines.temp] = lm.eunis2[[i]]$npp$residuals
}

# 8. analyse climate-residuals with random forests ------------------------------------------------

# create empty list to store results
rf.eunis2 = empty.list.for.results

# for each habitat, run random forest model for each CRP
for(i in 1:length(rf.eunis2)){
  
  dat.temp = filter(dat.scaled.eunis2, eunis_class2 == names(rf.eunis2)[i])
  
  #tuneRF(x  = dplyr::select(dat.temp, contains("CWM_pca")), y = dat.temp$resid_prop_sis_reflected, 
  #       ntreeTry=1000, stepFactor=1, improve=0.05)
  
  a = Sys.time()
  rf.eunis2[[i]]$prop_sis_reflected = 
    randomForest(resid_prop_sis_reflected ~ CWM_pca1 + CWM_pca2 + CWM_pca3 + CWM_pca4 + 
                   CWM_pca5 + CWM_pca6 + CWM_pca7 + CWM_pca8, 
                 data = dat.temp,
                 ntree = 2000, mtry = 2, type = "regression", importance = T, do.trace = 250)
  print(Sys.time() - a)
  
  a = Sys.time()
  rf.eunis2[[i]]$et = 
    randomForest(et ~ CWM_pca1 + CWM_pca2 + CWM_pca3 + CWM_pca4 + 
                   CWM_pca5 + CWM_pca6 + CWM_pca7 + CWM_pca8, 
                 data = dat.temp, 
                 ntree = 2000, mtry = 2, type = "regression", importance = T, do.trace = 250)
  print(Sys.time() - a)
  
  a = Sys.time()
  rf.eunis2[[i]]$npp = 
    randomForest(npp ~ CWM_pca1 + CWM_pca2 + CWM_pca3 + CWM_pca4 + 
                   CWM_pca5 + CWM_pca6 + CWM_pca7 + CWM_pca8, 
                 data = dat.temp, 
                 ntree = 2000, mtry = 2, type = "regression", importance = T, do.trace = 250)
  print(Sys.time() - a)
  print(paste0("done: ", names(rf.eunis2)[i]))
}

# 8. plot predicted variation in CFPs -------------------------------------

# create empty tables to store the predicted variation in CRP per habitat type nad CRP
data.to.plot.prop.sis.reflected = tibble(
  eunis_class2 = names(lm.eunis2), 
  eunis_class1 = gsub("_.*", "", names(lm.eunis2)), 
  cfp = rep("prop_sis_reflected", length(rf.eunis2)),
  lm.r2 = rep(NA, length(rf.eunis2)), 
  rf.r2 = rep(NA, length(rf.eunis2)), unexpl = rep(NA, length(rf.eunis2)))

data.to.plot.et = data.to.plot.prop.sis.reflected
data.to.plot.et$cfp = "et"

data.to.plot.npp = data.to.plot.prop.sis.reflected
data.to.plot.npp$cfp = "npp"

# for each habitat type, fill the table with the proportion of explained variation
for(i in 1:length(rf.eunis2)){
  
  # get climate R²
  data.to.plot.prop.sis.reflected$lm.r2[[i]] = 
    summary(lm.eunis2[[i]]$prop_sis_reflected)$adj.r.squared * 100
  data.to.plot.et$lm.r2[[i]] = summary(lm.eunis2[[i]]$et)$adj.r.squared * 100
  data.to.plot.npp$lm.r2[[i]] = summary(lm.eunis2[[i]]$npp)$adj.r.squared * 100
  
  # get trait R² (accounting for the variation that is already captured by climate)
  data.to.plot.prop.sis.reflected$rf.r2[[i]] = 
    (100 - data.to.plot.prop.sis.reflected$lm.r2[[i]]) *  max(rf.eunis2[[i]]$prop_sis_reflected$rsq)
  data.to.plot.et$rf.r2[[i]] = 
    (100 - data.to.plot.et$lm.r2[[i]]) * max(rf.eunis2[[i]]$et$rsq)
  data.to.plot.npp$rf.r2[[i]] = 
    (100 - data.to.plot.npp$lm.r2[[i]]) * max(rf.eunis2[[i]]$npp$rsq)
  
  # replace negative R² with zero
  if(data.to.plot.prop.sis.reflected$lm.r2[[i]] < 0){
    data.to.plot.prop.sis.reflected$lm.r2[[i]] = 0}
  if(data.to.plot.et$lm.r2[[i]] < 0){
    data.to.plot.et$lm.r2[[i]] = 0}
  if(data.to.plot.npp$lm.r2[[i]] < 0){
    data.to.plot.npp$lm.r2[[i]] = 0}
  
  if(data.to.plot.prop.sis.reflected$rf.r2[[i]] < 0){
    data.to.plot.prop.sis.reflected$rf.r2[[i]] = 0}
  if(data.to.plot.et$rf.r2[[i]] < 0){
    data.to.plot.et$rf.r2[[i]] = 0}
  if(data.to.plot.npp$rf.r2[[i]] < 0){
    data.to.plot.npp$rf.r2[[i]] = 0}
  
  # calculated unexplained variation
  data.to.plot.prop.sis.reflected$unexpl[[i]] = 
    100 - data.to.plot.prop.sis.reflected$lm.r2[i] - data.to.plot.prop.sis.reflected$rf.r2[i]
  data.to.plot.et$unexpl[[i]] = 
    100 - data.to.plot.et$lm.r2[i] - data.to.plot.et$rf.r2[i]
  data.to.plot.npp$unexpl[[i]] = 
    100 - data.to.plot.npp$lm.r2[i] - data.to.plot.npp$rf.r2[i]
}

# re-arrange data to long format
data.to.plot = bind_rows(data.to.plot.prop.sis.reflected, data.to.plot.et, data.to.plot.npp) %>%
  arrange(eunis_class1, eunis_class2) %>%
  pivot_longer(cols = c(4,5,6)) %>%
  mutate(eunis_class2 = gsub("_", " - ", eunis_class2)) %>%
  mutate(color.code = ifelse(name == "lm.r2", "Climate (lm)", ifelse(name == "rf.r2", eunis_class2, "Unexplained")))

# change order of color levels for correct plotting
data.to.plot = data.to.plot %>%
  mutate(color.code = 
           factor(color.code, levels = c(
             "Unexplained", "Climate (lm)", 
             "Forest - coniferous", "Forest - deciduous", "Forest - evergreen",
             "Shrubland - alpine", "Shrubland - heathland", "Shrubland - temperate",
             "Grassland - alpine", "Grassland - dry", "Grassland - mesic", "Grassland - wet")),
         eunis_class2 = factor(eunis_class2, levels = c(
           "Forest - coniferous", "Forest - deciduous", "Forest - evergreen",
           "Shrubland - alpine", "Shrubland - heathland", "Shrubland - temperate",
           "Grassland - alpine", "Grassland - dry", "Grassland - mesic", "Grassland - wet")))

# create plots
ggplot(data = data.to.plot %>%
         mutate(cfp = recode(cfp, 
                             "prop_sis_reflected" = "a) Reflected irradiation (%)<br>",
                             "et" = "b) Log evapotranspiration<br>(kg m<sup>-2</sup> year<sup>-1</sup>)",
                             "npp" = "c) Net primary productivity<br>(kg C m<sup>-2</sup> year<sup>-1</sup>)")) %>%
         mutate(eunis_class2 = factor(eunis_class2, levels = rev(levels(data.to.plot$eunis_class2))))) +
  geom_bar(aes(x = eunis_class2, y = value, fill = color.code), 
           color = "grey", alpha = 0.9, stat = "identity", linewidth = 0.01) + 
  geom_hline(yintercept = 0) + 
  facet_grid(. ~ cfp) + 
  coord_flip() + 
  scale_fill_manual(name = "", values = c("white", "#ededed", color.gradient)) +
  scale_x_discrete(position = "top") +
  xlab("") + ylab("Explained variation (%)") +
  theme_classic() +
  theme(panel.grid = element_blank(), strip.text = element_markdown(),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
        strip.background = element_blank())

# 9. plot variable importance ----------------------------------------------------

# create empty list to store importance values
rf.importance.eunis2 = empty.list.for.results

# for each habitat type, extract the proportion of explained variation for each CRP
for(i in 1:length(rf.eunis2)){
  
  rf.importance.eunis2[[names(rf.eunis2)[i]]]$prop_sis_reflected = 
    importance(rf.eunis2[[names(rf.eunis2)[i]]]$prop_sis_reflected)
  rf.importance.eunis2[[names(rf.eunis2)[i]]]$et = 
    importance(rf.eunis2[[names(rf.eunis2)[i]]]$et)
  rf.importance.eunis2[[names(rf.eunis2)[i]]]$npp = 
    importance(rf.eunis2[[names(rf.eunis2)[i]]]$npp)
}

# re-arrange variable importance tables
# for each habitat type
for(i in length(rf.eunis2)){
  # for each CRP
  for(j in 1:3){
    rf.importance.eunis2[[i]][[j]] = as.data.frame(rf.importance.eunis2[[i]][[j]]) %>%
      mutate(variable = rownames(rf.importance.eunis2[[i]][[j]]))
    names(rf.importance.eunis2[[i]][[j]]) = c("rf_inc_mse", "rf_inc_node_impurity", "trait")
    rownames(rf.importance.eunis2[[i]][[j]]) = NULL
    rf.importance.eunis2[[i]][[j]] = rf.importance.eunis2[[i]][[j]] %>%
      arrange(desc(rf_inc_mse))}}

# re-arrange variable importance tables into one table
data.to.plot = bind_rows(
  list("Forest_coniferous" = bind_rows(rf.importance.eunis2$Forest_coniferous, .id = "cfp"),
       "Forest_deciduous" = bind_rows(rf.importance.eunis2$Forest_deciduous, .id = "cfp"),
       "Forest_evergreen" = bind_rows(rf.importance.eunis2$Forest_evergreen, .id = "cfp"),
       "Shrubland_alpine" = bind_rows(rf.importance.eunis2$Shrubland_alpine, .id = "cfp"),
       "Shrubland_heathland" = bind_rows(rf.importance.eunis2$Shrubland_heathland, .id = "cfp"),
       "Shrubland_temperate" = bind_rows(rf.importance.eunis2$Shrubland_temperate, .id = "cfp"),
       "Grassland_alpine" = bind_rows(rf.importance.eunis2$Grassland_alpine, .id = "cfp"),
       "Grassland_dry" = bind_rows(rf.importance.eunis2$Grassland_dry, .id = "cfp"),
       "Grassland_mesic" = bind_rows(rf.importance.eunis2$Grassland_mesic, .id = "cfp"),
       "Grassland_wet" = bind_rows(rf.importance.eunis2$Grassland_wet, .id = "cfp")),
  .id = "habitat") %>%
  filter(trait %in% c("CWM_pca1", "CWM_pca2", "CWM_pca3", "CWM_pca4")) %>%
  mutate(trait = gsub("CWM_pca", "CWM PC", trait)) %>%
  mutate(eunis_class1 = gsub("_.*", "", habitat)) %>%
  group_by(habitat, cfp) %>%
  mutate(rank_mse = rank(-rf_inc_mse)) %>%
  ungroup() %>%
  mutate(habitat = gsub(".*_", "", habitat)) %>%
  arrange(trait) %>%
  mutate(trait = factor(trait, levels = rev(unique(trait)))) %>%
  unite(habitat_for_plotting, eunis_class1, habitat, sep = " - ", remove = F)

# calculate average rank-based correlation between variable importance ranks 
all.correlations = data.to.plot %>%
  distinct(eunis_class1, cfp) %>%
  mutate(r = NA)

for(i in 1:nrow(all.correlations)){
  
  dat.temp = filter(data.to.plot, 
                    eunis_class1 == all.correlations$eunis_class1[i] &
                      cfp == all.correlations$cfp[i])
  
  to.compare = expand_grid(habitat1 = unique(dat.temp$habitat), 
                           habitat2 = unique(dat.temp$habitat)) %>%
    filter(habitat1 != habitat2)
 
  cor.temp = c()
  
  for(j in 1:nrow(to.compare)){
    cor.temp = c(cor.temp, cor(
      dat.temp %>% filter(habitat == to.compare$habitat1[j]) %>% pull(rank_mse),
      dat.temp %>% filter(habitat == to.compare$habitat2[j]) %>% pull(rank_mse),
      method = "spearman"))}
  
  all.correlations$r[i] = round(mean(cor.temp), 2)
}

# replace names of the first two CWM trait axes
data.to.plot$trait = as.character(data.to.plot$trait)
data.to.plot$trait[data.to.plot$trait == "CWM PC1"] = "Leaf economics spectrum"
data.to.plot$trait[data.to.plot$trait == "CWM PC2"] = "Leaf mass"
data.to.plot$trait = factor(data.to.plot$trait, levels = c(
  "CWM PC4", "CWM PC3", "Leaf mass", "Leaf economics spectrum"))

# plot variable importance for proportion of reflected irradiation
ggplot(data.to.plot %>% 
         filter(cfp == "prop_sis_reflected") %>%
         left_join(all.correlations, by = c("eunis_class1", "cfp")) %>%
         unite(eunis_class1, eunis_class1, r, sep = ", r = ") %>%
         mutate(habitat_for_plotting = factor(habitat_for_plotting, levels = c(
           "Forest - coniferous", "Forest - deciduous", "Forest - evergreen",
           "Shrubland - alpine", "Shrubland - heathland", "Shrubland - temperate",
           "Grassland - alpine", "Grassland - dry", "Grassland - mesic", "Grassland - wet")),
           eunis_class1 = factor(eunis_class1, levels = unique(eunis_class1)))) +
  geom_point(aes(x = rf_inc_mse, y = trait, fill = habitat_for_plotting), 
             color = "black", pch = 21, alpha = 0.8, size = 6) +
  geom_text(aes(x = rf_inc_mse, y = trait, label = rank_mse), 
            size = 3,  show.legend = F, alpha = 0.8) +
  facet_grid(. ~ eunis_class1, scales = "free") +
  scale_fill_manual(name = "Habitat", values = color.gradient) +
  theme_bw() +
  ylab("") + xlab("% Increase in MSE") + 
  ggtitle("a) Reflected irradiation (%)") + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "right", plot.title = element_markdown(), strip.background = element_rect(fill = "white"))

# plot variable importance for evapotranspiration
ggplot(data.to.plot %>% 
         filter(cfp == "et") %>%
         left_join(all.correlations, by = c("eunis_class1", "cfp")) %>%
         unite(eunis_class1, eunis_class1, r, sep = ", r = ") %>%
         mutate(habitat_for_plotting = factor(habitat_for_plotting, levels = c(
           "Forest - coniferous", "Forest - deciduous", "Forest - evergreen",
           "Shrubland - alpine", "Shrubland - heathland", "Shrubland - temperate",
           "Grassland - alpine", "Grassland - dry", "Grassland - mesic", "Grassland - wet")),
           eunis_class1 = factor(eunis_class1, levels = unique(eunis_class1)))) +
  geom_point(aes(x = rf_inc_mse, y = trait, fill = habitat_for_plotting), 
             color = "black", pch = 21, alpha = 0.8, size = 6) +
  geom_text(aes(x = rf_inc_mse, y = trait, label = rank_mse), 
            size = 3,  show.legend = F, alpha = 0.8) +
  facet_grid(. ~ eunis_class1, scales = "free") +
  scale_fill_manual(name = "Habitat", values = color.gradient) +
  theme_bw() +
  ylab("") + xlab("% Increase in MSE") + 
  ggtitle("b) Log evapotranspiration (kg m<sup>-2</sup> year<sup>-1</sup>)") + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "right", plot.title = element_markdown(), strip.background = element_rect(fill = "white"))

# plot variable importance for net primary productivity
ggplot(data.to.plot %>% 
         filter(cfp == "npp") %>%
         left_join(all.correlations, by = c("eunis_class1", "cfp")) %>%
         unite(eunis_class1, eunis_class1, r, sep = ", r = ") %>%
         mutate(habitat_for_plotting = factor(habitat_for_plotting, levels = c(
           "Forest - coniferous", "Forest - deciduous", "Forest - evergreen",
           "Shrubland - alpine", "Shrubland - heathland", "Shrubland - temperate",
           "Grassland - alpine", "Grassland - dry", "Grassland - mesic", "Grassland - wet")),
           eunis_class1 = factor(eunis_class1, levels = unique(eunis_class1)))) +
  geom_point(aes(x = rf_inc_mse, y = trait, fill = habitat_for_plotting), 
             color = "black", pch = 21, alpha = 0.8, size = 6) +
  geom_text(aes(x = rf_inc_mse, y = trait, label = rank_mse), 
            size = 3,  show.legend = F, alpha = 0.8) +
  facet_grid(. ~ eunis_class1, scales = "free") +
  scale_fill_manual(name = "Habitat", values = color.gradient) +
  theme_bw() +
  ylab("") + xlab("% Increase in MSE") + 
  ggtitle("c) Net primary productivity (kg C m<sup>-2</sup> year<sup>-1</sup>)") + 
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "right", plot.title = element_markdown(), strip.background = element_rect(fill = "white"))

# 10. plot partial dependence plots --------------------------------------------

# create empty list to store partial dependence results
pd.list =  empty.list.for.results

# create empty table to store for all partial dependence analyses together
pd.results = tibble(value = 1,  
                    predicted = 1,
                    habitat = "test",
                    cfp = "test",
                    trait = "test")[0,]

# for each habitat type
for(habitat.temp in names(pd.list)){
  # for each CRP
  for(cfp.temp in names(pd.list[[habitat.temp]])){
    
    # get partial dependence relationships
    pd.results.temp = edarf::partial_dependence(
      fit = rf.eunis2[[habitat.temp]][[cfp.temp]],
      vars = traits.selected,
      n = c(100, 100),
      interaction = F,
      data = dat.scaled.eunis2 %>%
        filter(eunis_class2 == habitat.temp) %>%
        dplyr::select(all_of(c(cfp.temp, traits.selected))))  %>% # data must be the same as used for the model
      as.data.frame()
    
    # store results in data frame format
    pd.results.temp2 = data.frame(value = c(pd.results.temp[,1], pd.results.temp[,2], 
                                            pd.results.temp[,3], pd.results.temp[,4],
                                            pd.results.temp[,5], pd.results.temp[,6],
                                            pd.results.temp[,7], pd.results.temp[,8]),
                                  predicted = pd.results.temp[,ncol(pd.results.temp)],
                                  habitat = habitat.temp, 
                                  cfp = cfp.temp,
                                  trait = rep(names(pd.results.temp)[1:8], each = nrow(pd.results.temp))) %>%
      as_tibble() %>%
      drop_na()
    
    # add to to the overal partial dependence results table
    pd.results = bind_rows(pd.results, pd.results.temp2)
    
    # print current status
    print(paste0("done: ", habitat.temp,  "  ", cfp.temp))
  }
}

# rename habitat names
pd.results = pd.results %>%
  mutate(eunis_class1 = gsub("_.*", "", habitat),
         habitat = gsub("_", " - ", habitat)) %>%
  mutate(eunis_class1 = gsub(" -.*", "", eunis_class1))

# rename names of trait axes
to.replace = data.frame(
  "habitat" = rep(c("Forest", "Shrubland", "Grassland"), each = 8),
  "axis_old" = rep(paste("CWM_pca", 1:8, sep = ""), 3),
  "axis_new" = c(c("Leaf economics spectrum", "Leaf mass", "Leaf thickness", "Plant height",
                   "Stem density", "Seed mass", "Rooting depth", "Stem conduit diamter"),
                 c("Leaf economics spectrum", "Leaf mass", "Stem conduit diamter", "Plant height",
                   "Leaf thickness", "Specific root length", "Stem density", "Rooting depth"),
                 c("Leaf economics spectrum", "Leaf mass", "Leaf water content", "Plant height",
                   "Rooting depth", "Stem density", "Stem conduit diameter", "Stem conduit density")))

for(i in 1:nrow(to.replace)){
  pd.results$trait[pd.results$eunis_class1 == to.replace$habitat[i] &
                     pd.results$trait == to.replace$axis_old[i]] = to.replace$axis_new[i]}

# rescale x axis for plotting
 pd.results = pd.results %>%
   group_by(eunis_class1, cfp, trait) %>%
   mutate(value = scales::rescale(value, to = c(-1, 1))) %>%
   ungroup()

# plot partial dependence plots for proportion of reflected irradiation
pd.results %>%
  filter(cfp == "prop_sis_reflected" & trait %in% c("Leaf economics spectrum", 
                                                    "Leaf mass",
                                                    "Plant height")) %>%
  mutate(habitat = factor(habitat, levels = c(
    "Forest - coniferous", "Forest - deciduous", "Forest - evergreen",
    "Shrubland - alpine", "Shrubland - heathland", "Shrubland - temperate",
    "Grassland - alpine", "Grassland - dry", "Grassland - mesic", "Grassland - wet")),
    eunis_class1 = factor(eunis_class1, levels = c("Forest", "Shrubland", "Grassland"))) %>%
  ggplot() +
  geom_smooth(aes(y = predicted, x = value, color = habitat), 
              linewidth = 0.8, alpha = 0.8, fill = NA, method = "loess", span = 0.2) +
  scale_color_manual(name = "", values = color.gradient) +
  ggh4x::facet_grid2(eunis_class1 ~ trait, scales = "free_y", independent = "y") + 
  xlab("Trait axis scores low \u2192 high") + ylab("Climate-accounted residuals (predicted)") + 
  ggtitle("a) Reflected irradiation (%)") +
  theme_bw()  +
  theme(panel.grid = element_blank(), axis.line=element_line(), 
        strip.background.x = element_rect(fill = "#eaeaea"), strip.text.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.title = element_markdown(size = 20), legend.title = element_text(size=20), 
        legend.text = element_text(size=11))

# plot partial dependence plots for evapotranspiration
pd.results %>%
  filter(cfp == "et" & trait %in% c("Leaf economics spectrum", 
                                    "Leaf mass",
                                    "Plant height")) %>%
  mutate(habitat = factor(habitat, levels = c(
    "Forest - coniferous", "Forest - deciduous", "Forest - evergreen",
    "Shrubland - alpine", "Shrubland - heathland", "Shrubland - temperate",
    "Grassland - alpine", "Grassland - dry", "Grassland - mesic", "Grassland - wet")),
    eunis_class1 = factor(eunis_class1, levels = c("Forest", "Shrubland", "Grassland"))) %>%
  ggplot() +
  geom_smooth(aes(y = predicted, x = value, color = habitat), 
              linewidth = 0.8, alpha = 0.8, fill = NA, method = "loess", span = 0.2) +
  scale_color_manual(name = "", values = color.gradient) +
  ggh4x::facet_grid2(eunis_class1 ~ trait, scales = "free_y", independent = "y") + 
  xlab("Trait axis scores low \u2192 high") + ylab("Climate-accounted residuals (predicted)") + 
  ggtitle("b) Log evapotranspiration (kg m<sup>-2</sup> year<sup>-1</sup>)") +
  theme_bw()  +
  theme(panel.grid = element_blank(), axis.line=element_line(), 
        strip.background.x = element_rect(fill = "#eaeaea"), strip.text.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.title = element_markdown(size = 20), legend.title = element_text(size=20), 
        legend.text = element_text(size=11))

# plot partial dependence plots for net primary productivity
pd.results %>%
  filter(cfp == "npp" & trait %in% c("Leaf economics spectrum", 
                                     "Leaf mass",
                                     "Plant height")) %>%
  mutate(habitat = factor(habitat, levels = c(
    "Forest - coniferous", "Forest - deciduous", "Forest - evergreen",
    "Shrubland - alpine", "Shrubland - heathland", "Shrubland - temperate",
    "Grassland - alpine", "Grassland - dry", "Grassland - mesic", "Grassland - wet")),
    eunis_class1 = factor(eunis_class1, levels = c("Forest", "Shrubland", "Grassland"))) %>%
  ggplot() +
  geom_smooth(aes(y = predicted, x = value, color = habitat), 
              linewidth = 0.8, alpha = 0.8, fill = NA, method = "loess", span = 0.2) +
  scale_color_manual(name = "", values = color.gradient) +
  ggh4x::facet_grid2(eunis_class1 ~ trait, scales = "free_y", independent = "y") + 
  xlab("Trait axis scores low \u2192 high") + ylab("Climate-accounted residuals (predicted)") + 
  ggtitle("c) Net primary productivity (kg C m<sup>-2</sup> year<sup>-1</sup>)") +
  theme_bw()  +
  theme(panel.grid = element_blank(), axis.line=element_line(), 
        strip.background.x = element_rect(fill = "#eaeaea"), strip.text.y = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.title = element_markdown(size = 20), legend.title = element_text(size=20), 
        legend.text = element_text(size=11))