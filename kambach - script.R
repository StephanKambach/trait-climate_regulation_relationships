##########################################################
# Script to analyse effect on climate-regulating services
# script by: Dr. Stephan Kambach (stephan.kambach@gmail.com)
# initial date: 17.03.2022
########################################################

rm(list=ls())
gc()

# 1. load libraries -------------------------------------------------------
library(tidyverse)
library(tidyselect)
library(vroom)
library(lubridate)
library(FactoMineR)
library(ggforce)
library(corrplot)
library(dismo)
library(viridis)
library(PCAtest)
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

options(na.action = "na.fail")
setwd("")

# load additional functions
source("additional_functions.R")

# 2. load data ------------------------------------------------------------
dat.raw = read_delim("eva_data.csv")

# 3. transform data, calculate means per pixel and scale values ----------
# log transformation
dat.log = dat.raw %>%
  mutate(sp_richness = log(sp.richness),
         et = log(et),

         CWM_LDMC = log(CWM_LDMC), CWV_LDMC = log(CWV_LDMC),
         CWM_Leaf_delta_15N = CWM_Leaf_delta_15N, CWV_Leaf_delta_15N = CWV_Leaf_delta_15N, 
         CWM_LeafArea_leaf_undef = log(CWM_LeafArea_leaf_undef), CWV_LeafArea_leaf_undef = log(CWV_LeafArea_leaf_undef),
         CWM_LeafC_perdrymass = log(CWM_LeafC_perdrymass), CWV_LeafC_perdrymass = log(CWV_LeafC_perdrymass),
         CWM_LeafCN_ratio = log(CWM_LeafCN_ratio), CWV_LeafCN_ratio = log(CWV_LeafCN_ratio),
         CWM_LeafDryMass_single = log(CWM_LeafDryMass_single), CWV_LeafDryMass_single = log(CWV_LeafDryMass_single),
         CWM_Leaffreshmass = log(CWM_Leaffreshmass), CWV_Leaffreshmass = log(CWV_Leaffreshmass),
         CWM_LeafN = log(CWM_LeafN), CWV_LeafN = log(CWV_LeafN),
         CWM_LeafP = log(CWM_LeafP), CWV_LeafP = log(CWV_LeafP),
         CWM_LeafThickness = log(CWM_LeafThickness), CWV_LeafThickness = log(CWV_LeafThickness),
         CWM_LeafWaterCont = log(CWM_LeafWaterCont), CWV_LeafWaterCont = log(CWV_LeafWaterCont),
         CWM_PlantHeight = log(CWM_PlantHeight), CWV_PlantHeight = log(CWV_PlantHeight),
         CWM_RootingDepth = log(CWM_RootingDepth), CWV_RootingDepth = log(CWV_RootingDepth),
         CWM_SLA = log(CWM_SLA), CWV_SLA = log(CWV_SLA),
         CWM_SpecificRootLength = log(CWM_SpecificRootLength), CWV_SpecificRootLength = log(CWV_SpecificRootLength),
         CWM_Stem.cond.dens = log(CWM_Stem.cond.dens), CWV_Stem.cond.dens = log(CWV_Stem.cond.dens),
         CWM_StemConduitDiameter = log(CWM_StemConduitDiameter), CWV_StemConduitDiameter = log(CWV_StemConduitDiameter),
         CWM_StemDens = log(CWM_StemDens), CWV_StemDens = log(CWV_StemDens),
         CWM_StemDiam = log(CWM_StemDiam), CWV_StemDiam = log(CWV_StemDiam),
         
         eunis_class1 = factor(eunis_class1),
         eunis_class2 = factor(eunis_class2)) %>%
  
  drop_na(eunis_class2) 

# calculate mean values per pixel
dat.log = dat.log %>%
  group_by(eunis_class1, eunis_class2, 
           pixel_id, pixel_longitude, pixel_latitude) %>%
  summarise(n_plots = length(et), sp_richness = mean(sp.richness), cover_total_cwm = mean(cover_total_cwm),
            prop_sis_reflected = mean(prop_sis_reflected), et = mean(et), npp = mean(npp),
            
            CWM_LDMC = mean(CWM_LDMC), CWV_LDMC = mean(CWV_LDMC),
            CWM_Leaf_delta_15N = mean(CWM_Leaf_delta_15N), CWV_Leaf_delta_15N = mean(CWV_Leaf_delta_15N),
            CWM_LeafArea_leaf_undef = mean(CWM_LeafArea_leaf_undef), CWV_LeafArea_leaf_undef = mean(CWV_LeafArea_leaf_undef),
            CWM_LeafC_perdrymass = mean(CWM_LeafC_perdrymass), CWV_LeafC_perdrymass = mean(CWV_LeafC_perdrymass),
            CWM_LeafCN_ratio = mean(CWM_LeafCN_ratio), CWV_LeafCN_ratio = mean(CWV_LeafCN_ratio),
            CWM_LeafDryMass_single = mean(CWM_LeafDryMass_single), CWV_LeafDryMass_single = mean(CWV_LeafDryMass_single),
            CWM_Leaffreshmass = mean(CWM_Leaffreshmass), CWV_Leaffreshmass = mean(CWV_Leaffreshmass),
            CWM_LeafN = mean(CWM_LeafN), CWV_LeafN = mean(CWV_LeafN),
            CWM_LeafP = mean(CWM_LeafP), CWV_LeafP = mean(CWV_LeafP),
            CWM_LeafThickness = mean(CWM_LeafThickness), CWV_LeafThickness = mean(CWV_LeafThickness),
            CWM_LeafWaterCont = mean(CWM_LeafWaterCont), CWV_LeafWaterCont = mean(CWV_LeafWaterCont),
            CWM_PlantHeight = mean(CWM_PlantHeight), CWV_PlantHeight = mean(CWV_PlantHeight),
            CWM_RootingDepth = mean(CWM_RootingDepth), CWV_RootingDepth = mean(CWV_RootingDepth),
            CWM_SLA = mean(CWM_SLA), CWV_SLA = mean(CWV_SLA),
            CWM_SpecificRootLength = mean(CWM_SpecificRootLength), CWV_SpecificRootLength = mean(CWV_SpecificRootLength),
            CWM_Stem.cond.dens = mean(CWM_Stem.cond.dens), CWV_Stem.cond.dens = mean(CWV_Stem.cond.dens),
            CWM_StemConduitDiameter = mean(CWM_StemConduitDiameter), CWV_StemConduitDiameter = mean(CWV_StemConduitDiameter),
            CWM_StemDens = mean(CWM_StemDens), CWV_StemDens = mean(CWV_StemDens),
            CWM_StemDiam = mean(CWM_StemDiam), CWV_StemDiam = mean(CWV_StemDiam),
            
            chelsa_cmi = mean(chelsa_cmi), chelsa_gsl = mean(chelsa_gsl),
            chelsa_gsp = mean(chelsa_gsp), chelsa_gst = mean(chelsa_gst),
            chelsa_npp = mean(chelsa_npp), chelsa_pet = mean(chelsa_pet),
            chelsa_rsds = mean(chelsa_rsds), chelsa_sfcWind = mean(chelsa_sfcWind),
            chelsa_swb = mean(chelsa_swb), 
            bio_01 = mean(bio_01), bio_02 = mean(bio_02), bio_03 = mean(bio_03), 
            bio_04 = mean(bio_04), bio_05 = mean(bio_05), bio_06 = mean(bio_06),
            bio_07 = mean(bio_07), bio_08 = mean(bio_08), bio_09 = mean(bio_09), 
            bio_10 = mean(bio_10), bio_11 = mean(bio_11), bio_12 = mean(bio_12),
            bio_13 = mean(bio_13), bio_14 = mean(bio_14), bio_15 = mean(bio_15), 
            bio_16 = mean(bio_16), bio_17 = mean(bio_17), 
            bio_18 = mean(bio_18), bio_19 = mean(bio_19)) %>%
  ungroup()

# scale climate and trait variables within eunis2 habitats
dat.scaled.eunis2 = dat.log %>%
  group_by(eunis_class2) %>%
  mutate(CWM_LDMC = scale_this(CWM_LDMC), CWV_LDMC = scale_this(CWV_LDMC),
         CWM_Leaf_delta_15N = scale_this(CWM_Leaf_delta_15N), CWV_Leaf_delta_15N = scale_this(CWV_Leaf_delta_15N),
         CWM_LeafArea_leaf_undef = scale_this(CWM_LeafArea_leaf_undef), CWV_LeafArea_leaf_undef = scale_this(CWV_LeafArea_leaf_undef),
         CWM_LeafC_perdrymass = scale_this(CWM_LeafC_perdrymass), CWV_LeafC_perdrymass = scale_this(CWV_LeafC_perdrymass),
         CWM_LeafCN_ratio = scale_this(CWM_LeafCN_ratio), CWV_LeafCN_ratio = scale_this(CWV_LeafCN_ratio),
         CWM_LeafDryMass_single = scale_this(CWM_LeafDryMass_single), CWV_LeafDryMass_single = scale_this(CWV_LeafDryMass_single),
         CWM_Leaffreshmass = scale_this(CWM_Leaffreshmass), CWV_Leaffreshmass = scale_this(CWV_Leaffreshmass),
         CWM_LeafN = scale_this(CWM_LeafN), CWV_LeafN = scale_this(CWV_LeafN),
         CWM_LeafP = scale_this(CWM_LeafP), CWV_LeafP = scale_this(CWV_LeafP),
         CWM_LeafThickness = scale_this(CWM_LeafThickness), CWV_LeafThickness = scale_this(CWV_LeafThickness),
         CWM_LeafWaterCont = scale_this(CWM_LeafWaterCont), CWV_LeafWaterCont = scale_this(CWV_LeafWaterCont),
         CWM_PlantHeight = scale_this(CWM_PlantHeight), CWV_PlantHeight = scale_this(CWV_PlantHeight),
         CWM_RootingDepth = scale_this(CWM_RootingDepth), CWV_RootingDepth = scale_this(CWV_RootingDepth),
         CWM_SLA = scale_this(CWM_SLA), CWV_SLA = scale_this(CWV_SLA),
         CWM_SpecificRootLength = scale_this(CWM_SpecificRootLength), CWV_SpecificRootLength = scale_this(CWV_SpecificRootLength),
         CWM_Stem.cond.dens = scale_this(CWM_Stem.cond.dens), CWV_Stem.cond.dens = scale_this(CWV_Stem.cond.dens),
         CWM_StemConduitDiameter = scale_this(CWM_StemConduitDiameter), CWV_StemConduitDiameter = scale_this(CWV_StemConduitDiameter),
         CWM_StemDens = scale_this(CWM_StemDens), CWV_StemDens = scale_this(CWV_StemDens),
         CWM_StemDiam = scale_this(CWM_StemDiam), CWV_StemDiam = scale_this(CWV_StemDiam),
         
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
  ungroup()

color.gradient = c("#a6c6e1", "#2171b5", "#13436c", 
                   "#fed98e", "#fe9929", "#993404",
                   "#acdaba", "#45ac65", "#22723a", "#134121")

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

traits.selected = c("CWM_pca1", "CWM_pca2", "CWM_pca3", "CWM_pca4", 
                    "CWM_pca5", "CWM_pca6", "CWM_pca7", "CWM_pca8")

recode.trait.names = function(trait.vector.temp){
  trait.vector.temp[grepl("_LDMC", trait.vector.temp)] = "LDMC"
  trait.vector.temp[grepl("_Leaf_delta_15N", trait.vector.temp)] = "L15N"
  trait.vector.temp[grepl("_LeafArea_leaf_undef", trait.vector.temp)] = "LArea"
  trait.vector.temp[grepl("_LeafC_perdrymass", trait.vector.temp)] = "LCcont"
  trait.vector.temp[grepl("_LeafCN_ratio", trait.vector.temp)] = "LCNratio"
  trait.vector.temp[grepl("_LeafDryMass_single", trait.vector.temp)] = "LMdry"
  trait.vector.temp[grepl("_Leaffreshmass", trait.vector.temp)] = "LMfresh"
  trait.vector.temp[grepl("_LeafN", trait.vector.temp)] = "LNcont"
  trait.vector.temp[grepl("_LeafP", trait.vector.temp)] = "LPcont"
  trait.vector.temp[grepl("_LeafThickness", trait.vector.temp)] = "LThic"
  trait.vector.temp[grepl("_LeafWaterCont", trait.vector.temp)] = "LWcont"
  trait.vector.temp[grepl("_PlantHeight", trait.vector.temp)] = "PlantHeight"
  trait.vector.temp[grepl("_RootingDepth", trait.vector.temp)] = "RDepth"
  trait.vector.temp[grepl("_SLA", trait.vector.temp)] = "SLA"
  trait.vector.temp[grepl("_SpecificRootLength", trait.vector.temp)] = "SRL"
  trait.vector.temp[grepl("_Stem.cond.dens", trait.vector.temp)] = "SCdens"
  trait.vector.temp[grepl("_StemConduitDiameter", trait.vector.temp)] = "SCdiam"
  trait.vector.temp[grepl("_StemDens", trait.vector.temp)] = "Sdens"
  trait.vector.temp[grepl("_StemDiam", trait.vector.temp)] = "Sdiam"
  
  return(trait.vector.temp)
}

# 5. get trait PCA axes -----------------------------------------------
# show correlation
trait.cor = cor(dat.scaled.eunis2 %>%
                  dplyr::select(CWM_LDMC:CWM_StemDiam))
ggcorrplot::ggcorrplot(trait.cor, method = "square") + 
  scale_fill_stepsn(colours = c("blue", "white", "white", "red"))

# reduce variation with a PCA, using 5 axes
all.trait.pcas = list()
dat.scaled.eunis2$CWM_pca1 = as.numeric(NA)
dat.scaled.eunis2$CWM_pca2 = as.numeric(NA)
dat.scaled.eunis2$CWM_pca3 = as.numeric(NA)
dat.scaled.eunis2$CWM_pca4 = as.numeric(NA)
dat.scaled.eunis2$CWM_pca5 = as.numeric(NA)
dat.scaled.eunis2$CWM_pca6 = as.numeric(NA)
dat.scaled.eunis2$CWM_pca7 = as.numeric(NA)
dat.scaled.eunis2$CWM_pca8 = as.numeric(NA)

# calculate varimax rotated trait pcas with 8 axes
all.trait.pcas[["Forest"]] = make.a.nice.pca.plot.test(
  data.temp = dat.scaled.eunis2 %>% filter(eunis_class1 == "Forest") %>%
    dplyr::select(contains("CWM"), - contains("_pca"), - contains("cover_total")),
  nfactors.temp = 8, rotate.temp = "varimax", print.plot = T,
  recode.trait.names = recode.trait.names, switch.axes.temp = "Dim.5")

all.trait.pcas[["Shrubland"]] = make.a.nice.pca.plot.test(
  data.temp = dat.scaled.eunis2 %>% filter(eunis_class1 == "Shrubland") %>%
    dplyr::select(contains("CWM"), - contains("_pca"), - contains("cover_total")),
  nfactors.temp = 8, rotate.temp = "varimax", print.plot = T,
  recode.trait.names = recode.trait.names)

all.trait.pcas[["Grassland"]] = make.a.nice.pca.plot.test(
  data.temp = dat.scaled.eunis2 %>% filter(eunis_class1 == "Grassland") %>%
    dplyr::select(contains("CWM"), - contains("_pca"), - contains("cover_total")),
  nfactors.temp = 8, rotate.temp = "varimax", print.plot = T,
  recode.trait.names = recode.trait.names)

# plot PCAs
png(filename = "feedbacks_project2/results/pixel based analysis/trait PCA/trait_pca_forest.png",
    res = 300, height = 3000, width = 3000)
grid.arrange(all.trait.pcas[["Forest"]]$plots$pca12, 
             all.trait.pcas[["Forest"]]$plots$pca34,
             all.trait.pcas[["Forest"]]$plots$pca56, 
             all.trait.pcas[["Forest"]]$plots$pca78, nrow = 2)
graphics.off()

png(filename = "feedbacks_project2/results/pixel based analysis/trait PCA/trait_pca_grassland.png",
    res = 300, height = 3000, width = 3000)
grid.arrange(all.trait.pcas[["Grassland"]]$plots$pca12, 
             all.trait.pcas[["Grassland"]]$plots$pca34,
             all.trait.pcas[["Grassland"]]$plots$pca56, 
             all.trait.pcas[["Grassland"]]$plots$pca78, nrow = 2)
graphics.off()

png(filename = "feedbacks_project2/results/pixel based analysis/trait PCA/trait_pca_shrubland.png",
    res = 300, height = 3000, width = 3000)
grid.arrange(all.trait.pcas[["Shrubland"]]$plots$pca12, 
             all.trait.pcas[["Shrubland"]]$plots$pca34,
             all.trait.pcas[["Shrubland"]]$plots$pca56, 
             all.trait.pcas[["Shrubland"]]$plots$pca78, nrow = 2)
graphics.off()

# export factor loadings
data.to.export = all.trait.pcas[["Forest"]]$data_to_plot$var_temp_complete %>%
  mutate(habitat = "Forest") %>%
  bind_rows(all.trait.pcas[["Shrubland"]]$data_to_plot$var_temp_complete %>%
              mutate(habitat = "Shrubland")) %>%
  bind_rows(all.trait.pcas[["Grassland"]]$data_to_plot$var_temp_complete %>%
              mutate(habitat = "Grassland")) %>%
  mutate(Dim.1 = round(Dim.1, 2), Dim.2 = round(Dim.2, 2),
         Dim.3 = round(Dim.3, 2), Dim.4 = round(Dim.4, 2),
         Dim.5 = round(Dim.5, 2), Dim.6 = round(Dim.6, 2),
         Dim.7 = round(Dim.7, 2), Dim.8 = round(Dim.8, 2))
write.table(data.to.export, 
            "feedbacks_project2/results/pixel based analysis/trait PCA/all_pca_loadings.txt",
            sep = "\t", dec = ",", quote = F, row.names = F)

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

# 5. check relationship between bioclim and BIOCLIM+ variables -----------------
# check multivariate relationships
pca.climate = PCA(dplyr::select(dat.scaled.eunis2,
                                chelsa_cmi, chelsa_gsl, chelsa_gsp, chelsa_gst,
                                chelsa_npp, chelsa_pet, chelsa_rsds, chelsa_sfcWind, chelsa_swb,
                                bio_01, bio_02, bio_03, bio_04, bio_05, bio_06,
                                bio_07, bio_08, bio_09, bio_10, bio_11, bio_12,
                                bio_13, bio_14, bio_15, bio_16, bio_17, bio_18, bio_19))
cor.climate = cor(dplyr::select(dat.scaled.eunis2,
                                chelsa_cmi, chelsa_gsl, chelsa_gsp, chelsa_gst,
                                chelsa_npp, chelsa_pet, chelsa_rsds, chelsa_sfcWind, chelsa_swb,
                                bio_01, bio_02, bio_03, bio_04, bio_05, bio_06,
                                bio_07, bio_08, bio_09, bio_10, bio_11, bio_12,
                                bio_13, bio_14, bio_15, bio_16, bio_17, bio_18, bio_19))

corrplot(cor.climate)

# make a nice pca plot
data.to.plot = list(ind_temp = data.frame(pca.climate$ind$coord) %>% as_tibble(),
                    var_temp = data.frame(pca.climate$var$coord) %>% 
                      rownames_to_column(var = "climate_var") %>%
                      as_tibble() %>%
                      mutate(climate_var = gsub("_", ":", climate_var)) %>%
                      mutate(dataset = ifelse(grepl("chelsa", climate_var), "CHELSA bioclim+", "Bioclim")))
                    
plot1 = ggplot() +
  geom_hex(data = data.to.plot$ind_temp, 
           aes(x = Dim.1, y = Dim.2), bins = 100) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_segment(data = data.to.plot$var_temp, 
               aes(x = 0, y = 0, xend = Dim.1 * 14, yend = Dim.2 * 14),
               arrow = arrow(length = unit(0.8, "cm"))) +
  geom_label(data = data.to.plot$var_temp, 
             aes(x = Dim.1 * 15, y = Dim.2 * 15, label = climate_var, color = dataset), 
             alpha = 0.6) +
  scale_fill_gradientn(colours=c("#E8E8E8","#808080"),guide = "none", na.value = NA) +
  scale_color_manual(values=c("#808080","blue"), guide = "none") +
  xlab("PC1: 36.7% explained variation") +
  ylab("PC2: 25.7% explained variation") +
  theme_bw() +
  theme(panel.grid = element_blank())

# amount of explaine variation in bioclim variables with 6 bioclim+ variables
rda(dat.log %>% dplyr::select(bio_01, bio_02, bio_03, bio_04, bio_05, bio_06,
                              bio_07, bio_08, bio_09, bio_10, bio_11, bio_12,
                              bio_13, bio_14, bio_15, bio_16, bio_17, bio_18, bio_19) ~
      chelsa_cmi + chelsa_gsl + chelsa_gsp + chelsa_gst + chelsa_rsds + chelsa_sfcWind ,
    data = dat.log)
# 94% of variation

# check bioclim+ correlations
test.pca = PCA(dat.log %>% dplyr::select(
  chelsa_cmi, chelsa_gsl, chelsa_gsp, chelsa_gst, chelsa_rsds, chelsa_sfcWind))
test.cor = cor(dat.log %>% dplyr::select(
  chelsa_cmi, chelsa_gsl, chelsa_gsp, chelsa_gst, chelsa_rsds, chelsa_sfcWind))
corrplot(test.cor, method = "number")

# -> use  the following variables
# chelsa_cmi - climate moisture index
# chelsa_gsl - growthing season length
# chelsa_gsp - growing season precipitation
# chelsa_gst - growing season temperature
# chelsa_rsds - Surface downwelling shortwave radiation
# chelsa_sfcWind -  near-surface wind speed

# 8. get residuals after accounting for climate ---------------------------
# get residuals after climate and spatial dependencies
lm.eunis2 = empty.list.for.results
lm.expl.var.eunis2 = empty.list.for.results

dat.scaled.eunis2$resid_prop_sis_reflected = NA
dat.scaled.eunis2$resid_et = NA
dat.scaled.eunis2$resid_npp = NA

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
rf.eunis2 = empty.list.for.results

traits.selected = c("CWM_pca1", "CWM_pca2", "CWM_pca3", "CWM_pca4",
                    "CWM_pca5", "CWM_pca6", "CWM_pca7", "CWM_pca8")

formula.prop.sis.reflected = as.formula(
  paste0("resid_prop_sis_reflected ~ ", paste(traits.selected, collapse = " + ")))
formula.et = as.formula(
  paste0("resid_et ~ ", paste(traits.selected, collapse = " + ")))
formula.npp = as.formula(
  paste0("resid_npp ~ ", paste(traits.selected, collapse = " + ")))

# random forest models
for(i in 1:length(rf.eunis2)){
  
  dat.temp = filter(dat.scaled.eunis2, eunis_class2 == names(rf.eunis2)[i])
  
  #tuneRF(x  = dplyr::select(dat.temp, contains("CWM_pca")), y = dat.temp$resid_prop_sis_reflected, 
  #       ntreeTry=1000, stepFactor=1, improve=0.05)
  
  a = Sys.time()
  rf.eunis2[[i]]$prop_sis_reflected = 
    randomForest(formula.prop.sis.reflected, data = dat.temp,
                 ntree = 2000, mtry = 2, type = "regression", importance = T, do.trace = 250)
  print(Sys.time() - a)
  
  a = Sys.time()
  rf.eunis2[[i]]$et = 
    randomForest(formula.et, data = dat.temp, 
                 ntree = 2000, mtry = 2, type = "regression", importance = T, do.trace = 250)
  print(Sys.time() - a)
  
  a = Sys.time()
  rf.eunis2[[i]]$npp = 
    randomForest(formula.npp, data = dat.temp, 
                 ntree = 2000, mtry = 2, type = "regression", importance = T, do.trace = 250)
  print(Sys.time() - a)
  print(paste0("done: ", names(rf.eunis2)[i]))
}

# 8. plot predicted variation in CFPs -------------------------------------
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
  
  data.to.plot.prop.sis.reflected$unexpl[[i]] = 
    100 - data.to.plot.prop.sis.reflected$lm.r2[i] - data.to.plot.prop.sis.reflected$rf.r2[i]
  data.to.plot.et$unexpl[[i]] = 
    100 - data.to.plot.et$lm.r2[i] - data.to.plot.et$rf.r2[i]
  data.to.plot.npp$unexpl[[i]] = 
    100 - data.to.plot.npp$lm.r2[i] - data.to.plot.npp$rf.r2[i]
}

data.to.plot = bind_rows(data.to.plot.prop.sis.reflected, data.to.plot.et, data.to.plot.npp) %>%
  arrange(eunis_class1, eunis_class2) %>%
  pivot_longer(cols = c(4,5,6)) %>%
  mutate(eunis_class2 = gsub("_", " - ", eunis_class2)) %>%
  mutate(color.code = ifelse(name == "lm.r2", "Climate (lm)", ifelse(name == "rf.r2", eunis_class2, "Unexplained")))
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

# axes flipped
plot1 = ggplot(data = data.to.plot %>%
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
  theme(panel.grid = element_blank(), strip.text = element_markdown(), legend.position = "none",
        axis.line.y = element_blank(), axis.ticks.y = element_blank(), 
        strip.background = element_blank())

legend1 = get_legend(plot1 + 
                    scale_fill_manual(name = "", 
                                      breaks = c("Unexplained", "Climate (lm)"),
                                      values = c("white", "#ededed", color.gradient)) +
                    theme(legend.position = "bottom"))

grid.arrange(plot1, legend1, nrow = 2, heights = c(1,0.15))

# 8. get variable importance ----------------------------------------------------
rf.importance.eunis2 = empty.list.for.results

# rf.min.tree.depth.eunis2 = empty.list.for.results

for(i in 1:length(rf.eunis2)){
  
  rf.importance.eunis2[[names(rf.eunis2)[i]]]$prop_sis_reflected = 
    importance(rf.eunis2[[names(rf.eunis2)[i]]]$prop_sis_reflected)
  rf.importance.eunis2[[names(rf.eunis2)[i]]]$et = 
    importance(rf.eunis2[[names(rf.eunis2)[i]]]$et)
  rf.importance.eunis2[[names(rf.eunis2)[i]]]$npp = 
    importance(rf.eunis2[[names(rf.eunis2)[i]]]$npp)
  
  # a = Sys.time()
  # rf.min.tree.depth.eunis2[[names(rf.eunis2)[i]]]$prop_sis_reflected =
  #   measure_importance(rf.eunis2[[names(rf.eunis2)[i]]]$prop_sis_reflected, mean_sample  = "top_trees")
  # print(Sys.time() - a)
  # 
  # a = Sys.time()
  # rf.min.tree.depth.eunis2[[names(rf.eunis2)[i]]]$et =
  #   measure_importance(rf.eunis2[[names(rf.eunis2)[i]]]$et, mean_sample  = "top_trees")
  # print(Sys.time() - a)
  # 
  # a = Sys.time()
  # rf.min.tree.depth.eunis2[[names(rf.eunis2)[i]]]$npp =
  #   measure_importance(rf.eunis2[[names(rf.eunis2)[i]]]$npp, mean_sample  = "top_trees")
  # print(Sys.time() - a)
  
  print(paste0("done: ", names(rf.eunis2)[i]))
}

# 10. plot importance --------------------------------------------------------
# arrange importance tables
for(i in 1:10){
  for(j in 1:3){
    rf.importance.eunis2[[i]][[j]] = as.data.frame(rf.importance.eunis2[[i]][[j]]) %>%
      mutate(variable = rownames(rf.importance.eunis2[[i]][[j]]))
    names(rf.importance.eunis2[[i]][[j]]) = c("rf_inc_mse", "rf_inc_node_impurity", "trait")
    rownames(rf.importance.eunis2[[i]][[j]]) = NULL
    rf.importance.eunis2[[i]][[j]] = rf.importance.eunis2[[i]][[j]] %>%
      arrange(desc(rf_inc_mse))}}

# plot importance
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
  mutate(eunis_class1 = gsub("_.*", "", habitat)) %>%
  group_by(habitat, cfp) %>%
  mutate(rank_mse = rank(-rf_inc_mse)) %>%
  ungroup() %>%
  mutate(habitat = gsub(".*_", "", habitat)) %>%
  arrange(trait) %>%
  mutate(trait = factor(trait, levels = rev(unique(trait)))) %>%
  unite(habitat_for_plotting, eunis_class1, habitat, sep = " - ", remove = F)

# calculate average rank-based correlation between ranks 
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
data.to.plot = data.to.plot %>%
  mutate(trait = as.character(trait)) %>%
  mutate(trait = ifelse(eunis_class1 == "Forest" & trait == "CWM_pca1", "LES", trait),
         trait = ifelse(eunis_class1 == "Forest" & trait == "CWM_pca2", "Leaf size", trait),
         trait = ifelse(eunis_class1 == "Forest" & trait == "CWM_pca3", "Trait pc3/4", trait),
         trait = ifelse(eunis_class1 == "Forest" & trait == "CWM_pca4", "Plant height", trait),
         
         trait = ifelse(eunis_class1 == "Shrubland" & trait == "CWM_pca1", "LES", trait),
         trait = ifelse(eunis_class1 == "Shrubland" & trait == "CWM_pca2", "Leaf size", trait),
         trait = ifelse(eunis_class1 == "Shrubland" & trait == "CWM_pca3", "Plant height", trait),
         trait = ifelse(eunis_class1 == "Shrubland" & trait == "CWM_pca4", "Trait pc3/4", trait),
         
         trait = ifelse(eunis_class1 == "Grassland" & trait == "CWM_pca1", "LES", trait),
         trait = ifelse(eunis_class1 == "Grassland" & trait == "CWM_pca2", "Leaf size", trait),
         trait = ifelse(eunis_class1 == "Grassland" & trait == "CWM_pca3", "Plant height", trait),
         trait = ifelse(eunis_class1 == "Grassland" & trait == "CWM_pca4", "Trait pc3/4", trait)) %>%
  mutate(trait = factor(trait, levels = c("Trait pc3/4", "Plant height", "Leaf size", "LES")))

# check mean rank level of explained variation for all 4 traits
data.to.plot %>%
  group_by(trait) %>%
  summarise(mean_rank = mean(rank_mse))

# plotting
plot1 = ggplot(data.to.plot %>% 
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

plot2= ggplot(data.to.plot %>% 
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

plot3 = ggplot(data.to.plot %>% 
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

grid_arrange_shared_legend(plot1, plot2, plot3, ncol = 1, position = "right")

# 11. partial dependence plots --------------------------------------------
# get three best predictors
traits.selected = c("CWM_pca1", "CWM_pca2", "CWM_pca3", "CWM_pca4",
                    "CWM_pca5", "CWM_pca6", "CWM_pca7", "CWM_pca8")

pd.list =  empty.list.for.results

# create empty table to store results
pd.results = tibble(value = 1,  
                    predicted = 1,
                    habitat = "test",
                    cfp = "test",
                    trait = "test")[0,]

for(habitat.temp in names(pd.list)){ # per habitas
  for(cfp.temp in names(pd.list[[habitat.temp]])){ # per cfp
    
    # get partial results
    pd.results.temp = edarf::partial_dependence(
      fit = rf.eunis2[[habitat.temp]][[cfp.temp]],
      vars = traits.selected,
      n = c(100, 100),
      interaction = F,
      data = dat.scaled.eunis2 %>%
        filter(eunis_class2 == habitat.temp) %>%
        dplyr::select(all_of(c(cfp.temp, traits.selected))))  %>% # only the data for the model
      
      as.data.frame()
    
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
    
    # add to to the results table
    pd.results = bind_rows(pd.results, pd.results.temp2)
    
    # update
    print(paste0("done: ", habitat.temp,  "  ", cfp.temp))
  }
}

# 9. shorten fig 6 --------------------------------------------------------
# rename habitats
pd.results = pd.results %>%
  mutate(eunis_class1 = gsub("_.*", "", habitat),
         habitat = gsub("_", " - ", habitat)) %>%
  mutate(eunis_class1 = gsub(" -.*", "", eunis_class1))

# rename trait pca axes
pd.results = pd.results %>%
  mutate(trait = as.character(trait)) %>%
  mutate(trait = ifelse(eunis_class1 == "Forest" & trait == "CWM_pca1", "LES", trait),
         trait = ifelse(eunis_class1 == "Forest" & trait == "CWM_pca2", "Leaf size", trait),
         trait = ifelse(eunis_class1 == "Forest" & trait == "CWM_pca3", "Trait pc3/4", trait),
         trait = ifelse(eunis_class1 == "Forest" & trait == "CWM_pca4", "Plant height", trait),
         
         trait = ifelse(eunis_class1 == "Shrubland" & trait == "CWM_pca1", "LES", trait),
         trait = ifelse(eunis_class1 == "Shrubland" & trait == "CWM_pca2", "Leaf size", trait),
         trait = ifelse(eunis_class1 == "Shrubland" & trait == "CWM_pca3", "Plant height", trait),
         trait = ifelse(eunis_class1 == "Shrubland" & trait == "CWM_pca4", "Trait pc3/4", trait),
         
         trait = ifelse(eunis_class1 == "Grassland" & trait == "CWM_pca1", "LES", trait),
         trait = ifelse(eunis_class1 == "Grassland" & trait == "CWM_pca2", "Leaf size", trait),
         trait = ifelse(eunis_class1 == "Grassland" & trait == "CWM_pca3", "Plant height", trait),
         trait = ifelse(eunis_class1 == "Grassland" & trait == "CWM_pca4", "Trait pc3/4", trait)) %>%
  mutate(trait = factor(trait, levels = c("Trait pc3/4", "Plant height", "Leaf size", "LES")))

# rescale x and y-axis for plotting
 pd.results = pd.results %>%
   group_by(eunis_class1, cfp, trait) %>%
   #mutate(predicted = scales::rescale(predicted, to = c(-1, 1))) %>% 
   mutate(value = scales::rescale(value, to = c(-1, 1))) %>%
   ungroup()

# plotting
plot1 = pd.results %>%
  filter(cfp == "prop_sis_reflected" & trait %in% c("LES", "Leaf size", "Plant height", "Trait pc3/4")) %>%
  mutate(habitat = factor(habitat, levels = c(
    "Forest - coniferous", "Forest - deciduous", "Forest - evergreen",
    "Shrubland - alpine", "Shrubland - heathland", "Shrubland - temperate",
    "Grassland - alpine", "Grassland - dry", "Grassland - mesic", "Grassland - wet")),
    eunis_class1 = factor(eunis_class1, levels = c("Forest", "Shrubland", "Grassland")),
    trait = factor(trait, levels = c("LES", "Leaf size", "Plant height", "Trait pc3/4"))) %>%
  ggplot() +
  #geom_hline(yintercept = 0, linetype = "solid") + 
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

legend.save.for.later = get_legend(plot1)

plot2 = pd.results %>%
  filter(cfp == "et" & trait %in% c("LES", "Leaf size", "Plant height", "Trait pc3/4")) %>%
  mutate(habitat = factor(habitat, levels = c(
    "Forest - coniferous", "Forest - deciduous", "Forest - evergreen",
    "Shrubland - alpine", "Shrubland - heathland", "Shrubland - temperate",
    "Grassland - alpine", "Grassland - dry", "Grassland - mesic", "Grassland - wet")),
    eunis_class1 = factor(eunis_class1, levels = c("Forest", "Shrubland", "Grassland")),
    trait = factor(trait, levels = c("LES", "Leaf size", "Plant height", "Trait pc3/4"))) %>%
  ggplot() +
  #geom_hline(yintercept = 0, linetype = "solid") + 
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

plot3 = pd.results %>%
  filter(cfp == "npp" & trait %in% c("LES", "Leaf size", "Plant height", "Trait pc3/4")) %>%
  mutate(habitat = factor(habitat, levels = c(
    "Forest - coniferous", "Forest - deciduous", "Forest - evergreen",
    "Shrubland - alpine", "Shrubland - heathland", "Shrubland - temperate",
    "Grassland - alpine", "Grassland - dry", "Grassland - mesic", "Grassland - wet")),
    eunis_class1 = factor(eunis_class1, levels = c("Forest", "Shrubland", "Grassland")),
    trait = factor(trait, levels = c("LES", "Leaf size", "Plant height", "Trait pc3/4"))) %>%
  ggplot() +
  #geom_hline(yintercept = 0, linetype = "solid") + 
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

grid_arrange_shared_legend(plot1, plot2, plot3, position = "bottom", ncol = 1)
