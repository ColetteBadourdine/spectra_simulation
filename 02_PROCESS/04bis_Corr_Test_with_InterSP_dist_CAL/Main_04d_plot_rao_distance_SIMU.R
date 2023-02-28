# ==============================================================================
# PROGRAMMERS:
# Colette BADOURDINE <colette.badourdine@cirad.fr>
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Copyright 2021/06 Colette BADOURDINE
# ==============================================================================

################################################################################
# Always start a script with a clean environment
rm(list=ls(all=TRUE));gc()
# define working directory as the directory where the script is located
if (rstudioapi::isAvailable()){
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}
################################################################################

#import
source('../00_Libraries/Lib_ImageProcess.R')
source('../00_Libraries/Lib_variance.R')
source('../00_Libraries/Lib_diversity.R')
source('../00_Libraries/Lib_populations.R')
source('../00_Libraries/Lib_SpectralVarianceAnalysis.R')
library(doFuture)
library(readr)
library(zip)
library(ade4)
library(tools)
library(tidyverse)
library(progress)
library(ggplot2)
library(entropart)
library(dplyr)
library(ggpmisc)
library(ggpubr)


#--------------- define path for data and results ---------------
Path_Populations <- '../../03_RESULTS/03_Produce_POP/'
Path_reflectance <- '../../03_RESULTS/03_Produce_POP/Reflectance/'
Path_variance <- '../../03_RESULTS/04bis_Correlation_FullFeature_SIMU_VNIRCR/'

Path_Results <- '../../03_RESULTS/04_INTER_SPECIES_DISTANCE_CAL_SIMU'
dir.create(path = Path_Results, showWarnings = F, recursive = T)
path_figure <- file.path('../../03_RESULTS/FIGURES/')
dir.create(path = path_figure, showWarnings = F, recursive = T)


#--------------- define path for data and results ---------------
Path_Data <- '../../03_RESULTS/02_Produce_FilteredDatafiles'
File_Data <- 'TreeID_Filtered.zip'

#--------------- read data and discard crowns with too few pixels ---------------
# read data type
TreeIDZipFile <- file.path(Path_Data,File_Data)
fid <- unz(description = TreeIDZipFile, filename = zip_list(TreeIDZipFile)$filename,
           open = "rb", encoding = getOption("encoding"))
TreeID <- read_delim(file = fid,delim = ',')
close(fid)


############## SD = 1 ##############

ResPop = readRDS(file.path(Path_reflectance, "all_populations_simulated_spectra_sd1.rds"))
# compute spectral variance corresponding to a selection of PCs
nbPop <- ResPop$PopID

#--------------- read reflectance data  ---------------
ReflZipFile <- file.path(Path_reflectance,'Spectral_Simu_VNIRCR_sd1.zip')
fid <- unz(description = ReflZipFile, filename = zip::zip_list(ReflZipFile)$filename,
           open = "rb", encoding = getOption("encoding"))
ReflectanceData <- readr::read_delim(file = fid,delim = ',')
close(fid)
for (i in 3:ncol(ReflectanceData)){
  ReflectanceData[,i] <- scale(ReflectanceData[,i],center = T,scale = T)
}


AllVars <- seq(1,(ncol(ReflectanceData)-2))

#replace Spectralpop in ResPop to match the new reflectance data
for (pop in 1:nbPop){
  ResPop$Spectralpop[[pop]] = ReflectanceData[ReflectanceData$Population == pop, -1]
}

Simulated_reflectance_data = bind_rows(ResPop$Spectralpop, .id = "column_label")
Simulated_reflectance_data = Simulated_reflectance_data[, -1]

####JE CALCULE LES SPECTRES MOYENS PAR ESPECE ####
####JE LES UTILISE POUR CALCULER LES DISTANCES INTER-ESPECE####

# calcul des distances inter-espèce par population

ReflectanceData_mean <- compute_MeanSPectraSpecies(Simulated_reflectance_data, as.factor(Simulated_reflectance_data$SPID))
ReflectanceData_mean <- ReflectanceData_mean$mean_spectra_species
ReflectanceData_mean$SPID = as.factor(ReflectanceData_mean$SPID)
rownames(ReflectanceData_mean) = ReflectanceData_mean$SPID
ReflectanceData_mean = ReflectanceData_mean[, -1]

distance_inter_sp = dist(ReflectanceData_mean, method="euclidean")
distance_inter_sp_norm = distance_inter_sp/max(distance_inter_sp)

####JE CALCULE LA VARIANCE TOTALE ET LA VARIANCE INTER-ESPECE####


ResVar_all<- readRDS(file.path(Path_variance, "Var_DivTax_sd1.rds"))


####JE RECUPERE LES TABLES D'ABONDANCE DES 220 POPULATIONS PUIS JE CALCULE LES FREQUENCES####

Spectralpop_all = ResPop$Spectralpop

abundance_table_for_cal_pop <- as.data.frame(matrix(data = 0, nrow = ResPop$PopID, ncol = length(unique(TreeID$SPID))))
colnames(abundance_table_for_cal_pop) = unique(TreeID$SPID)

for (pop in 1:ResPop$PopID){
  abundance_temp <- as.data.frame(table(Spectralpop_all[[pop]]$SPID)/20)
  for (sp in abundance_temp[, 1]){
    abundance_table_for_cal_pop[pop, sp] = abundance_temp[abundance_temp$Var1 == sp, ]$Freq
  }
}


freq <- abundance_table_for_cal_pop/100

rao_all <- rao(freq, as.matrix(distance_inter_sp_norm))


####FIGURES####

data_all = cbind(ResVar_all, rao_all)


rplot_1 = ggplot(data = data_all, aes(x = rao_all, y = var_tot))+
  geom_point(aes(color = 'All pixels'))+
  stat_poly_line(data = data_all, aes(x = rao_all, y = var_tot, color = 'All pixels')) +
  stat_poly_eq(data = data_all, aes(x = rao_all, y = var_tot, color = "All pixels"), label.y = "top", label.x = "left", size = 7) +
  labs(x = 'Rao', y = "Total spectral Variance", color = "")+
  theme(legend.title = element_text(size=20),
        legend.text = element_text(size = 18),
        axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        strip.text.x = element_text(size = 18))+
  facet_grid(. ~"SD = 1")

rplot_2 = ggplot(data = data_all, aes(x = Simpson, y = var_tot))+
  geom_point(aes(color = 'All pixels'))+
  stat_poly_line(data = data_all, aes(x = Simpson, y = var_tot, color = 'All pixels')) +
  stat_poly_eq(data = data_all, aes(x = Simpson, y = var_tot, color = "All pixels"), label.y = "top", label.x = "left", size = 7) +
  labs(x = 'Simpson', y = "Total spectral Variance", color = "")+
  theme(legend.title = element_text(size=20),
        legend.text = element_text(size = 18),
        axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        strip.text.x = element_text(size = 18))+
  facet_grid(. ~"SD = 1")



rplot_3 = ggplot(data = data_all, aes(x = rao_all, y = var_sp))+
  geom_point(aes(color = 'All pixels'))+
  stat_poly_line(data = data_all, aes(x = rao_all, y = var_sp, color = 'All pixels')) +
  stat_poly_eq(data = data_all, aes(x = rao_all, y = var_sp, color = "All pixels"), label.y = "top", label.x = "left", size = 7) +
  labs(x = 'Rao', y = "Inter-species spectral Variance", color = "")+
  theme(legend.title = element_text(size=20),
        legend.text = element_text(size = 18),
        axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        strip.text.x = element_text(size = 18))+
  facet_grid(. ~"SD = 1")

rplot_4 = ggplot(data = data_all, aes(x = Simpson, y = var_sp))+
  geom_point(data = data_all, aes(x = Simpson, y = var_sp, color = 'All pixels'))+
  stat_poly_line(data = data_all, aes(x = Simpson, y = var_sp, color = 'All pixels')) +
  stat_poly_eq(data = data_all, aes(x = Simpson, y = var_sp, color = "All pixels"), label.y = "top", label.x = "left", size = 7) +
  labs(x = 'Simpson', y = "Inter-species spectral Variance", color = "")+
  theme(legend.title = element_text(size=20),
        legend.text = element_text(size = 18),
        axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        strip.text.x = element_text(size = 18))+
  facet_grid(. ~"SD = 1")


FullPlot_1 = ggarrange(rplot_1, rplot_2, rplot_3, rplot_4, common.legend = TRUE, legend = "bottom",
                       nrow=2, ncol = 2, labels = c("A", "B", "C", "D"))

filename = file.path(path_figure, "inter_species_distance_rao_sd1.png")
ggsave(filename, plot = FullPlot_1, device = "png", path = NULL,
       scale = 1, width = 10, height = 8, units = "in",
       dpi = 600)


############## SD = 0.5 ##############

ResPop = readRDS(file.path(Path_reflectance, "all_populations_simulated_spectra_sd05.rds"))
# compute spectral variance corresponding to a selection of PCs
nbPop <- ResPop$PopID

#--------------- read reflectance data  ---------------
ReflZipFile <- file.path(Path_reflectance,'Spectral_Simu_VNIRCR_sd05.zip')
fid <- unz(description = ReflZipFile, filename = zip::zip_list(ReflZipFile)$filename,
           open = "rb", encoding = getOption("encoding"))
ReflectanceData <- readr::read_delim(file = fid,delim = ',')
close(fid)
for (i in 3:ncol(ReflectanceData)){
  ReflectanceData[,i] <- scale(ReflectanceData[,i],center = T,scale = T)
}


AllVars <- seq(1,(ncol(ReflectanceData)-2))

#replace Spectralpop in ResPop to match the new reflectance data
for (pop in 1:nbPop){
  ResPop$Spectralpop[[pop]] = ReflectanceData[ReflectanceData$Population == pop, -1]
}


####JE CALCULE LES SPECTRES MOYENS PAR ESPECE POUR CHAQUE TABLEAU DE REFLECTANCE (DES 220 POP)####
####JE LES UTILISE POUR CALCULER LES DISTANCES INTER-ESPECE####

# calcul des distances inter-espèce par population
inter_sp_dist_all <- list()
ReflectanceData_mean_all <- list()

Spectralpop_all = ResPop$Spectralpop

for (pop in 1:length(Spectralpop_all)){
  ReflectanceData_mean <- compute_MeanSPectraSpecies(Spectralpop_all[[pop]], as.factor(Spectralpop_all[[pop]]$SPID))
  ReflectanceData_mean <- ReflectanceData_mean$mean_spectra_species
  ReflectanceData_mean$SPID = as.factor(ReflectanceData_mean$SPID)
  rownames(ReflectanceData_mean) = ReflectanceData_mean$SPID
  ReflectanceData_mean = ReflectanceData_mean[, -1]
  ReflectanceData_mean_all[[pop]] = ReflectanceData_mean

  distance_inter_sp = dist(ReflectanceData_mean, method="euclidean")
  inter_sp_dist_all[[pop]] = distance_inter_sp
}

#which distance is the max across all distance matrix
max_dist_pop = which.max((lapply(inter_sp_dist_all, function(x) x[which.max(x)])))
max_dist = max(inter_sp_dist_all[[max_dist_pop]])

inter_sp_dist_all = lapply(inter_sp_dist_all, function(x) x/max_dist)


####JE CALCULE LA VARIANCE TOTALE ET LA VARIANCE INTER-ESPECE####


ResVar_all<- readRDS(file.path(Path_variance, "Var_DivTax_sd05.rds"))


####JE RECUPERE LES TABLES D'ABONDANCE DES 220 POPULATIONS PUIS JE CALCULE LES FREQUENCES####

abundance_table_for_cal_pop <- as.data.frame(matrix(data = 0, nrow = ResPop$PopID, ncol = length(unique(TreeID$SPID))))
colnames(abundance_table_for_cal_pop) = unique(TreeID$SPID)

for (pop in 1:ResPop$PopID){
  abundance_temp <- as.data.frame(table(Spectralpop_all[[pop]]$SPID)/20)
  for (sp in abundance_temp[, 1]){
    abundance_table_for_cal_pop[pop, sp] = abundance_temp[abundance_temp$Var1 == sp, ]$Freq
  }
}


freq <- abundance_table_for_cal_pop/100


#overwrite rao
rao <- function(freq, distance_inter_sp){
  rao = 0
  for (i in colnames(distance_inter_sp)){
    for (j in colnames(distance_inter_sp)){
      rao = rao + (freq[, i] * freq[, j] * (distance_inter_sp[i, j]^2))
    }
  }
  return(rao)
}

####CALCUL DE RAO####

rao_all <- list()
for (pop in 1:nrow(freq)){
  print(pop)
  temp_rao = rao(freq[pop, ], as.matrix(inter_sp_dist_all[[pop]]))
  print(temp_rao)
  rao_all[[pop]] = temp_rao
}

####FIGURES####

data_all = cbind(ResVar_all, as.data.frame(do.call(rbind, rao_all)))
colnames(data_all)[7] = "rao_all"

rplot_1 = ggplot(data = data_all, aes(x = rao_all, y = var_tot))+
  geom_point(aes(color = 'All pixels'))+
  stat_poly_line(data = data_all, aes(x = rao_all, y = var_tot, color = 'All pixels')) +
  stat_poly_eq(data = data_all, aes(x = rao_all, y = var_tot, color = "All pixels"), label.y = "top", label.x = "left", size = 7) +
  labs(x = 'Rao', y = "Total spectral Variance", color = "")+
  theme(legend.title = element_text(size=20),
        legend.text = element_text(size = 18),
        axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        strip.text.x = element_text(size = 18))+
  facet_grid(. ~"SD = 0.5")

rplot_2 = ggplot(data = data_all, aes(x = Simpson, y = var_tot))+
  geom_point(aes(color = 'All pixels'))+
  stat_poly_line(data = data_all, aes(x = Simpson, y = var_tot, color = 'All pixels')) +
  stat_poly_eq(data = data_all, aes(x = Simpson, y = var_tot, color = "All pixels"), label.y = "top", label.x = "left", size = 7) +
  labs(x = 'Simpson', y = "Total spectral Variance", color = "")+
  theme(legend.title = element_text(size=20),
        legend.text = element_text(size = 18),
        axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        strip.text.x = element_text(size = 18))+
  facet_grid(. ~"SD = 0.5")



rplot_3 = ggplot(data = data_all, aes(x = rao_all, y = var_sp))+
  geom_point(aes(color = 'All pixels'))+
  stat_poly_line(data = data_all, aes(x = rao_all, y = var_sp, color = 'All pixels')) +
  stat_poly_eq(data = data_all, aes(x = rao_all, y = var_sp, color = "All pixels"), label.y = "top", label.x = "left", size = 7) +
  labs(x = 'Rao', y = "Inter-species spectral Variance", color = "")+
  theme(legend.title = element_text(size=20),
        legend.text = element_text(size = 18),
        axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        strip.text.x = element_text(size = 18))+
  facet_grid(. ~"SD = 0.5")

rplot_4 = ggplot(data = data_all, aes(x = Simpson, y = var_sp))+
  geom_point(data = data_all, aes(x = Simpson, y = var_sp, color = 'All pixels'))+
  stat_poly_line(data = data_all, aes(x = Simpson, y = var_sp, color = 'All pixels')) +
  stat_poly_eq(data = data_all, aes(x = Simpson, y = var_sp, color = "All pixels"), label.y = "top", label.x = "left", size = 7) +
  labs(x = 'Simpson', y = "Inter-species spectral Variance", color = "")+
  theme(legend.title = element_text(size=20),
        legend.text = element_text(size = 18),
        axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        strip.text.x = element_text(size = 18))+
  facet_grid(. ~"SD = 0.5")


FullPlot_2 = ggarrange(rplot_1, rplot_2, rplot_3, rplot_4, common.legend = TRUE, legend = "bottom",
                       nrow=2, ncol = 2, labels = c("A", "B", "C", "D"))

filename = file.path(path_figure, "inter_species_distance_rao_sd05.png")
ggsave(filename, plot = FullPlot_2, device = "png", path = NULL,
       scale = 1, width = 10, height = 8, units = "in",
       dpi = 600)


############## SD = 0 ##############

ResPop = readRDS(file.path(Path_reflectance, "all_populations_simulated_spectra_sd0.rds"))
# compute spectral variance corresponding to a selection of PCs
nbPop <- ResPop$PopID

#--------------- read reflectance data  ---------------
ReflZipFile <- file.path(Path_reflectance,'Spectral_Simu_VNIRCR_sd0.zip')
fid <- unz(description = ReflZipFile, filename = zip::zip_list(ReflZipFile)$filename,
           open = "rb", encoding = getOption("encoding"))
ReflectanceData <- readr::read_delim(file = fid,delim = ',')
close(fid)
# for (i in 3:ncol(ReflectanceData)){
#   ReflectanceData[,i] <- scale(ReflectanceData[,i],center = T,scale = T)
# }


AllVars <- seq(1,(ncol(ReflectanceData)-2))

#replace Spectralpop in ResPop to match the new reflectance data
for (pop in 1:nbPop){
  ResPop$Spectralpop[[pop]] = ReflectanceData[ReflectanceData$Population == pop, -1]
}


####JE CALCULE LES SPECTRES MOYENS PAR ESPECE POUR CHAQUE TABLEAU DE REFLECTANCE (DES 220 POP)####
####JE LES UTILISE POUR CALCULER LES DISTANCES INTER-ESPECE####

# calcul des distances inter-espèce par population
inter_sp_dist_all <- list()
ReflectanceData_mean_all <- list()

Spectralpop_all = ResPop$Spectralpop

for (pop in 1:length(Spectralpop_all)){
  ReflectanceData_mean <- compute_MeanSPectraSpecies(Spectralpop_all[[pop]], as.factor(Spectralpop_all[[pop]]$SPID))
  ReflectanceData_mean <- ReflectanceData_mean$mean_spectra_species
  ReflectanceData_mean$SPID = as.factor(ReflectanceData_mean$SPID)
  rownames(ReflectanceData_mean) = ReflectanceData_mean$SPID
  ReflectanceData_mean = ReflectanceData_mean[, -1]
  ReflectanceData_mean_all[[pop]] = ReflectanceData_mean

  distance_inter_sp = dist(ReflectanceData_mean, method="euclidean")
  inter_sp_dist_all[[pop]] = distance_inter_sp
}

#which distance is the max across all distance matrix
max_dist_pop = which.max((lapply(inter_sp_dist_all, function(x) x[which.max(x)])))
max_dist = max(inter_sp_dist_all[[max_dist_pop]])

inter_sp_dist_all = lapply(inter_sp_dist_all, function(x) x/max_dist)


####JE CALCULE LA VARIANCE TOTALE ET LA VARIANCE INTER-ESPECE####


ResVar_all<- readRDS(file.path(Path_variance, "Var_DivTax_sd0.rds"))


####JE RECUPERE LES TABLES D'ABONDANCE DES 220 POPULATIONS PUIS JE CALCULE LES FREQUENCES####

abundance_table_for_cal_pop <- as.data.frame(matrix(data = 0, nrow = ResPop$PopID, ncol = length(unique(TreeID$SPID))))
colnames(abundance_table_for_cal_pop) = unique(TreeID$SPID)

for (pop in 1:ResPop$PopID){
  abundance_temp <- as.data.frame(table(Spectralpop_all[[pop]]$SPID)/20)
  for (sp in abundance_temp[, 1]){
    abundance_table_for_cal_pop[pop, sp] = abundance_temp[abundance_temp$Var1 == sp, ]$Freq
  }
}


freq <- abundance_table_for_cal_pop/100


#overwrite rao
rao <- function(freq, distance_inter_sp){
  rao = 0
  for (i in colnames(distance_inter_sp)){
    for (j in colnames(distance_inter_sp)){
      rao = rao + (freq[, i] * freq[, j] * (distance_inter_sp[i, j]^2))
    }
  }
  return(rao)
}

####CALCUL DE RAO####

rao_all <- list()
for (pop in 1:nrow(freq)){
  print(pop)
  temp_rao = rao(freq[pop, ], as.matrix(inter_sp_dist_all[[pop]]))
  print(temp_rao)
  rao_all[[pop]] = temp_rao
}

####FIGURES####

data_all = cbind(ResVar_all, as.data.frame(do.call(rbind, rao_all)))
colnames(data_all)[7] = "rao_all"

rplot_1 = ggplot(data = data_all, aes(x = rao_all, y = var_tot))+
  geom_point(aes(color = 'All pixels'))+
  stat_poly_line(data = data_all, aes(x = rao_all, y = var_tot, color = 'All pixels')) +
  stat_poly_eq(data = data_all, aes(x = rao_all, y = var_tot, color = "All pixels"), label.y = "top", label.x = "left", size = 7) +
  labs(x = 'Rao', y = "Total spectral Variance", color = "")+
  theme(legend.title = element_text(size=20),
        legend.text = element_text(size = 18),
        axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        strip.text.x = element_text(size = 18))+
  facet_grid(. ~"SD = 0")

rplot_2 = ggplot(data = data_all, aes(x = Simpson, y = var_tot))+
  geom_point(aes(color = 'All pixels'))+
  stat_poly_line(data = data_all, aes(x = Simpson, y = var_tot, color = 'All pixels')) +
  stat_poly_eq(data = data_all, aes(x = Simpson, y = var_tot, color = "All pixels"), label.y = "top", label.x = "left", size = 7) +
  labs(x = 'Simpson', y = "Total spectral Variance", color = "")+
  theme(legend.title = element_text(size=20),
        legend.text = element_text(size = 18),
        axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        strip.text.x = element_text(size = 18))+
  facet_grid(. ~"SD = 0")



rplot_3 = ggplot(data = data_all, aes(x = rao_all, y = var_sp))+
  geom_point(aes(color = 'All pixels'))+
  stat_poly_line(data = data_all, aes(x = rao_all, y = var_sp, color = 'All pixels')) +
  stat_poly_eq(data = data_all, aes(x = rao_all, y = var_sp, color = "All pixels"), label.y = "top", label.x = "left", size = 7) +
  labs(x = 'Rao', y = "Inter-species spectral Variance", color = "")+
  theme(legend.title = element_text(size=20),
        legend.text = element_text(size = 18),
        axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        strip.text.x = element_text(size = 18))+
  facet_grid(. ~"SD = 0")

rplot_4 = ggplot(data = data_all, aes(x = Simpson, y = var_sp))+
  geom_point(data = data_all, aes(x = Simpson, y = var_sp, color = 'All pixels'))+
  stat_poly_line(data = data_all, aes(x = Simpson, y = var_sp, color = 'All pixels')) +
  stat_poly_eq(data = data_all, aes(x = Simpson, y = var_sp, color = "All pixels"), label.y = "top", label.x = "left", size = 7) +
  labs(x = 'Simpson', y = "Inter-species spectral Variance", color = "")+
  theme(legend.title = element_text(size=20),
        legend.text = element_text(size = 18),
        axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        strip.text.x = element_text(size = 18))+
  facet_grid(. ~"SD = 0.5")


FullPlot_3 = ggarrange(rplot_1, rplot_2, rplot_3, rplot_4, common.legend = TRUE, legend = "bottom",
                       nrow=2, ncol = 2, labels = c("A", "B", "C", "D"))

filename = file.path(path_figure, "inter_species_distance_rao_sd0.png")
ggsave(filename, plot = FullPlot_3, device = "png", path = NULL,
       scale = 1, width = 10, height = 8, units = "in",
       dpi = 600)


full_plot_tot = ggarrange(FullPlot_1, FullPlot_2, FullPlot_3, nrow = 1, ncol =3, common.legend = TRUE, legend = "bottom")

filename = file.path(path_figure, "inter_species_distance_rao_fullplot.png")
ggsave(filename, plot = full_plot_tot, device = "png", path = NULL,
       scale = 1, width = 10, height = 8, units = "in",
       dpi = 600)


#just one figure to see the gradient
ggplot(ResVar_all, aes(x=Richness, y=Simpson))+
  geom_point(aes(color = "Simpson"))+
  geom_point(aes(x=Richness, y = Shannon, color ="Shannon"))+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size=20),
        legend.text  = element_text(size=20),
        legend.title = element_text(size=20),
        aspect.ratio = 1)
