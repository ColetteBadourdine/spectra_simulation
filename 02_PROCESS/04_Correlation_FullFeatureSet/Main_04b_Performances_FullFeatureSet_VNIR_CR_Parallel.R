# ==============================================================================
# Description here
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

#--------------- import libraries ---------------
source('../00_Libraries/Lib_diversity.R')
source('../00_Libraries/Lib_populations.R')
source('../00_Libraries/Lib_ImageProcess.R')
source('../00_Libraries/Lib_variance.R')
library(zip)
library(readr)
library(dplyr)
library(tools)
library(doFuture)
library(ggpubr)
library(ggpmisc)

nbWorkers <- 36

#--------------- define path for data and results ---------------
Path_Populations <- '../../03_RESULTS/03_Produce_POP'
Path_reflectance <- '../../03_RESULTS/03_Produce_POP/Reflectance'
Path_Results <- file.path('../../03_RESULTS/04bis_Correlation_FullFeature_SIMU_VNIRCR')
dir.create(path = Path_Results, showWarnings = F, recursive = T)

Path_Figures <- file.path('../../03_RESULTS/FIGURES')
dir.create(path = Path_Figures, showWarnings = F, recursive = T)

BasePopulations_min = readRDS(file.path(Path_Populations,'BasePopulation_min.rds'))
BasePopulations_max = readRDS(file.path(Path_Populations,'BasePopulation_max.rds'))
Populations_mintomax = readRDS(file.path(Path_Populations,'Populations_mintomax.rds'))
Populations_maxtomin = readRDS(file.path(Path_Populations,'Populations_maxtomin.rds'))

# compute diversity metrics corresponding to
All_Populations_diversite = compute_diversity_index(BasePopulations_min, BasePopulations_max, Populations_mintomax, Populations_maxtomin,
                                                    name_column = c("BasePopulations_min", "BasePopulations_max", "Populations_mintomax",
                                                                    "Populations_maxtomin"))

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

ResVar<- lapply(X = ResPop$Spectralpop,FUN = spectral_variance_subset,
                vars = AllVars, stand='None')

ResVar <- as.data.frame(do.call(rbind,(lapply(ResVar,unlist)))[,1:2])
ResVar <- cbind(All_Populations_diversite, ResVar)

summary(lm(ResVar$var_tot~ResVar$Simpson))
summary(lm(ResVar$var_sp~ResVar$Simpson))


rplot1 = ggplot(data = ResVar, aes(x = Simpson, y = var_tot))+
  geom_point(aes(color = "Total spectral variance"))+
  stat_poly_line(aes(color = "Total spectral variance")) +
  stat_poly_eq(label.y = "top", label.x = "right", aes(color = "Total spectral variance"), size = 7) +
  geom_point(data = ResVar, aes(x = Simpson, y = var_sp, color = 'Inter-species spectral variance'))+
  stat_poly_line(data = ResVar, aes(x = Simpson, y = var_sp, color = 'Inter-species spectral variance')) +
  stat_poly_eq(data = ResVar, aes(x = Simpson, y = var_sp, color = "Inter-species spectral variance"), label.y = "bottom", label.x = "right", size = 7) +
  labs(x = 'Simpson', y = "Spectral Variance", color = "")+
  theme(legend.title = element_text(size=20),
        legend.text = element_text(size = 18),
        axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        strip.text.x = element_text(size = 18))+
  ylim(c(0, 250))+
  facet_grid(. ~"All spectral features - SD = 1")


saveRDS(ResVar, file = file.path(Path_Results, "Var_DivTax_sd1.rds"))

filename = file.path(Path_Figures, "Var_vs_Simpson_sd1.png")
ggsave(filename, plot = rplot1, device = "png", path = NULL,
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

#replace Spectralpop in ResPop to match the new reflectance data
for (pop in 1:nbPop){
  ResPop$Spectralpop[[pop]] = ReflectanceData[ReflectanceData$Population == pop, -1]
}

AllVars <- seq(1,(ncol(ReflectanceData)-2))

ResVar<- lapply(X = ResPop$Spectralpop,FUN = spectral_variance_subset,
                vars = AllVars, stand='None')

ResVar <- as.data.frame(do.call(rbind,(lapply(ResVar,unlist)))[,1:2])
ResVar <- cbind(All_Populations_diversite, ResVar)

summary(lm(ResVar$var_tot~ResVar$Simpson))
summary(lm(ResVar$var_sp~ResVar$Simpson))


rplot2 = ggplot(data = ResVar, aes(x = Simpson, y = var_tot))+
  geom_point(aes(color = "Total spectral variance"))+
  stat_poly_line(aes(color = "Total spectral variance")) +
  stat_poly_eq(label.y = "top", label.x = "right", aes(color = "Total spectral variance"), size = 7) +
  geom_point(data = ResVar, aes(x = Simpson, y = var_sp, color = 'Inter-species spectral variance'))+
  stat_poly_line(data = ResVar, aes(x = Simpson, y = var_sp, color = 'Inter-species spectral variance')) +
  stat_poly_eq(data = ResVar, aes(x = Simpson, y = var_sp, color = "Inter-species spectral variance"), label.y = "bottom", label.x = "right", size = 7) +
  labs(x = 'Simpson', y = "Spectral Variance", color = "")+
  theme(legend.title = element_text(size=20),
        legend.text = element_text(size = 18),
        axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        strip.text.x = element_text(size = 18))+
  facet_grid(. ~"All spectral features - SD = 0.5")


saveRDS(ResVar, file = file.path(Path_Results, "Var_DivTax_sd05.rds"))

filename = file.path(Path_Figures, "Var_vs_Simpson_sd05.png")
ggsave(filename, plot = rplot2, device = "png", path = NULL,
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

#replace Spectralpop in ResPop to match the new reflectance data
for (pop in 1:nbPop){
  ResPop$Spectralpop[[pop]] = as.data.frame(ReflectanceData[ReflectanceData$Population == pop, -1])
}


AllVars <- seq(1,(ncol(ReflectanceData)-2))

ResVar<- lapply(X = ResPop$Spectralpop,FUN = spectral_variance_subset,
                vars = AllVars, stand='None')

ResVar <- as.data.frame(do.call(rbind,(lapply(ResVar,unlist)))[,1:2])
ResVar <- cbind(All_Populations_diversite, ResVar)

summary(lm(ResVar$var_tot~ResVar$Simpson))
summary(lm(ResVar$var_sp~ResVar$Simpson))


rplot3 = ggplot(data = ResVar, aes(x = Simpson, y = var_tot))+
  geom_point(aes(color = "Total spectral variance"))+
  stat_poly_line(aes(color = "Total spectral variance")) +
  stat_poly_eq(label.y = "top", label.x = "right", aes(color = "Total spectral variance"), size = 7) +
  geom_point(data = ResVar, aes(x = Simpson, y = var_sp, color = 'Inter-species spectral variance'))+
  stat_poly_line(data = ResVar, aes(x = Simpson, y = var_sp, color = 'Inter-species spectral variance')) +
  stat_poly_eq(data = ResVar, aes(x = Simpson, y = var_sp, color = "Inter-species spectral variance"), label.y = "bottom", label.x = "right", size = 7) +
  labs(x = 'Simpson', y = "Spectral Variance", color = "")+
  theme(legend.title = element_text(size=20),
        legend.text = element_text(size = 18),
        axis.title = element_text(size=20),
        axis.text = element_text(size=18),
        strip.text.x = element_text(size = 18))+
  facet_grid(. ~"All spectral features - SD = 0")


saveRDS(ResVar, file = file.path(Path_Results, "Var_DivTax_sd0.rds"))

filename = file.path(Path_Figures, "Var_vs_Simpson_sd0.png")
ggsave(filename, plot = rplot3, device = "png", path = NULL,
       scale = 1, width = 10, height = 8, units = "in",
       dpi = 600)

legend = get_legend(rplot1)

fullplot = ggarrange(rplot1+theme(legend.position = "none"),
                     rplot2+theme(legend.position = "none"),
                     rplot3+theme(legend.position = "none"), legend, ncol = 2, nrow = 2)

filename = file.path(Path_Figures, "Var_vs_Simpson_fullplot.png")
ggsave(filename, plot = fullplot, device = "png", path = NULL,
       scale = 1, width = 10, height = 8, units = "in",
       dpi = 600)
