# ==============================================================================
# This script allows to generate artificial population based on inventory data.
# It also allows to compute diversity for the population and spectral diversity
# of each population.
# To compare spectral diversity and taxonomic diversity we generate graphics and
# use linear regression to to estimate the correlation.
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
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
################################################################################

#--------------- import libraries ---------------
source('../00_Libraries/Lib_populations.R')
source('../00_Libraries/Lib_diversity.R')
source('../00_Libraries/Lib_ImageProcess.R')

library(zip)
library(dplyr)
library(ade4)
library(readr)
library(tools)
library(progress)
library(ggplot2)
library(cowplot)
library(dplyr)
library(biodivMapR)


#--------------- define path for data and results ---------------
Path_Inventory <- '../../03_RESULTS/02_Produce_FilteredDatafiles'
Path_Populations <- '../../03_RESULTS/03_Produce_POP'
Path_Results <- '../../03_RESULTS/03_Produce_POP/Reflectance'
dir.create(path = Path_Results, showWarnings = F, recursive = T)

#--------------- read reflectance data  ---------------
ReflZipFile <- file.path(Path_Inventory,paste('ReflTrees_VNIR.zip',sep = ''))
fid <- unz(description = ReflZipFile, filename = zip::zip_list(ReflZipFile)$filename,
           open = "rb", encoding = getOption("encoding"))
ReflectanceData <- readr::read_delim(file = fid,delim = ',')
close(fid)

BasePopulations_min = readRDS(file.path(Path_Populations, 'BasePopulation_min.rds'))
BasePopulations_max = readRDS(file.path(Path_Populations, 'BasePopulation_max.rds'))
Populations_mintomax = readRDS(file.path(Path_Populations, 'Populations_mintomax.rds'))
Populations_maxtomin = readRDS(file.path(Path_Populations, 'Populations_maxtomin.rds'))

PopID <- 0
Spectralpop <- list()

############## SD = 1 ##############

sd_param = 1
#sd = 1 initial intra var sp
#sd = 0.5 half intra var sp
#sd = 0 no intra var sp

# minimum population
nbPop <- nrow(BasePopulations_min$pop)
ResPop <- GetSpectra_FromPops(Populations = BasePopulations_min$pop, ReflectanceData = ReflectanceData,
                              PopID=PopID, Spectralpop=Spectralpop, nbPop=nbPop, sd_param = sd_param)

# dilution of minimum population
pcDilution <- names(Populations_mintomax)
for (dilute in pcDilution){
  print(dilute)
  nbPop <- length(Populations_mintomax[[dilute]]$div$Simpson)
  ResPop <- GetSpectra_FromPops(Populations = Populations_mintomax[[dilute]]$pop, ReflectanceData = ReflectanceData,
                                PopID=ResPop$PopID, Spectralpop=ResPop$Spectralpop, nbPop=nbPop,  sd_param = sd_param)
}

# dilution of maximum population
pcDilution <- names(Populations_maxtomin)
for (dilute in pcDilution){
  print(dilute)
  nbPop <- length(Populations_maxtomin[[dilute]]$div$Simpson)
  ResPop <- GetSpectra_FromPops(Populations = Populations_maxtomin[[dilute]]$pop, ReflectanceData = ReflectanceData,
                                PopID=ResPop$PopID, Spectralpop=ResPop$Spectralpop, nbPop=nbPop,  sd_param = sd_param)
}

# maximum population
nbPop <- nrow(BasePopulations_max$pop)
ResPop <- GetSpectra_FromPops(Populations = BasePopulations_max$pop, ReflectanceData = ReflectanceData,
                              PopID=ResPop$PopID, Spectralpop=ResPop$Spectralpop, nbPop=nbPop,  sd_param = sd_param)

saveRDS(ResPop, file = file.path(Path_Results, "all_populations_simulated_spectra_sd1.rds"))



#apply CR

#ResPop = readRDS("../../03_RESULTS/03_Produce_POP/Reflectance/all_populations_simulated_spectra_sd1.rds")
Simulated_reflectance_data = bind_rows(ResPop$Spectralpop, .id = "column_label")
Spectral = list("Wavelength" = as.numeric(colnames(Simulated_reflectance_data[, -c(1,2)])))

Simu_Refl_CR_VNIR <- apply_continuum_removal(Spectral_Data= as.matrix(Simulated_reflectance_data[, -c(1,2)]),
                                        Spectral = Spectral, nbCPU = 6)

Simu_Refl_CR <- data.frame(Simu_Refl_CR_VNIR)
Simu_Refl_CR2 <- format(Simu_Refl_CR,digits=5,nsmall = 5)
names(Simu_Refl_CR2) <- as.character(Spectral$Wavelength[-c(1, length(Spectral$Wavelength))])
Simu_Refl_CR2 = cbind(Simulated_reflectance_data[, c(1,2)], Simu_Refl_CR2)
colnames(Simu_Refl_CR2)[1] = "Population"

fileName <- file.path(Path_Results,'Spectral_Simu_VNIRCR_sd1.csv')
write_csv(x = Simu_Refl_CR2,file = fileName,
          col_names = TRUE,append = FALSE)
ZipName <- paste(file_path_sans_ext(fileName),'.zip',sep = '')
zip(zipfile = ZipName,
    files =  fileName,mode = 'cherry-pick')
file.remove(fileName)


############## SD = 0.5 ##############

sd_param = 0.5
#sd = 1 initial intra var sp
#sd = 0.5 half intra var sp
#sd = 0 no intra var sp

# minimum population
nbPop <- nrow(BasePopulations_min$pop)
ResPop <- GetSpectra_FromPops(Populations = BasePopulations_min$pop, ReflectanceData = ReflectanceData,
                              PopID=PopID, Spectralpop=Spectralpop, nbPop=nbPop, sd_param = sd_param)

# dilution of minimum population
pcDilution <- names(Populations_mintomax)
for (dilute in pcDilution){
  print(dilute)
  nbPop <- length(Populations_mintomax[[dilute]]$div$Simpson)
  ResPop <- GetSpectra_FromPops(Populations = Populations_mintomax[[dilute]]$pop, ReflectanceData = ReflectanceData,
                                PopID=ResPop$PopID, Spectralpop=ResPop$Spectralpop, nbPop=nbPop,  sd_param = sd_param)
}

# dilution of maximum population
pcDilution <- names(Populations_maxtomin)
for (dilute in pcDilution){
  print(dilute)
  nbPop <- length(Populations_maxtomin[[dilute]]$div$Simpson)
  ResPop <- GetSpectra_FromPops(Populations = Populations_maxtomin[[dilute]]$pop, ReflectanceData = ReflectanceData,
                                PopID=ResPop$PopID, Spectralpop=ResPop$Spectralpop, nbPop=nbPop,  sd_param = sd_param)
}

# maximum population
nbPop <- nrow(BasePopulations_max$pop)
ResPop <- GetSpectra_FromPops(Populations = BasePopulations_max$pop, ReflectanceData = ReflectanceData,
                              PopID=ResPop$PopID, Spectralpop=ResPop$Spectralpop, nbPop=nbPop,  sd_param = sd_param)

saveRDS(ResPop, file = file.path(Path_Results, "all_populations_simulated_spectra_sd05.rds"))



#apply CR

#ResPop = readRDS("../../03_RESULTS/03_Produce_POP/Reflectance/all_populations_simulated_spectra_sd05.rds")
Simulated_reflectance_data = bind_rows(ResPop$Spectralpop, .id = "column_label")
Spectral = list("Wavelength" = as.numeric(colnames(Simulated_reflectance_data[, -c(1,2)])))

Simu_Refl_CR_VNIR <- apply_continuum_removal(Spectral_Data= as.matrix(Simulated_reflectance_data[, -c(1,2)]),
                                             Spectral = Spectral, nbCPU = 6)

Simu_Refl_CR <- data.frame(Simu_Refl_CR_VNIR)
Simu_Refl_CR2 <- format(Simu_Refl_CR,digits=5,nsmall = 5)
names(Simu_Refl_CR2) <- as.character(Spectral$Wavelength[-c(1, length(Spectral$Wavelength))])
Simu_Refl_CR2 = cbind(Simulated_reflectance_data[, c(1,2)], Simu_Refl_CR2)
colnames(Simu_Refl_CR2)[1] = "Population"

fileName <- file.path(Path_Results,'Spectral_Simu_VNIRCR_sd05.csv')
write_csv(x = Simu_Refl_CR2,file = fileName,
          col_names = TRUE,append = FALSE)
ZipName <- paste(file_path_sans_ext(fileName),'.zip',sep = '')
zip(zipfile = ZipName,
    files =  fileName,mode = 'cherry-pick')
file.remove(fileName)


############## SD = 0 ##############

sd_param = 0
#sd = 1 initial intra var sp
#sd = 0.5 half intra var sp
#sd = 0 no intra var sp

# minimum population
nbPop <- nrow(BasePopulations_min$pop)
ResPop <- GetSpectra_FromPops(Populations = BasePopulations_min$pop, ReflectanceData = ReflectanceData,
                              PopID=PopID, Spectralpop=Spectralpop, nbPop=nbPop, sd_param = sd_param)

# dilution of minimum population
pcDilution <- names(Populations_mintomax)
for (dilute in pcDilution){
  print(dilute)
  nbPop <- length(Populations_mintomax[[dilute]]$div$Simpson)
  ResPop <- GetSpectra_FromPops(Populations = Populations_mintomax[[dilute]]$pop, ReflectanceData = ReflectanceData,
                                PopID=ResPop$PopID, Spectralpop=ResPop$Spectralpop, nbPop=nbPop,  sd_param = sd_param)
}

# dilution of maximum population
pcDilution <- names(Populations_maxtomin)
for (dilute in pcDilution){
  print(dilute)
  nbPop <- length(Populations_maxtomin[[dilute]]$div$Simpson)
  ResPop <- GetSpectra_FromPops(Populations = Populations_maxtomin[[dilute]]$pop, ReflectanceData = ReflectanceData,
                                PopID=ResPop$PopID, Spectralpop=ResPop$Spectralpop, nbPop=nbPop,  sd_param = sd_param)
}

# maximum population
nbPop <- nrow(BasePopulations_max$pop)
ResPop <- GetSpectra_FromPops(Populations = BasePopulations_max$pop, ReflectanceData = ReflectanceData,
                              PopID=ResPop$PopID, Spectralpop=ResPop$Spectralpop, nbPop=nbPop,  sd_param = sd_param)

saveRDS(ResPop, file = file.path(Path_Results, "all_populations_simulated_spectra_sd0.rds"))



#apply CR

#ResPop = readRDS("../../03_RESULTS/03_Produce_POP/Reflectance/all_populations_simulated_spectra_sd0.rds")
Simulated_reflectance_data = bind_rows(ResPop$Spectralpop, .id = "column_label")
Spectral = list("Wavelength" = as.numeric(colnames(Simulated_reflectance_data[, -c(1,2)])))

Simu_Refl_CR_VNIR <- apply_continuum_removal(Spectral_Data= as.matrix(Simulated_reflectance_data[, -c(1,2)]),
                                             Spectral = Spectral, nbCPU = 6)

Simu_Refl_CR <- data.frame(Simu_Refl_CR_VNIR)
Simu_Refl_CR2 <- format(Simu_Refl_CR,digits=5,nsmall = 5)
names(Simu_Refl_CR2) <- as.character(Spectral$Wavelength[-c(1, length(Spectral$Wavelength))])
Simu_Refl_CR2 = cbind(Simulated_reflectance_data[, c(1,2)], Simu_Refl_CR2)
colnames(Simu_Refl_CR2)[1] = "Population"

fileName <- file.path(Path_Results,'Spectral_Simu_VNIRCR_sd0.csv')
write_csv(x = Simu_Refl_CR2,file = fileName,
          col_names = TRUE,append = FALSE)
ZipName <- paste(file_path_sans_ext(fileName),'.zip',sep = '')
zip(zipfile = ZipName,
    files =  fileName,mode = 'cherry-pick')
file.remove(fileName)
