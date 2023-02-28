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
library(ggdist)
library(ggpubr)



#--------------- define path for data and results ---------------
Path_Inventory <- '../../03_RESULTS/02_Produce_FilteredDatafiles'
Path_Spectra <- '../../03_RESULTS/03_Produce_POP/Reflectance'
Path_figures <- '../../03_RESULTS/FIGURES/verif_spectra'
dir.create(path = Path_figures, showWarnings = F, recursive = T)

#--------------- read reflectance data  ---------------
ReflZipFile <- file.path(Path_Inventory,paste('ReflTrees_VNIR.zip',sep = ''))
fid <- unz(description = ReflZipFile, filename = zip::zip_list(ReflZipFile)$filename,
           open = "rb", encoding = getOption("encoding"))
ReflectanceData <- readr::read_delim(file = fid,delim = ',')
close(fid)

############## SD = 1 ##############

ResPop = readRDS("../../03_RESULTS/03_Produce_POP/Reflectance/all_populations_simulated_spectra_sd1.rds")
Simulated_reflectance_data = bind_rows(ResPop$Spectralpop, .id = "column_label")
Simulated_reflectance_data = Simulated_reflectance_data[, -1]


make_plot <- function(sp_reflectance, facet){

  rplot = ggplot() +
    stat_lineribbon(data = sp_reflectance,
                    aes(x = Wavelength, y = value*100), alpha = 0.5, lwd=1) +
    guides(fill = guide_legend(reverse=TRUE))+
    scale_x_continuous(limits = c(400,900)) +
    scale_y_continuous(limits = c(0,100))+
    labs(y='Reflectance (%)') +
    facet_grid(facet$arg~unique(sp_reflectance$Species))

  return(rplot)

}

plot_list_simulated_VNIR = list()
plot_list_VNIR = list()

for (sp in unique(Simulated_reflectance_data$SPID)){
  print(sp)
  sp_reflectance_simulated = Simulated_reflectance_data[Simulated_reflectance_data$SPID == sp, ]
  sp_reflectance_simulated = as.data.frame(t(sp_reflectance_simulated[, -1]))
  sp_reflectance_simulated = cbind("Wavelength" = as.numeric(rownames(sp_reflectance_simulated)), sp_reflectance_simulated)

  sp_reflectance_simulated = sp_reflectance_simulated %>%
    tidyr::gather(Sample, value, -Wavelength)
  sp_reflectance_simulated = cbind(sp_reflectance_simulated, "Species" = sp)

  plot_list_simulated_VNIR[[sp]] = make_plot(sp_reflectance_simulated, facet = as.data.frame(cbind("arg" = "SIMULATED VNIR")))

  sp_reflectance = ReflectanceData[ReflectanceData$SPID == sp, ]
  sp_reflectance = as.data.frame(t(sp_reflectance[, -c(1:3)]))
  sp_reflectance = cbind("Wavelength" = as.numeric(rownames(sp_reflectance)), sp_reflectance)

  sp_reflectance = sp_reflectance %>%
    tidyr::gather(Sample, value, -Wavelength)
  sp_reflectance = cbind(sp_reflectance, "Species" = sp)

  plot_list_VNIR[[sp]] =  make_plot(sp_reflectance, facet = as.data.frame(cbind("arg" = "VNIR")))
}


#open simulated VNIR+CR

ReflZipFile <- file.path(Path_Spectra,'Spectral_Simu_VNIRCR_sd1.zip')
fid <- unz(description = ReflZipFile, filename = zip::zip_list(ReflZipFile)$filename,
           open = "rb", encoding = getOption("encoding"))
ReflectanceData_CR_simulated <- readr::read_delim(file = fid,delim = ',')
close(fid)

#apply CR to real reflectance data
Spectral <- list('Wavelength'= as.numeric(colnames(ReflectanceData[, -c(1:3)])))
Refl_CR <- apply_continuum_removal(Spectral_Data=as.matrix(ReflectanceData[, -c(1:3)]),
                                   Spectral = Spectral, nbCPU = 6)

Refl_CR <- data.frame(Refl_CR)
names(Refl_CR) <- as.character(Spectral$Wavelength[-c(1, length(Spectral$Wavelength))])
Refl_CR = cbind(ReflectanceData[, c(1:3)], Refl_CR)


plot_list_simulated_VNIRCR = list()
plot_list_VNIRCR = list()

for (sp in unique(Simulated_reflectance_data$SPID)){
  print(sp)
  sp_reflectance_simulated = ReflectanceData_CR_simulated[ReflectanceData_CR_simulated$SPID == sp, ]
  sp_reflectance_simulated = as.data.frame(t(sp_reflectance_simulated[, -c(1,2)]))
  sp_reflectance_simulated = cbind("Wavelength" = as.numeric(rownames(sp_reflectance_simulated)), sp_reflectance_simulated)

  sp_reflectance_simulated = sp_reflectance_simulated %>%
    tidyr::gather(Sample, value, -Wavelength)
  sp_reflectance_simulated = cbind(sp_reflectance_simulated, "Species" = sp)

  plot_list_simulated_VNIRCR[[sp]] = make_plot(sp_reflectance_simulated, facet = as.data.frame(cbind("arg" = "SIMULATED VNIR+CR")))

  sp_reflectance = Refl_CR[Refl_CR$SPID == sp, ]
  sp_reflectance = as.data.frame(t(sp_reflectance[, -c(1:3)]))
  sp_reflectance = cbind("Wavelength" = as.numeric(rownames(sp_reflectance)), sp_reflectance)

  sp_reflectance = sp_reflectance %>%
    tidyr::gather(Sample, value, -Wavelength)
  sp_reflectance = cbind(sp_reflectance, "Species" = sp)

  plot_list_VNIRCR[[sp]] =  make_plot(sp_reflectance, facet = as.data.frame(cbind("arg" = "VNIR+CR")))
}


#make_final_plot
for (sp in unique(Simulated_reflectance_data$SPID)){
  filename = file.path(Path_figures, paste(sp, ".png", sep=""))

  rplot = ggarrange(plot_list_simulated_VNIR[[sp]], plot_list_simulated_VNIRCR[[sp]],
                    plot_list_VNIR[[sp]], plot_list_VNIRCR[[sp]], common.legend = TRUE, legend = "bottom",
                    ncol = 2, nrow = 2)

  ggsave(filename, plot = rplot, device = "png", path = NULL,
         scale = 1, width = 10, height = 8, units = "in",
         dpi = 600)

  print(paste(filename, "is saved"))
}

