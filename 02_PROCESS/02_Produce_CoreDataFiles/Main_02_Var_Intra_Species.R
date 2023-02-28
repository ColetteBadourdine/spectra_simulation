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
source('../00_Libraries/Lib_variance.R')

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

#first need to apply CR + center-reduction
Spectral = list("Wavelength" = as.numeric(colnames(ReflectanceData[, -c(1:3)])))
ReflectanceData_CR = apply_continuum_removal(Spectral_Data= as.matrix(ReflectanceData[, -c(1:3)]),
                                             Spectral = Spectral, nbCPU = 6)

ReflectanceData_CR = as.data.frame(ReflectanceData_CR)
ReflectanceData_CR <- round(ReflectanceData_CR, 5)
names(ReflectanceData_CR) <- as.character(Spectral$Wavelength[-c(1, length(Spectral$Wavelength))])
ReflectanceData_CR = cbind(ReflectanceData[, c(1:3)], ReflectanceData_CR)

for (i in 3:ncol(ReflectanceData_CR)){
  ReflectanceData_CR[,i] <- scale(ReflectanceData_CR[,i],center = T,scale = T)
}

VarIntraSp = as.data.frame(matrix(nrow = length(unique(ReflectanceData$SPID)), ncol = 3))
colnames(VarIntraSp) = c("SPID", "var_intra", "NbPixels")
VarIntraSp$SPID = unique(ReflectanceData$SPID)

for (sp in unique(ReflectanceData$SPID)){
  VarIntraSp[VarIntraSp$SPID == sp, 2] = sum(as.vector(apply(ReflectanceData_CR[ReflectanceData_CR$SPID == sp, -c(1:3)],2,var)))
  VarIntraSp[VarIntraSp$SPID == sp, 3] = nrow(ReflectanceData_CR[ReflectanceData_CR$SPID == sp, ])
}

VarIntraSp$var_intra = round(VarIntraSp$var_intra, 4)

summary(VarIntraSp$var_intra)

### make histogram

rplot_histo = ggplot(data= VarIntraSp, aes(x=reorder(SPID, -var_intra), y=var_intra)) +
  geom_bar(stat="identity", color="black", fill="white")+
  geom_hline(yintercept=median(VarIntraSp$var_intra), color = "red")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Species")+
  ylab("Intra-specific variance")

filename = "../../03_RESULTS/FIGURES/intra_var_sp_bd_data.png"
ggsave(filename, plot = rplot_histo, device = "png", path = NULL,
       scale = 1, width = 10, height = 8, units = "in",
       dpi = 600)

rplot_varintra_vs_pixels = ggplot(VarIntraSp, aes(x=NbPixels, y = var_intra))+
  geom_point()

filename = "../../03_RESULTS/FIGURES/intra_var_vs_nbpixel_bd_data.png"
ggsave(filename, plot = rplot_varintra_vs_pixels, device = "png", path = NULL,
       scale = 1, width = 10, height = 8, units = "in",
       dpi = 600)
