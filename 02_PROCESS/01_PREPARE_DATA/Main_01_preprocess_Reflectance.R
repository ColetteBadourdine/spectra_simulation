# ==============================================================================
# This script aims at preprocessing reflectance data extracted from VSWIR imaging
# spectroscopy data
# ==============================================================================
# PROGRAMMERS:
# Colette BADOURDINE <colette.badourdine@cirad.fr>
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Copyright 2021/06 Colette BADOURDINE
# ==============================================================================
# This Library includes functions to pre-treat hyperspectral data VNIR and VSWIR
# These functions are based on the functions of the BiodivMapR package by Jb Féret
# and F De Boissieu
# ==============================================================================

################################################################################
# Always start a script with a clean environment
rm(list=ls(all=TRUE));gc()
# define working directory as the directory where the script is located
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
################################################################################
source('../00_Libraries/Lib_ImageProcess.R')

#--------------- import libraries ---------------
library(readr)
library(biodivMapR)
library(zip)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(stringr)
library(tools)
library(snow)
library(matrixStats)
library(future)
library(future.apply)

#--------------- define path for data and results ---------------
Path_Data <- '../../01_DATA'
Path_Results <- '../../03_RESULTS'
dir.create(path = Path_Results, showWarnings = F, recursive = T)
# number of CPU available for Continuum removal
nbCPU <- 6

#--------------- read and prepare data ---------------
# read reflectance data directly from zipfile
ReflZipFile <- file.path(Path_Data,'MosaPAR_SWIR_20160919_coefJB_extract3_entier.zip')
fid <- unz(description = ReflZipFile, filename = zip_list(ReflZipFile)$filename, open = "rb", encoding = getOption("encoding"))
InventoryData <- read_delim(file = fid,delim = ';')
close(fid)

# filter data based on OK_2016
SelPix <- which(InventoryData$OK_2016>=2)
InventoryData <- InventoryData[SelPix,]

# get data structured properly
SPID <- InventoryData$SPID                  # species name
PlotID <- InventoryData$NParcelle           # Plot ID
CrownID <- InventoryData$ID                 # CrownID
WL <- as.numeric(gsub(names(InventoryData)[-seq(1:5)],pattern = 'X',replacement = ''))
Refl <- as.matrix(InventoryData[,-seq(1:5)])

#--------------- filter pixels to discard shaded/cloudy/non-vegetated pixels ---------------
# - NIR thresholding: discard pixels with low reflectance in the NIR corresponding to shade
# - Blue thresholding: discard pixels with high reflectance in the Blue corresponding to clouds
# - NDVI thresholding: discard pixels with low NDVI corresponding to non-vegetated pixels
Blue <- 480
Red <- 680
NIR <- 835
Spectral_Bands <- c(Blue, Red, NIR)
Image_Bands <- get_image_bands(Spectral_Bands, WL)
NDVI <- c((Refl[,Image_Bands$ImBand[3]] - Refl[,Image_Bands$ImBand[2]]) /
            (Refl[,Image_Bands$ImBand[3]] + Refl[,Image_Bands$ImBand[2]]))

NDVI_Thresh <- 0.5
Blue_Thresh <- 500
NIR_Thresh <- 2000
# select pixels corresponding to filter criteria
SelPixels <- which(NDVI > NDVI_Thresh & Refl[,Image_Bands$ImBand[1]] < Blue_Thresh & Refl[,Image_Bands$ImBand[3]] > NIR_Thresh)
# exclude samples based on this filter
Refl_VSWIR <- Refl[SelPixels,]
colnames(Refl_VSWIR) <- NULL

SPID <- SPID[SelPixels]
PlotID <- PlotID[SelPixels]
CrownID <- CrownID[SelPixels]

#--------------- filter spectral data to discard noisy spectral domains ---------------
# 1- plot reflectance for a set of random samples
nbSamples <- 1000
Sel <- sample(x = seq(1,length(SelPixels)),nbSamples)
Rsel0 <- t(0.01*Refl_VSWIR[Sel,])
Rsel <- data.frame(cbind(WL,Rsel0))
names(Rsel) <- c('wavelength',paste('R#',seq(1,nbSamples),sep=''))
ResRefl <- Rsel %>%
  tidyr::gather(Sample, value, -wavelength)
# create dataframe to prepare for figure
ResRefl <- data.frame(ResRefl,'RT'='R')

rplot <- ggplot(ResRefl, aes() ) +
  geom_line(data = ResRefl,aes(x=wavelength, y=value, group=Sample, colour=Sample) ,size=0.5) +
  labs(x="Wavelength (nm)", y="Reflectance (%)") +
  theme(legend.title =element_blank(),legend.position ="none") +
  ylim(0,100) + xlim(400,2550)

filename <- file.path(Path_Results,'01_Raw_Reflectance.png')
ggsave(filename, plot = rplot, device = "png", path = NULL,
       scale = 1, width = 13, height = 10, units = "cm",
       dpi = 600)

#--------------- compute derivative to better identify noisy domains ---------------
# 1- plot derivative of reflectance for the same set of random samples
Rderiv0 <- compute_derivative_spectra(reflectance =Rsel0,
                                     spectralbands =WL)
Rderiv <- data.frame(cbind(WL[-length(WL)],Rderiv0))
names(Rderiv) <- c('wavelength',paste('R#',seq(1,nbSamples),sep=''))
ResRefl <- Rderiv %>%
  tidyr::gather(Sample, value, -wavelength)
# create dataframe to prepare for figure
ResRefl <- data.frame(ResRefl,'RT'='R')

rplot <- ggplot(ResRefl, aes() ) +
  geom_point(data = ResRefl,aes(x=wavelength, y=value, group=Sample, colour=Sample) ,size=0.5) +
  labs(x="Wavelength (nm)", y="delta(Reflectance) (%)") +
  theme(legend.title =element_blank(),legend.position ="none") +
  ylim(-2,2)

filename <- file.path(Path_Results,'02_Derivative_Reflectance.png')
ggsave(filename, plot = rplot, device = "png", path = NULL,
       scale = 1, width = 16, height = 8, units = "cm",
       dpi = 600)

rplot <- ggplot(ResRefl, aes() ) +
  geom_point(data = ResRefl,aes(x=wavelength, y=value, group=Sample, colour=Sample) ,size=0.5) +
  labs(x="Wavelength (nm)", y="delta(Reflectance) (%)") +
  theme(legend.title =element_blank(),legend.position ="none") +
  ylim(-2,2) + xlim(400,1100)

filename <- file.path(Path_Results,'02_Derivative_Reflectance_VNIR.png')
ggsave(filename, plot = rplot, device = "png", path = NULL,
       scale = 1, width = 16, height = 8, units = "cm",
       dpi = 600)

rplot <- ggplot(ResRefl, aes() ) +
  geom_point(data = ResRefl,aes(x=wavelength, y=value, group=Sample, colour=Sample) ,size=0.5) +
  labs(x="Wavelength (nm)", y="delta(Reflectance) (%)") +
  theme(legend.title =element_blank(),legend.position ="none") +
  ylim(-2,2) + xlim(900,2500)

filename <- file.path(Path_Results,'02_Derivative_Reflectance_SWIR.png')
ggsave(filename, plot = rplot, device = "png", path = NULL,
       scale = 1, width = 16, height = 8, units = "cm",
       dpi = 600)


# 2- identify spectral bands for which derivative is >X% beyond 800 nm
# and consider them as noisy
abs_Rderiv0 <- apply(abs(Rderiv0),MARGIN = 1,FUN = max)
abs_Rderiv <- data.frame(cbind(WL[-length(WL)],abs_Rderiv0))
names(abs_Rderiv) <- c('wavelength','MaxDerivative')
ResRefl <- abs_Rderiv %>%
  tidyr::gather(Sample, value, -wavelength)
# create dataframe to prepare for figure
ResRefl <- data.frame(ResRefl,'RT'='R')

rplot <- ggplot(ResRefl, aes() ) +
  geom_point(data = ResRefl,aes(x=wavelength, y=value, group=Sample, colour=Sample) ,size=0.5) +
  labs(x="Wavelength (nm)", y="delta(Reflectance) (%)") +
  theme(legend.title =element_blank(),legend.position ="none") +
  ylim(0,2) + xlim(400,2500)

filename <- file.path(Path_Results,'03_MaxAbs_Derivative_Reflectance.png')
ggsave(filename, plot = rplot, device = "png", path = NULL,
       scale = 1, width = 16, height = 8, units = "cm",
       dpi = 600)

Thresh_Noise <- 0.5
sup800 <- min(which(WL>800))
NoisyWL <- which(abs_Rderiv0>Thresh_Noise | abs_Rderiv0==0)
NoisyWL <- NoisyWL[which(NoisyWL>sup800)]
# add spectral domains corresponding to water absorption
WaterAbs <- which((WL>895 & WL<1005) |
                    (WL>1400 & WL<1500) |
                    (WL>1780 & WL<2000) |
                    (WL>2350 & WL<2600))

# combine noisy bands and water absorption bands
NoisyWL <- sort(unique(c(NoisyWL,WaterAbs)))
SelectedWL <- seq(1,length(WL))[-NoisyWL]

# 3- plot reflectance after removal of noisy bands
Rsel <- data.frame(cbind(WL[SelectedWL],Rsel0[SelectedWL,]))
names(Rsel) <- c('wavelength',paste('R#',seq(1,nbSamples),sep=''))
ResRefl <- Rsel %>%
  tidyr::gather(Sample, value, -wavelength)
# create dataframe to prepare for figure
ResRefl <- data.frame(ResRefl,'RT'='R')

rplot <- ggplot(ResRefl, aes() ) +
  geom_line(data = ResRefl,aes(x=wavelength, y=value, group=Sample, colour=Sample) ,size=0.5) +
  labs(x="Wavelength (nm)", y="Reflectance (%)") +
  theme(legend.title =element_blank(),legend.position ="none") +
  ylim(0,100) + xlim(400,2550)

filename <- file.path(Path_Results,'04_Reflectance_noNoise.png')
ggsave(filename, plot = rplot, device = "png", path = NULL,
       scale = 1, width = 13, height = 10, units = "cm",
       dpi = 600)

# 4- save reflectance data after removal of noisy bands
Refl_VSWIR_noNoise <- data.frame(Refl_VSWIR[,SelectedWL])
names(Refl_VSWIR_noNoise) <- format(WL[SelectedWL],digits=1,nsmall = 1)
fileName <- file.path(Path_Results,'01_PREPROCESS_REFL/Refl_Filter_NoisyBands_Removed.csv')
dir.create(dirname(fileName),showWarnings = FALSE)
write_csv(x = Refl_VSWIR_noNoise,file = fileName,
            col_names = TRUE,append = FALSE)
ZipName <- paste(file_path_sans_ext(fileName),'.zip',sep = '')
zip(zipfile = ZipName,
    files =  fileName,mode = 'cherry-pick')
file.remove(fileName)

#---------------          Perform continuum removal             ---------------
# 5- perform continuum removal on VNIR reflectance data and save file
WL_VNIR = WL[SelectedWL]
WL_VNIR = WL_VNIR[WL_VNIR<1000]

SelectedWL_VNIR = which(WL<1000)
SelectedWL_VNIR = which(SelectedWL_VNIR %in% SelectedWL)
Spectral_VNIR <- list('Wavelength'=WL_VNIR)
Refl_CR_VNIR <- apply_continuum_removal(Spectral_Data=0.01*Refl_VSWIR[,SelectedWL_VNIR],
                                        Spectral = Spectral_VNIR, nbCPU = nbCPU)
Refl_CR <- data.frame(Refl_CR_VNIR)
Refl_CR2 <- format(Refl_CR,digits=5,nsmall = 5)
names(Refl_CR2) <- format(WL[SelectedWL_VNIR[2:(length(SelectedWL_VNIR)-1)]],digits=1,nsmall = 1)
fileName <- file.path(Path_Results,'01_PREPROCESS_REFL/Refl_VNIR_ContinuumRemoval.csv')
write_csv(x = Refl_CR2,file = fileName,
          col_names = TRUE,append = FALSE)
ZipName <- paste(file_path_sans_ext(fileName),'.zip',sep = '')
zip(zipfile = ZipName,
    files =  fileName,mode = 'cherry-pick')
file.remove(fileName)

# 6- perform continuum removal on VSWIR reflectance data and save file
Spectral <- list('Wavelength'=WL[SelectedWL])
Refl_CR <- apply_continuum_removal(Spectral_Data=0.01*Refl_VSWIR[,SelectedWL],
                                   Spectral = Spectral, nbCPU = nbCPU)
Refl_CR <- data.frame(Refl_CR)
Refl_CR2 <- format(Refl_CR,digits=5,nsmall = 5)
names(Refl_CR2) <- format(WL[SelectedWL_VNIR[2:(length(SelectedWL_VNIR)-1)]],digits=1,nsmall = 1)
fileName <- file.path(Path_Results,'01_PREPROCESS_REFL/Refl_VSWIR_ContinuumRemoval.csv')
write_csv(x = Refl_CR2,file = fileName,
          col_names = TRUE,append = FALSE)
ZipName <- paste(file_path_sans_ext(fileName),'.zip',sep = '')
zip(zipfile = ZipName,
    files =  fileName,mode = 'cherry-pick')
file.remove(fileName)

TreeID <- data.frame(SPID=SPID,PlotID=PlotID,CrownID=CrownID)
fileName <- file.path(Path_Results,'01_PREPROCESS_REFL/TreeID.csv')
write_csv(x = TreeID,file = fileName,
          col_names = TRUE,append = FALSE)
