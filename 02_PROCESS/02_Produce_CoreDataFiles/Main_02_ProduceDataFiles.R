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
source('../00_Libraries/Lib_ImageProcess.R')
library(zip)
library(readr)

# data filtering: min number of pixels per spcies. all species with less pixels will be discarded from dataset
minPixels <- 125
# multiplying factorFactor
FactorRefl <- 10000


#--------------- define path for data and results ---------------
Path_Data <- '../../03_RESULTS/01_PREPROCESS_REFL'
File_Data <- list()
File_Data <- file.path(Path_Data,'Refl_Filter_NoisyBands_Removed.zip')
Path_Results <- '../../03_RESULTS/02_Produce_FilteredDatafiles'
dir.create(path = Path_Results, showWarnings = F, recursive = T)

#--------------- read data and discard crowns with too few pixels ---------------
# read reflectance data directly from zipfile
ReflZipFile <- File_Data
fid <- unz(description = ReflZipFile, filename = zip_list(ReflZipFile)$filename,
           open = "rb", encoding = getOption("encoding"))
Reflectance_VNIR <- read_delim(file = fid,delim = ',')
close(fid)

# read treeID data
TreeID <- read_delim(file = file.path(Path_Data,'TreeID.csv'),delim = ',')
Reflectance <- Reflectance_VNIR
colnames(TreeID) <- c("SPID", "NParcelle", "ID")
colnames_wl <- colnames(Reflectance)
ReflTrees <- cbind(TreeID, Reflectance)
#select vnir only
whichElim <- which(as.numeric(names(ReflTrees))>1000)
ReflTrees <- ReflTrees[,-whichElim]

# suppression de ttes les especes ayant moins de 125 px
sp_list <- as.data.frame(table(ReflTrees$SPID))
sp_list_minsize <- subset(sp_list, Freq >= minPixels)
ReflTrees_minsize <- subset(ReflTrees, SPID %in% sp_list_minsize$Var1)

TreeID_minsize <- subset(TreeID, SPID %in% sp_list_minsize$Var1)
fileName <- file.path(Path_Results,'TreeID_Filtered.csv')
Save_n_Zip(data2save = TreeID_minsize, fileName = fileName)

#save this VSWIR data table
ReflTrees_minsize[,-c(1,2,3)] <- ReflTrees_minsize[,-c(1,2,3)]/FactorRefl
fileName <- file.path(Path_Results,'ReflTrees_VNIR.csv')
Save_n_Zip(data2save = ReflTrees_minsize, fileName = fileName)

