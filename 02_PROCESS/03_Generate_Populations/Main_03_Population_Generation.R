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
# library(R.utils)
# R.utils::sourceDirectory("../00_Libraries/flashpca/")
# source('../00_Libraries/Lib_flashpca_jbf.R')
library(zip)
library(dplyr)
library(ade4)
library(readr)
library(tools)
library(progress)
library(ggplot2)
library(cowplot)


#--------------- define path for data and results ---------------
Path_Data <- '../../03_RESULTS/02_Produce_FilteredDatafiles'
File_Data <- 'TreeID_Filtered.zip'
Path_Results <- '../../03_RESULTS/03_Produce_POP'
dir.create(path = Path_Results, showWarnings = F, recursive = T)

#--------------- read data and discard crowns with too few pixels ---------------
# read data type
TreeIDZipFile <- file.path(Path_Data,File_Data)
fid <- unz(description = TreeIDZipFile, filename = zip_list(TreeIDZipFile)$filename,
           open = "rb", encoding = getOption("encoding"))
TreeID <- read_delim(file = fid,delim = ',')
close(fid)

# get ITC abundance for all species
abund_ <- sortAbundance(DataTable = TreeID)
speciesFilter <- abund_$speciesFilter

# generating population corresponding to low diversity and saving populations
BasePopulations_min <- low_diversity_mixture(speciesFilter = speciesFilter, nbSpeciesPop = 2, nbPop = 10)

fileName <- file.path(Path_Results,'BasePopulation_min.rds')
saveRDS(BasePopulations_min, file = fileName)

# generating population corresponding to high diversity and saving populations
BasePopulations_max <- high_diversity_mixture(speciesFilter = speciesFilter, nbSpeciesPop = 50, nbPop = 10)

fileName <- file.path(Path_Results, 'BasePopulation_max.rds')
saveRDS(BasePopulations_max, file = fileName)

# increase diversity of populations BasePopulationsMin
Populations_mintomax <- list()
for (percentage in seq(10, 100, by=10)){
  PCchar <- paste(as.character(percentage),'%',sep = '')
  Populations_mintomax[[PCchar]] = PopEnrichment(PopTable = BasePopulations_min$pop,
                                                 SpeciesName = unique(speciesFilter$SPID),
                                                 percentage = percentage)
  }
fileName <- file.path(Path_Results, 'Populations_mintomax.rds')
saveRDS(Populations_mintomax, file = fileName)

# decrease diversity of populations BasePopulationsMax
Populations_maxtomin <- list()
for (percentage in seq(10, 100, by=10)){
  PCchar <- paste(as.character(percentage),'%',sep = '')
  Populations_maxtomin[[PCchar]] = PopEnrichment(PopTable = BasePopulations_max$pop,
                                                 SpeciesName = unique(speciesFilter$SPID),
                                                 percentage = percentage)
  }
fileName <- file.path(Path_Results, 'Populations_maxtomin.rds')
saveRDS(Populations_maxtomin, file = fileName)


