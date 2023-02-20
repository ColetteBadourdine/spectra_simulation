# ==============================================================================
# Lib_diversity.R
# ==============================================================================
# PROGRAMMERS:
# Colette BADOURDINE <colette.badourdine@cirad.fr>
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Raphael Pelissier <raphael.pelissier@ird.fr>
# Copyright 2021/06 Colette BADOURDINE
# ==============================================================================
# This Library includes functions dedicated to create artificial populations
# and to modify the composition of these populations
# ==============================================================================
# library(ggplot2)
# library(tidyverse)
# library(gridExtra)
# library(cowplot)
# source('Lib_diversity.R')


#' generateSpectra
#'
#' @param DataTable initial and complete reflectance table
#' @param SpeciesName species we want to simulate spectra
#' @param NbSpectra number of spectra we want to generate
#' @param sd_param values of sd to choose intrasp variance
#'
#' @return spectra table
#' @export
generateSpectra <- function(DataTable, SpeciesName, NbSpectra, sd_param){
  tab_temp <- DataTable[DataTable$SPID == SpeciesName, ]

  L = chol(cov(tab_temp[, -c(1:3)]))
  nvars = dim(L)[1]

  r = t(L) %*% matrix(rnorm(nvars*NbSpectra, sd = sd_param), nrow = nvars, ncol = NbSpectra)
  spectra_simu = r + matrix(colMeans(tab_temp[, -c(1:3)]),nrow = nvars, ncol = NbSpectra)

  return(spectra_simu)
}


#' gets reflectance spectra from a population
#'
#' @param DataTable initial and complete reflectance table
#'
#' @return sorted abundance of the species.
#' @return speciesFilter correspondance of species id and crowns id
#' @export
GetSpectra_FromPops <- function(Populations, ReflectanceData, PopID, Spectralpop, nbPop, sd_param){
  for (i in 1:nbPop){
    Spectral_temp = c()
    abund = as.data.frame(table(t(Populations)[, i]))
    for (j in 1:nrow(abund)){
      sp = abund[j, 1]
      abund_sp = abund[j, 2]
      spectra_sp = generateSpectra(ReflectanceData, SpeciesName = sp, NbSpectra = (abund_sp * 20), sd_param = sd_param)
      spectra_sp = as.data.frame(t(spectra_sp))
      spectra_sp = cbind(rep(sp, (abund_sp * 20)), spectra_sp)
      Spectral_temp = rbind(Spectral_temp, spectra_sp)
    }
    colnames(Spectral_temp) = colnames(ReflectanceData[, -c(2,3)])
    PopID = PopID + 1
    Spectralpop[[PopID]] <- Spectral_temp
  }
  return(list('PopID'=PopID,'Spectralpop'=Spectralpop))
}



#' sortedAbundance
#'
#' @param DataTable initial and complete reflectance table
#'
#' @return sorted abundance of the species.
#' @return speciesFilter correspondance of species id and crowns id
#' @export
sortAbundance <- function(DataTable){
  SpeciesName <- unique(DataTable$SPID)
  # Filter data to get abundance of each species
  nbsamples <- length(DataTable$ID)

  speciesFilter <- unique(DataTable[c('SPID', 'ID')])

  abundance <- table(speciesFilter$SPID)
  sortedAbundance <- as.data.frame(sort(abundance,decreasing = TRUE))
  names(sortedAbundance) <- c('species','Frequency')
  list = list('sortedAbundance' = sortedAbundance, 'speciesFilter' = speciesFilter)
  return(list)
}


#' low_diversity_mixture
#' ##########################################################
#' create low diversity mixture
#' #########################################################
#' sort species by abundance and select most abundant ones with more than 50 crowns
#' in order to be able to produce populations of 100 species with as low as 2 species
#' among the most abundant
#'
#' create 10 populations of 50 crowns selected among 2 species
#'
#' @param speciesFilter matching species id and crowns id
#' @param nbSpeciesPop number of species for the populations
#' @param nbPop number of populations to generate
#'
#' @return BasePopulations_min nbPop populations (low diversity)
#' @export
low_diversity_mixture <- function(speciesFilter, nbSpeciesPop = 2, nbPop = 10){

  species_names <- unique(speciesFilter$SPID)
  pop <- div <- c()
  for(i in 1:nbPop){
    # sample 100 crowns from 2 randomly selected abundant species
    spe_list <- sample(species_names,nbSpeciesPop)
    abund_pop_temp = c()
    abund_pop_temp[1] = sample(49, 1)
    abund_pop_temp[2] = 50 - abund_pop_temp[1]
    fac_sp = as.factor(c(rep(spe_list[1], abund_pop_temp[1]),
                         rep(spe_list[2], abund_pop_temp[2])))
    diversite = divtaxo(fac_sp)
    div_temp <-  cbind(diversite$Simpson, diversite$Shannon, diversite$Richness)

    #generate spectra
    pop_temp <- t(as.data.frame(fac_sp))
    pop_temp <- cbind.data.frame(pop_temp)
    pop <- rbind(pop, pop_temp)
    div <- rbind(div, div_temp)
  }
  div = as.data.frame(div)
  colnames(div) = c("Simpson", "Shannon", "Richesse")
  rownames(div) = c(1:nbPop)
  rownames(pop) = c(1:nbPop)
  BasePopulations_min= list(pop=pop, div=div)
  return(BasePopulations_min)
}


#' high_diversity_mixture
#' ##########################################################
#' create high diversity mixture
#' create 10 populations
#' #########################################################
#' 55 species (one crown per species)
#'
#' @param speciesFilter matching species id and crowns id
#' @param nbSpeciesPop number of species for the populations
#' @param nbPop number of populations to generate
#'
#' @return BasePopulations_max nbPop populations (low diversity)
#' @export
high_diversity_mixture <- function(speciesFilter, nbSpeciesPop = 50, nbPop = 10){

  pop <- div <- c()
  SpeciesName <- unique(speciesFilter$SPID)
  for(i in 1:nbPop){
    # randomly select 100 species among all available
    spe_list <- sample(SpeciesName,nbSpeciesPop)
    fac_sp <- as.factor(spe_list)
    diversite <- divtaxo(fac_sp)
    div_temp <-  cbind(diversite$Simpson, diversite$Shannon, diversite$Richness)
    pop_temp <- t(as.data.frame(spe_list))
    pop_temp <- cbind.data.frame(pop_temp)
    pop <- rbind(pop, pop_temp)
    div <- rbind(div, div_temp)
  }
  div <- as.data.frame(div)
  colnames(div) <- c("Simpson", "Shannon", "Richesse")
  rownames(div) <- c(1:nbPop)
  rownames(pop) <- c(1:nbPop)
  BasePopulations_max = list(pop=pop, div=div)
  return(BasePopulations_max)
}



#' PopEnrichment
#' Fonction to enrich population
#'
#' @param PopTable dataframe.contains ID Crowns for a population
#' @param SpeciesName names of all the species in DataTable
#' @param percentage percentage of sample to remplace from PopTable
#' @return NewPopulations A new populations
#' @importFrom progress progress_bar
#' @export
PopEnrichment <- function(PopTable, SpeciesName, percentage){

  nbPop <- nrow(PopTable)
  nbInd <- ncol(PopTable)
  NewPop <- c()
  div <- c()
  # use progress bar
  pb <- progress_bar$new(
    format = "Producing artificial population [:bar] :percent in :elapsedfull",
    total = nbPop, clear = FALSE, width= 100)
  pc <- percentage/100
  for (pop in 1:nbPop){
    pb$tick()
    nind_to_delete <- nbInd * pc
    sel_to_delete <- sample(PopTable[pop,], nind_to_delete)
    newpop <- PopTable[pop,]
    newpop <- subset(newpop, select = -c(as.integer(colnames(sel_to_delete))))
    sel_to_add <- sample(SpeciesName, nind_to_delete, replace=FALSE)
    newpop <- cbind(newpop, t(sel_to_add))
    colnames(newpop) <- c(1:nbInd)
    fac_sp <- as.factor(newpop)
    diversite <- divtaxo(fac_sp)
    div_temp <- cbind(diversite$Simpson, diversite$Shannon, diversite$Richness)
    div <- rbind(div, div_temp)
    NewPop <- rbind(NewPop, newpop)
  }
  div <- as.data.frame(div)
  colnames(div) = c("Simpson", "Shannon", "Richesse")
  rownames(div) = c(1:nbPop)
  NewPopulations = list(pop=NewPop, div=div)
  return(NewPopulations)
}
