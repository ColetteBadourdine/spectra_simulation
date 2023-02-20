# ==============================================================================
# Lib_flashpca.R
# ==============================================================================
# PROGRAMMERS:
# Colette BADOURDINE <colette.badourdine@cirad.fr>
# Copyright 2021/02 Colette BADOURDINE
# ==============================================================================
# This Library includes functions to compute and analyse intra/inter specific
# variances
# ==============================================================================

#' weighted_intravarsp_pop
#'
#' @param DataTable dataframe. includes species names ('SPID') and crowns ID ('ID')
#' @param Intra_specific_variance. Obtained by running function 'compute_intra_specific_variance'
#' @param BasePopulations_min min populations obtained by running main_02_population_generation
#' @param BasePopulations_max max populations obtained by running main_02_population_generation
#' @param Populations_mintomax obtained by running main_02_population_generation
#' @param Populations_maxtomin obtained by running main_02_population_generation
#' @return mean_intravarsp_per_pop A data.frame containing mean intraspecific variance
#' weighted by abundance for each population
#' @export
weighted_intravarsp_pop <- function(DataTable, IntraSpecificVariance,
                                    Pop_min,Pop_mintomax,Pop_max, Pop_maxtomin){
  mean_intravarsp_per_pop = c()

  print("############### BasePopulationsMin ###############")

  abundance_table = apply(Pop_min$pop, 1, FUN=function(x) create_abundance_table(x, DataTable))
  for (i in 1:nrow(Pop_min$pop)){
    abundance = c()
    intraspvar = c()
    abundance_table[[i]] = as.data.frame(abundance_table[[i]])
    colnames(abundance_table[[i]]) = c('SPID', 'abundance')
    for (sp in abundance_table[[i]]$SPID){
      intraspvar = c(intraspvar, IntraSpecificVariance[rownames(IntraSpecificVariance)==sp,])
    }
    abundance = c(abundance,abundance_table[[i]]$abundance)
    intraspvar = as.data.frame(intraspvar)
    intraspvar = cbind(intraspvar, abundance)
    mean_intravarsp_per_pop = c(mean_intravarsp_per_pop, weighted.mean(intraspvar$intraspvar, intraspvar$abundance))
  }

  print("############### PopulationMintoMax ###############")

  for (i in 1:length(Pop_mintomax)){
    abundance_table = apply(Pop_mintomax[[i]]$pop, 1, FUN=function(x) create_abundance_table(x, DataTable))
    for (j in 1:nrow(Pop_mintomax[[i]]$pop)){
      abundance = c()
      intraspvar = c()
      abundance_table[[i]] = as.data.frame(abundance_table[[i]])
      colnames(abundance_table[[i]]) = c('SPID', 'abundance')
      for (sp in abundance_table[[i]]$SPID){
        intraspvar = c(intraspvar, IntraSpecificVariance[rownames(IntraSpecificVariance)==sp,])
      }
      abundance = c(abundance,abundance_table[[i]]$abundance)
      intraspvar = as.data.frame(intraspvar)
      intraspvar = cbind(intraspvar, abundance)
      mean_intravarsp_per_pop = c(mean_intravarsp_per_pop, weighted.mean(intraspvar$intraspvar, intraspvar$abundance))
    }
  }

  print("############### BasePopulationMax ###############")

  abundance_table = apply(Pop_max$pop, 1, FUN=function(x) create_abundance_table(x, DataTable))
  for (i in 1:nrow(Pop_max$pop)){
    abundance = c()
    intraspvar = c()
    abundance_table[[i]] = as.data.frame(abundance_table[[i]])
    colnames(abundance_table[[i]]) = c('SPID', 'abundance')
    for (sp in abundance_table[[i]]$SPID){
      intraspvar = c(intraspvar, IntraSpecificVariance[rownames(IntraSpecificVariance)==sp,])
    }
    abundance = c(abundance,abundance_table[[i]]$abundance)
    intraspvar = as.data.frame(intraspvar)
    intraspvar = cbind(intraspvar, abundance)
    mean_intravarsp_per_pop = c(mean_intravarsp_per_pop, weighted.mean(intraspvar$intraspvar, intraspvar$abundance))
  }

  print("############### PopulationMaxtoMin ###############")

  for (i in 1:length(Pop_maxtomin)){
    abundance_table = apply(Pop_maxtomin[[i]]$pop, 1, FUN=function(x) create_abundance_table(x, DataTable))
    for (j in 1:nrow(Pop_maxtomin[[i]]$pop)){
      abundance = c()
      intraspvar = c()
      abundance_table[[i]] = as.data.frame(abundance_table[[i]])
      colnames(abundance_table[[i]]) = c('SPID', 'abundance')
      for (sp in abundance_table[[i]]$SPID){
        intraspvar = c(intraspvar, IntraSpecificVariance[rownames(IntraSpecificVariance)==sp,])
      }
      abundance = c(abundance,abundance_table[[i]]$abundance)
      intraspvar = as.data.frame(intraspvar)
      intraspvar = cbind(intraspvar, abundance)
      mean_intravarsp_per_pop = c(mean_intravarsp_per_pop, weighted.mean(intraspvar$intraspvar, intraspvar$abundance))
    }
  }

  mean_intravarsp_per_pop = as.data.frame(mean_intravarsp_per_pop)
  return(mean_intravarsp_per_pop)
}


#' compute_MeanSPectraSpecies
#' mean spectra data.frame generation
#' generate au data.frame of the mean spectra per species and its transposition these two data.frame
#' can be use to calculate distance
#'
#' @param tab initial and complete reflectance table containing id crowns (ID)
#' species id (SPID) and parcelle (NParcelle)
#' @param idSP SPID as factor
#'
#' @return list containing mean_spectra_species and mean_spectra_species_t
#' @export
compute_MeanSPectraSpecies <- function(tab, idSP){
  mean_spectra_species = matrix(ncol = (ncol(tab)-3))
  colnames(mean_spectra_species) = colnames(tab[, -c(1:3)])
  mean_spectra_species = as.data.frame(mean_spectra_species)
  for (s in levels(idSP)){
    #print(s)
    mean_spectra_species = rbind(mean_spectra_species,t(colMeans(subset(tab, SPID==s)[,-c(1:3)])))
  }
  mean_spectra_species = mean_spectra_species[-1, ]
  mean_spectra_species = cbind(unique(levels(idSP)), mean_spectra_species)
  colnames(mean_spectra_species)[1] = 'SPID'
  mean_spectra_species_t = as.data.frame(t(mean_spectra_species))
  colnames(mean_spectra_species_t) = mean_spectra_species_t[1, ]
  mean_spectra_species_t = mean_spectra_species_t[-1, ]
  mean_spectra_species_t = as.data.frame(lapply(mean_spectra_species_t,as.numeric))
  return(list(mean_spectra_species = mean_spectra_species, mean_spectra_species_t = mean_spectra_species_t))
}


#' compute_MeanSPectraCrowns
#' mean spectra data.frame generation
#' generate au data.frame of the mean spectra per crown and its transposition these two data.frame
#' can be use to calculate distance
#'
#' @param tab initial and complete reflectance table containing id crowns (ID)
#' species id (SPID) and parcelle (NParcelle)
#' @param ITC ID as factor
#'
#' @return list containing mean_spectra_species and mean_spectra_species_t
#' @export
compute_MeanSPectraCrowns <- function(tab, ITC){
  mean_spectra_crowns = matrix(ncol = (ncol(tab)-3))
  colnames(mean_spectra_crowns) = colnames(tab[, -c(1:3)])
  mean_spectra_crowns = as.data.frame(mean_spectra_crowns)
  for (s in levels(ITC)){
    #print(s)
    mean_spectra_crowns = rbind(mean_spectra_crowns,t(colMeans(subset(tab, ID==s)[,-c(1:3)])))
  }
  mean_spectra_crowns = mean_spectra_crowns[-1, ]
  mean_spectra_crowns = cbind(unique(levels(ITC)), mean_spectra_crowns)
  colnames(mean_spectra_crowns)[1] = 'ID'
  mean_spectra_crowns_t = as.data.frame(t(mean_spectra_crowns))
  colnames(mean_spectra_crowns_t) = mean_spectra_crowns_t[1, ]
  mean_spectra_crowns_t = mean_spectra_crowns_t[-1, ]
  mean_spectra_crowns_t = as.data.frame(lapply(mean_spectra_crowns_t,as.numeric))
  return(list(mean_spectra_crowns = mean_spectra_crowns, mean_spectra_crowns_t = mean_spectra_crowns_t))
}
