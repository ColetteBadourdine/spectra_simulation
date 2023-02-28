# ==============================================================================
# Lib_diversity.R
# ==============================================================================
# PROGRAMMERS:
# Colette BADOURDINE <colette.badourdine@cirad.fr>
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Raphael Pelissier <raphael.pelissier@ird.fr>
# Copyright 2021/06 Colette BADOURDINE
# ==============================================================================
# This Library includes functions dedicated  to create abundance table from
# composition table also functions to compute diversity index
# ==============================================================================


#### library of functions

# library(ade4)
# library(vegan)
# library(ineq)
# library(rlist)

# source('Libraries/Lib_flashpca_jbf.R')
# source('Libraries/Lib_pretreatment.R')

#' divtaxo
#' compute taxonomic diversity indices such as richness shannon and simpson
#' @param facSp factor species names ('IDSP')
#'
#' @return list containing richness, simpson and shannon value
#' @export
divtaxo<-function(facSp){
  stopifnot(class(facSp)=="factor")
  Rich<-nlevels(facSp)
  n<-length(facSp)
  rabtab<-table(facSp)/n
  Simp<-1-sum(rabtab^2)
  Shan<-(-sum(rabtab*log(rabtab)))
  return(list(Richness=Rich,Simpson=round(Simp,2),Shannon=round(Shan,2)))
}

create_abundance_table <- function(PopTable, TAB){
  liste_sp = c()
  for (i in 1:length(PopTable)){
    cr_id = PopTable[i]
    liste_sp = c(liste_sp, subset(TAB, ID == cr_id)$SPID[1])
  }
  freq = as.data.frame(table(liste_sp))
  return(freq)
}

rao <- function(freq, distance_inter_sp){
  rao_pop = foreach(pop = 1:nrow(freq), .combine = 'c') %dopar% {
    rao = 0
    for (i in colnames(freq)){
      for (j in colnames(freq))
        rao = rao + (freq[pop, i] * freq[pop, j] * (distance_inter_sp[i, j]^2))
    }
    rao
  }
  return(rao_pop)
}
#' compute_diversity_index
#' Fonctions to compute results tables
#'
#' @param BasePopulations_min initial populations with low diversity (only 2 species)
#' @param BasePopulations_max initial populations with high diversity (100 species)
#' @param Populations_mintomax populations made from BasePopulations_min with an increasing diversity
#' @param Populations_maxtomin populations made from BasePopulations_max with a decreasing diversity
#' @param name_column names of the populations
#'
#' @return All_Populations_diversite dataframe with diversity index values for each populations
#' @export
compute_diversity_index <- function(BasePopulations_min, BasePopulations_max,
                                    Populations_mintomax, Populations_maxtomin,
                                    name_column = c("BasePopulations_min", "BasePopulations_max", "Populations_mintomax",
                                                    "Populations_maxtomin")){
  Simpson <- Shannon <- Richness <- Population <- c()
  nbPop <- length(BasePopulations_min$div$Simpson)
  for (i in 1:nbPop){
    Simpson = c(Simpson, BasePopulations_min$div$Simpson[i])
    Shannon = c(Shannon, BasePopulations_min$div$Shannon[i])
    Richness = c(Richness, BasePopulations_min$div$Richesse[i])
    Population = c(Population, name_column[1])
  }
  pcDilution <- names(Populations_mintomax)
  for (dilute in pcDilution){
    nbPop <- length(Populations_mintomax[[dilute]]$div$Simpson)
    for (j in 1:nbPop){
      Simpson <- c(Simpson, Populations_mintomax[[dilute]]$div$Simpson[j])
      Shannon <- c(Shannon, Populations_mintomax[[dilute]]$div$Shannon[j])
      Richness <- c(Richness, Populations_mintomax[[dilute]]$div$Richesse[j])
      Population <- c(Population, paste(name_column[3], "_", dilute, sep=''))
    }
  }
  pcDilution <- names(Populations_maxtomin)
  for (dilute in pcDilution){
    nbPop <- length(Populations_maxtomin[[dilute]]$div$Simpson)
    for (j in 1:nbPop){
      Simpson = c(Simpson, Populations_maxtomin[[dilute]]$div$Simpson[j])
      Shannon = c(Shannon, Populations_maxtomin[[dilute]]$div$Shannon[j])
      Richness = c(Richness, Populations_maxtomin[[dilute]]$div$Richesse[j])
      Population = c(Population, paste(name_column[4], "_", dilute, sep=''))
    }
  }
  nbPop <- length(BasePopulations_max$div$Simpson)
  for (i in 1:nbPop){
    Simpson <- c(Simpson, BasePopulations_max$div$Simpson[i])
    Shannon <- c(Shannon, BasePopulations_max$div$Shannon[i])
    Richness <- c(Richness, BasePopulations_max$div$Richesse[i])
    Population <- c(Population, name_column[2])
  }

  Simpson <- as.data.frame(Simpson)
  Shannon <- as.data.frame(Shannon)
  Richness <- as.data.frame(Richness)
  All_Populations_diversite <- cbind(Simpson, Shannon, Richness, Population)
  # All_Populations_diversite <- compute_evenness(All_Populations_diversite, BasePopulations_min,
  #                                               BasePopulations_max, Populations_mintomax, Populations_maxtomin)
  return(All_Populations_diversite)
}
#' compute_evenness
#'
#' @param All_Populations_diversite dataframe. contains Populations ID, Simpson, Shannon, Richness index
#' @param BasePopulations_min initial populations with low diversity (only 2 species)
#' @param BasePopulations_max initial populations with high diversity (100 species)
#' @param Populations_mintomax populations made from BasePopulations_min with an increasing diversity
#' @param Populations_maxtomin populations made from BasePopulations_max with a decreasing diversity
#'
#' @return All_Populations_diversite dataframe with Pielou and Gini values added
#' @export
# /!\ calcul de l'indice de Gini faux !!! A CORRIGER #
compute_evenness <- function(All_Populations_diversite, BasePopulations_min, BasePopulations_max,
                             Populations_mintomax, Populations_maxtomin){
  Pielou = c()
  #Gini = c()
  print("############ COMPUTATION OF PIELOU INDEX ############")
  for (i in 1:nrow(All_Populations_diversite)){
    Pielou = c(Pielou,
               # All_Populations_diversite[i, ]$Shannon/round(log(All_Populations_diversite[i, ]$Richness), 2))
               All_Populations_diversite[i, ]$Shannon/log(All_Populations_diversite[i, ]$Richness))
  }

  # print("############ COMPUTATION OF GINI INDEX ############")
  # for (i in 1:10){
  #   Gini = c(Gini, ineq(BasePopulations_min$pop[i,],type="Gini"))
  # }
  #
  # for (i in 1:10){
  #   for (j in 1:10){
  #     Gini = c(Gini, ineq(Populations_mintomax[[i]]$pop[j,],type="Gini"))
  #   }
  # }
  #
  # for (i in 1:10){
  #   Gini = c(Gini, ineq(BasePopulations_max$pop[i,],type="Gini"))
  # }
  #
  # for (i in 1:10){
  #   for (j in 1:10){
  #     Gini = c(Gini, ineq(Populations_maxtomin[[i]]$pop[j,],type="Gini"))
  #   }
  # }

  All_Populations_diversite = cbind(Pielou, All_Populations_diversite)
  #All_Populations_diversite = cbind(Pielou, Gini, All_Populations_diversite)
  return(All_Populations_diversite)
}


#' get_abundance
#' Fonction to get the abundance table of community
#'
#' @param PopTable dataframe.contains ID Crowns for a population
#' @param DataTable initial and complete reflectance table
#'
#' @return abundance_table abundance table for the population
#' @export
get_abundance <- function(PopTable, DataTable){
  nbPop =nrow(PopTable)
  nbInd = ncol(PopTable)
  abundance_table = list()
  for (pop in 1:nbPop){
    print("#############")
    print(pop)
    sp_list = c()
    for (ind in 1:nbInd){
      print(ind)
      crown = as.integer(PopTable[pop, ind])
      sp = unique(DataTable[DataTable$ID == crown, ]$SPID)
      print(sp)
      sp_list = c(sp_list, sp)
    }
    abundance_table[[pop]] = table(sp_list)
  }
  return(abundance_table)
}



