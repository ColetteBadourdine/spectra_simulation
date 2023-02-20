# ==============================================================================
# Lib_flashpca.R
# ==============================================================================
# PROGRAMMERS:
# Colette BADOURDINE <colette.badourdine@cirad.fr>
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Raphael Pelissier <raphael.pelissier@ird.fr>
# Copyright 2021/02 Colette BADOURDINE
# ==============================================================================
# This Library includes functions dedicated to pca and pcaiv using flashpca
# and functions decicated to simulation of population of individual tree crowns
# based on inventories
# ==============================================================================

# library(ade4)
# # library(flashpcaR)
# library(vegan)
# library(rlist)
# library(ggplot2)
# library(ggrepel)
# source('Libraries/Lib_diversity_jbf.R')



#' #' fpca
#' #'
#' #' @param mat dataframe. includes species names ('IDSP') and crowns ID ('ID')
#' #' @param stand character. allow to choose between scaled or centered data for the fpca
#' #' @param ndim numeric. number of dimension to be kept in fpca results
#' #'
#' #' @return res. results of pca and total variance
#' #' @export
#'
#' fpca<-function(mat,stand=c("scaled","centered"), ndim = 10) {
#'   mat<-as.matrix(mat)
#'   stand<-stand[1]
#'   if(stand=="centered") {
#'     res<-flashpca(as.matrix(tspc),ndim = ndim, stand="center",divisor="n1",do_loadings=TRUE)
#'     mat<-scale(mat,scale=FALSE)
#'   }
#'   else if (stand == "scaled") {
#'     res<-flashpca(as.matrix(tspc),ndim = ndim, stand="sd",divisor="n1",do_loadings=TRUE)
#'     mat<-scale(mat)
#'   }
#'   else
#'     stopifnot(stand%in%c("centered","scaled)"))
#'   res$var<-as.vector(apply(mat,2,var))
#'   #correlation des variables avec les axes
#'   res$co<-sweep(res$loadings, 2, sqrt(res$values), "*")
#'   #coordonnees des lignes
#'   res$li<-as.matrix(mat%*%as.matrix(res$loadings))
#'   class(res)<-"fpca"
#'   attr(res,"stand")<-stand
#'   res$call=match.call()
#'   return(res)
#' }
#'
#' #' plot.fpca
#' #'
#' #' @param fpca dataframe. includes species names ('IDSP') and crowns ID ('ID')
#' #' @param xax
#' #' @param yax numeric. number of dimension to be kept in fpca results
#' #' @param ptcol numeric. number of dimension to be kept in fpca results
#' #' @param stand character. allow to choose different parameters if the data are scaled or centered
#' #'
#' #' @return plot of the pca.
#' #' @export
#'
#' plot.fpca<-function(fpca,xax=1,yax=2,ptcol, stand) {
#'   stopifnot(stand%in%c("centered","scaled"))
#'   plot.new()
#'   par(mfrow=c(2,2))
#'   plot(x=range(fpca$li[,xax]),y=range(fpca$li[,yax]),type="n",ann=FALSE,frame.plot=FALSE,axes=FALSE)
#'   title(sub=paste(c(fpca$call),"\nstand =",attributes(fpca)$stand,"\nPC",xax,"-",yax,"\n\n","spectral diversity=",round(sum(fpca$var),2)))
#'   #text(min(fpca$li[,xax]),pos=1,labels=paste(c(fpca$call),"\nstand =",attributes(fpca)$stand,"\nPC",xax,"-",yax,"\n\n","spectral diversity=",round(sum(fpca$var),2)))
#'   barplot(fpca$pve)
#'   if (stand == "scaled"){
#'     s.corcircle(fpca$co,xax=xax,yax=yax)
#'   }
#'   else if (stand == "centered"){
#'     s.corcircle(fpca$co,xax=xax,yax=yax, fullcircle = FALSE)
#'   }
#'
#'   if(missing(ptcol)) ptcol<-1
#'   plot.default(fpca$li[,xax],fpca$li[,yax],asp=1,col=ptcol,cex=0.5,ann=FALSE)
#' }
#'
#' #' fpcaiv.fit
#' #'
#' #' @param mat dataframe. includes species names ('IDSP') and crowns ID ('ID')
#' #' @param fac factor. can be either crowns ID or species names
#' #' @param stand character. allow to choose between scaled or centered data for the fpca
#' #' @param ndim numeric. number of dimension to be kept in fpca results
#' #'
#' #' @return res. results of pcaiv and the ajusted table given a chosen factor
#' #' @export
#'
#' fpcaiv.fit<-function(mat,fac,stand=c("scaled","centered"), ndim = 10) {
#'   mat<-as.matrix(mat)
#'   stand<-stand[1]
#'   if(stand=="centered")
#'     mat<-scale(mat,scale=FALSE)
#'   else if (stand == "scaled")
#'     mat<-scale(mat)
#'   else
#'     stopifnot(stand%in%c("centered","scaled"))
#'   stopifnot(is.factor(fac))
#'   tmp<-lm(mat~fac)
#'   mat<-tmp$fitted.values
#'   res<-flashpca(mat,ndim = ndim, stand="none",divisor="n1",do_loadings=TRUE)
#'   res$var<-as.vector(apply(mat,2,var))
#'   #corrélation des varibales avec les axes
#'   res$co<-sweep(res$loadings, 2, sqrt(res$values), "*")
#'   #coordonnées des lignes
#'   res$li<-as.matrix(mat%*%as.matrix(res$loadings))
#'   class(res)<-"fpca"
#'   attr(res,"stand")<-stand
#'   res$call=match.call()
#'   return(res)
#' }

#' fpcaiv.res
#'
#' @param mat dataframe. includes species names ('IDSP') and crowns ID ('ID')
#' @param fac factor. can be either crowns ID or species names
#' @param stand character. allow to choose between scaled or centered data for the fpca
#' @param ndim numeric. number of dimension to be kept in fpca results
#'
#' @return res. results of pcaiv and the residual table of pcaiv given a chosen factor
#' @export

fpcaiv.res<-function(mat,fac,stand=c("scaled","centered"), ndim=10) {
  mat<-as.matrix(mat)
  stand<-stand[1]
  if(stand=="centered")
    mat<-scale(mat,scale=FALSE)
  else if (stand == "scaled")
    mat<-scale(mat)
  else
    stopifnot(stand%in%c("centered","scaled"))
  stopifnot(is.factor(fac))
  tmp<-lm(mat~fac)
  mat<-tmp$residuals
  res<-flashpca(mat,ndim = ndim, stand="none",divisor="n1",do_loadings=TRUE)
  res$var<-as.vector(apply(mat,2,var))
  #corrélation des varibales avec les axes
  res$co<-sweep(res$loadings, 2, sqrt(res$values), "*")
  #coordonnées des lignes
  res$li<-as.matrix(mat%*%as.matrix(res$loadings))
  class(res)<-"fpca"
  attr(res,"stand")<-stand
  res$call=match.call()
  return(res)
}
#'
#' #' fpcaiv.hier (hierarchical analysis)
#' #'
#' #' @param mat dataframe. includes species names ('IDSP') and crowns ID ('ID')
#' #' @param fac1 factor. can be either crowns ID or species names
#' #' @param fac2 factor. can be either crowns ID or species names but it should be different from fac1
#' #' @param stand character. allow to choose between scaled or centered data for the fpca
#' #' @param ndim numeric. number of dimension to be kept in fpca results
#' #'
#' #' @return res. results of pcaiv for two explanatory variables (hierachical analysis)
#' #'  and the ajusted table given two chosen factors
#' #' @export
#'
#' fpcaiv.hier<-function(mat,fac1,fac2,stand=c("scaled","centered"), ndim = 10) {
#'   mat<-as.matrix(mat)
#'   stand<-stand[1]
#'   if(stand=="centered")
#'     mat<-scale(mat,scale=FALSE)
#'   else if (stand == "scaled")
#'     mat<-scale(mat)
#'   else
#'     stopifnot(stand%in%c("centered","scaled"))
#'   stopifnot(is.factor(fac1))
#'   stopifnot(is.factor(fac2))
#'   tmp<-lm(mat~fac1)
#'   mat<-tmp$residuals
#'   tmp<-lm(mat~fac2)
#'   mat<-tmp$fitted.values
#'   res<-flashpca(mat,ndim = ndim, stand="none",divisor="n1",do_loadings=TRUE)
#'   res$var<-as.vector(apply(mat,2,var))
#'   #corrélation des varibales avec les axes
#'   res$co<-sweep(res$loadings, 2, sqrt(res$values), "*")
#'   #coordonnées des lignes
#'   res$li<-as.matrix(mat%*%as.matrix(res$loadings))
#'   class(res)<-"fpca"
#'   attr(res,"stand")<-stand
#'   res$call=match.call()
#'   return(res)
#' }
#'
#' #' fpcaiv.hier.res (hierarchical analysis)
#' #'
#' #' @param mat dataframe. includes species names ('IDSP') and crowns ID ('ID')
#' #' @param fac1 factor. can be either crowns ID or species names
#' #' @param fac2 factor. can be either crowns ID or species names but it should be different from fac1
#' #' @param stand character. allow to choose between scaled or centered data for the fpca
#' #' @param ndim numeric. number of dimension to be kept in fpca results
#' #'
#' #' @return res. results of pcaiv for two explanatory variables and the residual table
#' #' @export
#'
#' fpcaiv.hier.res<-function(mat,fac1,fac2,stand=c("scaled","centered"), ndim = 10) {
#'   mat<-as.matrix(mat)
#'   stand<-stand[1]
#'   if(stand=="centered")
#'     mat<-scale(mat,scale=FALSE)
#'   else if (stand == "scaled")
#'     mat<-scale(mat)
#'   else
#'     stopifnot(stand%in%c("centered","scaled"))
#'   stopifnot(is.factor(fac1))
#'   stopifnot(is.factor(fac2))
#'   tmp<-lm(mat~fac1)
#'   mat<-tmp$residuals
#'   tmp<-lm(mat~fac2)
#'   mat<-tmp$residuals
#'   res<-flashpca(mat,ndim = ndim, stand="none",divisor="n1",do_loadings=TRUE)
#'   res$var<-as.vector(apply(mat,2,var))
#'   #corrélation des varibales avec les axes
#'   res$co<-sweep(res$loadings, 2, sqrt(res$values), "*")
#'   #coordonnées des lignes
#'   res$li<-as.matrix(mat%*%as.matrix(res$loadings))
#'   class(res)<-"fpca"
#'   attr(res,"stand")<-stand
#'   res$call=match.call()
#'   return(res)
#' }


#' var_calc
#' calculate variance of each virtual community
#' estimation of spectral diversity
#'
#' @param mat dataframe.
#' @param fac_esp factor. species names ('SPID')
#' @param fac_crown factor. crowns ID ('ID')
#' @param stand character. allow to choose between scaled or centered data for the fpca
#'
#' @return res. list of variance totale, species variance, crown variance and residual variance
#' @export
#

var_calc <- function(stand = "scaled", mat, fac_esp,var_species = TRUE){
  stopifnot(is.factor(fac_esp))
  nbVars <- ncol(mat)
  if (length(mat)==length(fac_esp)){
    mat <- matrix(mat,ncol=1)
  }
  if(stand=="centered") {
    mat <- scale(mat,scale=FALSE)
  } else if (stand == "scaled"){
    mat <- scale(mat)
  } else {
    mat <- scale(mat,center = F,scale = F)
  }
  var_tot <- sum(as.vector(apply(mat,2,var)))
  if (var_species == TRUE){
    tmp <- lm(mat~fac_esp)
    mat2 <- tmp$fitted.values
    if (nbVars==1) {
      mat2 <- matrix(mat2,ncol=1)
    }
    var_species = sum(as.vector(apply(mat2,2,var)))
  }
  # tmp<-lm(mat~fac_crown)
  # mat<-tmp$fitted.values
  # var_cr_inter = round(sum(as.vector(apply(mat,2,var))),2)
  res = list(var_tot = var_tot, var_species = var_species)
  return(res)
}

var_calc_only_sp <- function(stand = "scaled", mat, fac_esp){
  if(stand=="centered") {
    mat<-scale(mat,scale=FALSE)
    var = sum(as.vector(apply(mat,2,var)))
  }
  else if (stand == "scaled"){
    mat<-scale(mat)
    var = ncol(mat)
  }
  else {
    mat <- scale(mat, center = F, scale = F)
    var = sum(as.vector(apply(mat,2,var)))
  }
  stopifnot(is.factor(fac_esp))
  tmp<-lm(mat~fac_esp)
  mat<-tmp$fitted.values
  var_sp = sum(as.vector(apply(mat,2,var)))
  mat<-tmp$residuals
  var_res = sum(as.vector(apply(mat,2,var)))
  res = list(var_tot=var, var_species=var_sp,var_residu=var_res)
  return(res)
}



#' var_calc
#' calculate coefficient àf variation of each virtual community
#' estimation of spectral diversity
#'
#' @param mat dataframe.
#' @return res. list of variance totale, species variance, crown variance and residual variance
#' @export
#

cv_calc <- function(mat, type=2){
  cv_tot <- sum(as.vector(apply(mat,type,function(x) sd(x)/mean(x))))/ncol(mat)
  return(cv_tot)
}


#' spectral_variance
#'
#' Fonctions to compute spectral variance
#' @param InputData table including SPID, ID and reflectance or PCA data
#' @param stand choice to centered or scaled variance
#'
#' @return list containing total variance , inter-species spectral variance,
#' residual spectral centered or scaled
#' @export
spectral_variance <- function(InputData, stand="centered"){

  res_spectral <- var_calc(stand=stand, mat = InputData[, -c(1:3)],
                           fac_esp = as.factor(InputData$SPID),var_species = T)

  if (stand == "scaled"){
    return(list("var_tot_scaled" = res_spectral$var_tot,
                "var_sp_scaled" = res_spectral$var_species))
  } else if (stand == "centered"){
    return(list("var_tot_centered" = res_spectral$var_tot,
                "var_sp_centered" = res_spectral$var_species))
  } else {
    return(list("var_tot" = res_spectral$var_tot,
                "var_sp" = res_spectral$var_species))
  }

}

#' spectral_variance_subset
#' computes variance for a subset of variables, assuming the 3 first columns
#' correspond to species ID info
#'
#' @param InputData dataframe. includes species names ('SPID') and crowns ID ('ID')
#' @param vars numeric. number of species (after 3 first columns)
#' @param stand character. allow to choose between scaled or centered data for the fpca, or none
#' @param var_crown boolean. should variance at crown scale be computed?
#' @param var_residu boolean. should residual variance be computed?
#'
#' @return res. list of variance totale, species variance, crown variance and residual variance
#' @export

spectral_variance_subset <- function(InputData, vars, stand="scaled"){

  Subset <- InputData[,c(1, (vars+1))]
  res_spectral <- var_calc(stand=stand, mat = as.data.frame(Subset[, -1]),
                           fac_esp = as.factor(InputData$SPID),var_species = T)

  if (stand == "scaled"){
    return(list("var_tot_scaled" = res_spectral$var_tot,
                "var_sp_scaled" = res_spectral$var_species))
  } else if (stand == "centered"){
    return(list("var_tot_centered" = res_spectral$var_tot,
                "var_sp_centered" = res_spectral$var_species))
  } else {
    return(list("var_tot" = res_spectral$var_tot,
                "var_sp" = res_spectral$var_species))
  }
}

#' cv_subset
#' computes coefficient of variation
#' for a subset of variables, assuming the 3 first columns
#' correspond to species ID info
#'
#' @param InputData dataframe. includes species names ('SPID') and crowns ID ('ID')
#' @param vars numeric. number of species (after 3 first columns)
#'
#' @return res. list of variance totale, species variance, crown variance and residual variance
#' @export

cv_subset <- function(InputData, vars){

  Subset <- InputData[,c(1,2,3,(vars+3))]
  res_spectral <- cv_calc(mat = Subset[, -c(1:3)])
  return(res_spectral)
}


#' XXXXXXX
#'
#' @param xxxxxx
#' @param xxxxxx
#' @param xxxxxx
#'
#' @return XXXXXX
#' @export

spectral_variance_for_mean_spectra_species <- function(InputData, vars, var_crown = T,
                                                       var_residu = T, stand="scaled"){
  Subset <- InputData[,c(1,(vars+1))]
  res_spectral <- var_calc(stand=stand, mat = as.data.frame(Subset[, -1]),
                           fac_esp = as.factor(InputData$SPID),
                           fac_crown =as.factor(InputData$ID),var_species = T,
                           var_crown = var_crown, var_residu = var_residu)

  if (stand == "scaled"){
    return(list("var_tot_scaled" = res_spectral$var_tot,
                "var_sp_scaled" = res_spectral$var_species,
                "var_cr_scaled" = res_spectral$var_crown,
                "var_res_scaled" = res_spectral$var_residu))
  } else if (stand == "centered"){
    return(list("var_tot_centered" = res_spectral$var_tot,
                "var_sp_centered" = res_spectral$var_species,
                "var_cr_centered" = res_spectral$var_crown,
                "var_res_centered" = res_spectral$var_residu))
  } else {
    return(list("var_tot" = res_spectral$var_tot,
                "var_sp" = res_spectral$var_species,
                "var_cr" = res_spectral$var_crown,
                "var_res" = res_spectral$var_residu))
  }
}




cbind.fill <- function(...){
  nm <- list(...)
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow))
  do.call(cbind, lapply(nm, function (x)
    rbind(x, matrix(, n-nrow(x), ncol(x)))))
}


#' This function computes rao entropy
#'
#' @param ReflTree dataframe. reflectance data with 3 first columns corresponding to species ID
#' @param mean_reflectance list. mean reflectance corresponding to each species in ReflTree
#' @param Pop double. date of acquisition
#'
#' @return raoVal numeric. rao entropy
#' @importFrom ecodist distance
#' @importFrom stats dist
#' @importFrom ade4 divc
#' @export

compute_rao <- function(ReflTree, mean_reflectance, Pop){
  rownames(mean_reflectance) <- mean_reflectance[, 1]
  mean_reflectance <- mean_reflectance[, -1]
  EuclideanDist <- ecodist::distance(mean_reflectance, method = "euclidean")
  EuclideanDist <- stats::dist(EuclideanDist)
  freq <- create_abundance_table(Pop, ReflTree)
  rownames_ <- freq[, 1]
  freq <- as.data.frame(freq[, -1])
  colnames(freq) <- NULL
  rownames(freq) <- rownames_
  raoVal <- ade4::divc(freq, EuclideanDist, scale = FALSE)$diversity
  return(raoVal)
}


#' compute_all_results
#' Fonctions to compute results tables
#'
#' @param BasePopulations_min initial populations with low diversity (only 2 species)
#' @param BasePopulations_max initial populations with high diversity (100 species)
#' @param Populations_mintomax populations made from BasePopulations_min with an increasing diversity
#' @param Populations_maxtomin populations made from BasePopulations_max with a decreasing diversity
#' @param All_Populations_diversite dataframe with diversity index values for each populations
#' @param DataTable reflectance table
#' @param name_column names of the populations
#'
#' @return All_Populations_results all populations results table with diversity and evenness
#' index and spectral variance values
#' @importFrom progress progress_bar
#' @export
compute_all_results <- function(BasePopulations_min, BasePopulations_max, Populations_mintomax,
                                Populations_maxtomin, All_Populations_diversite, DataTable,
                                name_column = c("BasePopulations_min",
                                                "BasePopulations_max",
                                                "Populations_mintomax",
                                                "Populations_maxtomin")){

  All_Populations_spectral_diversite <- list()
  Label <- c()
  # Minimum population
  nbPop <- nrow(BasePopulations_min$pop)
  pb <- progress_bar$new(
    format = "Generate populations with minimum diversity [:bar] :percent in :elapsedfull",
    total = nbPop, clear = FALSE, width= 100)
  # identifier for population
  PopID <- 0
  for (i in 1:nbPop){
    pb$tick()
    DataTable_temp = DataTable %>%
      group_by(ID) %>%
      sample_n(10)
    reflectance_table <- get_reflectance_table(BasePopulations_min$pop[i, ], DataTable_temp)
    SpectralVar <- spectral_variance(reflectance_table)
    raoVal <- compute_rao(ReflTree = DataTable_temp,
                          mean_reflectance = SpectralVar$mean_spectra,
                          Pop = BasePopulations_min$pop[i,])
    Label <- c(Label, "MIN")
    PopID <- PopID + 1
    All_Populations_spectral_diversite[[PopID]] <- data.frame("var_tot_centered" = SpectralVar$var_tot_centered,
                                                              "var_sp_centered" = SpectralVar$var_sp_centered,
                                                              "var_cr_centered" = SpectralVar$var_cr_centered,
                                                              "var_res_centered" = SpectralVar$var_res_centered,
                                                              "Rao" = raoVal,
                                                              "Population" = name_column[1])
  }


  # dilution of minimum population
  pcDilution <- names(Populations_mintomax)
  nbPop <- length(Populations_mintomax[[1]]$div$Simpson)
  pb <- progress_bar$new(
    format = "Dilute populations with minimum diversity [:bar] :percent in :elapsedfull",
    total = length(pcDilution)*nbPop, clear = FALSE, width= 100)
  # identifier for population
  for (dilute in pcDilution){
    nbPop <- length(Populations_mintomax[[dilute]]$div$Simpson)
    for (j in 1:nbPop){
      pb$tick()
      DataTable_temp = DataTable %>%
        group_by(ID) %>%
        sample_n(10)
      reflectance_table <- get_reflectance_table(Populations_mintomax[[dilute]]$pop[j,], DataTable_temp)
      SpectralVar <- spectral_variance(reflectance_table)
      raoVal <- compute_rao(ReflTree = DataTable_temp,
                            mean_reflectance = SpectralVar$mean_spectra,
                            Pop = Populations_mintomax[[dilute]]$pop[j,])
      Label <- c(Label, "MINTOMAX")
      PopID <- PopID + 1
      All_Populations_spectral_diversite[[PopID]] <- data.frame("var_tot_centered" = SpectralVar$var_tot_centered,
                                                                "var_sp_centered" = SpectralVar$var_sp_centered,
                                                                "var_cr_centered" = SpectralVar$var_cr_centered,
                                                                "var_res_centered" = SpectralVar$var_res_centered,
                                                                "Rao" = raoVal,
                                                                "Population" = paste(name_column[3], dilute, sep='_'))
    }
  }

  # dilution of maximum population
  pcDilution <- names(Populations_maxtomin)
  nbPop <- length(Populations_maxtomin[[1]]$div$Simpson)
  pb <- progress_bar$new(
    format = "Dilute populations with maximum diversity [:bar] :percent in :elapsedfull",
    total = length(pcDilution)*nbPop, clear = FALSE, width= 100)
  for (dilute in pcDilution){
    nbPop <- length(Populations_maxtomin[[dilute]]$div$Simpson)
    for (j in 1:nbPop){
      pb$tick()
      DataTable_temp = DataTable %>%
        group_by(ID) %>%
        sample_n(10)
      reflectance_table = get_reflectance_table(Populations_maxtomin[[dilute]]$pop[j,], DataTable_temp)
      SpectralVar <- spectral_variance(reflectance_table)
      raoVal <- compute_rao(ReflTree = DataTable_temp,
                            mean_reflectance = SpectralVar$mean_spectra,
                            Pop = Populations_maxtomin[[dilute]]$pop[j,])
      Label <- c(Label, "MAXTOMIN")
      PopID <- PopID + 1
      All_Populations_spectral_diversite[[PopID]] <- data.frame("var_tot_centered" = SpectralVar$var_tot_centered,
                                                                "var_sp_centered" = SpectralVar$var_sp_centered,
                                                                "var_cr_centered" = SpectralVar$var_cr_centered,
                                                                "var_res_centered" = SpectralVar$var_res_centered,
                                                                "Rao" = raoVal,
                                                                "Population" = paste(name_column[4], dilute, sep='_'))
    }
  }

  nbPop <- nrow(BasePopulations_max$pop)
  pb <- progress_bar$new(
    format = "Generate populations with maximum diversity [:bar] :percent in :elapsedfull",
    total = nbPop, clear = FALSE, width= 100)
  # identifier for population
  PopID <- 0
  for (i in 1:nbPop){
    pb$tick()
    DataTable_temp = DataTable %>%
      group_by(ID) %>%
      sample_n(10)
    reflectance_table <- get_reflectance_table(BasePopulations_max$pop[i, ], DataTable_temp)
    SpectralVar <- spectral_variance(reflectance_table)
    raoVal <- compute_rao(ReflTree = DataTable_temp,
                          mean_reflectance = SpectralVar$mean_spectra,
                          Pop = BasePopulations_max$pop[i,])
    Label <- c(Label, "MAX")
    PopID <- PopID + 1
    All_Populations_spectral_diversite[[PopID]] <- data.frame("var_tot_centered" = SpectralVar$var_tot_centered,
                                                              "var_sp_centered" = SpectralVar$var_sp_centered,
                                                              "var_cr_centered" = SpectralVar$var_cr_centered,
                                                              "var_res_centered" = SpectralVar$var_res_centered,
                                                              "Rao" = raoVal,
                                                              "Population" = name_column[2])
  }
  do.call(rbind,All_Populations_spectral_diversite)


  All_Populations_results = cbind(All_Populations_spectral_diversite, All_Populations_diversite)
  All_Populations_results = All_Populations_results[, -6]
  All_Populations_results$var_cr_centered = as.numeric(All_Populations_results$var_cr_centered)
  All_Populations_results$var_res_centered = as.numeric(All_Populations_results$var_res_centered)
  All_Populations_results$var_tot_centered = as.numeric(All_Populations_results$var_tot_centered)
  All_Populations_results$var_sp_centered = as.numeric(All_Populations_results$var_sp_centered)
  All_Populations_results$Rao = as.numeric(All_Populations_results$Rao)
  All_Populations_results$Population = as.factor(All_Populations_results$Population)


  numero_pop = c(1:220)
  All_Populations_results = cbind(All_Populations_results, numero_pop)
  All_Populations_results = cbind(All_Populations_results, Label)
  return(All_Populations_results)
}

compute_all_results_for_mean_spectra_species <- function(BasePopulations_min, BasePopulations_max,
                                                         Populations_mintomax, Populations_maxtomin,
                                                         All_Populations_diversite, DataTable,
                                                         mean_spectra_species,
                                                         name_column = c("BasePopulations_min",
                                                                         "BasePopulations_max",
                                                                         "Populations_mintomax",
                                                                         "Populations_maxtomin")){
  print("############### INITIALISATION ###############")
  All_Populations_spectral_diversite = c()
  DataTable_temp = DataTable %>%
    group_by(ID) %>%
    sample_n(10)
  reflectance_table = get_reflectance_table_for_mean_spectra_species(BasePopulations_min$pop[1, ],
                                                                     DataTable_temp, mean_spectra_species)
  temp = spectral_variance_for_mean_spectra_species(reflectance_table)
  All_Populations_spectral_diversite = cbind(All_Populations_spectral_diversite, temp$var_tot_centered,
                                             temp$var_sp_centered, temp$var_res_centered, name_column[1])
  All_Populations_spectral_diversite = as.data.frame(All_Populations_spectral_diversite)
  colnames(All_Populations_spectral_diversite) = c("var_tot_centered", "var_sp_centered",
                                                   "var_res_centered", "Population")
  Label = c("MIN")
  print("############### BasePopulationMin ###############")
  for (i in 2:10){
    print("SELECTION DE 10 PX PAR COURONNE")
    DataTable_temp = DataTable %>%
      group_by(ID) %>%
      sample_n(10)
    print("FIN SELECTION")
    reflectance_table = get_reflectance_table_for_mean_spectra_species(BasePopulations_min$pop[i, ],
                                                                       DataTable_temp, mean_spectra_species)
    temp = spectral_variance_for_mean_spectra_species(reflectance_table)
    All_Populations_spectral_diversite = rbind(All_Populations_spectral_diversite,
                                               c(temp$var_tot_centered, temp$var_sp_centered,
                                                 temp$var_res_centered, name_column[1])
    )
    Label = c(Label, "MIN")
  }
  print("############### PopulationMintoMax ###############")
  for (i in 1:10){
    for (j in 1:10){
      DataTable_temp = DataTable %>%
        group_by(ID) %>%
        sample_n(10)
      reflectance_table = get_reflectance_table_for_mean_spectra_species(Populations_mintomax[[i]]$pop[j,],
                                                                         DataTable_temp, mean_spectra_species)
      temp = spectral_variance_for_mean_spectra_species(reflectance_table)
      All_Populations_spectral_diversite = rbind(All_Populations_spectral_diversite,
                                                 c(temp$var_tot_centered, temp$var_sp_centered,
                                                   temp$var_res_centered,
                                                   paste(name_column[3], names(Populations_mintomax[i]), sep=''))
      )
      Label = c(Label, "MINTOMAX")
    }
  }
  print("############### BasePopulationMax ###############")
  for (i in 1:10){
    DataTable_temp = DataTable %>%
      group_by(ID) %>%
      sample_n(10)
    reflectance_table = get_reflectance_table_for_mean_spectra_species(BasePopulations_max$pop[i, ],
                                                                       DataTable_temp, mean_spectra_species)
    temp = spectral_variance_for_mean_spectra_species(reflectance_table)
    All_Populations_spectral_diversite = rbind(All_Populations_spectral_diversite,
                                               c(temp$var_tot_centered, temp$var_sp_centered,
                                                 temp$var_res_centered, name_column[2])
    )
    Label = c(Label, "MAX")
  }
  print("############### PopulationMaxtoMin ###############")
  for (i in 1:10){
    for (j in 1:10){
      DataTable_temp = DataTable %>%
        group_by(ID) %>%
        sample_n(10)
      reflectance_table = get_reflectance_table_for_mean_spectra_species(Populations_maxtomin[[i]]$pop[j,],
                                                                         DataTable_temp, mean_spectra_species)
      temp = spectral_variance_for_mean_spectra_species(reflectance_table)
      All_Populations_spectral_diversite = rbind(All_Populations_spectral_diversite,
                                                 c(temp$var_tot_centered, temp$var_sp_centered,
                                                   temp$var_res_centered,
                                                   paste(name_column[4], names(Populations_maxtomin[i]), sep=''))
      )
      Label = c(Label, "MAXTOMIN")
    }
  }

  All_Populations_results = cbind(All_Populations_spectral_diversite, All_Populations_diversite)
  All_Populations_results$var_res_centered = as.numeric(All_Populations_results$var_res_centered)
  All_Populations_results$var_tot_centered = as.numeric(All_Populations_results$var_tot_centered)
  All_Populations_results$var_sp_centered = as.numeric(All_Populations_results$var_sp_centered)
  All_Populations_results$Population = as.factor(All_Populations_results$Population)


  numero_pop = c(1:220)
  All_Populations_results = cbind(All_Populations_results, numero_pop)
  All_Populations_results = cbind(All_Populations_results, Label)
  return(All_Populations_results)
}





