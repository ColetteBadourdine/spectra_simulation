# ============================================================================= =
# Lib_FeatureSelection.R
# ============================================================================= =
# PROGRAMMERS:
# Jean-Baptiste FERET <jb.feret@teledetection.fr>
# Copyright 2021/06 Jean-Baptiste FERET
# ============================================================================= =
# This Library includes functions dedicated to perform feature selection for 
# tree species classification
# ============================================================================= =


#' This function samples calibration and validation data from a list containing species 
#' name and corresponding number of crowns
#'
#' @param CrownsPerSpecies list. contains species names and corresponding number of crowns
#' @param nbSpecies numeric. number of species randomly selected from CrownsPerSpecies
#' @param nbCal numeric. number of crowns included in calibration dataset
#' @param nbPixperCrownCAL numeric. number of pixels to be selected from each crown of the calibration dataset
#' @param nbVal numeric. number of crowns included in validation dataset. All crowns except CAL if FALSE
#' @param nbPixperCrownCAL numeric. number of pixels to be selected from each crown of the calibration dataset. All pixels if FALSE
#' @param verbose boolean. set TRUE to display selected species
#' 
#'
#' @return list. reflectance annd corresponding species for each pixel
#' @export

SampleCALVAL <- function(CrownsPerSpecies,nbSpecies,nbCal,
                         nbPixperCrownCAL,nbVal=FALSE,
                         nbPixperCrownVAL=FALSE,verbose=FALSE){
  SpeciesCAL <- SpeciesVAL <- ReflCAL <- ReflVAL <- c()
  spnum <- -1
  # random selection of nbSpecies among all species
  selSpecies <- sample(names(CrownsPerSpecies),size = nbSpecies)
  # for each species, select crowns and pixels to produce CAL and VAL datasets
  for (sp in selSpecies){
    # if nbVal = FALSE, then select all crowns which were not used for training
    if (nbVal==FALSE){
      nbValtmp <- CrownsPerSpecies[[sp]]-nbCal
    } else {
      nbValtmp <- nbVal
    }
    spnum <- spnum+1
    if (verbose==TRUE){
      print(sp)
    }
    # get a random selection of CAL and VAL crowns
    selCrowns <- sample(CrownsPerSpecies[[sp]],size = nbCal+nbValtmp)
    CrownCal <- selCrowns[1:nbCal]
    CrownVal <- selCrowns[(nbCal+1):(nbCal+nbValtmp)]
    # identify pixels corresponding to species
    pixSp <- which(TreeID$SPID==sp)
    # identify crown ID
    CrownID <- unique(TreeID$CrownID[pixSp])
    # get nbPixperCrown random pixels from each crown
    for (ccal in CrownID[CrownCal]){
      sel <- sample(which(TreeID$CrownID==ccal),nbPixperCrownCAL)
      SpeciesCAL <- c(SpeciesCAL,TreeID$SPID[sel])
      ReflCAL <- rbind(ReflCAL,Reflectance[sel,])
    }
    for (cval in CrownID[CrownVal]){
      if (nbPixperCrownVAL==FALSE){
        sel <- which(TreeID$CrownID==cval)
      } else {
        sel <- sample(which(TreeID$CrownID==cval),nbPixperCrown)
      }
      SpeciesVAL <- c(SpeciesVAL,TreeID$SPID[sel])
      ReflVAL <- rbind(ReflVAL,Reflectance[sel,])
    }
  }
  SpeciesCAL <- as.factor(SpeciesCAL)
  SpeciesVAL <- as.factor(SpeciesVAL)
  return(list('ReflCAL'=ReflCAL,'SpeciesCAL'=SpeciesCAL,'ReflVAL'=ReflVAL,'SpeciesVAL'=SpeciesVAL))
}

#' This function performs SVM classification on a training set and applies the classifier 
#' on a test set
#'
#' @param trainingSet dataframe. Must include a Y column corresponding to classes of the training samples
#' @param testX dataframe. test data contains same features as trainingSet, except Y
#' @param testY factor. classes corresponding to test samples
#' @param testCAL boolean. set to TRUE if the classifier is applied on training data in addition to test set
#' @param method character. which classification algorithm?
#'
#' @return list. confusion matrix and overall accuracy for test and training dataset: 
#' @importFrom liquidSVM mcSVM
#' @importFrom expandFunctions reset.warnings
#' @importFrom stats predict
#' @importFrom progress progress_bar
#' @importFrom graphics par
#' @importFrom stringr str_split
#' @import dplyr
#' @import ggplot2
#' @export

performClassif <- function(trainingSet,testX,testY,testCAL=FALSE,method='lda'){
  reset.warnings()
  if (method=='SVM'){
    model <- mcSVM(Y ~ ., trainingSet)
    if (!is.null(names(warnings()))>0){
      Msg <- names(warnings())
      reset.warnings()
      ValGamma <- str_split(string = Msg,pattern = 'gamma=')[[1]][2]
      ValLambda <- str_split(string = Msg,pattern = 'lambda=')[[1]][2]
      if (!is.na(as.numeric(ValGamma))){
        message('Adjusting Gamma accordingly')
        ValGamma <- as.numeric(ValGamma)
        if (method=='SVM'){
          model <- mcSVM(Y ~ ., trainingSet,max_gamma=ValGamma)
        }
      }
      if (!is.na(as.numeric(ValLambda))){
        message('Adjusting Lambda accordingly')
        ValLambda <- as.numeric(ValLambda)
        if (method=='SVM'){
          model <- mcSVM(Y ~ ., trainingSet,max_lambda=ValLambda)
        }
      }
    } 
  } else {
    model <- caret::train(Y~., trainingSet, method = method)
  }
  
  predY <- predict(model, testX)
  cm <- as.matrix(table(Actual = testY, Predicted = predY)) # create the confusion matrix
  accuracy <- sum(diag(cm)) / length(testY)
  if (testCAL==TRUE){
    trainY <- trainingSet$Y
    trainingSet$Y <- NULL
    predY <- predict(model, trainingSet)
    cmTrain <- as.matrix(table(Actual = trainY, Predicted = predY)) # create the confusion matrix
    accuracyTrain <- sum(diag(cmTrain)) / length(trainY)
  } else {
    cmTrain <- accuracyTrain <- NULL
  }
  return(list('OverallAcc'=accuracy,'ConfusionMat'=cm,'OverallAccTrain'=accuracyTrain,'ConfusionMatTrain'=cmTrain))
}

#' This function performs Sequential Feature Selection based on classification 
#' 
#' @param FeaturesToSelect numeric. number of features to select
#' @param ReflCAL dataframe. Reflectance of calbration dataset. Includes samples as rows and spectral bands as columns
#' @param SpeciesCAL factor. classes corresponding to calibration samples
#' @param ReflVAL dataframe. Reflectance of validation dataset. Includes samples as rows and spectral bands as columns
#' @param SpeciesVAL factor. classes corresponding to validation samples
#' 
#'
#' @return list. ResultsCAL = overall accuracy for calibration dataset at each step
#'               ResultsVAL = overall accuracy for validation dataset at each step
#'               SelectedBands = name of bands selected ranked following rank of identification
#' @importFrom progress progress_bar
#' @export

perform_SFS_Classif <- function(FeaturesToSelect,ReflCAL,SpeciesCAL,ReflVAL,SpeciesVAL,verbose=TRUE,method='lda'){
  
  # get name of features
  AllBands <- names(ReflCAL)
  # initialize feature selection
  SelectedBands <- c()
  nbSelectedBands <- 0
  Rcal_tmp <- ReflCAL[,FALSE]
  Rval_tmp <- ReflVAL[,FALSE]
  j <- 1
  res <- resTrain <- list()
  # while number of features < targeted number of features 
  while(nbSelectedBands<FeaturesToSelect){
    if (verbose==TRUE){
      message(paste('Selection of Feature #',j,' / ',FeaturesToSelect,sep = ''))
      pb <- progress_bar$new(
        format = "Selecting feature [:bar] :percent in :elapsedfull",
        total = length(AllBands), clear = FALSE, width= 100)
    }
    restmp <- resTraintmp <- c()
    # test all remaining bands 
    for (band in AllBands){
      if (verbose==TRUE){
        pb$tick()
      }
      seltmp <- which(names(ReflCAL)==band)
      dataCAL <- data.frame(cbind(Rcal_tmp,ReflCAL[,seltmp]),'Y' = as.factor(SpeciesCAL))
      dataVAL <- data.frame(cbind(Rval_tmp,ReflVAL[,seltmp]))
      res0 <- performClassif(trainingSet=dataCAL,testX=dataVAL,testY=SpeciesVAL,testCAL = TRUE,method=method)
      restmp <- c(restmp,res0$OverallAcc)
      resTraintmp <- c(resTraintmp,res0$OverallAccTrain)
    }
    res[[j]] <- restmp
    resTrain[[j]] <- resTraintmp
    selBand <- which(res[[j]] == max(res[[j]]))
    # in case several bands are identified
    selBand <- selBand[1]
    if (verbose==TRUE){
      print(c(as.numeric(AllBands[selBand]),max(res[[j]])))
    }
    # update datasets and features identified
    SelectedBands <- c(SelectedBands,AllBands[selBand])
    AllBands <- AllBands[-selBand]
    Rcal_tmp <- cbind(Rcal_tmp,ReflCAL[,selBand])
    Rval_tmp <- cbind(Rval_tmp,ReflVAL[,selBand])
    ReflCAL <- ReflCAL[,-selBand]
    ReflVAL <- ReflVAL[,-selBand]
    nbSelectedBands <- nbSelectedBands+1
    j <- j+1
  }
  return(list('ResultsCAL'=resTrain,'ResultsVAL'=res,'SelectedBands'=SelectedBands))
}

#' This function extracts optinal performances for each step of the SFS
#' 
#' @param ResSFS list. number of features to select
#'
#' @return 
#' @export

get_OA_per_Feature <- function(ResSFS){
  nbFeatures <- length(ResSFS)
  # get optimal performances for each feature added
  OA <- c()
  for (i in 1:nbFeatures){
    OA[i] <- max(ResSFS[[i]])
  }
  return(OA)
}

#' This function runs several instances of Sequential Feature Selection based on SVM classification in parallel
#' 
#' @param CALVAL list. CALVAL dataset
#' @param FeaturesToSelect numeric. number of features to select
#' @param Path_CALVAL character. path where to store results
#'
#' @return list. ResultsCAL = overall accuracy for calibration dataset at each step
#'               ResultsVAL = overall accuracy for validation dataset at each step
#'               SelectedBands = name of bands selected ranked following rank of identification
#' @export

perform_SFS_parallel <- function(CALVAL,FeaturesToSelect,Path_CALVAL,verbose=FALSE,method='lda'){
  
  #--------------- perform sequential feature selection based on SVM classification ---------------
  # perform classification with full spectral data
  dataCAL <- data.frame(CALVAL$ReflCAL,'Y' = CALVAL$SpeciesCAL)
  dataVAL <- data.frame(CALVAL$ReflVAL)
  resFull <- performClassif(trainingSet=dataCAL,testX=dataVAL,testY=CALVAL$SpeciesVAL,testCAL = TRUE,method=method)
  if (verbose ==TRUE){
    print(paste('Overall accuracy Training = ',resFull$OverallAccTrain,sep = ''))
    print(paste('Overall accuracy Test = ',resFull$OverallAcc,sep = ''))
  }
  # perform SFS on SVM classification
  ResultsSFS <- perform_SFS_Classif(FeaturesToSelect = FeaturesToSelect,
                                    ReflCAL = CALVAL$ReflCAL, SpeciesCAL = CALVAL$SpeciesCAL,
                                    ReflVAL = CALVAL$ReflVAL, SpeciesVAL = CALVAL$SpeciesVAL,
                                    method=method)
  listRes <- list('ResultsSFS'=ResultsSFS,'ResultsFULL',resFull)
  listOK <- list.files(Path_CALVAL)
  if (length(listOK)==0){
    savename <- file.path(Path_CALVAL,'CALVAL_1.RData')
    save(listRes,file = savename)
  } else {
    listFiles <- sort(as.numeric(gsub(pattern = '.RData',replacement = '',sapply(strsplit(listOK,split = '_'),tail,1))))
    lastFiles <- max(listFiles)
    NewFile <- lastFiles+1
    savename <- file.path(Path_CALVAL,paste('CALVAL_',NewFile,'.RData',sep = ''))
    save(listRes,file = savename)
  }
  return(listRes)
}

#' #' This function produces a plot to compare 
#' #' 
#' #' @param BandsToSelect numeric. number of features to select
#' #' @param ReflCAL dataframe. Reflectance of calbration dataset. Includes samples as rows and spectral bands as columns
#' #' @param SpeciesCAL factor. classes corresponding to calibration samples
#' #' @param ReflVAL dataframe. Reflectance of validation dataset. Includes samples as rows and spectral bands as columns
#' #' @param SpeciesVAL factor. classes corresponding to validation samples
#' #' 
#' #'
#' #' @return list. ResultsCAL = overall accuracy for calibration dataset at each step
#' #'               ResultsVAL = overall accuracy for validation dataset at each step
#' #'               SelectedBands = name of bands selected ranked following rank of identification
#' #' @export
#' 
#' # plot_resSFS <- function(OA_CAL_SFS,OA_CAL,OA_VAL_SFS,OA_VAL){
#' #   
#' #   
#' # }




#' This function performs feature ranking for multiclass classification using XGBoost
#' based on multi-class adaptation of xgboost found here:
#' https://rstudio-pubs-static.s3.amazonaws.com/233321_06bcdf2c8bc445dbb635740bb44f978b.html
#'
#' @param Xcal dataframe. Calibration data, containing values for each feature. the name of the features is expected to be wavelength of spectral band
#' @param Ycal factor. species corresponding to each sample in Xcal
#' @param Xval dataframe. Validation data, containing values for each feature. the name of the features is expected to be wavelength of spectral band
#' @param Yval factor. species corresponding to each sample in Xval
#' @param nthread numeric. number of threads for parallel computing
#'
#' @return list. ConfusionMatVAL = Confusion Matrix for validation data
#'               FeatureSelect = features ranked by importance
#'               FeatureImportance = importance criterion
#' @import xgboost
#' @export

MultiClass_XGBoost <- function(Xcal, Ycal, Xval, Yval,nthread=1){
  
  # set random seed
  set.seed(710)
  # assuming that name of features corresponds to spectral band, then delete space in feature name, and add WL at the beginning
  names(Xcal) <- paste('WL',sub(pattern = ' ',replacement = '',
                                x = names(Xcal)),sep = '')
  names(Xval) <- paste('WL',sub(pattern = ' ',replacement = '',
                                x = names(Xval)),sep = '')
  # Prepare training data set
  train_label <- as.numeric(Ycal)-1
  train_data <- as.matrix(Xcal)
  train_matrix <- xgb.DMatrix(data = train_data, label = train_label)
  
  # Prepare test data set
  test_label <- as.numeric(Yval)-1
  test_data <- as.matrix(Xval)
  test_matrix <- xgb.DMatrix(data = test_data, label = test_label)
  
  numberOfClasses <- length(unique(train_label))
  xgb_params <- list("objective" = "multi:softprob",
                     "eval_metric" = "mlogloss",
                     "num_class" = numberOfClasses,
                     "nthread" = nthread)
  
  nround    <- 100 # number of XGBoost rounds
  # # Fit cv.nfold * cv.nround XGB models and save OOF predictions
  # cv.nfold  <- 5
  # cv_model <- xgb.cv(params = xgb_params,
  #                    data = train_matrix, 
  #                    nrounds = nround,
  #                    nfold = cv.nfold,
  #                    verbose = FALSE,
  #                    prediction = TRUE)
  
  # OOF_prediction <- data.frame(cv_model$pred) %>%
  #   mutate(max_prob = max.col(., ties.method = "last"),
  #          label = train_label + 1)
  
  # confusionMatrix(factor(OOF_prediction$max_prob),
  #                 factor(OOF_prediction$label),
  #                 mode = "everything")
  
  bst_model <- xgb.train(params = xgb_params,
                         data = train_matrix,
                         nrounds = nround)
  
  # Predict hold-out test set
  test_pred <- predict(bst_model, newdata = test_matrix)
  test_prediction <- matrix(test_pred, nrow = numberOfClasses,
                            ncol=length(test_pred)/numberOfClasses) %>%
    t() %>%
    data.frame() %>%
    mutate(label = test_label + 1,
           max_prob = max.col(., "last"))
  # confusion matrix of test set
  ConfusionMatVAL <- confusionMatrix(factor(test_prediction$max_prob),
                                     factor(test_prediction$label),
                                     mode = "everything")
  # get the feature real names
  names <- (sub(pattern = 'WL',replacement = '',x = colnames(Xcal)))
  # compute feature importance matrix
  importance_matrix = xgb.importance(feature_names = names, model = bst_model)
  # # plot
  # gp <- xgb.ggplot.importance(importance_matrix[1:40,])
  # print(gp)
  
  # get feature #ID and importance
  FeatureSelect <- as.numeric(sub(pattern = 'WL',replacement = '',x = importance_matrix$Feature))
  FeatureImportance <- importance_matrix$Gain
  
  res <- list('ConfusionMatVAL'=ConfusionMatVAL,'FeatureSelect'=FeatureSelect,'FeatureImportance'=FeatureImportance)
  return(res)
  
}











MultiClass_XGBoost_CB <- function(Xcal, Ycal, Xval, Yval){
  
  # set random seed
  set.seed(717)
  # assuming that name of features corresponds to spectral band, then delete space in feature name, and add WL at the beginning
  names(Xcal) <- paste('WL',sub(pattern = ' ',replacement = '',
                                x = names(Xcal)),sep = '')
  names(Xval) <- paste('WL',sub(pattern = ' ',replacement = '',
                                x = names(Xval)),sep = '')
  # Prepare training data set
  train_label <- as.numeric(Ycal)-1
  train_data <- as.matrix(Xcal)
  train_matrix <- xgb.DMatrix(data = train_data, label = train_label)
  
  # Prepare test data set
  test_label <- as.numeric(Yval)-1
  test_data <- as.matrix(Xval)
  test_matrix <- xgb.DMatrix(data = test_data, label = test_label)
  
  numberOfClasses <- length(unique(train_label))
  xgb_params <- list("objective" = "multi:softprob",
                     "eval_metric" = "mlogloss",
                     "num_class" = numberOfClasses)
  
  nround    <- 50 # number of XGBoost rounds
  # # Fit cv.nfold * cv.nround XGB models and save OOF predictions
  # cv.nfold  <- 5
  # cv_model <- xgb.cv(params = xgb_params,
  #                    data = train_matrix, 
  #                    nrounds = nround,
  #                    nfold = cv.nfold,
  #                    verbose = FALSE,
  #                    prediction = TRUE)
  
  # OOF_prediction <- data.frame(cv_model$pred) %>%
  #   mutate(max_prob = max.col(., ties.method = "last"),
  #          label = train_label + 1)
  
  # confusionMatrix(factor(OOF_prediction$max_prob),
  #                 factor(OOF_prediction$label),
  #                 mode = "everything")
  
  bst_model <- xgb.train(params = xgb_params,
                         data = train_matrix,
                         nrounds = nround)
  
  # Predict hold-out test set
  test_pred <- predict(bst_model, newdata = test_matrix)
  test_prediction <- matrix(test_pred, nrow = numberOfClasses,
                            ncol=length(test_pred)/numberOfClasses) %>%
    t() %>%
    data.frame() %>%
    mutate(label = test_label + 1,
           max_prob = max.col(., "last"))
  # confusion matrix of test set
  ConfusionMatVAL <- confusionMatrix(factor(test_prediction$max_prob),
                                     factor(test_prediction$label))
  # get the feature real names
  names <- (sub(pattern = 'WL',replacement = '',x = colnames(Xcal)))
  # compute feature importance matrix
  importance_matrix = xgb.importance(feature_names = names, model = bst_model)
  # # plot
  # gp <- xgb.ggplot.importance(importance_matrix[1:40,])
  # print(gp)
  
  # get feature #ID and importance
  FeatureSelect <- as.numeric(sub(pattern = 'WL',replacement = '',x = importance_matrix$Feature))
  FeatureImportance <- importance_matrix$Gain
  
  res <- list('ConfusionMatVAL'=ConfusionMatVAL,'FeatureSelect'=FeatureSelect,'FeatureImportance'=FeatureImportance)
  return(res)
  
}

#' This function plots results corresponding to overall accuracy for the classification over a set of calval datasets

plot_OA_CALVAL_FullSpectrum_vs_FeatureSelect <- function(FullSpectrum,
                                                         FeatureSelect,
                                                         filename,Labs){
  resFull_vect_CAL <- resFull_vect_VAL <- resXGB_vect_CAL <- resXGB_vect_VAL <- c()
  for (dbnum in 1:length(FullSpectrum)){
    resFull_vect_CAL <- c(resFull_vect_CAL,FullSpectrum[[dbnum]]$OverallAccTrain)
    resFull_vect_VAL <- c(resFull_vect_VAL,FullSpectrum[[dbnum]]$OverallAcc)
    resXGB_vect_CAL <- c(resXGB_vect_CAL,FeatureSelect[[dbnum]]$OverallAccTrain)
    resXGB_vect_VAL <- c(resXGB_vect_VAL,FeatureSelect[[dbnum]]$OverallAcc)
  }
  
  df <- data.frame(x=c(resFull_vect_CAL,resFull_vect_VAL), y=c(resXGB_vect_CAL,resXGB_vect_VAL), 
                   col=gl(2,length(resFull_vect_CAL),labels=c('CAL','VAL')))
  library(ggplot2)
  p1 <- ggplot(data = df, aes(x = x,y = y,color=col))  + 
    geom_point() + theme(aspect.ratio=1) +
    ggtitle('Overall accuracy') +
    xlim(0, 1) + ylim(0, 1) +
    labs(x=Labs[1],y=Labs[2]) +
    theme(plot.title = element_text(size=14,hjust = 0.5),
          legend.position="bottom",legend.title = element_text(color = "white", size = 10),legend.text = element_text(size = 11),
          axis.text = element_text(size=15),
          axis.title.x = element_text(size=16, face="bold"),
          axis.title.y = element_text(size=16, face="bold")) +
    guides(fill=guide_legend(nrow = 2),size = FALSE)
  
  # Add 1:1 line
  p1 <- p1 + geom_abline(slope = 1, intercept = 0,linetype='dashed',size=1.25)
  ggsave(filename, plot = last_plot(), device = "png", path = NULL,
         scale = 1, width = 16, height = 14, units = "cm",
         dpi = 600)
}

plot_OA_CALVAL_FullSpectrum_vs_FeatureSelect_CB <- function(FullSpectrum,
                                                         FeatureSelect,
                                                         filename,Labs){
  resFull_vect_CAL <- resFull_vect_VAL <- resXGB_vect_CAL <- resXGB_vect_VAL <- c()
  resFull_vect_CAL <- c(resFull_vect_CAL,FullSpectrum$OverallAccTrain)
  resFull_vect_VAL <- c(resFull_vect_VAL,FullSpectrum$OverallAcc)
  resXGB_vect_CAL <- c(resXGB_vect_CAL,FeatureSelect$OverallAccTrain)
  resXGB_vect_VAL <- c(resXGB_vect_VAL,FeatureSelect$OverallAcc)
  
  df <- data.frame(x=c(resFull_vect_CAL,resFull_vect_VAL), y=c(resXGB_vect_CAL,resXGB_vect_VAL), 
                   col=gl(2,length(resFull_vect_CAL),labels=c('CAL','VAL')))
  library(ggplot2)
  p1 <- ggplot(data = df, aes(x = x,y = y,color=col))  + 
    geom_point() + theme(aspect.ratio=1) +
    ggtitle('Overall accuracy') +
    xlim(0, 1) + ylim(0, 1) +
    labs(x=Labs[1],y=Labs[2]) +
    theme(plot.title = element_text(size=14,hjust = 0.5),
          legend.position="bottom",legend.title = element_text(color = "white", size = 10),legend.text = element_text(size = 11),
          axis.text = element_text(size=15),
          axis.title.x = element_text(size=16, face="bold"),
          axis.title.y = element_text(size=16, face="bold")) +
    guides(fill=guide_legend(nrow = 2),size = FALSE)
  
  # Add 1:1 line
  p1 <- p1 + geom_abline(slope = 1, intercept = 0,linetype='dashed',size=1.25)
  ggsave(filename, plot = last_plot(), device = "png", path = NULL,
         scale = 1, width = 16, height = 14, units = "cm",
         dpi = 600)
}
