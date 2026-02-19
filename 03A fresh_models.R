setwd("C:/Users/Shan Kothari/Dropbox/PostdocProjects/PosadaSpectra/")

library(spectrolab)
library(caret)
library(pls)

############################################
## useful functions

## root mean squared deviation
RMSD<-function(measured,predicted){
  not.na<-which(!is.na(measured) & !is.na(predicted))
  return(sqrt(sum((measured-predicted)^2,na.rm=T)/(length(not.na)-1)))
}

## percent RMSD (based on data quantiles)
## set min and max to 0 and 1 for range as denominator
## or to 0.25 and 0.75 for IQR as denominator
percentRMSD<-function(measured,predicted,min,max,na.rm=T){
  RMSD_data<-RMSD(measured,predicted)
  range<-unname(quantile(measured,probs=max,na.rm=na.rm)-quantile(measured,probs=min,na.rm=na.rm))
  return(RMSD_data/range)
}

###############################################
## read in data
## define parameters for all repeated nested cross-validation

fresh_spectra<-readRDS("ProcessedData/fresh_spectra_and_traits.rds")
fresh_spectra_sub<-fresh_spectra[!meta(fresh_spectra)$leaf_type %in% c("needleleaf","scaleleaf"),]

repeats <- 50    # repeats of nested CV
outer_folds <- 5       # number of outer folds
max_comps <- 20    # maximum PLS components to test

## build a function to do the repeated nested CV
plsr_rnCV<-function(repeats, outer_folds, max_comps, yvar, xmat){
  
  if(length(yvar)!=nrow(xmat)){
    stop("different number of samples in predictor matrix and response")
  }
  
  # to store results from every single outer fold across all repeats
  all_r2 <- matrix(NA, nrow = repeats, ncol = outer_folds)
  all_rmse <- matrix(NA, nrow = repeats, ncol = outer_folds)
  all_ncomp <- matrix(NA, nrow = repeats, ncol = outer_folds)
  # matrix to store predictions for each sample
  pred_matrix <- matrix(NA, nrow = nrow(xmat), ncol = repeats)
  
  # outside loop that gets repeated
  for(r in 1:repeats) {
    
    # create new random folds for this repetition
    folds_outer <- sample(cut(seq(1, nrow(xmat)), breaks = outer_folds, labels = FALSE))
    
    for(i in 1:outer_folds) {
      # the ith fold is reserved for testing
      test_idx <- which(folds_outer == i)
      train_xmat <- xmat[-test_idx, ]
      test_xmat  <- xmat[test_idx, ]
      
      train_yvar <- yvar[-test_idx]
      test_yvar <- yvar[test_idx]
      
      # plsr built-in CV acts as our inner loop
      inner_model <- plsr(train_yvar ~ train_xmat,
                          ncomp = max_comps, 
                          validation = "CV", 
                          segments = 10)
      
      # select n comp from inner loop with one-sigma rule
      # to use for prediction in outer loop
      best_ncomp <- selectNcomp(inner_model, method = "onesigma", plot = FALSE)
      if(best_ncomp == 0) best_ncomp <- 1 ## minimum of one component
      
      # outer loop evaluation
      preds <- predict(inner_model, newdata = test_xmat, ncomp = best_ncomp)
      # save predictions for these test data in this repeat
      pred_matrix[test_idx, r] <- preds
      
      # calculate R-squared (or RMSEP)
      # there are NAs in the measured values so we just drop those
      all_r2[r,i] <- cor(preds, test_yvar,
                         use="complete.obs")^2
      all_rmse[r,i] <- RMSD(preds, test_yvar)
      all_ncomp[r,i] <- best_ncomp
    }
    
    # track progress
    message(sprintf("Completed Repetition %d/%d", r, repeats))
    
  }
  
  return(list(pred_matrix=pred_matrix,
              all_r2=all_r2,
              all_rmse=all_rmse,
              all_ncomp=all_ncomp))
}

#########################################
## apply to individual traits

## Asat
Asat_preds<-plsr_rnCV(repeats = repeats,
                      outer_folds = outer_folds,
                      max_comps = max_comps,
                      yvar=meta(fresh_spectra_sub)$Asat,
                      xmat=as.matrix(fresh_spectra_sub))

Asat_plot_df <- data.frame(
  measured = meta(fresh_spectra_sub)$Asat,
  pred_mean = rowMeans(Asat_preds$pred_matrix),
  pred_sd   = apply(Asat_preds$pred_matrix, 1, sd)
)

Asat_lims <- c(min(c(Asat_plot_df$measured,
                     Asat_plot_df$pred_mean-Asat_plot_df$pred_sd),na.rm=T),
               max(c(Asat_plot_df$measured,
                     Asat_plot_df$pred_mean+Asat_plot_df$pred_sd),na.rm=T))

ggplot(Asat_plot_df, aes(y = measured, x = pred_mean)) +
  theme_bw()+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_errorbar(aes(xmin = pred_mean - pred_sd, 
                    xmax = pred_mean + pred_sd), 
                alpha = 0.3, color = "darkslategray") +
  theme(text = element_text(size=20))+
  coord_cartesian(xlim=Asat_lims,ylim=Asat_lims)+
  labs(y = "Measured Asat",
       x = "Predicted Asat")

## ETR
ETR_preds<-plsr_rnCV(repeats = repeats,
                     outer_folds = outer_folds,
                     max_comps = max_comps,
                     yvar=meta(fresh_spectra_sub)$ETR,
                     xmat=as.matrix(fresh_spectra_sub))

ETR_plot_df <- data.frame(
  measured = meta(fresh_spectra_sub)$ETR,
  pred_mean = rowMeans(ETR_preds$pred_matrix),
  pred_sd   = apply(ETR_preds$pred_matrix, 1, sd)
)

ETR_lims <- c(min(c(ETR_plot_df$measured,
                    ETR_plot_df$pred_mean-ETR_plot_df$pred_sd),na.rm=T),
              max(c(ETR_plot_df$measured,
                    ETR_plot_df$pred_mean+ETR_plot_df$pred_sd),na.rm=T))

ggplot(ETR_plot_df, aes(y = measured, x = pred_mean)) +
  theme_bw()+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_errorbar(aes(xmin = pred_mean - pred_sd, 
                    xmax = pred_mean + pred_sd), 
                alpha = 0.3, color = "darkslategray") +
  theme(text = element_text(size=20))+
  coord_cartesian(xlim=ETR_lims,ylim=ETR_lims)+
  labs(y = "Measured ETR",
       x = "Predicted ETR")
