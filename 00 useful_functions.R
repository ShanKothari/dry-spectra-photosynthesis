########################################
## define functions

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

## repeated nested CV
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
      # here we divide by the full range
      all_perrmse[r,i] <- percentRMSD(preds, test_yvar, min=0, max=1)
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