setwd("C:/Users/Shan Kothari/Dropbox/PostdocProjects/PosadaSpectra/")

library(spectrolab)
library(caret)
library(pls)

fresh_spectra<-readRDS("ProcessedData/fresh_spectra_and_traits.rds")
## unit-vector normalization
fresh_spectra_norm<-normalize(fresh_spectra)

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

###########################################################################
## non-normalized, drop needleleaves

fresh_spectra<-fresh_spectra[!meta(fresh_spectra)$leaf_type=="needleleaf",]

## split training and testing data
train_sample <- createDataPartition(
  y = meta(fresh_spectra)$Species,
  p = .6,
  list = FALSE
)

spectra_train<-fresh_spectra[train_sample]
spectra_test<-fresh_spectra[-train_sample]

## LDMC
LDMC_model<-plsr(meta(spectra_train)$LDMC~as.matrix(spectra_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)

ncomp_LDMC <- selectNcomp(LDMC_model, method = "onesigma", plot = FALSE)

LDMC_val_pred<-data.frame(sample=meta(spectra_test)$sample_id,
                          val_pred=predict(LDMC_model,newdata=as.matrix(spectra_test),ncomp=ncomp_LDMC)[,,1],
                          measured=meta(spectra_test)$LDMC)

ggplot(LDMC_val_pred,aes(y=measured,x=val_pred))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0.20,0.5),ylim=c(0.20,0.5))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting LDMC")

## EWT
EWT_model<-plsr(meta(spectra_train)$EWT~as.matrix(spectra_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)

ncomp_EWT <- selectNcomp(EWT_model, method = "onesigma", plot = FALSE)

EWT_val_pred<-data.frame(sample=meta(spectra_test)$sample_id,
                          val_pred=predict(EWT_model,newdata=as.matrix(spectra_test),ncomp=ncomp_EWT)[,,1],
                          measured=meta(spectra_test)$EWT)

ggplot(EWT_val_pred,aes(y=measured,x=val_pred))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0,0.3),ylim=c(0.,0.3))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting EWT")

## LMA
LMA_model<-plsr(meta(spectra_train)$LMA~as.matrix(spectra_train),
                ncomp=30,method = "oscorespls",validation="CV",segments=10)

ncomp_LMA <- selectNcomp(LMA_model, method = "onesigma", plot = FALSE)

LMA_val_pred<-data.frame(sample=meta(spectra_test)$sample_id,
                         val_pred=predict(LMA_model,newdata=as.matrix(spectra_test),ncomp=ncomp_LMA)[,,1],
                         measured=meta(spectra_test)$LMA)

ggplot(LMA_val_pred,aes(y=measured,x=val_pred))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting LMA")

## Asat
Asat_model<-plsr(meta(spectra_train)$Asat~as.matrix(spectra_train),
                ncomp=30,method = "oscorespls",validation="CV",segments=10)

ncomp_Asat <- selectNcomp(Asat_model, method = "onesigma", plot = FALSE)

Asat_val_pred<-data.frame(sample=meta(spectra_test)$sample_id,
                         val_pred=predict(Asat_model,newdata=as.matrix(spectra_test),ncomp=ncomp_Asat)[,,1],
                         measured=meta(spectra_test)$Asat)

mylims <- range(with(Asat_val_pred, c(measured, val_pred)), na.rm = T)
ggplot(Asat_val_pred,aes(y=measured,x=val_pred))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured Asat",x="Predicted Asat")+
  ggtitle("Fresh leaves")+
  coord_cartesian(xlim = mylims, ylim = mylims)

## Vcmax25
Vcmax25_model<-plsr(meta(spectra_train)$Vcmax25~as.matrix(spectra_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)

ncomp_Vcmax25 <- selectNcomp(Vcmax25_model, method = "onesigma", plot = FALSE)

Vcmax25_val_pred<-data.frame(sample=meta(spectra_test)$sample_id,
                          val_pred=predict(Vcmax25_model,newdata=as.matrix(spectra_test),ncomp=ncomp_Vcmax25)[,,1],
                          measured=meta(spectra_test)$Vcmax25)

mylims <- range(with(Vcmax25_val_pred, c(measured, val_pred)), na.rm = T)
ggplot(Vcmax25_val_pred,aes(y=measured,x=val_pred))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting Vcmax25")+
  coord_cartesian(xlim = mylims, ylim = mylims)

## ETR
ETR_model<-plsr(meta(spectra_train)$ETR~as.matrix(spectra_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)

ncomp_ETR <- selectNcomp(ETR_model, method = "onesigma", plot = FALSE)

ETR_val_pred<-data.frame(sample=meta(spectra_test)$sample_id,
                          val_pred=predict(ETR_model,newdata=as.matrix(spectra_test),ncomp=ncomp_ETR)[,,1],
                          measured=meta(spectra_test)$ETR)

mylims <- range(with(ETR_val_pred, c(measured, val_pred)), na.rm = T)
ggplot(ETR_val_pred,aes(y=measured,x=val_pred))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured ETR",x="Predicted ETR")+
  ggtitle("Fresh leaves")+
  coord_cartesian(xlim = mylims, ylim = mylims)


#########################
## normalized, keep needleleaves

## split training and testing data
train_sample <- createDataPartition(
  y = meta(fresh_spectra_norm)$Species,
  p = .6,
  list = FALSE
)

spectra_train<-fresh_spectra_norm[train_sample]
spectra_test<-fresh_spectra_norm[-train_sample]

## LDMC
LDMC_model<-plsr(meta(spectra_train)$LDMC~as.matrix(spectra_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)

ncomp_LDMC <- selectNcomp(LDMC_model, method = "onesigma", plot = FALSE)

LDMC_val_pred<-data.frame(sample=meta(spectra_test)$sample_id,
                          val_pred=predict(LDMC_model,newdata=as.matrix(spectra_test),ncomp=ncomp_LDMC)[,,1],
                          measured=meta(spectra_test)$LDMC)

ggplot(LDMC_val_pred,aes(y=measured,x=val_pred))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  coord_cartesian(xlim=c(0.20,0.5),ylim=c(0.20,0.5))+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting LDMC")

## LMA
LMA_model<-plsr(meta(spectra_train)$LMA~as.matrix(spectra_train),
                ncomp=30,method = "oscorespls",validation="CV",segments=10)

ncomp_LMA <- selectNcomp(LMA_model, method = "onesigma", plot = FALSE)

LMA_val_pred<-data.frame(sample=meta(spectra_test)$sample_id,
                         val_pred=predict(LMA_model,newdata=as.matrix(spectra_test),ncomp=ncomp_LMA)[,,1],
                         measured=meta(spectra_test)$LMA)

ggplot(LMA_val_pred,aes(y=measured,x=val_pred))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting LMA")

## Asat
Asat_model<-plsr(meta(spectra_train)$Asat~as.matrix(spectra_train),
                ncomp=30,method = "oscorespls",validation="CV",segments=10)

ncomp_Asat <- selectNcomp(Asat_model, method = "onesigma", plot = FALSE)

Asat_val_pred<-data.frame(sample=meta(spectra_test)$sample_id,
                         val_pred=predict(Asat_model,newdata=as.matrix(spectra_test),ncomp=ncomp_Asat)[,,1],
                         measured=meta(spectra_test)$Asat)

ggplot(Asat_val_pred,aes(y=measured,x=val_pred))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting Asat")

## ETR
ETR_model<-plsr(meta(spectra_train)$ETR~as.matrix(spectra_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)

ncomp_ETR <- selectNcomp(ETR_model, method = "onesigma", plot = FALSE)

ETR_val_pred<-data.frame(sample=meta(spectra_test)$sample_id,
                          val_pred=predict(ETR_model,newdata=as.matrix(spectra_test),ncomp=ncomp_ETR)[,,1],
                          measured=meta(spectra_test)$ETR)

ggplot(ETR_val_pred,aes(y=measured,x=val_pred))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting ETR")

###############################################
## trying repeated nested cross-validation

repeats <- 50    # repeats of nested CV
outer_folds <- 5       # number of outer folds
max_comps <- 20    # maximum PLS components to test
nsamp <- dim(fresh_spectra)[1]

# to store results from every single outer fold across all repeats
all_r2 <- matrix(NA, nrow = repeats, ncol = outer_folds)
all_rmse <- matrix(NA, nrow = repeats, ncol = outer_folds)
all_ncomp <- matrix(NA, nrow = repeats, ncol = outer_folds)
# matrix to store predictions for each sample
pred_matrix <- matrix(NA, nrow = nsamp, ncol = repeats)

# outside loop that gets repeated
for(r in 1:repeats) {
  
  # create new random folds for this repetition
  folds_outer <- sample(cut(seq(1, nsamp), breaks = outer_folds, labels = FALSE))
  
  for(i in 1:outer_folds) {
    # the ith fold is reserved for testing
    test_idx <- which(folds_outer == i)
    train_data <- fresh_spectra[-test_idx, ]
    test_data  <- fresh_spectra[test_idx, ]
    
    # plsr built-in CV acts as our inner loop
    inner_model <- plsr(meta(train_data)$Asat ~ as.matrix(train_data),
                        ncomp = max_comps, 
                        validation = "CV", 
                        segments = 10)
    
    # select n comp from inner loop with one-sigma rule
    # to use for prediction in outer loop
    best_ncomp <- selectNcomp(inner_model, method = "onesigma", plot = FALSE)
    if(best_ncomp == 0) best_ncomp <- 1 ## minimum of one component
    
    # outer loop evaluation
    preds <- predict(inner_model, newdata = as.matrix(test_data), ncomp = best_ncomp)
    # save predictions for these test data in this repeat
    pred_matrix[test_idx, r] <- preds
    
    # calculate R-squared (or RMSEP)
    # there are NAs in the measured values so we just drop those
    all_r2[r,i] <- cor(preds, meta(test_data)$Asat,
                           use="complete.obs")^2
    all_rmse[r,i] <- RMSD(preds,meta(test_data)$Asat)
    all_ncomp[r,i] <- best_ncomp
  }
  message(sprintf("Completed Repetition %d/%d", r, repeats))
}

# 4. Final Performance Summary
cat("\n--- Final Unbiased Evaluation ---\n")
cat(sprintf("Median R2: %.3f\n", median(all_r2)))
cat(sprintf("95%% quantiles: %.3f - %.3f\n", 
            quantile(all_r2,probs = 0.025), 
            quantile(all_r2,probs = 0.975)))

plot_df <- data.frame(
  measured = meta(fresh_spectra)$Asat,
  pred_mean = rowMeans(pred_matrix),
  pred_sd   = apply(pred_matrix, 1, sd)
)

ggplot(plot_df, aes(y = measured, x = pred_mean)) +
  geom_errorbar(aes(xmin = pred_mean - pred_sd, 
                    xmax = pred_mean + pred_sd), 
                alpha = 0.3, color = "gray") +
  geom_point(size = 2, color="blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size=2) +
  theme_minimal() +
  labs(y = "Measured Asat",
       x = "Predicted Asat")
  # annotate("text", x = min(y), y = max(y), 
  #          label = paste("Avg R2 =", round(mean(cor(pred_matrix, y)^2), 3)),
  #          hjust = 0, vjust = 1, fontface = "italic")