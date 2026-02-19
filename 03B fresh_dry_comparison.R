setwd("C:/Users/Shan Kothari/Dropbox/PostdocProjects/PosadaSpectra/")
source("Scripts/dry-spectra-photosynthesis/00 useful_functions.R")

library(spectrolab)
library(caret)
library(pls)

###############################################
## read in data
## define parameters for all repeated nested cross-validation

fresh_spectra<-readRDS("ProcessedData/fresh_spectra_and_traits.rds")
fresh_spectra_sub<-fresh_spectra[!meta(fresh_spectra)$leaf_type %in% c("needleleaf","scaleleaf"),]

repeats <- 50    # repeats of nested CV
outer_folds <- 5       # number of outer folds
max_comps <- 20    # maximum PLS components to test

#########################################
## apply to individual traits

## Asat
Asat_fresh_preds<-plsr_rnCV(repeats = repeats,
                            outer_folds = outer_folds,
                            max_comps = max_comps,
                            yvar=meta(fresh_spectra_sub)$Asat,
                            xmat=as.matrix(fresh_spectra_sub))

Asat_fresh_plot_df <- data.frame(
  measured = meta(fresh_spectra_sub)$Asat,
  pred_mean = rowMeans(Asat_fresh_preds$pred_matrix),
  pred_sd   = apply(Asat_fresh_preds$pred_matrix, 1, sd)
)

Asat_fresh_lims <- define_lims(Asat_fresh_plot_df)

ggplot(Asat_fresh_plot_df, aes(y = measured, x = pred_mean)) +
  theme_bw()+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_errorbar(aes(xmin = pred_mean - pred_sd, 
                    xmax = pred_mean + pred_sd), 
                alpha = 0.3, color = "darkslategray") +
  theme(text = element_text(size=20))+
  coord_cartesian(xlim=Asat_fresh_lims,ylim=Asat_fresh_lims)+
  labs(y = "Measured Asat",
       x = "Predicted Asat")

## ETR
ETR_fresh_preds<-plsr_rnCV(repeats = repeats,
                           outer_folds = outer_folds,
                           max_comps = max_comps,
                           yvar=meta(fresh_spectra_sub)$ETR,
                           xmat=as.matrix(fresh_spectra_sub))

ETR_fresh_plot_df <- data.frame(
  measured = meta(fresh_spectra_sub)$ETR,
  pred_mean = rowMeans(ETR_fresh_preds$pred_matrix),
  pred_sd   = apply(ETR_fresh_preds$pred_matrix, 1, sd)
)

ETR_fresh_lims <- define_lims(ETR_fresh_plot_df)

ggplot(ETR_fresh_plot_df, aes(y = measured, x = pred_mean)) +
  theme_bw()+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_errorbar(aes(xmin = pred_mean - pred_sd, 
                    xmax = pred_mean + pred_sd), 
                alpha = 0.3, color = "darkslategray") +
  theme(text = element_text(size=20))+
  coord_cartesian(xlim=ETR_fresh_lims,ylim=ETR_fresh_lims)+
  labs(y = "Measured ETR",
       x = "Predicted ETR")

## Rd
Rd_fresh_preds<-plsr_rnCV(repeats = repeats,
                          outer_folds = outer_folds,
                          max_comps = max_comps,
                          yvar=meta(fresh_spectra_sub)$Rd,
                          xmat=as.matrix(fresh_spectra_sub))

Rd_fresh_plot_df <- data.frame(
  measured = meta(fresh_spectra_sub)$Rd,
  pred_mean = rowMeans(Rd_fresh_preds$pred_matrix),
  pred_sd   = apply(Rd_fresh_preds$pred_matrix, 1, sd)
)

Rd_fresh_lims <- define_lims(Rd_fresh_plot_df)

ggplot(Rd_fresh_plot_df, aes(y = measured, x = pred_mean)) +
  theme_bw()+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_errorbar(aes(xmin = pred_mean - pred_sd, 
                    xmax = pred_mean + pred_sd), 
                alpha = 0.3, color = "darkslategray") +
  theme(text = element_text(size=20))+
  coord_cartesian(xlim=Rd_fresh_lims,ylim=Rd_fresh_lims)+
  labs(y = "Measured Rd",
       x = "Predicted Rd")
