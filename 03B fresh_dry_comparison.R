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

dry_spectra<-readRDS("ProcessedData/dry_spectra_and_traits.rds")
dry_spectra_sub<-dry_spectra[!meta(dry_spectra)$leaf_type %in% c("needleleaf","scaleleaf"),]

repeats <- 50    # repeats of nested CV
outer_folds <- 5       # number of outer folds
max_comps <- 20    # maximum PLS components to test

#########################################
## apply to individual traits - Asat first

## fresh
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

## dry
Asat_dry_preds<-plsr_rnCV(repeats = repeats,
                          outer_folds = outer_folds,
                          max_comps = max_comps,
                          yvar=meta(dry_spectra_sub)$Asat,
                          xmat=as.matrix(dry_spectra_sub))

Asat_dry_plot_df <- data.frame(
  measured = meta(dry_spectra_sub)$Asat,
  pred_mean = rowMeans(Asat_dry_preds$pred_matrix),
  pred_sd   = apply(Asat_dry_preds$pred_matrix, 1, sd)
)

Asat_lims <- define_lims_comparison(Asat_dry_plot_df,Asat_fresh_plot_df)

Asat_fresh_plot <- ggplot(Asat_fresh_plot_df,
                          aes(y = measured, x = pred_mean)) +
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

Asat_dry_plot <- ggplot(Asat_dry_plot_df,
                        aes(y = measured, x = pred_mean)) +
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

#########################################
## Vcmax25

## fresh
Vcmax25_fresh_preds<-plsr_rnCV(repeats = repeats,
                            outer_folds = outer_folds,
                            max_comps = max_comps,
                            yvar=meta(fresh_spectra_sub)$Vcmax25,
                            xmat=as.matrix(fresh_spectra_sub))

Vcmax25_fresh_plot_df <- data.frame(
  measured = meta(fresh_spectra_sub)$Vcmax25,
  pred_mean = rowMeans(Vcmax25_fresh_preds$pred_matrix),
  pred_sd   = apply(Vcmax25_fresh_preds$pred_matrix, 1, sd)
)

## dry
Vcmax25_dry_preds<-plsr_rnCV(repeats = repeats,
                             outer_folds = outer_folds,
                             max_comps = max_comps,
                             yvar=meta(dry_spectra_sub)$Vcmax25,
                             xmat=as.matrix(dry_spectra_sub))

Vcmax25_dry_plot_df <- data.frame(
  measured = meta(dry_spectra_sub)$Vcmax25,
  pred_mean = rowMeans(Vcmax25_dry_preds$pred_matrix),
  pred_sd   = apply(Vcmax25_dry_preds$pred_matrix, 1, sd)
)

Vcmax25_lims <- define_lims_comparison(Vcmax25_dry_plot_df,Vcmax25_fresh_plot_df)

Vcmax25_fresh_plot <- ggplot(Vcmax25_fresh_plot_df,
                             aes(y = measured, x = pred_mean)) +
  theme_bw()+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_errorbar(aes(xmin = pred_mean - pred_sd, 
                    xmax = pred_mean + pred_sd), 
                alpha = 0.3, color = "darkslategray") +
  theme(text = element_text(size=20))+
  coord_cartesian(xlim=Vcmax25_lims,ylim=Vcmax25_lims)+
  labs(y = "Measured Vcmax25",
       x = "Predicted Vcmax25")

Vcmax25_dry_plot <- ggplot(Vcmax25_dry_plot_df,
                           aes(y = measured, x = pred_mean)) +
  theme_bw()+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_errorbar(aes(xmin = pred_mean - pred_sd, 
                    xmax = pred_mean + pred_sd), 
                alpha = 0.3, color = "darkslategray") +
  theme(text = element_text(size=20))+
  coord_cartesian(xlim=Vcmax25_lims,ylim=Vcmax25_lims)+
  labs(y = "Measured Vcmax25",
       x = "Predicted Vcmax25")

#######################################
## ETR

## fresh
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

## dry
ETR_dry_preds<-plsr_rnCV(repeats = repeats,
                         outer_folds = outer_folds,
                         max_comps = max_comps,
                         yvar=meta(dry_spectra_sub)$ETR,
                         xmat=as.matrix(dry_spectra_sub))

ETR_dry_plot_df <- data.frame(
  measured = meta(dry_spectra_sub)$ETR,
  pred_mean = rowMeans(ETR_dry_preds$pred_matrix),
  pred_sd   = apply(ETR_dry_preds$pred_matrix, 1, sd)
)

## plot
ETR_lims <- define_lims_comparison(ETR_dry_plot_df,ETR_fresh_plot_df)

ETR_fresh_plot <- ggplot(ETR_fresh_plot_df,
                         aes(y = measured, x = pred_mean)) +
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

ETR_dry_plot <- ggplot(ETR_dry_plot_df,
                       aes(y = measured, x = pred_mean)) +
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

##################################
## Rd

## fresh
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

## dry
Rd_dry_preds<-plsr_rnCV(repeats = repeats,
                          outer_folds = outer_folds,
                          max_comps = max_comps,
                          yvar=meta(dry_spectra_sub)$Rd,
                          xmat=as.matrix(dry_spectra_sub))

Rd_dry_plot_df <- data.frame(
  measured = meta(dry_spectra_sub)$Rd,
  pred_mean = rowMeans(Rd_dry_preds$pred_matrix),
  pred_sd   = apply(Rd_dry_preds$pred_matrix, 1, sd)
)

Rd_lims <- define_lims_comparison(Rd_dry_plot_df,Rd_fresh_plot_df)

Rd_fresh_lims <- define_lims(Rd_fresh_plot_df)

Rd_fresh_plot <- ggplot(Rd_fresh_plot_df,
                        aes(y = measured, x = pred_mean)) +
  theme_bw()+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_errorbar(aes(xmin = pred_mean - pred_sd, 
                    xmax = pred_mean + pred_sd), 
                alpha = 0.3, color = "darkslategray") +
  theme(text = element_text(size=20))+
  coord_cartesian(xlim=Rd_lims,ylim=Rd_lims)+
  labs(y = "Measured Rd",
       x = "Predicted Rd")

Rd_dry_plot <- ggplot(Rd_dry_plot_df,
                      aes(y = measured, x = pred_mean)) +
  theme_bw()+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_errorbar(aes(xmin = pred_mean - pred_sd, 
                    xmax = pred_mean + pred_sd), 
                alpha = 0.3, color = "darkslategray") +
  theme(text = element_text(size=20))+
  coord_cartesian(xlim=Rd_lims,ylim=Rd_lims)+
  labs(y = "Measured Rd",
       x = "Predicted Rd")
