setwd("C:/Users/Shan Kothari/Dropbox/PostdocProjects/PosadaSpectra/")
source("Scripts/dry-spectra-photosynthesis/00 useful_functions.R")

library(spectrolab)
library(caret)
library(pls)

##################################################################
## read data and set parameters

dry_spectra<-readRDS("ProcessedData/dry_spectra_and_traits.rds")

repeats <- 50    # repeats of nested CV
outer_folds <- 5       # number of outer folds
max_comps <- 20    # maximum PLS components to test

######################################
## apply repeated nested CV to each trait

## Asat
Asat_preds<-plsr_rnCV(repeats = repeats,
                      outer_folds = outer_folds,
                      max_comps = max_comps,
                      yvar=meta(dry_spectra)$Asat,
                      xmat=as.matrix(dry_spectra))

Asat_plot_df <- data.frame(
  measured = meta(dry_spectra)$Asat,
  pred_mean = rowMeans(Asat_preds$pred_matrix),
  pred_sd   = apply(Asat_preds$pred_matrix, 1, sd)
)

Asat_lims <- define_lims(Asat_plot_df)

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
                      yvar=meta(dry_spectra)$ETR,
                      xmat=as.matrix(dry_spectra))

ETR_plot_df <- data.frame(
  measured = meta(dry_spectra)$ETR,
  pred_mean = rowMeans(ETR_preds$pred_matrix),
  pred_sd   = apply(ETR_preds$pred_matrix, 1, sd)
)

ETR_lims <- define_lims(ETR_plot_df)

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

## Rd
Rd_preds<-plsr_rnCV(repeats = repeats,
                     outer_folds = outer_folds,
                     max_comps = max_comps,
                     yvar=meta(dry_spectra)$Rd,
                     xmat=as.matrix(dry_spectra))

Rd_plot_df <- data.frame(
  measured = meta(dry_spectra)$Rd,
  pred_mean = rowMeans(Rd_preds$pred_matrix),
  pred_sd   = apply(Rd_preds$pred_matrix, 1, sd)
)

Rd_lims <- define_lims(Rd_plot_df)

ggplot(Rd_plot_df, aes(y = measured, x = pred_mean)) +
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
