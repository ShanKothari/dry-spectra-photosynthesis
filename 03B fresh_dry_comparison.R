setwd("C:/Users/Shan Kothari/Dropbox/PostdocProjects/PosadaSpectra/")
source("Scripts/dry-spectra-photosynthesis/00 useful_functions.R")

library(spectrolab)
library(caret)
library(pls)
library(ggplot2)
library(patchwork)

###############################################
## read in data
## define parameters for all repeated nested cross-validation

# only broadleaf
# only 400-900 nm (need to trim dry)
fresh_spectra<-readRDS("ProcessedData/fresh_spectra_and_traits.rds")
fresh_spectra_sub<-fresh_spectra[!meta(fresh_spectra)$leaf_type %in% c("needleleaf","scaleleaf"),]

dry_spectra<-readRDS("ProcessedData/dry_spectra_and_traits.rds")
dry_spectra_sub<-dry_spectra[!meta(dry_spectra)$leaf_type %in% c("needleleaf","scaleleaf"),400:900]

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
  labs(y = expression(paste("Measured ",italic(A[sat])," (",mu,"mol ",CO[2]," ",m^-2," ",s^-1,")")),
       x = expression(paste("Predicted ",italic(A[sat])," (",mu,"mol ",CO[2]," ",m^-2," ",s^-1,")")))
  ggtitle("Fresh")

Asat_dry_plot <- ggplot(Asat_dry_plot_df,
                        aes(y = measured, x = pred_mean)) +
  theme_bw()+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_errorbar(aes(xmin = pred_mean - pred_sd, 
                    xmax = pred_mean + pred_sd), 
                alpha = 0.3, color = "darkslategray") +
  theme(text = element_text(size=20),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())+
  coord_cartesian(xlim=Asat_lims,ylim=Asat_lims)+
  labs(y = expression(paste("Measured ",italic(A[sat])," (",mu,"mol ",CO[2]," ",m^-2," ",s^-1,")")),
       x = expression(paste("Predicted ",italic(A[sat])," (",mu,"mol ",CO[2]," ",m^-2," ",s^-1,")")))
  ggtitle("Dried")

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
  labs(y = expression(paste("Measured ",italic(V[cmax25])," (",mu,"mol ",CO[2]," ",m^-2," ",s^-1,")")),
       x = expression(paste("Predicted ",italic(V[cmax25])," (",mu,"mol ",CO[2]," ",m^-2," ",s^-1,")")))

Vcmax25_dry_plot <- ggplot(Vcmax25_dry_plot_df,
                           aes(y = measured, x = pred_mean)) +
  theme_bw()+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_errorbar(aes(xmin = pred_mean - pred_sd, 
                    xmax = pred_mean + pred_sd), 
                alpha = 0.3, color = "darkslategray") +
  theme(text = element_text(size=20),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())+
  coord_cartesian(xlim=Vcmax25_lims,ylim=Vcmax25_lims)+
  labs(y = expression(paste("Measured ",italic(V[cmax25])," (",mu,"mol ",CO[2]," ",m^-2," ",s^-1,")")),
       x = expression(paste("Predicted ",italic(V[cmax25])," (",mu,"mol ",CO[2]," ",m^-2," ",s^-1,")")))

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
  labs(y = expression(paste("Measured ETR (",mu,"mol ",m^-2," ",s^-1,")")),
       x = expression(paste("Predicted ETR (",mu,"mol ",m^-2," ",s^-1,")")))

ETR_dry_plot <- ggplot(ETR_dry_plot_df,
                       aes(y = measured, x = pred_mean)) +
  theme_bw()+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_errorbar(aes(xmin = pred_mean - pred_sd, 
                    xmax = pred_mean + pred_sd), 
                alpha = 0.3, color = "darkslategray") +
  theme(text = element_text(size=20),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())+
  coord_cartesian(xlim=ETR_lims,ylim=ETR_lims)+
  labs(y = expression(paste("Measured ETR (",mu,"mol ",m^-2," ",s^-1,")")),
       x = expression(paste("Predicted ETR (",mu,"mol ",m^-2," ",s^-1,")")))

#########################################
## Rd25

## fresh
Rd25_fresh_preds<-plsr_rnCV(repeats = repeats,
                            outer_folds = outer_folds,
                            max_comps = max_comps,
                            yvar=meta(fresh_spectra_sub)$Rd25,
                            xmat=as.matrix(fresh_spectra_sub))

Rd25_fresh_plot_df <- data.frame(
  measured = meta(fresh_spectra_sub)$Rd25,
  pred_mean = rowMeans(Rd25_fresh_preds$pred_matrix),
  pred_sd   = apply(Rd25_fresh_preds$pred_matrix, 1, sd)
)

## dry
Rd25_dry_preds<-plsr_rnCV(repeats = repeats,
                          outer_folds = outer_folds,
                          max_comps = max_comps,
                          yvar=meta(dry_spectra_sub)$Rd25,
                          xmat=as.matrix(dry_spectra_sub))

Rd25_dry_plot_df <- data.frame(
  measured = meta(dry_spectra_sub)$Rd25,
  pred_mean = rowMeans(Rd25_dry_preds$pred_matrix),
  pred_sd   = apply(Rd25_dry_preds$pred_matrix, 1, sd)
)

Rd25_lims <- define_lims_comparison(Rd25_dry_plot_df,Rd25_fresh_plot_df)

Rd25_fresh_plot <- ggplot(Rd25_fresh_plot_df,
                          aes(y = measured, x = pred_mean)) +
  theme_bw()+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_errorbar(aes(xmin = pred_mean - pred_sd, 
                    xmax = pred_mean + pred_sd), 
                alpha = 0.3, color = "darkslategray") +
  theme(text = element_text(size=20))+
  coord_cartesian(xlim=Rd25_lims,ylim=Rd25_lims)+
  labs(y = expression(paste("Measured ",italic(R[d25])," (",mu,"mol ",CO[2]," ",m^-2," ",s^-1,")")),
       x = expression(paste("Predicted ",italic(R[d25])," (",mu,"mol ",CO[2]," ",m^-2," ",s^-1,")")))

Rd25_dry_plot <- ggplot(Rd25_dry_plot_df,
                        aes(y = measured, x = pred_mean)) +
  theme_bw()+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_errorbar(aes(xmin = pred_mean - pred_sd, 
                    xmax = pred_mean + pred_sd), 
                alpha = 0.3, color = "darkslategray") +
  theme(text = element_text(size=20),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())+
  coord_cartesian(xlim=Rd25_lims,ylim=Rd25_lims)+
  labs(y = expression(paste("Measured ",italic(R[d25])," (",mu,"mol ",CO[2]," ",m^-2," ",s^-1,")")),
       x = expression(paste("Predicted ",italic(R[d25])," (",mu,"mol ",CO[2]," ",m^-2," ",s^-1,")")))

#########################################
## LMA

## fresh
LMA_fresh_preds<-plsr_rnCV(repeats = repeats,
                           outer_folds = outer_folds,
                           max_comps = max_comps,
                           yvar=meta(fresh_spectra_sub)$LMA,
                           xmat=as.matrix(fresh_spectra_sub))

LMA_fresh_plot_df <- data.frame(
  measured = meta(fresh_spectra_sub)$LMA,
  pred_mean = rowMeans(LMA_fresh_preds$pred_matrix),
  pred_sd   = apply(LMA_fresh_preds$pred_matrix, 1, sd)
)

## dry
LMA_dry_preds<-plsr_rnCV(repeats = repeats,
                         outer_folds = outer_folds,
                         max_comps = max_comps,
                         yvar=meta(dry_spectra_sub)$LMA,
                         xmat=as.matrix(dry_spectra_sub))

LMA_dry_plot_df <- data.frame(
  measured = meta(dry_spectra_sub)$LMA,
  pred_mean = rowMeans(LMA_dry_preds$pred_matrix),
  pred_sd   = apply(LMA_dry_preds$pred_matrix, 1, sd)
)

LMA_lims <- define_lims_comparison(LMA_dry_plot_df,LMA_fresh_plot_df)

LMA_fresh_plot <- ggplot(LMA_fresh_plot_df,
                         aes(y = measured, x = pred_mean)) +
  theme_bw()+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_errorbar(aes(xmin = pred_mean - pred_sd, 
                    xmax = pred_mean + pred_sd), 
                alpha = 0.3, color = "darkslategray") +
  theme(text = element_text(size=20))+
  coord_cartesian(xlim=LMA_lims,ylim=LMA_lims)+
  labs(y=expression("Measured LMA (kg m"^-2*")"),
       x=expression("Predicted LMA (kg m"^-2*")"))+
  ggtitle("Fresh")

LMA_dry_plot <- ggplot(LMA_dry_plot_df,
                       aes(y = measured, x = pred_mean)) +
  theme_bw()+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_errorbar(aes(xmin = pred_mean - pred_sd, 
                    xmax = pred_mean + pred_sd), 
                alpha = 0.3, color = "darkslategray") +
  theme(text = element_text(size=20),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())+
  coord_cartesian(xlim=LMA_lims,ylim=LMA_lims)+
  labs(y=expression("Measured LMA (kg m"^-2*")"),
       x=expression("Predicted LMA (kg m"^-2*")"))+
  ggtitle("Dried")


#########################################
## EWT

## fresh
EWT_fresh_preds<-plsr_rnCV(repeats = repeats,
                           outer_folds = outer_folds,
                           max_comps = max_comps,
                           yvar=meta(fresh_spectra_sub)$EWT,
                           xmat=as.matrix(fresh_spectra_sub))

EWT_fresh_plot_df <- data.frame(
  measured = meta(fresh_spectra_sub)$EWT,
  pred_mean = rowMeans(EWT_fresh_preds$pred_matrix),
  pred_sd   = apply(EWT_fresh_preds$pred_matrix, 1, sd)
)

## dry
EWT_dry_preds<-plsr_rnCV(repeats = repeats,
                         outer_folds = outer_folds,
                         max_comps = max_comps,
                         yvar=meta(dry_spectra_sub)$EWT,
                         xmat=as.matrix(dry_spectra_sub))

EWT_dry_plot_df <- data.frame(
  measured = meta(dry_spectra_sub)$EWT,
  pred_mean = rowMeans(EWT_dry_preds$pred_matrix),
  pred_sd   = apply(EWT_dry_preds$pred_matrix, 1, sd)
)

EWT_lims <- define_lims_comparison(EWT_dry_plot_df,EWT_fresh_plot_df)

EWT_fresh_plot <- ggplot(EWT_fresh_plot_df,
                         aes(y = measured, x = pred_mean)) +
  theme_bw()+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_errorbar(aes(xmin = pred_mean - pred_sd, 
                    xmax = pred_mean + pred_sd), 
                alpha = 0.3, color = "darkslategray") +
  theme(text = element_text(size=20))+
  coord_cartesian(xlim=EWT_lims,ylim=EWT_lims)+
  labs(y = "Measured EWT (mm)",
       x = "Predicted EWT (mm)")

EWT_dry_plot <- ggplot(EWT_dry_plot_df,
                       aes(y = measured, x = pred_mean)) +
  theme_bw()+
  geom_point(size=2)+
  geom_smooth(method="lm",se=F)+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  geom_errorbar(aes(xmin = pred_mean - pred_sd, 
                    xmax = pred_mean + pred_sd), 
                alpha = 0.3, color = "darkslategray") +
  theme(text = element_text(size=20),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())+
  coord_cartesian(xlim=EWT_lims,ylim=EWT_lims)+
  labs(y = "Measured EWT (mm)",
       x = "Predicted EWT (mm)")

##########################################
## plotting together

pdf("Results/FigXX.pdf",height=20,width=11)
(Asat_fresh_plot/Vcmax25_fresh_plot/ETR_fresh_plot/Rd25_fresh_plot)|
  (Asat_dry_plot/Vcmax25_dry_plot/ETR_dry_plot/Rd25_dry_plot) +
  plot_layout(guides="collect") & theme(legend.position = "right")
dev.off()

pdf("Results/FigXX2.pdf",height=10,width=11)
(LMA_fresh_plot/EWT_fresh_plot)|
  (LMA_dry_plot/EWT_dry_plot) +
  plot_layout(guides="collect") & theme(legend.position = "right")
dev.off()

############################################
## output tables of results


