setwd("C:/Users/Shan Kothari/Dropbox/PostdocProjects/PosadaSpectra/")

library(spectrolab)
library(caret)
library(pls)

all_spectra_avg<-readRDS("ProcessedData/all_spectra_avg.rds")
all_spectra_norm_avg<-readRDS("ProcessedData/all_spectra_norm_avg.rds")

###########################################################################
## non-normalized, drop needleleaves

all_spectra_avg<-all_spectra_avg[!meta(all_spectra_avg)$leaf_type=="needleleaf",]

## split training and testing data
train_sample <- createDataPartition(
  y = meta(all_spectra_avg)$Species,
  p = .6,
  list = FALSE
)

spectra_train<-all_spectra_avg[train_sample]
spectra_test<-all_spectra_avg[-train_sample]

## LDMC
LDMC_model<-plsr(meta(spectra_train)$LDMC~as.matrix(spectra_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)

ncomp_LDMC <- selectNcomp(LDMC_model, method = "onesigma", plot = FALSE)

LDMC_val_pred<-data.frame(sample=meta(spectra_test)$sample_name,
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

LMA_val_pred<-data.frame(sample=meta(spectra_test)$sample_name,
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

Asat_val_pred<-data.frame(sample=meta(spectra_test)$sample_name,
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

#########################
## normalized, keep conifers and pubescent species

## split training and testing data
train_sample <- createDataPartition(
  y = meta(all_spectra_norm_avg)$Species,
  p = .6,
  list = FALSE
)

spectra_train<-all_spectra_norm_avg[train_sample]
spectra_test<-all_spectra_norm_avg[-train_sample]

## LDMC
LDMC_model<-plsr(meta(spectra_train)$LDMC~as.matrix(spectra_train),
                 ncomp=30,method = "oscorespls",validation="CV",segments=10)

ncomp_LDMC <- selectNcomp(LDMC_model, method = "onesigma", plot = FALSE)

LDMC_val_pred<-data.frame(sample=meta(spectra_test)$sample_name,
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

LMA_val_pred<-data.frame(sample=meta(spectra_test)$sample_name,
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

Asat_val_pred<-data.frame(sample=meta(spectra_test)$sample_name,
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
