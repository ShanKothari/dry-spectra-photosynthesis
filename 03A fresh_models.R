setwd("C:/Users/Shan Kothari/Dropbox/PostdocProjects/PosadaSpectra/")

library(spectrolab)
library(caret)
library(pls)

fresh_spectra<-readRDS("ProcessedData/fresh_spectra_and_traits.rds")
## unit-vector normalization
fresh_spectra_norm<-normalize(fresh_spectra)

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

ggplot(Asat_val_pred,aes(y=measured,x=val_pred))+
  geom_point(size=2)+geom_smooth(method="lm",se=F)+
  theme_bw()+
  geom_abline(slope=1,intercept=0,linetype="dashed",size=2)+
  theme(text = element_text(size=20),
        legend.position = c(0.8, 0.2))+
  labs(y="Measured",x="Predicted")+
  ggtitle("Predicting Asat")

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
