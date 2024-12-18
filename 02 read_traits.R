setwd("C:/Users/Shan Kothari/Dropbox/PostdocProjects/PosadaSpectra/")

library(spectrolab)
library(dplyr)

########################################
## trait and summary data

## read information from summary data
data_summary<-read.csv("Traits/data_summary.csv")
data_summary$Species<-trimws(data_summary$Species, which = "both")
## get rid of the extra Picea M8
data_summary<-data_summary[!duplicated(data_summary$ID),]

## binary classification as sun/shade
data_summary$Type[data_summary$ID=="ACSA1 Real Int"]<-"Shade"
data_summary$Type[data_summary$ID=="ACSA2 Real Int"]<-"Sun"
data_summary$Type[data_summary$Type=="Int/Shade"]<-"Shade"

area_data<-read.csv("Traits/leaf_area_11_19.csv")
## get rid of the extra Picea M8
area_data<-area_data[!duplicated(area_data$SpectraName),]

area_match_fresh<-match(data_summary$SpectraName,area_data$SpectraName)
## sq in to sq cm
data_summary$area<-area_data$Area..in2.[area_match_fresh]*2.54^2
## g per sq cm to g per sq m
data_summary$LMA<-data_summary$Dry.weight..g./data_summary$area*10000

## classifications of species into groups
needleleaf<-c("Abies balsamea","Larix laricina","Larix decidua",
              "Picea abies","Picea glauca","Picea omorika",
              "Picea rubens","Pinus strobus")

scaleleaf<-"Thuja occidentalis"
data_summary$leaf_type<-"broadleaf"
data_summary$leaf_type[data_summary$Species %in% needleleaf]<-"needleleaf"
data_summary$leaf_type[data_summary$Species %in% scaleleaf]<-"scaleleaf"

pubescent<-c("Rhododendron degronianum","Rhododendron makinoi",
             "Rhododendron maximum","Rhododendron Paul Bosley")
data_summary$pubescent<-FALSE
data_summary$pubescent[data_summary$Species %in% pubescent]<-TRUE

evergreen<-c("Abies balsamea","Picea abies","Picea glauca","Picea omorika",
             "Picea rubens","Pinus strobus","Thuja occidentalis",
             "Rhododendron degronianum","Rhododendron makinoi",
             "Rhododendron maximum","Rhododendron Paul Bosley")
data_summary$leaf_habit<-"deciduous"
data_summary$leaf_habit[data_summary$Species %in% evergreen]<-"evergreen"

data_summary$pathway<-"C3"
data_summary$pathway[data_summary$Species=="Zea mays"]<-"C4"

# write.csv(data_summary,"ProcessedData/data_summary_processed.csv",row.names=F)

#############################################
## gas exchange and chlorophyll fluorescence data

## these have Juan's manual corrections
## e.g. for leaf area, for empty chamber measured
## already done

Asat<-read.csv("Traits/Asat_sheet.csv")
lightFvFm<-read.csv("Traits/lightFvFm_sheet.csv")
darkFvFm<-read.csv("Traits/darkFvFm_sheet.csv")
Rd<-read.csv("Traits/Rd_sheet.csv")

## note: lightFvFm has values halved for "Betula alleg. shade 2"
## environmental variables for some reason, but that's inconsequential

## filter out data where the proper corrections
## couldn't be done for Asat and Rd
bad_corrections<-c("Lala sun2")
Asat$Amax[which(Asat$Leaf.ID %in% bad_corrections)]<-NA
Asat$Corrected.gs[which(Asat$Leaf.ID %in% bad_corrections)]<-NA
# this sample already has NA for Rd

## filter out data marked as having
## unreliable Rd values
bad_Rd<-c("Abiesbalsamea sun/inte 1","Larix decidua? Sun",
          "Picea glauca ? Sun","Picea glauca shade 1",
          "Picea ormonica shade 1","Picea plot m4 sun/int",
          "Picea plotm8shade","Picea rubens shade 1",
          "Piom1","Rhodo degronianum shade 1",
          "Rhodo makinoi shade 2","Abies balsamea sun2")
Rd$Rd..dark.resp.[which(Rd$Leaf.ID %in% bad_Rd)]<-NA

## remove 'PAR->0' measurements, which are less
## reliable than lamp off measurements
## these are always the second one of their label
Rd$Rd..dark.resp.[which(duplicated(Rd$Leaf.ID))]<-NA

## replace erroneous value with better one from raw data
darkFvFm$Fv.Fm[which(darkFvFm$Leaf.ID=="Picea rubens shade 1")]<-0.682972273

## combine datasets
Asat_lightFvFm<-match(Asat$Leaf.ID,lightFvFm$Leaf.ID)
Asat$Fv..Fm.<-lightFvFm$Fv..Fm.[Asat_lightFvFm]
Asat$qP<-lightFvFm$qP[Asat_lightFvFm]
Asat$qN<-lightFvFm$qN[Asat_lightFvFm]
Asat$NPQ<-lightFvFm$NPQ[Asat_lightFvFm]
Asat$ETR<-lightFvFm$ETR[Asat_lightFvFm]

Asat_darkFvFm<-match(Asat$Leaf.ID,darkFvFm$Leaf.ID)
Asat$Fv.Fm<-darkFvFm$Fv.Fm[Asat_darkFvFm]

Asat_Rd<-match(Asat$Leaf.ID,Rd$Leaf.ID)
Asat$Rd..dark.resp.<-Rd$Rd..dark.resp.[Asat_Rd]
Asat$Tleaf_Rdark<-Rd$Average.of.Tleaf[Asat_Rd]

## attach proper names
Asat$SampleID<-data_summary$ID[match(Asat$Leaf.ID,data_summary$GEName)]

Asat<-rename(Asat,
             A=Amax,
             CO2S=Average.of.CO2S,
             Patm=Average.of.Press,
             Qin=Average.of.PARi,
             RHs=Average.of.RH_S,
             Tleaf=Average.of.Tleaf,
             gsw=Corrected.gs,
             Rdark=Rd..dark.resp.)

###############################
## read Julien's Farquhar model fits

load("ProcessedData/LamourAnalyses/2_Fitted_ACi_data.Rdata")
# what about temperature-corrected Rd?

#############################################
## attach to spectral data

fresh_spectra<-readRDS("ProcessedData/fresh_spectra_processed.rds")
dry_spectra<-readRDS("ProcessedData/dry_spectra_processed.rds")

## note: dry samples were scanned with labels "PIST 2 sun"
## and "PIST2 sun". currently we use only the latter
## they are presumably the same sample, just mistakenly
## measured twice
dry_spectra<-dry_spectra[-which(names(dry_spectra)=="PIST 2 sun")]

summary_match_fresh<-match(meta(fresh_spectra)$sample_id,data_summary$SpectraName)
meta(fresh_spectra)$SampleID<-data_summary$ID[summary_match_fresh]
meta(fresh_spectra)$Species<-data_summary$Species[summary_match_fresh]
meta(fresh_spectra)$LDMC<-data_summary$LDMC[summary_match_fresh]
meta(fresh_spectra)$LMA<-data_summary$LMA[summary_match_fresh]
meta(fresh_spectra)$leaf_type<-data_summary$leaf_type[summary_match_fresh]
meta(fresh_spectra)$pubescent<-data_summary$pubescent[summary_match_fresh]
meta(fresh_spectra)$leaf_habit<-data_summary$leaf_habit[summary_match_fresh]
meta(fresh_spectra)$Photosynthetic_pathway<-data_summary$pathway[summary_match_fresh]

spectra_Asat_match_fresh<-match(meta(fresh_spectra)$SampleID,Asat$SampleID)
meta(fresh_spectra)$Asat<-Asat$A[spectra_Asat_match_fresh]
meta(fresh_spectra)$Rd<-Asat$Rdark[spectra_Asat_match_fresh]
meta(fresh_spectra)$ETR<-Asat$ETR[spectra_Asat_match_fresh]
meta(fresh_spectra)$darkFvFm<-Asat$Fv.Fm[spectra_Asat_match_fresh]

bilan_match_fresh<-match(meta(fresh_spectra)$SampleID,Bilan$SampleID)
meta(fresh_spectra)$Vcmax25<-Bilan$Vcmax25[bilan_match_fresh]

## attach to dry spectra...
summary_match_dry<-match(meta(dry_spectra)$sample_id,data_summary$DrySpectraName)
meta(dry_spectra)$SampleID<-data_summary$ID[summary_match_dry]
meta(dry_spectra)$Species<-data_summary$Species[summary_match_dry]
meta(dry_spectra)$LDMC<-data_summary$LDMC[summary_match_dry]
meta(dry_spectra)$LMA<-data_summary$LMA[summary_match_dry]
meta(dry_spectra)$leaf_type<-data_summary$leaf_type[summary_match_dry]
meta(dry_spectra)$pubescent<-data_summary$pubescent[summary_match_dry]
meta(dry_spectra)$leaf_habit<-data_summary$leaf_habit[summary_match_dry]
meta(dry_spectra)$Photosynthetic_pathway<-data_summary$pathway[summary_match_dry]

spectra_Asat_match_dry<-match(meta(dry_spectra)$SampleID,Asat$SampleID)
meta(dry_spectra)$Asat<-Asat$A[spectra_Asat_match_dry]
meta(dry_spectra)$Rd<-Asat$Rdark[spectra_Asat_match_dry]
meta(dry_spectra)$ETR<-Asat$ETR[spectra_Asat_match_dry]
meta(dry_spectra)$darkFvFm<-Asat$Fv.Fm[spectra_Asat_match_dry]

bilan_match_dry<-match(meta(dry_spectra)$SampleID,Bilan$SampleID)
meta(dry_spectra)$Vcmax25<-Bilan$Vcmax25[bilan_match_dry]

## write processed files
saveRDS(fresh_spectra,"ProcessedData/fresh_spectra_and_traits.rds")
saveRDS(dry_spectra,"ProcessedData/dry_spectra_and_traits.rds")