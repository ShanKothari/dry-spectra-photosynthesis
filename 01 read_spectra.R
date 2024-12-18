setwd("C:/Users/Shan Kothari/Dropbox/PostdocProjects/PosadaSpectra/")

library(spectrolab)
library(reshape2)
library(dplyr)

###############################################
## fresh spectra

## function to read the lines of the spectra files
## and only output reflectance values
## resampled to 1 nm resolution

read.spectrum<-function(lines){
  blanks<-which(lines=="")
  lines_sub<-lines[(blanks[2]+1):(blanks[3]-1)]
  lines_split<-strsplit(lines_sub,split=",")
  
  ## if any lines have more or fewer than three elements
  ## spit out an error
  if(sum(unlist(lapply(lines_split,length))!=3)!=0){
    stop("check data file format")
  }
  
  lines_split_df<-do.call(rbind.data.frame, lines_split[-1])
  colnames(lines_split_df)<-c("Wavelength","Raw Data","Reflectance")
  
  lines_split_df$Reflectance<-as.numeric(gsub("%","",lines_split_df$Reflectance))
  lines_split_df$Wavelength<-as.numeric(lines_split_df$Wavelength)
  
  ## trim very noisy ends
  ## later we will trim to 400-900 nm
  lines_split_df<-subset(lines_split_df,Wavelength>340 & Wavelength<1010)
  reflectance_resample<-prospectr::resample(X = lines_split_df$Reflectance,
                                            wav = lines_split_df$Wavelength,
                                            new.wav = 350:1000)
  
  return(reflectance_resample)
}

## get paths of all spectra files
all_spectra_files<-list.files(path="RawFreshSpectra",
                              pattern=".csv",
                              recursive = T,
                              full.names = T)

## remove those we don't care about
all_spectra_files<-all_spectra_files[-which(grepl("_Calculations",all_spectra_files))]
all_spectra_files<-all_spectra_files[-which(grepl("Compilation",all_spectra_files))]
all_spectra_files<-all_spectra_files[-which(grepl("White_Standard",all_spectra_files))]
all_spectra_files<-all_spectra_files[-which(grepl("Black_CardBoard",all_spectra_files))]
all_spectra_files<-all_spectra_files[-which(grepl("Empty_Chamber",all_spectra_files))]
all_spectra_files<-all_spectra_files[-which(grepl("Test",all_spectra_files))]
all_spectra_files<-all_spectra_files[-which(grepl("Second",all_spectra_files))]

## get information about samples from file paths
all_spectra_files_split<-strsplit(all_spectra_files,split="/")
sample_ids<-unlist(lapply(all_spectra_files_split,
                            function(x) x[[2]]))
measurement_names<-unlist(lapply(all_spectra_files_split,
                                 function(x) x[[3]]))

## read data and pull out reflectances
all_lines<-lapply(all_spectra_files,readLines)
all_reflectance<-lapply(all_lines,read.spectrum)

## turn into a spectra object
all_reflectance_df<-do.call(rbind.data.frame,all_reflectance)
all_spectra<-spectra(all_reflectance_df,
                     bands = 350:1000,
                     names = measurement_names)

meta(all_spectra)$sample_id<-sample_ids
meta(all_spectra)$dummy<-NA

## filtering obviously bad data
all_spectra<-all_spectra[-which(all_spectra[,470]>9)]

## trim to 400-900
## this is a bit conservative, just that on
## visual inspection to 900-950 nm area looks
## a bit funky for some data
all_spectra<-all_spectra[,400:900]

## averaging by sample
all_spectra_avg<-aggregate(all_spectra,
                           by=meta(all_spectra)$sample_id,
                           FUN=try_keep_txt(mean))

saveRDS(all_spectra_avg,"ProcessedData/fresh_spectra_processed.rds")
rm(list=ls())

################################################
## dry spectra

# metadata
dry_meta<-read.csv("RawDrySpectra/metadata.csv")

# spectra
dry_spectra_dirs<-list.dirs(path="RawDrySpectra/")
dry_spectra_list<-lapply(dry_spectra_dirs[-1],read_spectra,format="sig")
dry_spectra_raw<-Reduce(spectrolab::combine,dry_spectra_list)
names(dry_spectra_raw)<-gsub(pattern=".sig",replacement="",names(dry_spectra_raw))

dry_meta$file_name<-paste0(dry_meta$basename,"_",
                           formatC(dry_meta$file_number,
                                   width = 4,
                                   format = "d",
                                   flag = "0"))

## cleaning, resampling, and matching sensors

sensor.ends<-NULL
# find indices of last wavelength for first two sensors
for(i in 1:(ncol(dry_spectra_raw)-1)){if(bands(dry_spectra_raw)[i]>bands(dry_spectra_raw)[i+1]) sensor.ends<-c(sensor.ends,i)}

spectra_s1_bands<-bands(dry_spectra_raw)[1:sensor.ends[1]]
spectra_s2_bands<-bands(dry_spectra_raw)[(sensor.ends[1]+1):sensor.ends[2]]
spectra_s3_bands<-bands(dry_spectra_raw)[(sensor.ends[2]+1):ncol(dry_spectra_raw)]

spectra_s1<-dry_spectra_raw[,spectra_s1_bands]
spectra_s2<-dry_spectra_raw[,spectra_s2_bands]
spectra_s3<-dry_spectra_raw[,spectra_s3_bands]

wvl_begin<-350
wvl_end<-2500

spectra_s1_resamp<-as.matrix(resample(spectra_s1,new_bands = wvl_begin:1010))
spectra_s2_resamp<-as.matrix(resample(spectra_s2,new_bands = 1000:1910))
spectra_s3_resamp<-as.matrix(resample(spectra_s3,new_bands = 1900:wvl_end))
spectra_resamp<-do.call(cbind,args = list(spectra_s1_resamp,
                                          spectra_s2_resamp,
                                          spectra_s3_resamp))

spectra_resamp_long<-melt(spectra_resamp,id.vars="sample_id")
colnames(spectra_resamp_long)<-c("sample_id","wavelength","reflectance")

spectra_resamp_long$wvl_id<-paste(spectra_resamp_long$sample_id,spectra_resamp_long$wavelength,sep="_")
dup_ids_ref<-spectra_resamp_long$wvl_id[duplicated(spectra_resamp_long$wvl_id)]
spectra_resamp_long_no_dups<-spectra_resamp_long[-which(spectra_resamp_long$wvl_id %in% dup_ids_ref),]
wvl_range<-wvl_begin:wvl_end

## this function is for linear interpolation over the sensor overlap region
interpolate <- function(x) {
  wvls <- x$wavelength
  ref <- x$reflectance
  new_refs <- approx(wvls, ref, xout = wvl_range)$y
  tmp <- data_frame(wavelength = wvl_range, reflectance = new_refs)
  return(tmp)
}

## apply linear interpolation step over sensor overlap
spectra_resamp_long_cleaned <-spectra_resamp_long_no_dups%>%
  group_by(sample_id) %>%
  do(interpolate(.))

spectra_cleaned<-reshape(data.frame(spectra_resamp_long_cleaned),
                         idvar="sample_id",
                         timevar="wavelength",
                         direction="wide")

spectra_cleaned<-spectra(spectra_cleaned[,-1],
                         bands=wvl_range,
                         names=spectra_cleaned[,1])
dry_spectra<-match_sensors(spectra_cleaned,splice_at = 1005,interpolate_wvl = 10)

## filtering, averaging, and correcting spectra

meta_match<-match(names(dry_spectra),dry_meta$file_name)
meta(dry_spectra)<-dry_meta[meta_match,]

# remove material that's not leaf/filter paper or bad quality
dry_spectra<-dry_spectra[!(meta(dry_spectra)$quality=="bad"),]
dry_spectra<-dry_spectra[!(meta(dry_spectra)$spectra_type %in% c("white","black","test","mistake")),]

# separate out filter paper and average
filter_spectra<-dry_spectra[(meta(dry_spectra)$spectra_type=="only filter paper"),]
filter_spectra_avg<-colMeans(as.data.frame(filter_spectra,metadata = F)[,-1])
dry_spectra<-dry_spectra[!(meta(dry_spectra)$spectra_type=="only filter paper"),]

# split 'good' leaf (and leaf + filter paper) spectra by specimen 
dry_df<-as.data.frame(dry_spectra)
dry_df_split<-split(dry_df, f = dry_df$sample_id)

## which leaves have only target spectra (no filter paper)? (big leaves)
spectra_types<-unlist(lapply(dry_df_split,function(x) length(unique(x$spectra_type))))
big_leaves<-names(dry_df_split)[which(spectra_types<2)]

big_spectra<-dry_df[dry_df$sample_id %in% big_leaves,]
# average!
big_spectra_mean<-aggregate(big_spectra,
                            by=list(big_spectra$sample_id),
                            FUN=try_keep_txt(mean))

dry_spectra_big<-spectra(big_spectra_mean[,which(colnames(big_spectra_mean) %in% wvl_range)],
                         names=big_spectra_mean$sample_id,
                         bands=wvl_range,
                         meta = data.frame(sample_id = big_spectra_mean$sample_id))

## perform corrections for the remainder (small leaves)

small_leaves<-dry_df_split[!(names(dry_df_split) %in% big_leaves)]
# further processing here
small_leaves_avg<-lapply(small_leaves,
                         function(sample) aggregate(sample,
                                                    by=list(sample$spectra_type),
                                                    FUN=try_keep_txt(mean)))

# assumption about spectralon reflectance at 400 nm
spectralon_ref<-0.99

small_leaves_corr<-lapply(small_leaves_avg,
                               function(sample){
                                 wvl_cols<-which(colnames(sample) %in% wvl_range)
                                   
                                 with_paper<-sample[sample$spectra_type=="array filter paper",wvl_cols]
                                 without_paper<-sample[sample$spectra_type=="target",wvl_cols]
                                 
                                 wvl_400<-which(colnames(with_paper)==400)
                                 gap_fraction<-(with_paper[wvl_400]-without_paper[wvl_400])*spectralon_ref/filter_spectra_avg[wvl_400]
                                 corr_spec<-without_paper*spectralon_ref/(1-as.numeric(gap_fraction))
                                 
                                 sample_id<-sample$sample_id[[1]]
                                 return(list(gap_fraction,sample_id,corr_spec))
                               })

small_leaves_corr_list<-lapply(small_leaves_corr,function(x) x[[3]])
small_leaves_corr_df<-do.call(rbind.data.frame,small_leaves_corr_list)

dry_spectra_small<-spectra(small_leaves_corr_df,
                           bands = wvl_range,
                           names = names(small_leaves_corr))

meta(dry_spectra_small)$sample_id<-unlist(lapply(small_leaves_corr,
                                                    function(x) x[[2]]))
meta(dry_spectra_small)$gap_fraction<-unlist(lapply(small_leaves_corr,
                                                    function(x) x[[1]]))

dry_spectra_processed<-spectrolab::combine(dry_spectra_big,dry_spectra_small)

saveRDS(dry_spectra_processed,"ProcessedData/dry_spectra_processed.rds")
