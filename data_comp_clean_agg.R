############################################################################
#~~     Incorporation of Sam Bashevkin's LTMRdata data compilation package
#~~     to convert dataset into cluster analysis format.
#~~     Data are aggregated into 3 "data groups": 
#~~         1) open water, year-round surveys (Bay Study MWT and UCD OTR)
#~~         2) open water, fall-only surveys (Bay Study MWT, FMWT, and UCD OTR)
#~~         3) demersal, year-round surveys (Bay Study OTR and UCD OTR)
#~~
#~~     Script initially created by JW Gaeta on Apr 9, 2020
#~~     Script last modified by JW Gaeta on October 28, 2020
############################################################################

# rm(list=ls())
# graphics.off()
# cat("\f")

############################################################################
#~~     STEP 1: Upload data and subset data based on minimum years surveyed,
#~~             minimum annual detection, and minimum length thresholds
############################################################################

library(tidyr)
# install.packages("devtools")
# library("devtools")
# devtools::install_github("sbashevkin/LTMRdata")
require(LTMRdata)
data("Species")
print(Species, n=dim(Species)[1])

if(fall_only == TRUE){
  threshold = 33 # minimum number of years a station is surveyed to be included
} else {
  threshold = 38 # minimum number of years a station is surveyed to be included
}

check_exists=function(x) tryCatch(if(class(x) == 'logical') 1 else 1, error=function(e) 0) 

if(check_exists(catch_cutoff_threshold)==0){
  catch_cutoff_threshold = 1 #   observation minimum per month
}
min_length_threshold_40=TRUE # applying minimum length threshold of 40mm

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~    UC Davis Suisun Marsh Data
data(Suisun)
Suisun$yr=format(x = Suisun$Date, format="%Y")
ucd_sta = unique(Suisun$Station)

ucd_sta_above_thresh = c()
for(i in 1:length(ucd_sta)){
  sub = subset(Suisun, Suisun$Station==ucd_sta[i])
  sub_yr = unique(sub$yr)
  if(length(sub_yr)>=threshold){
    ucd_sta_above_thresh = c(ucd_sta_above_thresh, ucd_sta[i])
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~    CDFW Bay Study Data
data(Baystudy)
Baystudy$yr=format(x = Baystudy$Date, format="%Y")
bs_sta = unique(Baystudy$Station)
Baystudy$station_int = as.integer(Baystudy$Station)
Baystudy = Baystudy[order(Baystudy$yr, Baystudy$station_int),]
unique(Baystudy$station_int)

# Find "no_fish_caught" and add zeros to dataset
Baystudy_wide=Baystudy %>%
  pivot_wider(names_from = yr,
              values_from = Tow_volume,
              values_fn = list(Tow_volume = sum),
              id_cols = c(station_int))
Baystudy_wide=data.frame(Baystudy_wide)
Baystudy_wide=Baystudy_wide[order(Baystudy_wide$station_int),]
head(Baystudy_wide)
Baystudy_wide_mat=as.matrix(Baystudy_wide[,2:dim(Baystudy_wide)[2]])
rownames(Baystudy_wide_mat) = Baystudy_wide[,1]
Baystudy_wide_mat = t(Baystudy_wide_mat)
Baystudy_wide_mat[!is.na(Baystudy_wide_mat)]=1

na0s_Baystudy_wide_mat = Baystudy_wide_mat
na0s_Baystudy_wide_mat[is.na(na0s_Baystudy_wide_mat)]=0


bs_sta_above_thresh = c()
for(i in 1:length(bs_sta)){
  sub = subset(Baystudy, Baystudy$Station==bs_sta[i])
  sub_yr = unique(sub$yr)
  if(length(sub_yr)>=threshold){
    bs_sta_above_thresh = c(bs_sta_above_thresh, bs_sta[i])
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~    CDFW Fall Midwater Trawl Study Data
data(FMWT)
FMWT$yr=as.integer(format(x = FMWT$Date, format="%Y"))
FMWT = subset(FMWT, FMWT$yr>1979 & FMWT$yr<2018)
FMWT$Survey = as.integer(format(FMWT$Date, format="%m"))
FMWT = subset(FMWT, FMWT$Survey %in% c(9,10,11))

FMWT$station_int = as.integer(FMWT$Station)
FMWT = FMWT[order(FMWT$yr, FMWT$station_int),]
unique(FMWT$station_int)
FMWT_wide=FMWT %>%
  pivot_wider(names_from = yr,
              values_from = Tow_volume,
              values_fn = list(Tow_volume = sum),
              id_cols = c(station_int))
FMWT_wide=data.frame(FMWT_wide)
FMWT_wide=FMWT_wide[order(FMWT_wide$station_int),]

FMWT_wide_mat=as.matrix(FMWT_wide[,2:39])
rownames(FMWT_wide_mat) = FMWT_wide[,1]
FMWT_wide_mat = t(FMWT_wide_mat)
FMWT_wide_mat[!is.na(FMWT_wide_mat)]=1

# Find "no_fish_caught" and add zeros to dataset
na0s_FMWT_wide_mat = FMWT_wide_mat
na0s_FMWT_wide_mat[is.na(na0s_FMWT_wide_mat)]=0

fmwt_sta_above_thresh = names(which(colSums(na0s_FMWT_wide_mat)>=threshold))


############################################################################
#~~     STEP 2: Aggregate study-specific datasets based on data groups, 
#~~             restrucutre for Principal Tensor Analysis, and eliminate 
#~~             non-fish species and unidentified fishes.
#~~       Data groups are as follows:
#~~         1) open water, year-round surveys (Bay Study MWT and UCD OTR)
#~~         2) open water, fall-only surveys (Bay Study MWT, FMWT, and UCD OTR)
#~~         3) demersal, year-round surveys (Bay Study OTR and UCD OTR)
############################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~    develop study-specific station IDs

dat<- LTMRpilot(convert_lengths = TRUE,
                remove_unconverted_lengths = FALSE)
str(dat)
unique(dat$Source)
dat$lme = rep("ucd", times = dim(dat)[1])
dat$lme[which(dat$Source=="Bay Study")] = "bs"
dat$lme[which(dat$Source=="FMWT")] = "fmwt"
unique(dat$lme)
dat$sta_lme = paste(dat$lme, dat$Station, sep="_")
str(Species)
dat = merge(dat, Species[,c(2,7)])


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#~~ Note: length conversions were only available for a subset of species.
#~~     Refer to the LTMR report for additional details
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(min_length_threshold_40==TRUE){
  dat = dat[which(dat$Length>=40),]
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Add survey number to non-numbered surveys (i.e., UCD)
dat$Survey=as.integer(format(dat$Date, "%m"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~   calculate CPUE (NOTE: for UCD Otter Trawl conversion from area into volume: 
#~~                   multiply by 1.5m)


if(calc_tow_volume==TRUE){
  dat$cpue = c(dat$Count/(dat$Tow_area*1.5))
} else {
  dat$cpue = c(dat$Count/(dat$Tow_area))
}

dat$cpue[which(dat$Method=="Midwater trawl")] = 
  dat$Count[which(dat$Method=="Midwater trawl")]/dat$Tow_volume[which(dat$Method=="Midwater trawl")]


fmwt_index = which(dat$Station %in% fmwt_sta_above_thresh & dat$Source=="FMWT")
bs_index = which(dat$Station %in% bs_sta_above_thresh & dat$Source=="Bay Study")
ucd_index = which(dat$Station %in% ucd_sta_above_thresh & dat$Source=="Suisun")


red = dat[c(fmwt_index,bs_index,ucd_index), ]
red$yr = format(red$Date, format="%Y")
red$Method[which(red$Method=="Otter trawl")] = "Otter Trawl"
red2=red[!is.na(red$cpue),]

lat_long_tab = aggregate(cbind(Longitude, Latitude) ~ sta_lme, data = red2, FUN = median)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~     Pivot Wider
cpue=red2 %>%
  pivot_wider(names_from = CommonName,
              values_from = cpue,
            values_fn = list(cpue = sum),
            id_cols = c(Method, Source, lme, sta_lme, Station, Survey, yr),
            values_fill = list(cpue=0))
cpue=data.frame(cpue)
unique(cpue$Method)

count=red2 %>%
  pivot_wider(names_from = CommonName,
              values_from = Count,
              values_fn = list(Count = sum),
              id_cols = c(Method, Source, lme, sta_lme, Station, Survey, yr),
              values_fill = list(Count=0))
count=data.frame(count)
unique(count$Method)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~     Break into all season long and fall midwater taxa

all_sea = cpue[which(cpue$Method=="Midwater trawl" & cpue$lme=="bs" |
                   cpue$Method=="Otter Trawl" & cpue$lme=="ucd" ),]
demersal = cpue[which(cpue$Method=="Otter Trawl" & cpue$lme=="bs" |
                        cpue$Method=="Otter Trawl" & cpue$lme=="ucd" ),]
fall = cpue[which(cpue$Method=="Midwater trawl" &
                     cpue$lme=="bs" & cpue$Survey %in% c(9,10,11) |
                     cpue$lme=="fmwt" & cpue$Survey %in% c(9,10,11) |
                     cpue$lme=="ucd" & cpue$Survey %in% c(9,10,11)),]

count_all_sea = count[which(count$Method=="Midwater trawl" & count$lme=="bs" |
                       count$Method=="Otter Trawl" & count$lme=="ucd" ),]
count_demersal = count[which(count$Method=="Otter Trawl" & count$lme=="bs" |
                        count$Method=="Otter Trawl" & count$lme=="ucd" ),]
count_fall = count[which(count$Method=="Midwater trawl" &
                    count$lme=="bs" & count$Survey %in% c(9,10,11) |
                    count$lme=="fmwt" & count$Survey %in% c(9,10,11) |
                    count$lme=="ucd" & count$Survey %in% c(9,10,11)),]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~  eliminate non-fish or unidentified fishes from the dataset
taxa_names = colnames(all_sea)
cut_taxa_index = c()
cut_taxa_index = c(cut_taxa_index, grep("Clam", taxa_names))
cut_taxa_index = c(cut_taxa_index, grep("Shrimp", taxa_names))
cut_taxa_index = c(cut_taxa_index, grep("Maeotias", taxa_names))
cut_taxa_index = c(cut_taxa_index, grep("Prawn", taxa_names))
cut_taxa_index = c(cut_taxa_index, grep("Sea.Gooseberry", taxa_names))
cut_taxa_index = c(cut_taxa_index, grep("Jellyfish", taxa_names))
cut_taxa_index = c(cut_taxa_index, grep("UnID", taxa_names))
cut_taxa_index = c(cut_taxa_index, grep("Unid", taxa_names))
cut_taxa_index = c(cut_taxa_index, grep("Crab", taxa_names))
cut_taxa_index = c(cut_taxa_index, grep("No.Catch", taxa_names))
cut_taxa_index = c(cut_taxa_index, grep("No.Trawl", taxa_names))
cut_taxa_index_sorted = unique(sort(cut_taxa_index))

all_sea = all_sea[,-c(cut_taxa_index_sorted)]
count_all_sea = count_all_sea[,-c(cut_taxa_index_sorted)]

taxa_count = c()
for(i in 0:12){
  taxa_count=c(taxa_count, length(which(colSums(count_all_sea[,-c(1:7)])>i*38)))
}


#~~~~~
taxa_fall_names = colnames(fall)
cut_taxa_fall_index = c()
cut_taxa_fall_index = c(cut_taxa_fall_index, grep("Clam", taxa_fall_names))
cut_taxa_fall_index = c(cut_taxa_fall_index, grep("Shrimp", taxa_fall_names))
cut_taxa_fall_index = c(cut_taxa_fall_index, grep("Maeotias", taxa_fall_names))
cut_taxa_fall_index = c(cut_taxa_fall_index, grep("Prawn", taxa_fall_names))
cut_taxa_fall_index = c(cut_taxa_fall_index, grep("Sea.Gooseberry", taxa_fall_names))
cut_taxa_fall_index = c(cut_taxa_fall_index, grep("Jellyfish", taxa_fall_names))
cut_taxa_fall_index = c(cut_taxa_fall_index, grep("UnID", taxa_fall_names))
cut_taxa_fall_index = c(cut_taxa_fall_index, grep("Unid", taxa_fall_names))
cut_taxa_fall_index = c(cut_taxa_fall_index, grep("Crab", taxa_fall_names))
cut_taxa_fall_index = c(cut_taxa_fall_index, grep("No.", taxa_fall_names))
cut_taxa_fall_index_sorted = unique(sort(cut_taxa_fall_index))

fall = fall[,-c(cut_taxa_fall_index_sorted)]
count_fall = count_fall[,-c(cut_taxa_fall_index_sorted)]

taxa_fall_count = c()
for(i in 0:12){
  taxa_fall_count=c(taxa_fall_count, length(which(colSums(count_fall[,-c(1:7)])>i*38)))
}

#~~~~~
taxa_demersal_names = colnames(demersal)
cut_taxa_demersal_index = c()
cut_taxa_demersal_index = c(cut_taxa_demersal_index, grep("Clam", taxa_demersal_names))
cut_taxa_demersal_index = c(cut_taxa_demersal_index, grep("Shrimp", taxa_demersal_names))
cut_taxa_demersal_index = c(cut_taxa_demersal_index, grep("Maeotias", taxa_demersal_names))
cut_taxa_demersal_index = c(cut_taxa_demersal_index, grep("Prawn", taxa_demersal_names))
cut_taxa_demersal_index = c(cut_taxa_demersal_index, grep("Sea.Gooseberry", taxa_demersal_names))
cut_taxa_demersal_index = c(cut_taxa_demersal_index, grep("Jellyfish", taxa_demersal_names))
cut_taxa_demersal_index = c(cut_taxa_demersal_index, grep("UnID", taxa_demersal_names))
cut_taxa_demersal_index = c(cut_taxa_demersal_index, grep("Unid", taxa_demersal_names))
cut_taxa_demersal_index = c(cut_taxa_demersal_index, grep("Crab", taxa_demersal_names))
cut_taxa_demersal_index = c(cut_taxa_demersal_index, grep("No.", taxa_demersal_names))
cut_taxa_demersal_index_sorted = unique(sort(cut_taxa_demersal_index))

demersal = demersal[,-c(cut_taxa_demersal_index_sorted)]
count_demersal = count_demersal[,-c(cut_taxa_demersal_index_sorted)]

taxa_demersal_count = c()
for(i in 0:12){
  taxa_demersal_count=c(taxa_demersal_count, length(which(colSums(count_demersal[,-c(1:7)])>i*38)))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~  eliminate taxa below threshold cpue

#sort(colSums(all_sea[,-c(1:7)]), decreasing = TRUE)
all_cut_index=which(colSums(count_all_sea[,-c(1:7)])<c(catch_cutoff_threshold*38*12))
all_sea = all_sea[,-c(all_cut_index+7)]
sort(colSums(all_sea[,-c(1:7)]), decreasing=TRUE)

#sort(colSums(fall[,-c(1:7)]), decreasing = TRUE)
fall_cut_index=which(colSums(count_fall[,-c(1:7)])<c(catch_cutoff_threshold*33*4))
fall = fall[,-c(fall_cut_index+7)]
sort(colSums(fall[,-c(1:7)]), decreasing=TRUE)

#sort(colSums(fall[,-c(1:7)]), decreasing = TRUE)
demersal_cut_index=which(colSums(count_demersal[,-c(1:7)])<c(catch_cutoff_threshold*38*12))
demersal = demersal[,-c(demersal_cut_index+7)]
sort(colSums(demersal[,-c(1:7)]), decreasing=TRUE)





