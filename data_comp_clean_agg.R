############################################################################
#~~     Incorporation of Sam Bashevkin's LTMRdata data compilation package
#~~     to convert dataset into cluster analysis format.
#~~     Data are aggregated into 2 data types: 
#~~         1) open water year-round surveys (Bay Study MWT and UCD OTR)
#~~         2) open water Fall-only surveys (Bay Study MWT, FMWT, and UCD OTR)
#~~
#~~     Script initially created by JW Gaeta on Apr 9, 2020
#~~     Script last modified by JW Gaeta on October 27, 2020
############################################################################

# rm(list=ls())
# graphics.off()
# cat("\f")
# threshold = 38


library(tidyr)
# install.packages("devtools")
# library("devtools")
# devtools::install_github("sbashevkin/LTMRdata")
require(LTMRdata)
data("Species")
print(Species, n=dim(Species)[1])
############################################################################
#~~     Data upload and reduction to sites those at or above threshold
############################################################################

#calc_tow_volume=TRUE
if(fall_only == TRUE){
  threshold = 33 # minimum number of years a station is surveyed to be included
} else {
  threshold = 38 # minimum number of years a station is surveyed to be included
}

catch_cutoff_threshold = 1 # observation minimum per month
min_length_threshold_40=TRUE

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~    UC Davis Suisun Marsh Data
data(Suisun)
#str(Suisun)
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
#str(Baystudy)
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



#bs_sta_above_thresh = names(which(colSums(na0s_Baystudy_wide_mat)>=threshold))

bs_sta_above_thresh = c()
for(i in 1:length(bs_sta)){
  sub = subset(Baystudy, Baystudy$Station==bs_sta[i])
  sub_yr = unique(sub$yr)
  if(length(sub_yr)>=threshold){
    bs_sta_above_thresh = c(bs_sta_above_thresh, bs_sta[i])
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~    CDFW Bay Study Data
data(FMWT)
#str(FMWT)
FMWT$yr=as.integer(format(x = FMWT$Date, format="%Y"))
FMWT = subset(FMWT, FMWT$yr>1979 & FMWT$yr<2018)
FMWT$Survey = as.integer(format(FMWT$Date, format="%m"))
FMWT = subset(FMWT, FMWT$Survey %in% c(9,10,11))

# FMWT_na = FMWT[is.na(FMWT$Tow_volume),]
# FMWT_na$Survey = as.integer(format(FMWT_na$Date, format="%m"))
# FMWT_na = subset(FMWT_na, FMWT_na$Survey %in% c(9,10,11))
# table(FMWT_na$Survey)
# table(FMWT_na$Station)
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
head(FMWT_wide)
FMWT_wide_mat=as.matrix(FMWT_wide[,2:39])
rownames(FMWT_wide_mat) = FMWT_wide[,1]
FMWT_wide_mat = t(FMWT_wide_mat)
FMWT_wide_mat[!is.na(FMWT_wide_mat)]=1

#FMWT_wide_mat[,125:130]

na0s_FMWT_wide_mat = FMWT_wide_mat
na0s_FMWT_wide_mat[is.na(na0s_FMWT_wide_mat)]=0

fmwt_sta_above_thresh = names(which(colSums(na0s_FMWT_wide_mat)>=threshold))



############################################################################
#~~     Generate Aggregated Datasets
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
#~~ test of removing all individuals <40mm NOT taking into account length conversion

if(min_length_threshold_40==TRUE){
  dat = dat[which(dat$Length>=40),]
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Add survey number to non-numbered surveys (UCD)
dat$Survey=as.integer(format(dat$Date, "%m"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~   calculate CPUE (NOTE: 1.5m conversion into tow volume)
#UCD otter net  = 4.5mx1.5m mouth *70% opening * tow time/60 minutes *4000m/hr 

if(calc_tow_volume==TRUE){
  dat$cpue = c(dat$Count/(dat$Tow_area*1.5))
} else {
  dat$cpue = c(dat$Count/(dat$Tow_area))
}

dat$cpue[which(dat$Method=="Midwater trawl")] = 
  dat$Count[which(dat$Method=="Midwater trawl")]/dat$Tow_volume[which(dat$Method=="Midwater trawl")]



fmwt_index = which(dat$Station %in% fmwt_sta_above_thresh & dat$Source=="FMWT")
#sort(unique(dat$Station[fmwt_index]))
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




#~~~~~~~~~~~~~~~~~~~~~~~~~~~
#plot

setwd("/Users/jeremewgaeta/Files/DFW/Review_2019_2020/Cluster/integrated_data_analyses/")

# pdf(file = "Taxa_thresholds.pdf", width=6, height=10, paper = "special")
# quartz(width=6, height=10)
par(mfrow=c(3,1), mar=c(4,4,6,4)+0.1)
plot(taxa_count ~ c(0:12), type='n', las=1,
     ylab="", xlab="")
abline(v=seq(0,25,by=1), col="gray80")
abline(h=seq(0,110,by=10), col="gray80")
lines(taxa_count ~ c(0:12), lwd=2.5)
axis(side = 3, at = seq(from=0, to=12, by=1), labels=paste(0:12, "/", "12", sep=""))
axis(side = 4, at = taxa_count[1]*c(1,0.9,0.8,0.7,0.6, 0.5,0.4,0.3,0.2),
     labels=c(1,0.9,0.8,0.7,0.6, 0.5,0.4,0.3,0.2)*100, las=1)
mtext(text = "Minimum Detection Count (average per survey)",
      side = 3, line=2.25, cex=0.75)
mtext(text = "Minimum Detection Count (average per year)",
      side = 1, line=2.25, cex=0.75)
mtext(text = "Percent of Total Species",
      side = 4, line=2.25, cex=0.75)
mtext(text = "Species cpue",
      side = 2, line=2.25, cex=0.75)
mtext(text = "Year-round open water surveys (BS MWT; UCD OTR)",
      side = 3, line=3.75, font=2, cex=1)

plot(taxa_demersal_count ~ c(0:12), type='n', las=1,
     ylab="", xlab="")
abline(v=seq(0,25,by=1), col="gray80")
abline(h=seq(0,130,by=10), col="gray80")
lines(taxa_demersal_count ~ c(0:12), lwd=2.5)
axis(side = 3, at = seq(from=0, to=12, by=1), labels=paste(0:12, "/", "12", sep=""))
axis(side = 4, at = taxa_demersal_count[1]*c(1,0.9,0.8,0.7,0.6, 0.5,0.4,0.3,0.2),
     labels=c(1,0.9,0.8,0.7,0.6, 0.5,0.4,0.3,0.2)*100, las=1)
mtext(text = "Minimum Detection Count (average per survey)",
      side = 3, line=2.25, cex=0.75)
mtext(text = "Minimum Detection Count (average per year)",
      side = 1, line=2.25, cex=0.75)
mtext(text = "Percent of Total Species",
      side = 4, line=2.25, cex=0.75)
mtext(text = "Species cpue",
      side = 2, line=2.25, cex=0.75)
mtext(text = "Year-round demersal surveys (BS OTR; UCD OTR)",
      side = 3, line=3.75, font=2, cex=1)

plot(taxa_fall_count ~ c(0:12), type='n', las=1,
     ylab="", xlab="")
abline(v=seq(0,25,by=1), col="gray80")
abline(h=seq(0,130,by=10), col="gray80")
lines(taxa_fall_count ~ c(0:12), lwd=2.5)
axis(side = 3, at = seq(from=0, to=12, by=1), labels=paste(0:12, "/", "4", sep=""))
axis(side = 4, at = taxa_fall_count[1]*c(1,0.9,0.8,0.7,0.6, 0.5,0.4,0.3,0.2),
     labels=c(1,0.9,0.8,0.7,0.6, 0.5,0.4,0.3,0.2)*100, las=1)
mtext(text = "Minimum Detection Count (average per survey)",
      side = 3, line=2.25, cex=0.75)
mtext(text = "Minimum Detection Count (average per year)",
      side = 1, line=2.25, cex=0.75)
mtext(text = "Percent of Total Species",
      side = 4, line=2.25, cex=0.75)
mtext(text = "Species cpue",
      side = 2, line=2.25, cex=0.75)
mtext(text = "Fall-only open water surveys (BS MWT; FMWT; UCD OTR)",
      side = 3, line=3.75, font=2, cex=1)
# dev.off()




