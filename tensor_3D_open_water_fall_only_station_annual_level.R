
#############################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~    Cluster analysis of open water, fall-only data
#~    Script created on March 19, 2020 by JW Gaeta
#~    Script last modified on November 10, 2020 by JW Gaeta
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#############################################################################

# rm(list=ls())
# graphics.off()
# cat("\f")

library(vegan)
library(dendextend)
library(ade4)
library(PTAk)
library(cluster)
library(sfsmisc)
library(devtools)

############################################################################
#~~     STEP 1: set data aggregation and cleaning parameters, then run data 
#~~             compilation & processing script
############################################################################

calc_tow_volume=TRUE
fall_only = TRUE

source_url("https://raw.githubusercontent.com/JGaetaFish/LTMR_2020_Fish_Community/main/data_comp_clean_agg.R")

fall = fall[which(fall$yr>1984 & fall$yr<2018),]
fall = fall[order(fall$sta_lme, fall$yr),]

############################################################################
#~~     STEP 2: Convert surveys into annual timestep:
#~~                 Fall   = Sept-Nov
#~~         NOTE: too many survey-stations combinations were missing to run 
#~~               the analysis at the survey-station level
#~~               (see the report for details)
############################################################################

fall2=fall %>%
  pivot_longer(cols = !Method:yr,
               names_to = "CommonName",
               values_to="cpue")

relative_sum = function(x){sum(x)/length(x)}

fall3=fall2 %>%
  pivot_wider(names_from = CommonName,
              values_from = cpue,
              values_fn = list(cpue = relative_sum),
              id_cols = c(yr, lme, sta_lme),
              values_fill = list(cpue=0))
fall3=as.data.frame(fall3)

############################################################################
#~~     STEP 3: Add zeros for seasons-stations flagged as "No fish caught"
#~~             i.e., instances of zero detection NOT skipped surveys
############################################################################

######
#       NOTE: the following resulted in the addition of:
#             164 FMWT station-years without fish detected
#               5 BS station-years without fish detected
#               0 UCD station-years without fish detected
#       The addition of these zero-catch station-years increased the fall3
#       dataset by 3.65 %


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~ Fall Midwater Trawl "no_fish_caught"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

data(FMWT)
FMWT$Survey = as.integer(format(FMWT$Date, format="%m"))
FMWT$yr = as.integer(format(FMWT$Date, format="%Y"))
FMWT = FMWT[which(FMWT$yr>1984 & FMWT$yr<2018),]

fallMWT = subset(FMWT, FMWT$Method=="Midwater trawl" & FMWT$Survey %in% c(9,10,11))
fmwt_no_fish = fallMWT[which(fallMWT$Length_NA_flag=="No fish caught" & fallMWT$Station %in% fmwt_sta_above_thresh),]

#~~~~~~~~~~~
# Aggregate annually; compare total surveys per year per station with
#       total surveys per year per station without fish captured
#~~~~~~~~~~~

#~~ Aggregate to survey level
fmwt_survey_station_no_fish = aggregate(SampleID ~ Source + yr + Survey + Station, 
                                      data = fmwt_no_fish, FUN = unique)
fmwt_survey_station = aggregate(SampleID ~ Source + yr + Survey + Station, 
                              data = fallMWT, FUN = unique)

fmwt_survey_station_no_fish$event_count = lengths(fmwt_survey_station_no_fish$SampleID)
fmwt_survey_station$event_count = lengths(fmwt_survey_station$SampleID)

#~~ Aggregate to annual resolution
fmwt_survey_year_station_no_fish = aggregate(event_count ~ Source + yr + Station, 
                                           data = fmwt_survey_station_no_fish, FUN = sum)
fmwt_survey_year_station = aggregate(event_count ~ Source + yr + Station, 
                                        data = fmwt_survey_station, FUN = sum)

colnames(fmwt_survey_year_station)[length(names(fmwt_survey_year_station))] = 
  "total_surveys_per_season"
colnames(fmwt_survey_year_station_no_fish)[length(names(fmwt_survey_year_station_no_fish))] = 
  "surveys_without_fish"

fmwt_no_fish_merged = merge(fmwt_survey_year_station, fmwt_survey_year_station_no_fish)
fmwt_zero_catch_index = which(fmwt_no_fish_merged$surveys_without_fish ==
                              fmwt_no_fish_merged$total_surveys_per_season)
fmwt_add_zeros = fmwt_no_fish_merged[fmwt_zero_catch_index,]

# Add events without fish detected to the original "sea" dataset adding 0s for all fishes
dim(fmwt_add_zeros)[1]/dim(fmwt_survey_year_station)[1]*100
#~~  note 164 station-years were 0 detections; representing 4.34% of total FMWT station-years
fmwt_no_fish_caught = data.frame(yr = fmwt_add_zeros$yr,
                               lme = rep("fmwt", dim(fmwt_add_zeros)[1]),
                               sta_lme = paste(rep("fmwt", dim(fmwt_add_zeros)[1]),
                                               fmwt_add_zeros$Station, sep="_"))

for(i in 4:length(names(fall3))){
  fmwt_no_fish_caught = cbind(fmwt_no_fish_caught,rep(0, times=dim(fmwt_no_fish_caught)[1]))
  colnames(fmwt_no_fish_caught)[i]=names(fall3)[i]
}
fall3 = rbind(fall3,fmwt_no_fish_caught)
fall3 = fall3[order(fall3$sta_lme, fall3$yr),]

fall3[which(fall3$yr==2011 & fall3$sta_lme=="fmwt_902"),]

# Check for duplicate study-station-years
paste(fall3$yr, fall3$sta_lme, sep="_")[duplicated(paste(fall3$yr, fall3$sta_lme, sep="_"))]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~ Bay Study "no_fish_caught"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load Bay Study original dataset and add season
data(Baystudy)
Baystudy$Survey = as.integer(format(Baystudy$Date, format="%m"))
Baystudy$yr = as.integer(format(Baystudy$Date, format="%Y"))
Baystudy_red = Baystudy[which(Baystudy$Survey %in% c(9,10,11)),]

# Identify survey events without fish detected
BayStudyMWT = subset(Baystudy_red, Baystudy_red$Method=="Midwater trawl")
bs_no_fish = BayStudyMWT[which(BayStudyMWT$Length_NA_flag=="No fish caught" &
                                 BayStudyMWT$Station %in% bs_sta_above_thresh),]

#~~~~~~~~~~~
# Aggregate annually; compare total surveys per year per station with
#       total surveys per year per station without fish captured
#~~~~~~~~~~~
bs_survey_station_no_fish = aggregate(SampleID ~ Source + yr + Survey + Station, 
                                           data = bs_no_fish, FUN = unique)
bs_survey_station = aggregate(SampleID ~ Source + yr + Survey + Station, 
                                   data = BayStudyMWT, FUN = unique)

bs_survey_station_no_fish$event_count = lengths(bs_survey_station_no_fish$SampleID)
bs_survey_station$event_count = lengths(bs_survey_station$SampleID)

#~~ Aggregate to annual resolution
bs_survey_year_station_no_fish = aggregate(event_count ~ Source + yr + Station, 
                                      data = bs_survey_station_no_fish, FUN = sum)
bs_survey_year_station_list = aggregate(event_count ~ Source + yr + Station, 
                              data = bs_survey_station, FUN = sum)

colnames(bs_survey_year_station_list)[length(names(bs_survey_year_station_list))] = 
  "total_surveys_per_season"
colnames(bs_survey_year_station_no_fish)[length(names(bs_survey_year_station_no_fish))] = 
  "surveys_without_fish"

bs_no_fish_merged = merge(bs_survey_year_station_list, bs_survey_year_station_no_fish)
bs_zero_catch_index = which(bs_no_fish_merged$surveys_without_fish ==
                              bs_no_fish_merged$total_surveys_per_season)
bs_add_zeros = bs_no_fish_merged[bs_zero_catch_index,]

# Add events without fish detected to the original "sea" dataset adding 0s for all fishes
dim(bs_add_zeros)[1]/dim(bs_survey_year_station_list)[1]*100
#~~  note 5 station-years were 0 detections; representing 0.29% of total Bay Study station-years
bs_no_fish_caught = data.frame(yr = bs_add_zeros$yr,
                               lme = rep("bs", dim(bs_add_zeros)[1]),
                               sta_lme = paste(rep("bs", dim(bs_add_zeros)[1]),
                                               bs_add_zeros$Station, sep="_"))

for(i in 4:length(names(fall3))){
  bs_no_fish_caught = cbind(bs_no_fish_caught,rep(0, times=dim(bs_no_fish_caught)[1]))
  colnames(bs_no_fish_caught)[i]=names(fall3)[i]
}
fall3 = rbind(fall3,bs_no_fish_caught)
fall3 = fall3[order(fall3$sta_lme, fall3$yr),]

fall3[which(fall3$yr==2000 & fall3$sta_lme=="bs_213"),]

# Check for duplicate study-station-years
paste(fall3$yr, fall3$sta_lme, sep="_")[duplicated(paste(fall3$yr, fall3$sta_lme, sep="_"))]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~ Suisun Marsh "no_fish_caught"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load Suisun Marsh Study original dataset and add season
data(Suisun)
Suisun$Survey = as.integer(format(Suisun$Date, format="%m"))
Suisun$yr = as.integer(format(Suisun$Date, format="%Y"))
Suisun_red = Suisun[which(Suisun$Survey %in% c(9,10,11)),]

# Identify survey events without fish detected
SuisunOTR = subset(Suisun_red, Suisun_red$Method=="Otter trawl")
ucd_no_fish = SuisunOTR[which(SuisunOTR$Length_NA_flag=="No fish caught" &
                                SuisunOTR$Station %in% ucd_sta_above_thresh),]

# Aggregate annually; compare total surveys per year per station with
#       total surveys per season without fish captured
ucd_survey_station_no_fish = aggregate(SampleID ~ Source + yr + Survey + Station, 
                                       data = ucd_no_fish, FUN = unique)
ucd_survey_station = aggregate(SampleID ~ Source + yr + Survey + Station, 
                               data = SuisunOTR, FUN = unique)

ucd_survey_station_no_fish$event_count = lengths(ucd_survey_station_no_fish$SampleID)
ucd_survey_station$event_count = lengths(ucd_survey_station$SampleID)

#~~ Aggregate to annual resolution
ucd_survey_year_station_no_fish = aggregate(event_count ~ Source + yr + Station, 
                                            data = ucd_survey_station_no_fish, FUN = sum)
ucd_survey_year_station_list = aggregate(event_count ~ Source + yr + Station, 
                                         data = ucd_survey_station, FUN = sum)

colnames(ucd_survey_year_station_list)[length(names(ucd_survey_year_station_list))] = 
  "total_surveys_per_season"
colnames(ucd_survey_year_station_no_fish)[length(names(ucd_survey_year_station_no_fish))] = 
  "surveys_without_fish"

ucd_no_fish_merged = merge(ucd_survey_year_station_list, ucd_survey_year_station_no_fish)
ucd_zero_catch_index = which(ucd_no_fish_merged$surveys_without_fish ==
                               ucd_no_fish_merged$total_surveys_per_season)
ucd_add_zeros = ucd_no_fish_merged[ucd_zero_catch_index,]

# Add events without fish detected to the original "fall3" dataset adding 0s for all fishes
#~~  note 0 station-years were 0 detections; representing 0% of total Suisun Marsh Study station-years


############################################################################
#~~     STEP 5: Scale and restructure data into an array framework for 
#~~             Principal Tensor Analysis
############################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Standardize by z-scoring per taxa per study-gear

fall3 = fall3[order(fall3$sta_lme, fall3$yr),]

sub2 = fall3[,-c(1:3)]
sub2 = (sub2)^(1/3)
sub_scale=apply(sub2, MARGIN = 2,
                FUN = function(x){
                  (x-mean(x))/sd(x)
                }
)
sub_scale[is.na(sub_scale)]=0
fall3_scaled=cbind(fall3[,c(1:3)], sub_scale)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Convert data into array format
#~  [season_timeseries_number, taxa, station]

fall3_stations = sort(unique(fall3_scaled$sta_lme))
year_seq=min(fall3_scaled$yr):max(fall3_scaled$yr)
taxa=colnames(fall3_scaled[,-c(1:3)])

fall3_scaled_arr=array(data = NA, dim = c(length(year_seq),
                                        length(taxa),
                                        length(fall3_stations)),
                     dimnames = list("ts_fall3"=year_seq,
                                     "taxa"=taxa,
                                     "station"=fall3_stations)
)

dim(fall3_scaled_arr)
for(i in 1:length(fall3_stations)){
  sub = fall3_scaled[which(fall3_scaled$sta_lme==fall3_stations[i]),]
  #print(sub)
  sub2 = sub[,-c(1:3)]
  row_order=intersect(year_seq , sub$yr)
  for(x in 1:dim(sub2)[2]){
    fall3_scaled_arr[row_order,x,i] = sub2[,x]
  }
  
}

na_list=list()
for(i in 1:dim(fall3_scaled_arr)[3]){
  names(fall3_scaled_arr[is.na(fall3_scaled_arr[,2,i]),1,i])
  na_list[i]=list(names(fall3_scaled_arr[is.na(fall3_scaled_arr[,2,i]),1,i]))
}

fall3_scaled_arr[,,146]
fall3_stations[127]
fall3[which(fall3$sta_lme=="fmwt_908"),]

# [fall3son_timeseries_number, taxa, station]
# fall3_scaled_arr[,,1]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ identify survey count timeseries without all stations sampled

vec = c()
na_count = c()
for(i in 1:dim(fall3_scaled_arr)[3]){
  sub = fall3_scaled_arr[,,i]
  sub_na_element = c(1:dim(fall3_scaled_arr)[1])[is.na(sub[,1])]
  vec=c(vec, sub_na_element)
  na_count = c(na_count, length(sub_na_element))
}
na_count
vec=na.omit(vec)
unique(vec)


fall3_scaled_arr[1,1,]
fall3_scaled_arr[,,1]

missing_survey_years = sort(unique(vec))
rownames(fall3_scaled_arr[,,1])[missing_survey_years]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~  Write reduced fall3sonal array (cutting 47 fall3sons or 152) (-18.2%):
dim(fall3_scaled_arr)
red_fall3_arr = fall3_scaled_arr[-missing_survey_years,,]

length(missing_survey_years)/dim(fall3_scaled_arr)[1]

red_fall3_arr[1,1,]; length(red_fall3_arr[1,1,])
red_fall3_arr[1,,1]; length(red_fall3_arr[1,,1])
red_fall3_arr[,1,1]; length(red_fall3_arr[,1,1])
red_fall3_arr[,,40]
red_fall3_arr[8,,40:45]

####################################################################
#~  Principal tensor analysis 'a la' Frelat etal 2017 PlosOne
##################################################################
rm(list=c("pta"))
pta<-PTA3(red_fall3_arr, nbPT = 3, nbPT2 = 3, minpct = 0.1)
summary.PTAk(pta,testvar = 0)
par(mfrow=c(1,1), mar=c(4.5,4.5,1.5,1.5)+0.1)
plot(pta, scree=TRUE)
keep = c(1,3,21)

pta[[1]]$n
pta[[2]]$n
pta[[3]]$n

#Create the scree plot
out <- !substr(pta[[3]]$vsnam, 1, 1) == "*"
gct<-(pta[[3]]$pct*pta[[3]]$ssX/pta[[3]]$ssX[1])[out]
par(mfrow=c(1,1), mar=c(4.5,4.5,1.5,1.5)+0.1)
barplot(sort(gct, decreasing = TRUE), xlab="PT",
        ylab="Percentage of variance", las=1); box(which="plot")

#Top 3 PT = 1, 11, 3
par(mfrow=c(1,1))

# fall3_summary = aggregate(yr ~ fall3son_number + fall3son_count + fall3son_year, FUN=unique, data=fall3)
# fall3_summary$yr=substr(as.character(fall3_summary$yr), start=3, stop=4)
# fall3_summary = subset(fall3_summary, fall3_summary$fall3son_count %in% as.numeric(pta[[1]]$n))

col_fall3_year = c("darkblue", "darkgreen", "darkred", "tan4")

col_fall3_year_start = c("#00D3EB", "#00FF33", "#CDA8A8", "#FFBD04")
col_fall3_year_end = c("#001D6C", "#095500", "#640000", "#372900")

#~~~~~~~~
# Year/fall3son plots

unique_fall3_year_plot = 1:4
# plot(pta[[1]]$v[1,] , pta[[1]]$v[11,], type='n', las=1,
#      xlab  = "Primary Tensor 1 (6.64% Global)",
#      ylab  = "Primary Tensor 3 (2.41% Global)")
# abline(h=0, col=gray(0.2,0.2))
# abline(v=0, col=gray(0.2,0.2))
# for(s in 1:length(unique_fall3_year_plot)){
#   sub_dat = pta[[1]]$v[,which(fall3_summary[,1]==unique_fall3_year_plot[s])]
#   color_range <- colorRampPalette( c(col_fall3_year_start[s], col_fall3_year_end[s]) ) 
#   color_spectra <- color_range(dim(sub_dat)[1])
#   
#   points(sub_dat[1,], sub_dat[11,], col=color_spectra, pch=20)
# }
# 
# 
# plot(pta[[1]]$v[1,] , pta[[1]]$v[11,], type='n', las=1,
#      xlab  = "Primary Tensor 1 (6.64% Global)",
#      ylab  = "Primary Tensor 3 (2.41% Global)")
# abline(h=0, col=gray(0.2,0.2))
# abline(v=0, col=gray(0.2,0.2))
# text(pta[[1]]$v[1,] , pta[[1]]$v[11,],
#      col=col_fall3_year[fall3_summary[,1]],labels =fall3_summary[,4],
#      adj=c(0.5,0.5), cex=0.5, font=2)
# 
# plot(pta[[1]]$v[1,] , pta[[1]]$v[11,], type='n', las=1)
# unique_fall3_year_plot=sort(unique(fall3_summary[,1]))
# for(s in 1:length(unique_fall3_year_plot)){
#   sub_dat = pta[[1]]$v[,which(fall3_summary[,1]==unique_fall3_year_plot[s])]
#   for(t in 1:c(dim(sub_dat)[2]-1)){
#     arrows(x0 = sub_dat[1,t], y0 = sub_dat[11,t],
#            x1 = sub_dat[1,t+1], y1 = sub_dat[11,t+1],
#            length = 0.1, col = col_fall3_year[s])
#   }
# }

plot(pta[[1]]$v[keep[1],] , pta[[1]]$v[keep[2],], type='n', las=1,
     xlab  = "Primary Tensor 1 (6.64% Global)",
     ylab  = "Primary Tensor 3 (2.41% Global)")
abline(h=0, col=gray(0.2,0.2))
abline(v=0, col=gray(0.2,0.2))
text(pta[[1]]$v[keep[1],] , pta[[1]]$v[keep[2],],
     labels =pta[[1]]$n,
     adj=c(0.5,0.5), cex=0.5,font=2)

plot(pta[[2]]$v[keep[1],] , pta[[2]]$v[keep[2],], type='n', las=1,
     xlab  = "Primary Tensor 1 (6.64% Global)",
     ylab  = "Primary Tensor 3 (2.41% Global)")
abline(h=0, col=gray(0.2,0.2))
abline(v=0, col=gray(0.2,0.2))
text(pta[[2]]$v[keep[1],] , pta[[2]]$v[keep[2],],
     labels =pta[[2]]$n,
     adj=c(0.5,0.5), cex=0.5,font=2)


plot(pta[[3]]$v[keep[1],] , pta[[3]]$v[keep[2],], type='n', las=1,
     xlab  = "Primary Tensor 1 (6.64% Global)",
     ylab  = "Primary Tensor 3 (2.41% Global)")
abline(h=0, col=gray(0.2,0.2))
abline(v=0, col=gray(0.2,0.2))
text(pta[[3]]$v[keep[1],] , pta[[3]]$v[keep[2],],
     labels =pta[[3]]$n,
     adj=c(0.5,0.5), cex=0.5,font=2)


# 
# unique_fall3_year_plot=sort(unique(fall3_summary[,1]))
# par(mfrow=c(2,2))
# for(s in 1:length(unique_fall3_year_plot)){
#   index=which(fall3_summary[,1]==unique_fall3_year_plot[s])
#   sub_dat = pta[[1]]$v[,index]
#   plot(x = as.numeric(as.character(fall3_summary[index,3])), y=colSums(sub_dat), type='l', las=1)
#   abline(h=0, col=gray(0.2,0.2))
# }


######################################################################
#       EVALUATE LUSTER FIT AND COMPARE METHODS
######################################################################

#Create the matrix with the projection of species on the 4 PT
coo<-t(pta[[2]]$v[c(keep),])

rownames(coo) = pta[[2]]$n
PT_associations = c("PT1: Spatial", "PT2: Spatial", 
                    "PT3: Temporal (Annual)")
labkeep <- paste0(PT_associations, "\n(",pta[[3]]$vsnam[keep], " - ", round((100 * (pta[[3]]$d[keep])^2)/pta[[3]]$ssX[1],1), "%)")





###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####  Cophenetic Correlation: compares clustering methods   ####
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Numerical Ecology in R p.63: cophenetic dist is dist where 2 objects become members of same group
# correlation is between original dist matrix and cophenetic distance matrix
# higher value is good - implies the cophenetic dist matrix reproduces original dist matrix
# above .7 is good
#1. Compute the distance between species
dist1=dist(coo, method = "euclidean")

#2. Build a trees with various linkage


##  Ward_2 clustering
hclust.ward2<- hclust(dist1, method = "ward.D2")
ward2.coph <- cophenetic(hclust.ward2)
cophcorr.ward2 <- cor(dist1, ward2.coph)

##  Single linkage
hclust.single <- hclust(dist1, method="single")
single.coph <- cophenetic(hclust.single)
cophcorr.single <- cor(dist1, single.coph)

##  Complete linkage
hclust.complete <- hclust(dist1, method="complete")
complete.coph <- cophenetic(hclust.complete)
cophcorr.complete <- cor(dist1, complete.coph)

##  Average clustering
hclust.avg <- hclust(dist1, method="average")
avg.coph <- cophenetic(hclust.avg)
cophcorr.avg <- cor(dist1, avg.coph)

#####

data.frame(cophcorr.ward2, cophcorr.single,
           cophcorr.complete, cophcorr.avg)


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####   Gower (1983) distances: to compare clustering methods   ####
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# p65 Numerical Ecology in R: computed as the sum of squared differences between original and cophenetic distances
# p413 Numerical Ecology: measure only has relative value for clustering results from same distance matrix

## The clustering method that produces the smallest distance is best

gow.dist.single <- sum((dist1 - single.coph)^2)
gow.dist.complete <- sum((dist1 - complete.coph)^2)
gow.dist.average <- sum((dist1 - avg.coph)^2)
gow.dist.ward2 <- sum((dist1 - ward2.coph)^2)

###
#####  Combine coph corr and Gower dist in a dataframe   #####
###

fit.metrics = data.frame(clustmethod = c("single", "complete", "avg", "ward2"),
                         cophcorr = c(cophcorr.single, cophcorr.complete,
                                      cophcorr.avg, cophcorr.ward2),
                         gowdist = c(gow.dist.single, gow.dist.complete,
                                     gow.dist.average, gow.dist.ward2))
fit.metrics

###
####  Shepard-like diagrams: to compare cluster methods    ####
###
# These plot pairwise (original) distance against hclust cophenetic distance


par(mfrow=c(2, 2),mar=c(4.5,4.5,3,1.5)+0.1, oma=c(0,0,0,0))
plot(dist1, single.coph, xlab="Chord distance", 
     ylab="Cophenetic distance", asp=1, xlim=c(0, sqrt(1)), ylim=c(0, sqrt(0.5)), 
     main=c("Single linkage",paste("Cophenetic correlation ",
                                   round(cophcorr.single, 3))),
     pch=19, col=gray(0.25,0.25))
abline(0,1)
lines(lowess(dist1, single.coph), col="red")

plot(dist1, complete.coph, xlab="Chord distance", 
     ylab="Cophenetic distance", asp=1, xlim=c(0, sqrt(2)), ylim=c(0, sqrt(2)),
     main=c("Complete linkage", paste("Cophenetic correlation ",
                                      round(cophcorr.complete, 3))),
     pch=19, col=gray(0.25,0.25))
abline(0,1)
lines(lowess(dist1, complete.coph), col="red")

plot(dist1, avg.coph, xlab="Chord distance", 
     ylab="Cophenetic distance", asp=1, xlim=c(0,sqrt(0.5)), ylim=c(0,sqrt(0.5)), 
     main=c("UPGMA (Average)", paste("Cophenetic correlation ",
                                     round(cophcorr.avg, 3))),
     pch=19, col=gray(0.25,0.25))
abline(0,1)
lines(lowess(dist1, avg.coph), col="red")

plot(dist1, ward2.coph, xlab="Chord distance", 
     ylab="Cophenetic distance", asp=1, xlim=c(0,sqrt(2)), ylim=c(0,sqrt(2)),
     main=c("Ward clustering", paste("Cophenetic correlation ",
                                     round(cophcorr.ward2, 3))),
     pch=19, col=gray(0.25,0.25))
abline(0,1)
lines(lowess(dist1, ward2.coph), col="red")

clustObj = list("hclust.single" = hclust.single, 
                "hclust.complete" = hclust.complete,
                "hclust.avg" = hclust.avg,
                "hclust.ward2" = hclust.ward2,
                "hclust.fit" = fit.metrics)
clustObj

#############################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Optimal number of clusters (k) analysis 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#############################################################################

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Fusion Levels (pg 66 in Numerical Ecology)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

names(hclust.avg)
par(mfrow=c(1,2), mar=c(4.5,4.5,1.5,1.5)+0.1)
plot(hclust.avg$height, nrow(coo):2, type="S", 
     main="Fusion levels - Average", ylab="number of clusters",
     xlab="h (node height)", col="gray")
text(hclust.avg$height, nrow(coo):2, nrow(coo):2,
     col="red", cex=0.8)

plot(hclust.single$height, nrow(coo):2, type="S", 
     main="Fusion levels - Single", ylab="number of clusters",
     xlab="h (node height)", col="gray")
text(hclust.single$height, nrow(coo):2, nrow(coo):2,
     col="red", cex=0.8)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Silhouette Widths (pg 70 in Numerical Ecology)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# create empty vector to write silhouette values
asw = numeric(nrow(coo))

# write values
for(k in 2:(nrow(coo)-1)) {
  sil = silhouette(cutree(hclust.avg, k=k), dist1)
  asw[k] = summary(sil)$avg.width
}

#best (larggest Silhouette width)
k.best.silhouette = which.max(asw)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Mantel statistic (Pearson) (pg 71 in Numerical Ecology)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This approach compares the original distance matrix to the binary matricies computed from the dendrogram at various levels of k
grpdist = function(X) {
  require(cluster)
  gr = as.data.frame(as.factor(X))
  distgr = daisy(gr, "gower")
  distgr
}

kt = data.frame(k=1:nrow(coo), r=0) 

for(i in 2:nrow((coo)-1)) {
  gr = cutree(hclust.avg, i)
  distgr = grpdist(gr)
  mt = cor(dist1, distgr, method = "pearson")
  kt[i,2] = mt
}

kt = na.omit(kt)
k.best.mantel = which.max(kt$r)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Optimal number of clusters (k) Plot 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("/Users/jeremewgaeta/Files/DFW/Review_2019_2020/Cluster/integrated_data_analyses/open_water_fall/")

# pdf(file = "K_n_OpenWater_fall_cluster_eval.pdf", width = 11, height=5, paper = "special")
par(mfrow=c(1,2), mar=c(4.5,4.5,1.5,1.5)+0.1)
plot(1:nrow(coo), asw, 
     main="", xlab="Number of Clusters",
     ylab="Average Silhouette Width", type='h', lend=1,
     las=1, lwd=2, ylim=c(0, max(asw)*1.05), yaxs='i')
lines(y = c(0,max(asw)), x=c(k.best.silhouette, k.best.silhouette),
      col="red3", lwd=4, lend=1)
axis(1, k.best.silhouette,  col.axis="red3", col.ticks = "red3", font=2)
axis(1, k.best.silhouette, "(optimum k)", col="red3",
     font=2, col.axis="red3", line = 0.75, tick = FALSE)

plot(kt$k, kt$r, 
     main="", xlab="Number of Clusters",
     ylab=bquote("Pearson's"~italic(rho)), type='h', lend=1,
     las=1, lwd=2, ylim=c(0, max(kt$r)*1.05), yaxs='i')
lines(y = c(0,max(kt$r)), x=c(k.best.mantel, k.best.mantel),
      col="red3", lwd=4, lend=1)
axis(1, k.best.mantel,  col.axis="red3", col.ticks = "red3", font=2)
axis(1, k.best.mantel, "(optimum k)", col="red3",
     font=2, col.axis="red3", line = 0.75, tick = FALSE)
# dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Silhouette plot 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
red_dend_ave <- as.dendrogram(hclust.avg, hang=-1)
red_dend_ave <- rotate(red_dend_ave, 1:length(hclust.avg$labels))

par(mfrow=c(1,1), mar=c(10,4.5,1.5,1.5)+0.1)
#plot(red_dend_ave)
#rect.hclust(hclust.avg, k=9, border=gray(0.2,0.5))

optim_k = 7
cutg = cutree(hclust.avg, k = optim_k)
sil = silhouette(cutg, dist1)
silo = sortSilhouette(sil)
silo_names = data.frame(code=row.names(coo)[attr(silo, "iOrd")])
#silo_names_ave=as.character(merge(silo_names, merged_names, sort=FALSE)[,2])
row.names(silo)=silo_names$code
silo_cluster_label = (silo[,1])

plot(silo)
clust_avg_wd=summary(silo)$clus.avg.widths
total_avg_wd=summary(silo)$avg.width

break_binary = c(0.5)
for(i in 1:(length(silo_cluster_label)-1)){
  if(silo_cluster_label[i]!=silo_cluster_label[i+1]){
    break_binary[i+1] = break_binary[i]+1.5
  } else {
    break_binary[i+1] = break_binary[i]+1.1
  }
}

data.frame(silo_cluster_label[2:length(break_binary)],
           break_binary[2:length(break_binary)],
           brk=diff(break_binary))

sil_space=c(0.1)
for(i in 1:(length(silo_cluster_label)-1)){
  if(silo_cluster_label[i]!=silo_cluster_label[i+1]){
    sil_space[i] = 0.5
  } else {
    sil_space[i] = 0.1
  }
}

# pdf(file = "OpenWater_fall_Silhouette.pdf", width = 7, height=7, paper = "special")
# quartz(width=7, height=7)
par(mfrow=c(1,1), mar=c(4.5,10,2,2)+0.1)
barplot(silo[,3], horiz=TRUE, las=1, xlim=c(range(silo[,3])+
                                              c(-0.01, diff(range(silo[,3]))*0.15)),
        space=c(0.1,sil_space), border=NA, col=NA,
        names.arg=FALSE, xlab="Silhouette Width", ylim=c(min(break_binary),max(break_binary)))
box(which="plot")
abline(v=0)
axis(side = 2, at=break_binary, labels = row.names(silo),
     las=1, tick=F, cex.axis=0.8,font=2)

poly_group=c(which(as.character(diff(break_binary))=="1.5"), length(break_binary))
binary_poly_break = c(break_binary, 100)

for(i in seq(from=1, to=(length(poly_group)), by=2)){
  poly_value1=mean(c(binary_poly_break[poly_group[i]],
                     binary_poly_break[poly_group[i]+1]))
  poly_value2=mean(c(binary_poly_break[poly_group[i+1]],
                     binary_poly_break[poly_group[i+1]+1]))
  polygon(x=c(-2,2,2,-2,-2),
          y=c(poly_value1, poly_value1, poly_value2, poly_value2, poly_value1)+0.1,
          col=gray(0.6, 0.2), border=gray(0.6, 0.2))
}
for(i in 1:length(clust_avg_wd)){
  if(clust_avg_wd[i]!=0){
    x_loc = c(diff(range(silo[,3]))*0.14)+max(silo[,3])
    val=formatC(signif(clust_avg_wd[i],digits=2), digits=2,format="fg", flag="#")
    y_loc = mean(break_binary[which(silo[,1]==i)]) 
    text(x = x_loc, y=y_loc,
         labels = bquote(italic(bar(s))*'='*.(val)),
         adj=c(1,0.5), cex=0.85)
  }
}
par(new=TRUE)
barplot(silo[,3], horiz=TRUE, las=1, xlim=c(range(silo[,3])+
                                              c(-0.01, diff(range(silo[,3]))*0.15)),
        space=c(0.1,sil_space), border="black", names.arg=FALSE,
        xlab="", ylim=c(min(break_binary),max(break_binary)),
        axes=FALSE)
mtext(text = paste("Average Silhouette Width = ",
                   formatC(signif(total_avg_wd,digits=2),
                           digits=2,format="fg", flag="#")),
      font=2, line = 0.25)
# dev.off()

#############################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Plot final dendrogram
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#############################################################################

taxa_dendro = as.dendrogram(hclust.avg)
taxaclusters <- cutree(taxa_dendro, k=optim_k)[order.dendrogram(taxa_dendro)]
taxacol_id = c(hcl.colors(n=9, palette = "Viridis")[c(1,3,7)],
               hcl.colors(n=9, palette = "Heat")[c(5)],
               hcl.colors(n=10, palette = "RdBu")[c(10,2)],
               hcl.colors(n=7, palette = "Dark 2")[c(3,5)],
               hcl.colors(n=10, palette = "BrBG")[c(2,10)],
               hcl.colors(n=10, palette = "PiYG")[c(2,10)])


taxa_dendro = taxa_dendro %>%
  branches_attr_by_clusters(taxaclusters, values = taxacol_id) %>%
  set("branches_lwd", 2) %>%
  set("labels_cex", 0.8)

nclust = optim_k
par(mar=c(1,3,1,1))
par(mfrow=c(1,1), mar=c(10,1.5,1.5,1.5)+0.1)
plot(taxa_dendro, xaxt='n', axes=FALSE)
rect.hclust(hclust.avg, k=nclust, border=gray(0.2,0.2))

clust_3D <- as.factor(cutree(taxa_dendro, k=nclust))
dendro_df = data.frame(taxa = unlist(partition_leaves(taxa_dendro)[1]),
                       dendr_order=1:(length(unlist(partition_leaves(taxa_dendro)[1]))))
taxa_cluster_df=data.frame( taxa = row.names(coo), cluster = clust_3D)
taxa_cluster_df=merge(taxa_cluster_df ,dendro_df)
taxa_cluster_df=taxa_cluster_df[order(taxa_cluster_df$dendr_order),]
clust_col_order=taxa_cluster_df$cluster[!duplicated(taxa_cluster_df$cluster)]

taxa_cluster_df = taxa_cluster_df[order(taxa_cluster_df$cluster),]
taxa_dendro_black=as.dendrogram(hclust.avg)
taxa_dendro_black = taxa_dendro_black %>%
  set("branches_lwd", 1.5) %>%
  set("labels_cex", 0.9) 

# pdf(file = "PTA_rect_dendo.pdf", width = 12, height=7, paper = "special")
# windows(width = 7, height=9)
# quartz(width = 12, height=7)
par(mar=c(1,3,1,1))
par(mfrow=c(1,1), mar=c(12,1.5,1.5,1.5)+0.1)
plot(taxa_dendro_black, xaxt='n', axes=FALSE, lwd=5)
rect.dendrogram(taxa_dendro_black, k=nclust,
                border=hcl.colors(n=nclust+1,
                                  palette = "Viridis")[clust_col_order],
                lwd=2)
# dev.off()

#############################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Plot PTA Biplot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#############################################################################
s.label(coo, xax=1, yax=2, clabel=0.5)
s.label(coo, xax=1, yax=3, clabel=0.5)
taxa_cluster_df
taxa_cluster_names = c("YFG & STRBAS", 
                       "Fresh: Starry+",
                       "DELSME & Shads",
                       "Fresh: Sucker+",
                       "Top & Jacksmelt, Herr.",
                       "Marine: Sardine+",
                       "Longfin")

# pdf(file = "PTA_fall_ordination.pdf", width = 10, height=4.5, paper = "special")
# windows(width = 7, height=9)
# quartz(width = 10, height=4.5)
par(mfrow=c(1,2), mar=c(4,4,1.5,1.5)+0.1)
s.class(coo, fac = clust_3D, xax=1, yax=2,
        col=hcl.colors(n=nclust+1, palette = "Viridis"),
        clabel = 0.75, label=taxa_cluster_names,
        addaxes = TRUE, cgrid = 0)
text(((par('usr')[2]-par('usr')[1])*0.05)+par('usr')[1], 0, labkeep[2], srt=90, xpd=NA, cex=1)
text(0,((par('usr')[4]-par('usr')[3])*0.05)+par('usr')[3], labkeep[1], xpd=NA, cex=1)
s.class(coo, fac = clust_3D, xax=1, yax=3,
        col=hcl.colors(n=nclust+1, palette = "Viridis"),
        clabel = 0.75, label=taxa_cluster_names,
        addaxes = TRUE, cgrid = 0)
text(((par('usr')[2]-par('usr')[1])*0.05)+par('usr')[1], 0,
     labkeep[3], srt=90, xpd=NA, cex=1)
text(0,((par('usr')[4]-par('usr')[3])*0.05)+par('usr')[3],
     labkeep[1], xpd=NA, cex=1)
# dev.off()

#####################################################
#~~~ Visualize Principal tensors
#####################################################

#~~~~~  PT1
pta1_time=pta[[1]]$v[keep[1],]
pta1_region=(pta[[3]]$v[keep[1],])
pta1_mat=matrix(data = NA, nrow=length(pta[[3]]$n), ncol=length(pta[[1]]$n))
colnames(pta1_mat)=pta[[1]]$n
rownames(pta1_mat)=pta[[3]]$n
for(i in 1:length(pta[[3]]$n)){
  pta1_mat[i,] = pta1_region[i]*pta1_time
}

par(mfrow=c(1,1), mar=c(4,4,1,1))
image(t(pta1_mat), col=hcl.colors(9, "BrBg", rev = FALSE), axes=F)
axis(side = 1, at=seq(0,1, length.out=length(pta[[1]]$n)),
     labels=pta[[1]]$n)
axis(side = 2, at=seq(0,1, length.out=length(pta[[3]]$n)),
     labels=c(pta[[3]]$n), las=1, cex.axis=0.85)
box(which="plot")

#~~~~~  PT11
pta2_time=pta[[1]]$v[keep[2],]
pta2_region=(pta[[3]]$v[keep[2],])
pta2_mat=matrix(data = NA, nrow=length(pta[[3]]$n), ncol=length(pta[[1]]$n))
colnames(pta2_mat)=pta[[1]]$n
rownames(pta2_mat)=pta[[3]]$n
for(i in 1:length(pta[[3]]$n)){
  pta2_mat[i,] = pta2_region[i]*pta2_time
}

par(mfrow=c(1,1), mar=c(4,4,1,1))
image(t(pta2_mat), col=hcl.colors(9, "BrBg", rev = FALSE), axes=F)
axis(side = 1, at=seq(0,1, length.out=length(pta[[1]]$n)), labels=pta[[1]]$n)
axis(side = 2, at=seq(0,1, length.out=length(pta[[3]]$n)),
     labels=c(pta[[3]]$n), las=1, cex.axis=0.85)
box(which="plot")

#~~~~~  PT3
pta3_time=pta[[1]]$v[keep[3],]
pta3_region=(pta[[3]]$v[keep[3],])
pta3_mat=matrix(data = NA, nrow=length(pta[[3]]$n), ncol=length(pta[[1]]$n))
colnames(pta3_mat)=pta[[1]]$n
rownames(pta3_mat)=pta[[3]]$n
for(i in 1:length(pta[[3]]$n)){
  pta3_mat[i,] = pta3_region[i]*pta3_time
}

par(mfrow=c(1,1), mar=c(4,4,1,1))
image(t(pta3_mat), col=hcl.colors(9, "BrBg", rev = FALSE), axes=F)
axis(side = 1, at=seq(0,1, length.out=length(pta[[1]]$n)), labels=pta[[1]]$n)
axis(side = 2, at=seq(0,1, length.out=length(pta[[3]]$n)),
     labels=c(pta[[3]]$n), las=1, cex.axis=0.85)
box(which="plot")

plot(pta[[1]]$v[keep[3],] ~ as.numeric(pta[[1]]$n),
     pch=19, cex=2, type='l',xaxt='n', las=1)
axis(side = 1, at = seq(from=1, to=152, by=4), labels=1980:2017)
abline(h=0, col=gray(0.5,0.3))
points(pta[[1]]$v[keep[3],] ~ as.numeric(pta[[1]]$n), 
       col=col_fall3_year[fall3_summary[,1]], pch=19)

#~~~~~~~~~~~~~~~~~~
color.bar <- function(lut, min, max=-min, nticks=11,
                      ticks=seq(min, max, len=nticks),
                      title='') {
  scale = (length(lut))/(max-min)
  plot(c(min,max), c(0,1), type='n', bty='n',
       xaxt='n', xlab='', yaxt='n', ylab='',
       main=title)
  for (i in 1:(length(lut))) {
    x = (i-1)/scale + min
    rect(0,x, 1, x+1/scale, col=lut[i], border=NA)
  }
}

PT_axis = paste0("PT", 1:3, ": ", round((100 * (pta[[3]]$d[keep])^2)/pta[[3]]$ssX[1],1), "%")

# pdf(file = "PTA_fall_heatgrids.pdf", width = 12, height=6.5, paper = "special")
# windows(width = 7, height=9)
# quartz(width = 12, height=6)
par(mfrow=c(1,3), mar=c(4,5.5,2.5,2)+0.1, oma=c(0,0.5,0,7.5))
image(t(pta1_mat), col=hcl.colors(150, "BrBg", rev = FALSE), axes=F)
axis(side = 1, at=seq(0,1, length.out=length(pta[[1]]$n)), labels=pta[[1]]$n)
axis(side = 2, at=seq(0,1, length.out=length(pta[[3]]$n)),
     labels=NA, las=1, cex.axis=0.65)
axis(side = 2, at=seq(0,1, length.out=length(pta[[3]]$n))[c(TRUE, FALSE)], labels=c(pta[[3]]$n[c(TRUE, FALSE)]), las=1, cex.axis=0.65, tick = FALSE)
axis(side = 4, at=seq(0,1, length.out=length(pta[[3]]$n)),
     labels=NA, las=1, cex.axis=0.65)
axis(side = 4, at=seq(0,1, length.out=length(pta[[3]]$n))[c(FALSE, TRUE)], labels=c(pta[[3]]$n[c(FALSE,TRUE)]), las=1, cex.axis=0.65, tick = FALSE)
mtext(text = PT_axis[1], font=2, line = 0.5)
mtext(text = "Time Step", side = 1, line = 2.25, cex=0.85)
box(which="plot")

image(t(pta2_mat), col=hcl.colors(150, "BrBg", rev = FALSE), axes=F)
axis(side = 1, at=seq(0,1, length.out=length(pta[[1]]$n)), labels=pta[[1]]$n)
axis(side = 2, at=seq(0,1, length.out=length(pta[[3]]$n)),
     labels=NA, las=1, cex.axis=0.65)
axis(side = 2, at=seq(0,1, length.out=length(pta[[3]]$n))[c(TRUE, FALSE)], labels=c(pta[[3]]$n[c(TRUE, FALSE)]), las=1, cex.axis=0.65, tick = FALSE)
axis(side = 4, at=seq(0,1, length.out=length(pta[[3]]$n)),
     labels=NA, las=1, cex.axis=0.65)
axis(side = 4, at=seq(0,1, length.out=length(pta[[3]]$n))[c(FALSE, TRUE)], labels=c(pta[[3]]$n[c(FALSE,TRUE)]), las=1, cex.axis=0.65, tick = FALSE)
box(which="plot")
mtext(text = PT_axis[2], font=2, line = 0.5)
mtext(text = "Time Step", side = 1, line = 2.25, cex=0.85)
box(which="plot")

image(t(pta3_mat), col=hcl.colors(150, "BrBg", rev = FALSE), axes=F)
axis(side = 1, at=seq(0,1, length.out=length(pta[[1]]$n)), labels=pta[[1]]$n)
axis(side = 2, at=seq(0,1, length.out=length(pta[[3]]$n)),
     labels=NA, las=1, cex.axis=0.65)
axis(side = 2, at=seq(0,1, length.out=length(pta[[3]]$n))[c(TRUE, FALSE)], labels=c(pta[[3]]$n[c(TRUE, FALSE)]), las=1, cex.axis=0.65, tick = FALSE)
axis(side = 4, at=seq(0,1, length.out=length(pta[[3]]$n)),
     labels=NA, las=1, cex.axis=0.65)
axis(side = 4, at=seq(0,1, length.out=length(pta[[3]]$n))[c(FALSE, TRUE)], labels=c(pta[[3]]$n[c(FALSE,TRUE)]), las=1, cex.axis=0.65, tick = FALSE)
box(which="plot")
mtext(text = PT_axis[3], font=2, line = 0.5)
mtext(text = "Time Step", side = 1, line = 2.25, cex=0.85)
box(which="plot")

par(fig=c(0,1,0,1), mar=c(4.5,4.5,4.5,1.5)+0.1, oma=c(0,81,0,0), new=TRUE)
color.bar(hcl.colors(150, "BrBg", rev = FALSE), min=0, max=1)
mtext(text = bquote(atop(bold(Anomaly),bold(Gradient))),
      side = 3,line = 0.25, adj=0.5, las=0)
# dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# pdf(file = "PT3_fall_temporal.pdf", width = 12, height=7, paper = "special")
# windows(width = 7, height=9)
# quartz(width = 12, height=7)
par(mfrow=c(1,2), mar=c(4,4.5,2.5,4.5)+0.1)
image(t(pta3_mat), col=hcl.colors(150, "BrBg", rev = FALSE), axes=F)
axis(side = 1, at=seq(0,1, length.out=length(pta[[1]]$n)), labels=pta[[1]]$n)
axis(side = 2, at=seq(0,1, length.out=length(pta[[3]]$n)),
     labels=NA, las=1, cex.axis=0.5)
axis(side = 2, at=seq(0,1, length.out=length(pta[[3]]$n))[c(TRUE, FALSE)], labels=c(pta[[3]]$n[c(TRUE, FALSE)]), las=1, cex.axis=0.5, tick = FALSE)
axis(side = 4, at=seq(0,1, length.out=length(pta[[3]]$n)),
     labels=NA, las=1, cex.axis=0.5)
axis(side = 4, at=seq(0,1, length.out=length(pta[[3]]$n))[c(FALSE, TRUE)], labels=c(pta[[3]]$n[c(FALSE,TRUE)]), las=1, cex.axis=0.5, tick = FALSE)
box(which="plot")
mtext(text = PT_axis[3], font=2, line = 0.5)
mtext(text = "Time Step", side = 1, line = 2.25, cex=0.85)
box(which="plot")

par(mar=c(4,6,2.5,1.5)+0.1)
plot(pta[[1]]$v[keep[3],] ~ as.numeric(pta[[1]]$n),
     pch=19, cex=2, type='l',xaxt='n', las=1,
     ylab="Relative Anomaly",
     xlab="")
mtext(text = "Time Step", side = 1, line = 2.25, cex=0.85)
axis(side = 1, at=as.numeric(pta[[1]]$n), las=1, cex=0.7)
abline(h=0, col=gray(0.5,0.3))
mtext(text = PT_axis[3], font=2, line = 0.5)
# dev.off()

