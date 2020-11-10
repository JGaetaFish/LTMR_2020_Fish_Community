
#############################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~    Cluster analysis of demersal, year-round data
#~    Script created on March 19, 2020 by JW Gaeta
#~    Script last modified on October 30, 2020 by JW Gaeta
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

calc_tow_volume=FALSE
fall_only=FALSE

source_url("https://raw.githubusercontent.com/JGaetaFish/LTMR_2020_Fish_Community/main/data_comp_clean_agg.R")


demersal = demersal[which(demersal$yr>1979 & demersal$yr<2018),]
demersal = demersal[order(demersal$lme, demersal$yr, demersal$Survey),]

############################################################################
#~~     STEP 2: Convert surveys into seasons:
#~~                 Spring = March-May
#~~                 Summer = June-Aug
#~~                 Fall   = Sept-Nov
#~~                 Winter = Dec-Feb
#~~         NOTE: too many survey-stations combinations were missing to run 
#~~               the analysis at the survey-station level
#~~               (see the report for details)
############################################################################

season_number = rep(1, times = length(demersal$Method))
season_number[which(demersal$Survey %in% c(3,4,5))]=2
season_number[which(demersal$Survey %in% c(6,7,8))]=3
season_number[which(demersal$Survey %in% c(9,10,11))]=4
demersal = data.frame(data.frame(season_number), demersal)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# assign December into the following year's winter
season_year = demersal$yr
season_year[which(demersal$Survey==12)]=
  as.character(as.numeric(demersal$yr[which(demersal$Survey==12)])+1)
demersal = data.frame(data.frame(season_year), demersal)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~  Create unique season-year variables
season_yr = paste(demersal$season_year, demersal$season_number, sep="_")
demersal = data.frame(data.frame(season_yr), demersal)
demersal = demersal[-which(as.character(demersal$season_yr)=="2018_1"),]

demersal$season_yr = factor(demersal$season_yr)
table(demersal$season_yr)


############################################################################
#~~     STEP 3: Aggregate data into season-station level
############################################################################

demersal2=demersal %>%
  pivot_longer(cols = !season_yr:yr,
               names_to = "CommonName",
               values_to="cpue")

demersal_no_0 = demersal2[which(demersal2$cpue>0),]

relative_sum = function(x){sum(x)/length(x)}

demersal3=demersal2 %>%
  pivot_wider(names_from = CommonName,
              values_from = cpue,
              values_fn = list(cpue = relative_sum),
              id_cols = c(season_yr, season_number,
                          season_year, Method, Source, lme,
                          sta_lme, Station),
              values_fill = list(cpue=0))
demersal3=as.data.frame(demersal3)

############################################################################
#~~     STEP 4: Add zeros for seasons-stations flagged as "No fish caught"
#~~             i.e., instances of zero detection NOT skipped surveys
############################################################################

# NOTE: added 27 season-stations with zero catch, which is 0.35% of season-stations
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~ Bay Study "no_fish_caught"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load Bay Study original dataset and add season
data(Baystudy)
Baystudy$Survey = as.integer(format(Baystudy$Date, format="%m"))
Baystudy$yr = as.integer(format(Baystudy$Date, format="%Y"))
bs_season_number = rep(1, times = length(Baystudy$Method))
bs_season_number[which(Baystudy$Survey %in% c(3,4,5))]=2
bs_season_number[which(Baystudy$Survey %in% c(6,7,8))]=3
bs_season_number[which(Baystudy$Survey %in% c(9,10,11))]=4
Baystudy = data.frame(data.frame(season_number=bs_season_number), Baystudy)

bs_season_year = Baystudy$yr
bs_season_year[which(Baystudy$Survey==12)]=
  as.character(as.numeric(Baystudy$yr[which(Baystudy$Survey==12)])+1)
Baystudy$season_year=bs_season_year

bs_season_yr = paste(Baystudy$season_year, Baystudy$season_number, sep="_")
Baystudy = data.frame(data.frame(season_yr=bs_season_yr), Baystudy)

# Identify survey events without fish detected
BayStudyOTR = subset(Baystudy, Baystudy$Method=="Otter trawl")
bs_no_fish = BayStudyOTR[which(BayStudyOTR$Length_NA_flag=="No fish caught" &
                                 BayStudyOTR$Station %in% bs_sta_above_thresh),]

# Aggregate seasonally; compare total surveys per season per station with
#       total surveys per season without fish captured
bs_survey_year_station_no_fish = aggregate(yr ~ Source + season_number  + season_yr + Survey 
                                           + season_year + Station, data = bs_no_fish, FUN = length)
bs_survey_year_station = aggregate(yr ~ Source + season_number  + season_yr + Survey 
                                   + season_year + Station, data = BayStudyOTR, FUN = length)

bs_surveys_per_season = aggregate(Survey ~ Source + season_number  + season_yr +  
                                    season_year + Station, data = bs_survey_year_station, FUN = length)
colnames(bs_surveys_per_season)[length(names(bs_surveys_per_season))] = "total_surveys_per_season"

bs_no_fish_surveys_per_season = aggregate(Survey ~ Source + season_number  + season_yr +  
                                            season_year + Station,
                                          data = bs_survey_year_station_no_fish, FUN = length)
colnames(bs_no_fish_surveys_per_season)[length(names(bs_no_fish_surveys_per_season))] = 
  "surveys_without_fish"

bs_no_fish_merged = merge(bs_surveys_per_season, bs_no_fish_surveys_per_season)

bs_zero_catch_index = which(bs_no_fish_merged$surveys_without_fish ==
                              bs_no_fish_merged$total_surveys_per_season)
bs_add_zeros = bs_no_fish_merged[bs_zero_catch_index,]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~     FLAG: Code gets angry about surveys with two dates
#~~
#~  BS 102 2017 Summer, BS 317 2017 Summer, and BS 428 1981 Fall have catch 
#~  data in the LTMRdata

#check on error with 2017_3 bs_102; remove from "no fish caught dataset
demersal3[which(demersal3$season_yr=="2017_3" & demersal3$Station== 102),]
bs_add_zeros[which(bs_add_zeros$season_yr=="2017_3" & bs_add_zeros$Station== 102),]
bs_add_zeros = bs_add_zeros[-which(bs_add_zeros$season_yr=="2017_3" & bs_add_zeros$Station== 102),]

#check on error with 2017_3 bs_317; remove from "no fish caught dataset
demersal3[which(demersal3$season_yr=="2017_3" & demersal3$Station== 317),]
bs_add_zeros[which(bs_add_zeros$season_yr=="2017_3" & bs_add_zeros$Station== 317),]
bs_add_zeros = bs_add_zeros[-which(bs_add_zeros$season_yr=="2017_3" & bs_add_zeros$Station== 317),]

#check on error with 1981_4 bs_428; remove from "no fish caught dataset
demersal3[which(demersal3$season_yr=="1981_4" & demersal3$Station == 428),]
bs_add_zeros[which(bs_add_zeros$season_yr=="1981_4" & bs_add_zeros$Station == 428),]
bs_add_zeros = bs_add_zeros[-which(bs_add_zeros$season_yr=="1981_4" & bs_add_zeros$Station == 428),]

#~~     END FLAG:
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Add events without fish detected to the original "sea" dataset adding 0s for all fishes
#~~  note 26 BS season-stations were 0 detections; representing 0.50% of total Bay Study season-stations
bs_no_fish_caught = data.frame(season_yr = bs_add_zeros$season_yr,
                               season_number = bs_add_zeros$season_number, 
                               season_year = bs_add_zeros$season_year,
                               Method = rep("Otter trawl", dim(bs_add_zeros)[1]),
                               Source = bs_add_zeros$Source,
                               lme = rep("bs", dim(bs_add_zeros)[1]),
                               sta_lme = paste(rep("bs", dim(bs_add_zeros)[1]),
                                               bs_add_zeros$Station, sep="_"),
                               Station = bs_add_zeros$Station)

for(i in 9:length(names(demersal3))){
  bs_no_fish_caught = cbind(bs_no_fish_caught,rep(0, times=dim(bs_no_fish_caught)[1]))
  colnames(bs_no_fish_caught)[i]=names(demersal3)[i]
}
demersal3 = rbind(demersal3,bs_no_fish_caught)
demersal3 = demersal3[order(demersal3$sta_lme, demersal3$season_year),]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~ UC Davis Suisun Marsh Study "no_fish_caught"
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data(Suisun)
Suisun$Survey = as.integer(format(Suisun$Date, format="%m"))
Suisun$yr = as.integer(format(Suisun$Date, format="%Y"))
ucd_season_number = rep(1, times = length(Suisun$Method))
ucd_season_number[which(Suisun$Survey %in% c(3,4,5))]=2
ucd_season_number[which(Suisun$Survey %in% c(6,7,8))]=3
ucd_season_number[which(Suisun$Survey %in% c(9,10,11))]=4
Suisun = data.frame(data.frame(season_number=ucd_season_number), Suisun)

ucd_season_year = Suisun$yr
ucd_season_year[which(Suisun$Survey==12)]=
  as.character(as.numeric(Suisun$yr[which(Suisun$Survey==12)])+1)
Suisun$season_year=ucd_season_year

ucd_season_yr = paste(Suisun$season_year, Suisun$season_number, sep="_")
Suisun = data.frame(data.frame(season_yr=ucd_season_yr), Suisun)


# Identify survey events without fish detected

ucd_no_fish = Suisun[which(Suisun$Length_NA_flag=="No fish caught" &
                             Suisun$Station %in% ucd_sta_above_thresh),]

# Aggregate seasonally; compare total surveys per season per station with
#       total surveys per season without fish captured
ucd_survey_year_station_no_fish = aggregate(yr ~ Source + season_number  + season_yr + Survey 
                                            + season_year + Station, data = ucd_no_fish, FUN = length)
ucd_survey_year_station = aggregate(yr ~ Source + season_number  + season_yr + Survey 
                                    + season_year + Station, data = Suisun, FUN = length)

ucd_surveys_per_season = aggregate(Survey ~ Source + season_number  + season_yr +  
                                     season_year + Station, data = ucd_survey_year_station, FUN = length)
colnames(ucd_surveys_per_season)[length(names(ucd_surveys_per_season))] = "total_surveys_per_season"

ucd_no_fish_surveys_per_season = aggregate(Survey ~ Source + season_number  + season_yr +  
                                             season_year + Station,
                                           data = ucd_survey_year_station_no_fish, FUN = length)
colnames(ucd_no_fish_surveys_per_season)[length(names(ucd_no_fish_surveys_per_season))] = 
  "surveys_without_fish"

ucd_no_fish_merged = merge(ucd_surveys_per_season, ucd_no_fish_surveys_per_season)

ucd_zero_catch_index = which(ucd_no_fish_merged$surveys_without_fish == ucd_no_fish_merged$total_surveys_per_season)
ucd_add_zeros = ucd_no_fish_merged[ucd_zero_catch_index,]

# Add events without fish detected to the original "sea" dataset adding 0s for all fishes
#~~  note 1 season-stations had 0 detections; representing 0.04% of total season-stations
ucd_no_fish_caught = data.frame(season_yr = ucd_add_zeros$season_yr,
                                season_number = ucd_add_zeros$season_number, 
                                season_year = ucd_add_zeros$season_year,
                                Method = rep("Otter trawl", dim(ucd_add_zeros)[1]),
                                Source = ucd_add_zeros$Source,
                                lme = rep("ucd", dim(ucd_add_zeros)[1]),
                                sta_lme = paste(rep("ucd", dim(ucd_add_zeros)[1]),
                                                ucd_add_zeros$Station, sep="_"),
                                Station = ucd_add_zeros$Station)

for(i in 9:length(names(demersal3))){
  ucd_no_fish_caught = cbind(ucd_no_fish_caught,rep(0, times=dim(ucd_no_fish_caught)[1]))
  colnames(ucd_no_fish_caught)[i]=names(demersal3)[i]
}
demersal3 = rbind(demersal3,ucd_no_fish_caught)
demersal3 = demersal3[order(demersal3$sta_lme, demersal3$season_year),]


############################################################################
#~~     STEP 5: Restructure data into an array framework for 
#~~             Principal Tensor Analysis
############################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ create unique variables for season-station spatiotemporal resolution
season_count_df = data.frame(season_count = c(1:length(unique(demersal3$season_yr))),
                             season_yr=sort(unique(demersal3$season_yr)))
demersal3 = merge(season_count_df, demersal3, all=TRUE)

test = paste(demersal3$season_yr, demersal3$sta_lme, sep="_")
test[duplicated(test)]


sea_lme = unique(demersal3$lme)

demersal_scaled = demersal3
sea_lme = unique(demersal3$lme)




for(i in 1:length(sea_lme)){
  sub = demersal3[which(demersal3$lme==sea_lme[i]),]
  sub2 = sub[,-c(1:9)]
  sub2 = (sub2)^(1/3) #~ Increase normality by cube-root transformation
  sub_scale=apply(sub2, MARGIN = 2,
                  FUN = function(x){
                    (x-mean(x))/sd(x) #~ Standardize by z-scoring per taxa per study-gear
                  }
  )
  sub_scale[is.na(sub_scale)]=0
  demersal_scaled[which(demersal3$lme==sea_lme[i]),-c(1:9)]=sub_scale
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Convert data into array format
#~    Array dimensions: [season_timeseries_number, taxa, station]

demersal_stations = sort(unique(demersal3$sta_lme))
season_seq=1:max(demersal3$season_count)
taxa=colnames(demersal3[,-c(1:9)])

demersal_scaled_arr=array(data = NA, dim = c(length(season_seq),
                                        length(taxa),
                                        length(demersal_stations)),
                     dimnames = list("ts_sea"=season_seq,
                                     "taxa"=taxa,
                                     "station"=demersal_stations)
)
dim(demersal_scaled_arr)

for(i in 1:length(demersal_stations)){
  sub = demersal_scaled[which(demersal_scaled$sta_lme==demersal_stations[i]),]
  sub2 = sub[,-c(1:9)]
  row_order=intersect(season_seq , sub$season_count)
  for(x in 1:dim(sub2)[2]){
    demersal_scaled_arr[row_order,x,i] = sub2[,x]
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ identify survey count timeseries without all stations sampled

vec = c()
na_count = c()
for(i in 1:dim(demersal_scaled_arr)[3]){
  sub = demersal_scaled_arr[,,i]
  sub_na_element = c(1:dim(demersal_scaled_arr)[1])[is.na(sub[,1])]
  vec=c(vec, sub_na_element)
  na_count = c(na_count, length(sub_na_element))
}
na_count
vec=na.omit(vec)
unique(vec)

missing_survey_seasons = sort(unique(vec))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~  Write reduced seasonal array (cutting 47 seasons or 152) (-31%):
dim(demersal_scaled_arr)
red_demersal_array = demersal_scaled_arr[-missing_survey_seasons,,]

length(missing_survey_seasons)/dim(demersal_scaled_arr)[1]

############################################################################
#~~     STEP 5: Perform Principal tensor analysis 
#~~             Methods based on:
#~~               Leibovici. 2010. Journal of Statistical Software
#~~               ichocki et al. 2015. IEEE Signal Processing Magazine 
#~~               Frelat et al. 2017. PlosOne
############################################################################

rm(list=c("pta"))
pta<-PTA3(red_demersal_array, nbPT = 3, nbPT2 = 3, minpct = 0.1)
summary.PTAk(pta,testvar = 0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Use Scree plot to identify dominant pricipal tensors
out <- !substr(pta[[3]]$vsnam, 1, 1) == "*"
gct<-(pta[[3]]$pct*pta[[3]]$ssX/pta[[3]]$ssX[1])[out]
par(mfrow=c(1,1), mar=c(4.5,4.5,1.5,1.5)+0.1)
barplot(sort(gct, decreasing = TRUE), xlab="PT",
        ylab="Percentage of variance", las=1); box(which="plot")
keep = c(1,3,11,4)

demersal_summary = aggregate(season_year ~ season_number + season_count + season_year, FUN=unique, data=demersal3)

demersal_summary$yr=substr(as.character(demersal_summary$season_year), start=3, stop=4)
demersal_summary = subset(demersal_summary, demersal_summary$season_count %in% as.numeric(pta[[1]]$n))

#Create the matrix with the projection of species on the PTs
coo<-t(pta[[2]]$v[c(keep),])

rownames(coo) = pta[[2]]$n
PT_associations = c("PT1: Spatial (region)", "PT2: Spatial (sub-region)", 
                    "PT3: Temporal (season & year)",
                    "PT4: Spatial (sub-region)")
labkeep <- paste0(PT_associations, "\n(",pta[[3]]$vsnam[keep], " - ", round((100 * (pta[[3]]$d[keep])^2)/pta[[3]]$ssX[1],1), "%)")

############################################################################
#~~     STEP 6: Evaluate & select optimal clusting algorithm
#~~             Methods based on Sofaer et al. 2019. Diversity and Distributions
############################################################################

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####  Cophenetic Correlation: compares clustering methods   ####
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####   Gower (1983) distances: to compare clustering methods   ####
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

############################################################################
#~~     STEP 7: Determine optimal number of clusters based on Silhouette width
#~~             Methods based on Sofaer et al. 2019. Diversity and Distributions
#~~             Refer to Rousseeuw 1987. J. of Computational and Applied Math.
############################################################################

# create empty vector to write silhouette values
asw = numeric(nrow(coo))

# write values
for(k in 2:(nrow(coo)-1)) {
  sil = silhouette(cutree(hclust.avg, k=k), dist1)
  asw[k] = summary(sil)$avg.width
}

#best (largest Silhouette width)
k.best.silhouette = which.max(asw)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~   Alternative approach to find fewer clusters based on local maxima 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#find local maxima with 'm' points on either side
find_peaks <- function (x, m = 2){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}

asw_maxima=find_peaks(asw)

# pdf(file = "K_n_demersal_all_season_cluster_eval.pdf", width = 11, height=5, paper = "special")
par(mfrow=c(1,2), mar=c(4.5,4.5,1.5,1.5)+0.1)
plot(1:nrow(coo), asw, 
     main="", xlab="Number of Clusters",
     ylab="Average Silhouette Width", type='h', lend=1,
     las=1, lwd=2, ylim=c(0, max(asw)*1.05), yaxs='i')
lines(y = c(0,max(asw)), x=c(k.best.silhouette, k.best.silhouette),
      col="red3", lwd=4, lend=1)
axis(1, k.best.silhouette,  col.axis="red3", col.ticks = "red3",
     font=2, labels=FALSE)
axis(1, k.best.silhouette,  col.axis="red3", col.ticks = "red3", font=2,
     line = -0.7)
axis(1, k.best.silhouette+1, "(optimum k)", col="red3",
     font=2, col.axis="red3", line = 0.75, tick = FALSE)
lines(y = c(0,asw[asw_maxima[2]]), x=c(asw_maxima[2], asw_maxima[2]),
      col="darkgreen", lwd=4, lend=1)
axis(1, asw_maxima[2],  col.axis="darkgreen", col.ticks = "darkgreen",
     font=2, labels=FALSE)
axis(1, asw_maxima[2],  col.axis="darkgreen", col.ticks = "darkgreen",
     font=2, line=-0.7)
axis(1, asw_maxima[2]-1, "(alternative k)", col="darkgreen",
     font=2, col.axis="darkgreen", line = 0.75, tick = FALSE)
# dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Silhouette plot 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
red_dend_ave <- as.dendrogram(hclust.avg, hang=-1)
red_dend_ave <- rotate(red_dend_ave, 1:length(hclust.avg$labels))

optim_k = 12
cutg = cutree(hclust.avg, k = optim_k)
sil = silhouette(cutg, dist1)
silo = sortSilhouette(sil)
silo_names = data.frame(code=row.names(coo)[attr(silo, "iOrd")])
row.names(silo)=silo_names$code

par(mfrow=c(1,1), mar=c(6,8,1.5,1.5)+0.1, oma=c(0,0,0,0))
plot(silo, max.strlen = 25, cex.names = 0.5, main="", nmax.lab = 50)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Dendrogram plot 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
taxa_dendro = as.dendrogram(hclust.avg)
nclust = optim_k

clust_3D <- as.factor(cutree(taxa_dendro, k=nclust))
dendro_df = data.frame(taxa = unlist(partition_leaves(taxa_dendro)[1]),
                       dendr_order=1:(length(unlist(partition_leaves(taxa_dendro)[1]))))
taxa_cluster_df=data.frame( taxa = row.names(coo), cluster = clust_3D)
taxa_cluster_df=merge(taxa_cluster_df ,dendro_df)
taxa_cluster_df=taxa_cluster_df[order(taxa_cluster_df$dendr_order),]
clust_col_order=taxa_cluster_df$cluster[!duplicated(taxa_cluster_df$cluster)]
taxa_dendro_black=as.dendrogram(hclust.avg)
taxa_dendro_black = taxa_dendro_black %>%
  set("branches_lwd", 1.5) %>%
  set("labels_cex", 0.9) 

par(mfrow=c(1,1), mar=c(12,1.5,1.5,1.5)+0.1)
plot(taxa_dendro_black, xaxt='n', axes=FALSE, lwd=5)
rect.dendrogram(taxa_dendro_black, k=nclust,
                border=hcl.colors(n=nclust+1,
                                  palette = "Viridis")[clust_col_order],
                lwd=2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~   Biplot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

taxa_cluster_names = c("Sturgeon; DS+", 
                       "Silverside+",
                       "Wht. Catfish",
                       "Gunnel+",
                       "Tomcod+", #5
                       "Sole+",
                       "Splittail",
                       "C. Gobies",
                       "Halibut & Tonguefish",
                       "Starry",
                       "Longfin",
                       "Sh. Gobies")

par(mfrow=c(1,3), mar=c(4,4,1.5,1.5)+0.1)
s.class(coo, fac = clust_3D, xax=1, yax=2,
        col=hcl.colors(n=nclust+1, palette = "Viridis"),
        clabel = 1, label=taxa_cluster_names,
        addaxes = TRUE, cgrid = 0)
text(((par('usr')[2]-par('usr')[1])*0.05)+par('usr')[1], 0, labkeep[2], srt=90, xpd=NA, cex=1)
text(0,((par('usr')[4]-par('usr')[3])*0.05)+par('usr')[3], labkeep[1], xpd=NA, cex=1)
s.class(coo, fac = clust_3D, xax=1, yax=3,
        col=hcl.colors(n=nclust+1, palette = "Viridis"),
        clabel = 1, label=taxa_cluster_names,
        addaxes = TRUE, cgrid = 0)
text(((par('usr')[2]-par('usr')[1])*0.05)+par('usr')[1], 0,
     labkeep[3], srt=90, xpd=NA, cex=1)
text(0,((par('usr')[4]-par('usr')[3])*0.05)+par('usr')[3],
     labkeep[1], xpd=NA, cex=1)
s.class(coo, fac = clust_3D, xax=1, yax=4,
        col=hcl.colors(n=nclust+1, palette = "Viridis"),
        clabel = 1, label=taxa_cluster_names,
        addaxes = TRUE, cgrid = 0)
text(((par('usr')[2]-par('usr')[1])*0.05)+par('usr')[1], 0,
     labkeep[4], srt=90, xpd=NA, cex=1)
text(0,((par('usr')[4]-par('usr')[3])*0.05)+par('usr')[3],
     labkeep[1], xpd=NA, cex=1)

#~~~~~~~~~~~~~~~~~~~~
par(mfrow=c(1,2), mar=c(4,4,2,2)+0.1, oma=c(0,0,0,0))
s.class(coo, fac = clust_3D, xax=1, yax=3,
        col=hcl.colors(n=nclust+1, palette = "Viridis"),
        clabel = 1, label=taxa_cluster_names,
        addaxes = TRUE, cgrid = 0)
text(((par('usr')[2]-par('usr')[1])*0.05)+par('usr')[1], 0,
     labkeep[3], srt=90, xpd=NA, cex=1)
text(0,((par('usr')[4]-par('usr')[3])*0.05)+par('usr')[3],
     labkeep[1], xpd=NA, cex=1)

s.label(coo, xax=1, yax=3,clabel = 0.5)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~ Visualize Principal tensors
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
col_season = c("darkblue", "darkgreen", "darkred", "tan4")

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
image(t(pta1_mat), col=hcl.colors(50, "BrBg", rev = FALSE), axes=F)
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
image(t(pta2_mat), col=hcl.colors(50, "BrBg", rev = FALSE), axes=F)
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
image(t(pta3_mat), col=hcl.colors(50, "BrBg", rev = FALSE), axes=F)
axis(side = 1, at=seq(0,1, length.out=length(pta[[1]]$n)), labels=pta[[1]]$n)
axis(side = 2, at=seq(0,1, length.out=length(pta[[3]]$n)),
     labels=c(pta[[3]]$n), las=1, cex.axis=0.85)
box(which="plot")

plot(pta[[1]]$v[keep[3],] ~ as.numeric(pta[[1]]$n),
     pch=19, cex=2, type='l',xaxt='n', las=1)
axis(side = 1, at = seq(from=1, to=152, by=4), labels=1980:2017)
abline(h=0, col=gray(0.5,0.3))
points(pta[[1]]$v[keep[3],] ~ as.numeric(pta[[1]]$n), 
       col=col_season[demersal_summary[,1]], pch=19)
legend("topleft", legend = c("Spring", "Summer", "Fall", "Winter"),
       pch=19, col = col_season[c(2:4,1)], bty='n', inset=0.025)

#~~~~~  PT4
pta4_time=pta[[1]]$v[keep[4],]
pta4_region=(pta[[3]]$v[keep[4],])
pta4_mat=matrix(data = NA, nrow=length(pta[[3]]$n), ncol=length(pta[[1]]$n))
colnames(pta4_mat)=pta[[1]]$n
rownames(pta4_mat)=pta[[3]]$n
for(i in 1:length(pta[[3]]$n)){
  pta4_mat[i,] = pta4_region[i]*pta4_time
}

par(mfrow=c(1,1), mar=c(4,4,1,1))
image(t(pta4_mat), col=hcl.colors(50, "BrBg", rev = FALSE), axes=F)
axis(side = 1, at=seq(0,1, length.out=length(pta[[1]]$n)), labels=pta[[1]]$n)
axis(side = 2, at=seq(0,1, length.out=length(pta[[3]]$n)),
     labels=c(pta[[3]]$n), las=1, cex.axis=0.85)
box(which="plot")
#~~~~~~~~~~~~~~~~~~
