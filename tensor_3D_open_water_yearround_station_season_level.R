
#############################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~    Cluster analysis of open water, year-round data
#~    Script created on March 19, 2020 by JW Gaeta
#~    Script last modified on October 28, 2020 by JW Gaeta
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#############################################################################


# rm(list=ls())
# graphics.off()
# cat("\f")

library(vegan)
library(dendextend)
library(circlize)
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
fall_only=FALSE

source_url("https://raw.githubusercontent.com/JGaetaFish/LTMR_2020_Fish_Community/main/data_comp_clean_agg.R")

all_sea = all_sea[which(all_sea$yr>1979 & all_sea$yr<2018),]
all_sea = all_sea[order(all_sea$lme, all_sea$yr, all_sea$Survey),]

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Add December into the subsequent year's winter
season_year = all_sea$yr
season_year[which(all_sea$Survey==12)]=
  as.character(as.numeric(all_sea$yr[which(all_sea$Survey==12)])+1)
all_sea = data.frame(data.frame(season_year), all_sea)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~  Classify seasons
season_number = rep(1, times = length(all_sea$Method))
season_number[which(all_sea$Survey %in% c(3,4,5))]=2
season_number[which(all_sea$Survey %in% c(6,7,8))]=3
season_number[which(all_sea$Survey %in% c(9,10,11))]=4
all_sea = data.frame(data.frame(season_number), all_sea)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~  Create unique season-year variables
season_yr = paste(all_sea$season_year, all_sea$season_number, sep="_")
all_sea = data.frame(data.frame(season_yr), all_sea)
all_sea = all_sea[-which(as.character(all_sea$season_yr)=="2018_1"),]

all_sea$season_yr = factor(all_sea$season_yr)
table(all_sea$season_yr)

season_count_df=data.frame(season_count = c(1:length(unique(all_sea$season_yr))),
           season_yr=sort(unique(all_sea$season_yr)))

all_sea=merge(season_count_df, all_sea, all=TRUE)


############################################################################
#~~     STEP 3: Aggregate data into season-station level
############################################################################

all_sea2=all_sea %>%
  pivot_longer(cols = !season_yr:yr,
              names_to = "CommonName",
              values_to="cpue")

all_sea_no_0 = all_sea2[which(all_sea2$cpue>0),]

study_obs_diff=all_sea_no_0 %>%
  pivot_wider(names_from = lme,
              values_from = cpue,
              values_fn = list(cpue = length),
              id_cols = c(CommonName),
              values_fill = list(cpue=0))
study_obs_diff = data.frame(study_obs_diff)

study_obs_diff$perc_bs = round((study_obs_diff[,2]/rowSums(study_obs_diff[,2:3]))*100, digits=1)
study_obs_diff$perc_ucd = round((study_obs_diff[,3]/rowSums(study_obs_diff[,2:3]))*100, digits=1)
study_obs_diff$CommonName = gsub("[.]", " ", study_obs_diff$CommonName)

relative_sum = function(x){sum(x)/length(x)}

sea=all_sea2 %>%
  pivot_wider(names_from = CommonName,
              values_from = cpue,
              values_fn = list(cpue = relative_sum),
              id_cols = c(season_yr, season_count, season_number,
                          season_year, Method, Source, lme,
                          sta_lme, Station),
              values_fill = list(cpue=0))
sea=as.data.frame(sea)

############################################################################
#~~     STEP 4: Restructure data into an array framework for 
#~~             Principal Tensor Analysis
############################################################################

sea_lme = unique(sea$lme)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Increase normality by cube-root transformation
sub2 = sea[,-c(1:9)]
sub2 = (sub2)^(1/3)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Standardize by z-scoring per taxa per study-gear
sub_scale=apply(sub2, MARGIN = 2,
                FUN = function(x){
                  (x-mean(x))/sd(x)
                }
)
sub_scale[is.na(sub_scale)]=0
sea_scaled=cbind(sea[,c(1:9)], sub_scale)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Convert data into array format
#~    Array dimensions: [season_timeseries_number, taxa, station]

sea_stations = sort(unique(sea$sta_lme))
season_seq=1:max(sea$season_count)
taxa=colnames(sea[,-c(1:9)])

sea_scaled_arr=array(data = NA, dim = c(length(season_seq),
                         length(taxa),
                         length(sea_stations)),
      dimnames = list("ts_sea"=season_seq,
                      "taxa"=taxa,
                      "station"=sea_stations)
      )
dim(sea_scaled_arr)
for(i in 1:length(sea_stations)){
  sub = sea_scaled[which(sea_scaled$sta_lme==sea_stations[i]),]
  sub2 = sub[,-c(1:9)]
  row_order=intersect(season_seq , sub$season_count)
  for(x in 1:dim(sub2)[2]){
    sea_scaled_arr[row_order,x,i] = sub2[,x]
  }
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ identify survey count timeseries without all stations sampled

vec = c()
na_count = c()
for(i in 1:dim(sea_scaled_arr)[3]){
  sub = sea_scaled_arr[,,i]
  sub_na_element = c(1:dim(sea_scaled_arr)[1])[is.na(sub[,1])]
  vec=c(vec, sub_na_element)
  na_count = c(na_count, length(sub_na_element))
}
na_count
vec=na.omit(vec)
unique(vec)

missing_survey_seasons = sort(unique(vec))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~  Write reduced seasonal array
dim(sea_scaled_arr)
red_sea_arr = sea_scaled_arr[-missing_survey_seasons,,]

length(missing_survey_seasons)/dim(sea_scaled_arr)[1]

############################################################################
#~~     STEP 5: Perform Principal tensor analysis 
#~~             Methods based on:
#~~               Leibovici. 2010. Journal of Statistical Software
#~~               ichocki et al. 2015. IEEE Signal Processing Magazine 
#~~               Frelat et al. 2017. PlosOne
############################################################################

rm(list=c("pta"))
pta<-PTA3(red_sea_arr, nbPT = 3, nbPT2 = 3, minpct = 0.1)
summary.PTAk(pta,testvar = 0)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Perform Scree plot evaluation to identify dominant pricipal tensors
out <- !substr(pta[[3]]$vsnam, 1, 1) == "*"
gct<-(pta[[3]]$pct*pta[[3]]$ssX/pta[[3]]$ssX[1])[out]
par(mfrow=c(1,1), mar=c(4.5,4.5,1.5,1.5)+0.1)
barplot(sort(gct, decreasing = TRUE), xlab="PT",
        ylab="Percentage of variance", las=1); box(which="plot")
keep = c(1,3,11)


sea_summary = aggregate(season_year ~ season_number + season_count + season_year, FUN=unique, data=sea)
sea_summary$yr=substr(as.character(sea_summary$season_year), start=3, stop=4)
sea_summary = subset(sea_summary, sea_summary$season_count %in% as.numeric(pta[[1]]$n))

#Create the matrix with the projection of species on the 4 PT
coo<-t(pta[[2]]$v[c(keep),])

rownames(coo) = pta[[2]]$n
PT_associations = c("PT1: Spatial (region)", "PT2: Spatial (sub-region)", 
                    "PT3: Temporal (season & year)")
labkeep <- paste0(PT_associations, "\n(",pta[[3]]$vsnam[keep], " - ", round((100 * (pta[[3]]$d[keep])^2)/pta[[3]]$ssX[1],1), "%)")



############################################################################
#~~     STEP 6: Evaluate & select optimal clusting algorithm
#~~             Methods based on Sofaer et al. 2019. Diversity and Distributions
############################################################################

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####  Cophenetic Correlation: compares clustering methods   ####
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#1. Compute the euclidean distance between species
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
find_peaks <- function (x, m = 3){
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
par(mfrow=c(1,1), mar=c(4.5,4.5,1.5,1.5)+0.1)
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
axis(1, k.best.silhouette+1.5, "(optimum k)", col="red3",
     font=2, col.axis="red3", line = 0.75, tick = FALSE)
lines(y = c(0,asw[asw_maxima[1]]), x=c(asw_maxima[1], asw_maxima[1]),
      col="darkgreen", lwd=4, lend=1)
axis(1, asw_maxima[1],  col.axis="darkgreen", col.ticks = "darkgreen",
     font=2, labels=FALSE)
axis(1, asw_maxima[1],  col.axis="darkgreen", col.ticks = "darkgreen",
     font=2, line=-0.7)
axis(1, asw_maxima[1]-1.5, "(alternative k)", col="darkgreen",
     font=2, col.axis="darkgreen", line = 0.75, tick = FALSE)
# dev.off()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Silhouette plot 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
red_dend_ave <- as.dendrogram(hclust.avg, hang=-1)
red_dend_ave <- rotate(red_dend_ave, 1:length(hclust.avg$labels))

optim_k = 7
cutg = cutree(hclust.avg, k = optim_k)
sil = silhouette(cutg, dist1)
silo = sortSilhouette(sil)
silo_names = data.frame(code=row.names(coo)[attr(silo, "iOrd")])
row.names(silo)=silo_names$code

par(mfrow=c(1,1), mar=c(6,8,1.5,1.5)+0.1)
plot(silo, max.strlen = 25, cex.names = 0.6, main="")

red_dend_ave <- as.dendrogram(hclust.avg, hang=-1)
red_dend_ave <- rotate(red_dend_ave, 1:length(hclust.avg$labels))

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

taxa_cluster_names = c("Freshwater", 
                       "Shads",
                       "Marine: Sardine+",
                       "Anchovy & Jacksmelt",
                       "Marine: Croaker+",
                       "ESA Smelts",
                       "Staghorn+")

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

#####################################################
#~~~ Visualize Principal tensors
#####################################################
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
       col=col_season[sea_summary[,1]], pch=19)
legend("bottomleft", legend = c("Spring", "Summer", "Fall", "Winter"),
       pch=19, col = col_season[c(2:4,1)], bty='n', inset=0.025)

