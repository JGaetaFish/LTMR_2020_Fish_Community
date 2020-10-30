
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

season_count_df=data.frame(season_count = c(1:length(unique(demersal$season_yr))),
                           season_yr=sort(unique(demersal$season_yr)))

demersal=merge(season_count_df, demersal, all=TRUE)

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
              id_cols = c(season_yr, season_count, season_number,
                          season_year, Method, Source, lme,
                          sta_lme, Station),
              values_fill = list(cpue=0))
demersal3=as.data.frame(demersal3)

##################################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Add stations-years with no catch of any fishes

#~~~~~~~ Bay Study "no_fish_caught"
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
#Baystudy = Baystudy[-which(as.character(Baystudy$season_yr)=="2018_1"),]

BayStudyOTR = subset(Baystudy, Baystudy$Method=="Otter trawl")
bs_no_fish = BayStudyOTR[which(BayStudyOTR$Length_NA_flag=="No fish caught" & BayStudyOTR$Station %in% bs_sta_above_thresh),]
bs_no_fish_agg=aggregate(Temp_surf ~ Source + season_yr + season_year + season_number + Station, 
                         data=bs_no_fish, FUN = length)
test=bs_no_fish_agg[which(bs_no_fish_agg$Station=="101"),]
test[order(test$season_yr),]
BayStudyOTR[which(BayStudyOTR$Station=="101" & BayStudyOTR$season_year==2009),]

bs_no_fish_agg = bs_no_fish_agg[which(bs_no_fish_agg$season_year>1984 & bs_no_fish_agg$season_year<2018),]
bs_no_fish_agg$lme = rep("bs", times=dim(bs_no_fish_agg)[1])
bs_no_fish_agg$sta_lme = paste(bs_no_fish_agg$lme, bs_no_fish_agg$Station, sep="_")
bs_no_fish_agg$sta_lme_yr = paste(bs_no_fish_agg$sta_lme, bs_no_fish_agg$season_year, sep="_")

open_fish_caught = paste(demersal3$sta_lme, demersal3$season_year, sep="_")
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#NOTE: all year-seasons have data for BS-OTR
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
bs_no_fish_agg = bs_no_fish_agg[-which(bs_no_fish_agg$sta_lme_yr %in% unique(open_fish_caught)),]



# bs_no_fish_caught = data.frame(season_year = bs_no_fish_agg$season_year, lme = bs_no_fish_agg$lme,
#                                sta_lme = bs_no_fish_agg$sta_lme)
# 
# for(i in 4:length(names(fall3))){
#   bs_no_fish_caught = cbind(bs_no_fish_caught,rep(0, times=dim(bs_no_fish_caught)[1]))
#   colnames(bs_no_fish_caught)[i]=names(fall3)[i]
# }
# demersal3 = rbind(demersal3,bs_no_fish_caught)
# 
# demersal3 = demersal3[order(demersal3$sta_lme, fall3$yr),]

#~~~~~~~ UCD "no_fish_caught"
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
#Suisun = Suisun[-which(as.character(Suisun$season_yr)=="2018_1"),]

ucdOTR = subset(Suisun, Suisun$Method=="Otter trawl")
ucd_no_fish = ucdOTR[which(ucdOTR$Length_NA_flag=="No fish caught" & ucdOTR$Station %in% ucd_sta_above_thresh),]
ucd_no_fish_agg=aggregate(Temp_surf ~ Source + season_yr + season_year + season_number + Station, 
                          data=ucd_no_fish, FUN = length)
test=ucd_no_fish_agg[which(ucd_no_fish_agg$Station=="BY1"),]
test[order(test$season_yr),]
ucdOTR[which(ucdOTR$Station=="BY1" & ucdOTR$season_year==1992),]

ucd_no_fish_agg = ucd_no_fish_agg[which(ucd_no_fish_agg$season_year>1984 & ucd_no_fish_agg$season_year<2018),]
ucd_no_fish_agg$lme = rep("ucd", times=dim(ucd_no_fish_agg)[1])
ucd_no_fish_agg$sta_lme = paste(ucd_no_fish_agg$lme, ucd_no_fish_agg$Station, sep="_")
ucd_no_fish_agg$sta_lme_yr = paste(ucd_no_fish_agg$sta_lme, ucd_no_fish_agg$season_year, sep="_")

open_fish_caught = paste(demersal3$sta_lme, demersal3$season_year, sep="_")

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#NOTE: all years have data for UCD
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ucd_no_fish_agg = ucd_no_fish_agg[-which(ucd_no_fish_agg$sta_lme_yr %in% unique(open_fish_caught)),]

##################################################################################
#~~   Convert data into an array for PTA analysis
##################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Standardize by z-scoring per taxa per study-gear

# sea_lme = unique(demersal3$lme)
# 
# sub2 = demersal3[,-c(1:9)]
# sub2 = (sub2)^(1/3)
# sub_scale=apply(sub2, MARGIN = 2,
#                 FUN = function(x){
#                   (x-mean(x))/sd(x)
#                 }
# )
# sub_scale[is.na(sub_scale)]=0
# demersal_scaled=cbind(demersal3[,c(1:9)], sub_scale)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Standardize by z-scoring per taxa per study-gear
# iqr_score = function(x){
#   (x-median(x))/1.5*(IQR(x))
# }

demersal_scaled = demersal3
sea_lme = unique(demersal3$lme)

for(i in 1:length(sea_lme)){
  sub = demersal3[which(demersal3$lme==sea_lme[i]),]
  sub2 = sub[,-c(1:9)]
  sub2 = (sub2)^(1/3)
  sub_scale=apply(sub2, MARGIN = 2,
                  FUN = function(x){
                    (x-mean(x))/sd(x)
                  }
  )
  sub_scale[is.na(sub_scale)]=0
  demersal_scaled[which(demersal3$lme==sea_lme[i]),-c(1:9)]=sub_scale
}

names(demersal_scaled)
##############
par(mfrow=c(2,3), mar=c(4.5, 2, 1,1)+0.1)
hist(demersal3[,"Yellowfin.Goby"][which(demersal3[,"Yellowfin.Goby"]>0)], breaks=30)
hist((demersal3[,"Yellowfin.Goby"][which(demersal3[,"Yellowfin.Goby"]>0)])^(1/3), breaks=30)
hist(demersal_scaled[,"Yellowfin.Goby"][which(demersal3[,"Yellowfin.Goby"]>0)], breaks=30)
hist(demersal3[,"American.shad"][which(demersal3[,"American.shad"]>0)], breaks=30)
hist((demersal3[,"American.shad"][which(demersal3[,"American.shad"]>0)])^(1/3), breaks=30)
hist(demersal_scaled[,"American.shad"][which(demersal3[,"American.shad"]>0)], breaks=30)

# pdf(file = "Count_v_CPUE.pdf", width=10, height=7, paper = "special")
# quartz(height=5, width=8)
par(mfrow=c(2,3), mar=c(4.5, 5, 1,1)+0.1, oma=c(0,0,0,0))
boxplot((count_demersal[,"American.shad"]^(1/3))[which(count_demersal[,"American.shad"]>0)] ~
          count_demersal[,"lme"][which(count_demersal[,"American.shad"]>0)],
        las=1, ylab="cube-root(raw count)", xlab="American shad",
        pch=20,col=gray(0.2,0.2), outcol=gray(0.2,0.2))
boxplot((demersal3[,"American.shad"]^(1/3))[which(demersal3[,"American.shad"]>0)] ~
          demersal3[,"lme"][which(demersal3[,"American.shad"]>0)],
        las=1, ylab="cube-root(raw CPUE)", xlab="American shad",
        pch=20,col=gray(0.2,0.2),outcol=gray(0.2,0.2))
boxplot(demersal_scaled[,"American.shad"][which(demersal3[,"American.shad"]>0)] ~ 
          demersal_scaled[,"lme"][which(demersal3[,"American.shad"]>0)],
        las=1, ylab="z-scored cube-root CPUE", xlab="American shad",
        pch=20,col=gray(0.2,0.2), outcol=gray(0.2,0.2))
boxplot((count_demersal[,"Yellowfin.Goby"]^(1/3))[which(count_demersal[,"Yellowfin.Goby"]>0)] ~
          count_demersal[,"lme"][which(count_demersal[,"Yellowfin.Goby"]>0)],
        las=1, ylab="cube-root(raw count)", xlab="Yellowfin Goby",
        pch=20,col=gray(0.2,0.2), outcol=gray(0.2,0.2))
boxplot((demersal3[,"Yellowfin.Goby"]^(1/3))[which(demersal3[,"Yellowfin.Goby"]>0)] ~
          demersal3[,"lme"][which(demersal3[,"Yellowfin.Goby"]>0)],
        las=1, ylab="cube-root(raw CPUE)", xlab="Yellowfin Goby",
        pch=20,col=gray(0.2,0.2), outcol=gray(0.2,0.2))
boxplot(demersal_scaled[,"Yellowfin.Goby"][which(demersal3[,"Yellowfin.Goby"]>0)] ~
          demersal_scaled[,"lme"][which(demersal3[,"Yellowfin.Goby"]>0)],
        las=1, ylab="z-scored cube-root CPUE", xlab="Yellowfin Goby",
        pch=20,col=gray(0.2,0.2), outcol=gray(0.2,0.2))
# dev.off()

###################



par(mfrow=c(1,2), mar=c(8, 4.5, 1,1)+0.1)
boxplot(demersal3[,-c(1:9)], las=2, cex.axis=0.65,
        pch=20, outcol=gray(0.1,0.1), yaxt='n')
axis(2, las=1, cex=0.8)
boxplot(demersal_scaled[,-c(1:9)], las=2,cex.axis=0.65,
        pch=20, outcol=gray(0.1, 0.1), yaxt='n')
axis(2, las=1, cex=0.8)

par(mfrow=c(1,2), mar=c(8, 4.5, 1,1)+0.1)
boxplot(demersal3[,-c(1:9)], las=2,
        ylim=c(quantile(unlist(demersal3[,-c(1:9)]), probs = c(0,0.975))),
        ylab = "min to 97.5%-tile", cex.axis=0.65,
        pch=20, outcol=gray(0.1,0.1), yaxt='n')
axis(2, las=1, cex=0.8)
boxplot(demersal_scaled[,-c(1:9)], las=2,
        ylim=c(quantile(unlist(demersal_scaled[,-c(1:9)]), probs = c(0,0.975))),
        ylab = "min to 97.5%-tile", cex.axis=0.65,
        pch=20, outcol=gray(0.1,0.1), yaxt='n')
axis(2, las=1, cex=0.8)

raw_xlims=range(demersal3[,"American.shad"])
cuberoot_xlims=range((demersal3[,"American.shad"])^(1/3))
scaled_xlims=range(demersal_scaled[,"American.shad"])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~    Kernal Density visualizations and overlap calculation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~ RAW
lower_raw <- min(c(demersal3[,"American.shad"], demersal3[,"Yellowfin.Goby"]))
upper_raw <- max(c(demersal3[,"American.shad"], demersal3[,"Yellowfin.Goby"]))
dNAC_raw <- density(demersal3[,"American.shad"], from=lower_raw, to=upper_raw,
                    adjust=2, cut=TRUE)
dYFG_raw <- density(demersal3[,"Yellowfin.Goby"], from=lower_raw, to=upper_raw,
                    adjust=2, cut=TRUE)
df_raw <- data.frame(x=dNAC_raw$x, a=dNAC_raw$y, b=dYFG_raw$y)
# calculate intersection densities
df_raw$w <- pmin(df_raw$a, df_raw$b)
# integrate areas under curves
total_raw <- integrate.xy(df_raw$x, df_raw$a) + integrate.xy(df_raw$x, df_raw$b)
intersection_raw <- integrate.xy(df_raw$x, df_raw$w)
# compute overlap coefficient
overlap_raw <- 2 * intersection_raw / total_raw

#~~ CUBE-ROOT
lower_cubeRoot <- min(c((demersal3[,"American.shad"])^(1/3), (demersal3[,"Yellowfin.Goby"])^(1/3)))
upper_cubeRoot <- max(c((demersal3[,"American.shad"])^(1/3), (demersal3[,"Yellowfin.Goby"])^(1/3)))
dNAC_cubeRoot <- density((demersal3[,"American.shad"])^(1/3), from=lower_cubeRoot, to=upper_cubeRoot,
                         adjust=2, cut=TRUE)
dYFG_cubeRoot <- density((demersal3[,"Yellowfin.Goby"])^(1/3), from=lower_cubeRoot, to=upper_cubeRoot,
                         adjust=2, cut=TRUE)
df_cubeRoot <- data.frame(x=dNAC_cubeRoot$x, a=dNAC_cubeRoot$y, b=dYFG_cubeRoot$y)
# calculate intersection densities
df_cubeRoot$w <- pmin(df_cubeRoot$a, df_cubeRoot$b)
# integrate areas under curves
total_cubeRoot <- integrate.xy(df_cubeRoot$x, df_cubeRoot$a) + integrate.xy(df_cubeRoot$x, df_cubeRoot$b)
intersection_cubeRoot <- integrate.xy(df_cubeRoot$x, df_cubeRoot$w)
# compute overlap coefficient
overlap_cubeRoot <- 2 * intersection_cubeRoot / total_cubeRoot

#~~ SCALED
lower_scaled <- min(c(demersal_scaled[,"American.shad"], demersal_scaled[,"Yellowfin.Goby"]))
upper_scaled <- max(c(demersal_scaled[,"American.shad"], demersal_scaled[,"Yellowfin.Goby"]))
dNAC_scaled <- density(demersal_scaled[,"American.shad"], from=lower_scaled, to=upper_scaled,
                       adjust=1, cut=TRUE)
dYFG_scaled <- density(demersal_scaled[,"Yellowfin.Goby"], from=lower_scaled, to=upper_scaled,
                       adjust=1, cut=TRUE)
df_scaled <- data.frame(x=dNAC_scaled$x, a=dNAC_scaled$y, b=dYFG_scaled$y)
# calculate intersection densities
df_scaled$w <- pmin(df_scaled$a, df_scaled$b)
# integrate areas under curves
total_scaled <- integrate.xy(df_scaled$x, df_scaled$a) + integrate.xy(df_scaled$x, df_scaled$b)
intersection_scaled <- integrate.xy(df_scaled$x, df_scaled$w)
# compute overlap coefficient
overlap_scaled <- 2 * intersection_scaled / total_scaled



#~~~~~~~~~~~~
#~ Plot

# pdf(file = "Standardization_density_plot_demersal.pdf", width=10, height=4, paper = "special")
# quartz(height=4, width=10)
par(mfrow=c(1,3), mar=c(4, 4, 2.25,1.5)+0.1, oma=c(0.5,2.5,0.5,0.5))
plot(range(dNAC_raw$x, dYFG_raw$x), range(dNAC_raw$y, dYFG_raw$y), las=1, type='n',
     ylab="", xlab="")
polygon(x = c(dYFG_raw$x[1], dYFG_raw$x[512], rev(dYFG_raw$x), dYFG_raw$x[1]), 
        y = c(0, 0, rev(dYFG_raw$y), dYFG_raw$y[1]), 
        col=gray(0.2,0.6), lwd=0.1)
polygon(x = c(dNAC_raw$x[1], dNAC_raw$x[512], rev(dNAC_raw$x), dNAC_raw$x[1]), 
        y = c(0, 0, rev(dNAC_raw$y), dNAC_raw$y[1]), 
        col=rgb(0.1,0.5,0.85, maxColorValue = 1, alpha=0.5), lwd=0.1)
mtext(text = "Density (a metric of relative probability)", side = 2,line = 4, cex=0.85)
mtext(text = "CPUE", side = 1,line = 2.5, cex=0.85)
text(((par('usr')[2]-par('usr')[1])*0.5)+par('usr')[1],
     ((par('usr')[4]-par('usr')[3])*0.95)+par('usr')[3],
     cex=1, pos=1,
     paste("American shad range =",
           round(min(demersal3[,"American.shad"]), digits=2), "-",
           round(max(demersal3[,"American.shad"]), digits=2), sep=" "))
text(((par('usr')[2]-par('usr')[1])*0.5)+par('usr')[1],
     ((par('usr')[4]-par('usr')[3])*0.9)+par('usr')[3],
     cex=1, pos=1,
     paste("Yellowfin Goby range =",
           round(min(demersal3[,"Yellowfin.Goby"]), digits=2), "-",
           round(max(demersal3[,"Yellowfin.Goby"]), digits=2), sep=" "))
text(((par('usr')[2]-par('usr')[1])*0.5)+par('usr')[1],
     ((par('usr')[4]-par('usr')[3])*0.85)+par('usr')[3],
     cex=1, pos=1,
     paste("Distribution Overlap = ~",
           round(overlap_raw, digits=2)*100, "%", sep=""))

plot(range(dNAC_cubeRoot$x, dYFG_cubeRoot$x), range(dNAC_cubeRoot$y, dYFG_cubeRoot$y), las=1, type='n',
     ylab="", xlab="")
polygon(x = c(dYFG_cubeRoot$x[1], dYFG_cubeRoot$x[512], rev(dYFG_cubeRoot$x), dYFG_cubeRoot$x[1]), 
        y = c(0, 0, rev(dYFG_cubeRoot$y), dYFG_cubeRoot$y[1]), 
        col=gray(0.2,0.6), lwd=0.1)
polygon(x = c(dNAC_cubeRoot$x[1], dNAC_cubeRoot$x[512], rev(dNAC_cubeRoot$x), dNAC_cubeRoot$x[1]), 
        y = c(0, 0, rev(dNAC_cubeRoot$y), dNAC_cubeRoot$y[1]), 
        col=rgb(0.1,0.5,0.85, maxColorValue = 1, alpha=0.5), lwd=0.1)
mtext(text = bquote(sqrt(CPUE, 3)), side = 1,line = 2.5, cex=0.85)
text(((par('usr')[2]-par('usr')[1])*0.5)+par('usr')[1],
     ((par('usr')[4]-par('usr')[3])*0.95)+par('usr')[3],
     cex=1, pos=1,
     paste("American shad range =",
           round(min((demersal3[,"American.shad"])^(1/3)), digits=2), "-",
           round(max((demersal3[,"American.shad"])^(1/3)), digits=2), sep=" "))
text(((par('usr')[2]-par('usr')[1])*0.5)+par('usr')[1],
     ((par('usr')[4]-par('usr')[3])*0.9)+par('usr')[3],
     cex=1, pos=1,
     paste("Yellowfin Goby range =",
           round(min((demersal3[,"Yellowfin.Goby"])^(1/3)), digits=2), "-",
           round(max((demersal3[,"Yellowfin.Goby"])^(1/3)), digits=2), sep=" "))
text(((par('usr')[2]-par('usr')[1])*0.5)+par('usr')[1],
     ((par('usr')[4]-par('usr')[3])*0.85)+par('usr')[3],
     cex=1, pos=1,
     paste("Distribution Overlap = ~",
           round(overlap_cubeRoot, digits=2)*100, "%", sep=""))


plot(range(dNAC_scaled$x, dYFG_scaled$x), range(dNAC_scaled$y, dYFG_scaled$y), las=1, type='n',
     ylab="", xlab="")
polygon(x = c(dYFG_scaled$x[1], dYFG_scaled$x[512], rev(dYFG_scaled$x), dYFG_scaled$x[1]), 
        y = c(0, 0, rev(dYFG_scaled$y), dYFG_scaled$y[1]), 
        col=gray(0.2,0.6), lwd=0.1)
polygon(x = c(dNAC_scaled$x[1], dNAC_scaled$x[512], rev(dNAC_scaled$x), dNAC_scaled$x[1]), 
        y = c(0, 0, rev(dNAC_scaled$y), dNAC_scaled$y[1]), 
        col=rgb(0.1,0.5,0.85, maxColorValue = 1, alpha=0.5), lwd=0.1)
mtext(text = bquote(z-scored~sqrt(CPUE, 3)), side = 1,line = 2.5, cex=0.85)
text(((par('usr')[2]-par('usr')[1])*0.5)+par('usr')[1],
     ((par('usr')[4]-par('usr')[3])*0.95)+par('usr')[3],
     cex=1, pos=1,
     paste("American shad range =",
           round(min(demersal_scaled[,"American.shad"]), digits=2), "-",
           round(max(demersal_scaled[,"American.shad"]), digits=2), sep=" "))
text(((par('usr')[2]-par('usr')[1])*0.5)+par('usr')[1],
     ((par('usr')[4]-par('usr')[3])*0.9)+par('usr')[3],
     cex=1, pos=1,
     paste("Yellowfin Goby range =",
           round(min(demersal_scaled[,"Yellowfin.Goby"]), digits=2), "-",
           round(max(demersal_scaled[,"Yellowfin.Goby"]), digits=2), sep=" "))
text(((par('usr')[2]-par('usr')[1])*0.5)+par('usr')[1],
     ((par('usr')[4]-par('usr')[3])*0.85)+par('usr')[3],
     cex=1, pos=1,
     paste("Distribution Overlap = ~",
           round(overlap_scaled, digits=2)*100, "%", sep=""))

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top", legend = c("American shad", "Yellowfin Goby") , pch=22,
       pt.bg = c(rgb(0.1,0.5,0.85, maxColorValue = 1, alpha=0.5), gray(0.2,0.6))
       , xpd=TRUE, inset=c(0, 0), bty='n', horiz=TRUE, pt.cex=3, cex=1.25) 
# dev.off()

range(demersal3[,-c(1:9)])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ Convert data into array format
#~  [season_timeseries_number, taxa, station]

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

# [season_timeseries_number, taxa, station]
# demersal_scaled_arr[,,1]

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

red_demersal_array[1,1,]; length(red_demersal_array[1,1,])
red_demersal_array[1,,1]; length(red_demersal_array[1,,1])
red_demersal_array[,1,1]; length(red_demersal_array[,1,1])

####################################################################
#~  Principal tensor analysis 'a la' Frelat etal 2017 PlosOne
##################################################################
rm(list=c("pta"))
pta<-PTA3(red_demersal_array, nbPT = 3, nbPT2 = 3, minpct = 0.1)
summary.PTAk(pta,testvar = 0)
par(mfrow=c(1,1), mar=c(4.5,4.5,1.5,1.5+0.1), oma=c(0,0,0,0))
plot(pta, scree=TRUE)
keep = c(1,3,11,4)

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

demersal_summary = aggregate(season_year ~ season_number + season_count + season_year, FUN=unique, data=demersal3)

demersal_summary$yr=substr(as.character(demersal_summary$season_year), start=3, stop=4)
demersal_summary = subset(demersal_summary, demersal_summary$season_count %in% as.numeric(pta[[1]]$n))

col_season = c("darkblue", "darkgreen", "darkred", "tan4")

col_season_start = c("#00D3EB", "#00FF33", "#CDA8A8", "#FFBD04")
col_season_end = c("#001D6C", "#095500", "#640000", "#372900")

#~~~~~~~~
# Year/season plots

unique_season_plot = 1:4
plot(pta[[1]]$v[1,] , pta[[1]]$v[11,], type='n', las=1,
     xlab  = "Primary Tensor 1 (6.64% Global)",
     ylab  = "Primary Tensor 3 (2.41% Global)")
abline(h=0, col=gray(0.2,0.2))
abline(v=0, col=gray(0.2,0.2))
for(s in 1:length(unique_season_plot)){
  sub_dat = pta[[1]]$v[,which(demersal_summary[,1]==unique_season_plot[s])]
  color_range <- colorRampPalette( c(col_season_start[s], col_season_end[s]) ) 
  color_spectra <- color_range(dim(sub_dat)[1])
  
  points(sub_dat[1,], sub_dat[11,], col=color_spectra, pch=20)
}


plot(pta[[1]]$v[1,] , pta[[1]]$v[11,], type='n', las=1,
     xlab  = "Primary Tensor 1 (6.64% Global)",
     ylab  = "Primary Tensor 3 (2.41% Global)")
abline(h=0, col=gray(0.2,0.2))
abline(v=0, col=gray(0.2,0.2))
text(pta[[1]]$v[1,] , pta[[1]]$v[11,],
     col=col_season[demersal_summary[,1]],labels =demersal_summary[,4],
     adj=c(0.5,0.5), cex=0.5, font=2)

plot(pta[[1]]$v[1,] , pta[[1]]$v[11,], type='n', las=1)
unique_season_plot=sort(unique(demersal_summary[,1]))
for(s in 1:length(unique_season_plot)){
  sub_dat = pta[[1]]$v[,which(demersal_summary[,1]==unique_season_plot[s])]
  for(t in 1:c(dim(sub_dat)[2]-1)){
    arrows(x0 = sub_dat[1,t], y0 = sub_dat[11,t],
           x1 = sub_dat[1,t+1], y1 = sub_dat[11,t+1],
           length = 0.1, col = col_season[s])
  }
}

plot(pta[[2]]$v[1,] , pta[[2]]$v[11,], type='n', las=1,
     xlab  = "Primary Tensor 1 (6.64% Global)",
     ylab  = "Primary Tensor 3 (2.41% Global)")
abline(h=0, col=gray(0.2,0.2))
abline(v=0, col=gray(0.2,0.2))
text(pta[[2]]$v[1,] , pta[[2]]$v[11,],
     labels =pta[[2]]$n,
     adj=c(0.5,0.5), cex=0.5,font=2)


plot(pta[[3]]$v[1,] , pta[[3]]$v[11,], type='n', las=1,
     xlab  = "Primary Tensor 1 (6.64% Global)",
     ylab  = "Primary Tensor 3 (2.41% Global)")
abline(h=0, col=gray(0.2,0.2))
abline(v=0, col=gray(0.2,0.2))
text(pta[[3]]$v[1,] , pta[[3]]$v[11,],
     labels =pta[[3]]$n,
     adj=c(0.5,0.5), cex=0.5,font=2)



unique_season_plot=sort(unique(demersal_summary[,1]))
par(mfrow=c(2,2))
for(s in 1:length(unique_season_plot)){
  index=which(demersal_summary[,1]==unique_season_plot[s])
  sub_dat = pta[[1]]$v[,index]
  plot(x = as.numeric(as.character(demersal_summary[index,3])), y=colSums(sub_dat), type='l', las=1)
  abline(h=0, col=gray(0.2,0.2))
}


######################################################################
#       EVALUATE LUSTER FIT AND COMPARE METHODS
######################################################################

#Create the matrix with the projection of species on the 4 PT
coo<-t(pta[[2]]$v[c(keep),])

rownames(coo) = pta[[2]]$n
PT_associations = c("PT1: Spatial (region)", "PT2: Spatial (sub-region)", 
                    "PT3: Temporal (season & year)",
                    "PT4: Spatial (sub-region)")
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


setwd("/Users/jeremewgaeta/Files/DFW/Review_2019_2020/Cluster/integrated_data_analyses/demersal_all_seasons/")


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

plot(kt$k, kt$r, 
     main="", xlab="Number of Clusters",
     ylab=bquote("Pearson's"~italic(rho)), type='h', lend=1,
     las=1, lwd=2, ylim=c(0, max(kt$r)*1.05), yaxs='i')
lines(y = c(0,max(kt$r)), x=c(k.best.mantel, k.best.mantel),
      col="red3", lwd=4, lend=1)
axis(1, k.best.mantel,  col.axis="red3", col.ticks = "red3",
     font=2, labels=FALSE)
axis(1, k.best.mantel,  col.axis="red3", col.ticks = "red3", font=2,line = -0.7)
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

optim_k = 11
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

setwd("/Users/jeremewgaeta/Files/DFW/Review_2019_2020/Cluster/integrated_data_analyses/demersal_all_seasons/")
# pdf(file = "Demersal_all_season_Silhouette.pdf", width = 7, height=9, paper = "special")
# quartz(width=7, height=9)
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
par(new=TRUE, mar=c(4.5,10,2,2)+0.1)
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

# pdf(file = "Demersal_PTA_rect_dendo.pdf", width = 12, height=7, paper = "special")
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
taxa_cluster_df
taxa_cluster_names = c("Sturgeon+", 
                       "Silverside+",
                       "Wht. Catfish",
                       "Sanddab+",
                       "Sole+", #5
                       "Splittail+",
                       "C. Gobies",
                       "PSS & Longfin",
                       "Halibut & Tonguefish",
                       "Shokihaze",#10
                       "Shimofuri")

# pdf(file = "PTA_demersals_ordination.pdf", width = 12, height=4.5, paper = "special")
# windows(width = 7, height=9)
# quartz(width = 12, height=4.5)
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
# dev.off()


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

#~~~~~~~~~~~~~~~~~~~~
par(mfrow=c(1,2), mar=c(4,4,2,2)+0.1, oma=c(0,0,0,0))
s.class(coo, fac = clust_3D, xax=1, yax=2,
        col=hcl.colors(n=nclust+1, palette = "Viridis"),
        clabel = 1, label=taxa_cluster_names,
        addaxes = TRUE, cgrid = 0)
text(((par('usr')[2]-par('usr')[1])*0.05)+par('usr')[1], 0,
     labkeep[3], srt=90, xpd=NA, cex=1)
text(0,((par('usr')[4]-par('usr')[3])*0.05)+par('usr')[3],
     labkeep[1], xpd=NA, cex=1)

s.label(coo, xax=1, yax=2,clabel = 0.5)

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
image(t(pta4_mat), col=hcl.colors(9, "BrBg", rev = FALSE), axes=F)
axis(side = 1, at=seq(0,1, length.out=length(pta[[1]]$n)), labels=pta[[1]]$n)
axis(side = 2, at=seq(0,1, length.out=length(pta[[3]]$n)),
     labels=c(pta[[3]]$n), las=1, cex.axis=0.85)
box(which="plot")
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

PT_axis = paste0("PT", 1:4, ": ", round((100 * (pta[[3]]$d[keep])^2)/pta[[3]]$ssX[1],1), "%")

# pdf(file = "Demersal_PTA_allseason_heatgrids.pdf", width = 12, height=5, paper = "special")
# windows(width = 7, height=9)
# quartz(width = 12, height=5)
par(mfrow=c(1,4), mar=c(4,4.5,2.5,0.75)+0.1, oma=c(0,0.5,0,7))
image(t(pta1_mat), col=hcl.colors(150, "BrBg", rev = FALSE), axes=F)
axis(side = 1, at=seq(0,1, length.out=length(pta[[1]]$n)), labels=pta[[1]]$n)
axis(side = 2, at=seq(0,1, length.out=length(pta[[3]]$n)), labels=c(pta[[3]]$n),
     las=1, cex.axis=0.7)
mtext(text = PT_axis[1], font=2, line = 0.5)
mtext(text = "Time Step", side = 1, line = 2.25, cex=0.85)
box(which="plot")

image(t(pta2_mat), col=hcl.colors(150, "BrBg", rev = FALSE), axes=F)
axis(side = 1, at=seq(0,1, length.out=length(pta[[1]]$n)), labels=pta[[1]]$n)
axis(side = 2, at=seq(0,1, length.out=length(pta[[3]]$n)), labels=c(pta[[3]]$n),
     las=1, cex.axis=0.7)
box(which="plot")
mtext(text = PT_axis[2], font=2, line = 0.5)
mtext(text = "Time Step", side = 1, line = 2.25, cex=0.85)
box(which="plot")

image(t(pta3_mat), col=hcl.colors(150, "BrBg", rev = FALSE), axes=F)
axis(side = 1, at=seq(0,1, length.out=length(pta[[1]]$n)), labels=pta[[1]]$n)
axis(side = 2, at=seq(0,1, length.out=length(pta[[3]]$n)), labels=c(pta[[3]]$n),
     las=1, cex.axis=0.7)
box(which="plot")
mtext(text = PT_axis[3], font=2, line = 0.5)
mtext(text = "Time Step", side = 1, line = 2.25, cex=0.85)
box(which="plot")

image(t(pta4_mat), col=hcl.colors(150, "BrBg", rev = FALSE), axes=F)
axis(side = 1, at=seq(0,1, length.out=length(pta[[1]]$n)), labels=pta[[1]]$n)
axis(side = 2, at=seq(0,1, length.out=length(pta[[3]]$n)), labels=c(pta[[3]]$n),
     las=1, cex.axis=0.7)
box(which="plot")
mtext(text = PT_axis[4], font=2, line = 0.5)
mtext(text = "Time Step", side = 1, line = 2.25, cex=0.85)
box(which="plot")


par(fig=c(0,1,0,1), mar=c(4.5,4.5,4.5,1.5)+0.1, oma=c(0,81,0,0), new=TRUE)
color.bar(hcl.colors(150, "BrBg", rev = FALSE), min=0, max=1)
mtext(text = bquote(atop(bold(Anomaly),bold(Gradient))),
      side = 3,line = 0.25, adj=0.5, las=0)
# dev.off()




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hex_col = c("#ABD9E9", "#D7191C", "#FDAE61", "#2C7BB6")

# pdf(file = "Demersal_PT3_temporal.pdf", width = 12, height=6, paper = "special")
# windows(width = 7, height=9)
# quartz(width = 12, height=6)
par(mfrow=c(1,2), mar=c(3.75,5,3.25,1.5)+0.1)
image(t(pta3_mat), col=hcl.colors(150, "BrBg", rev = FALSE), axes=F)
axis(side = 1, at=seq(0,1, length.out=length(pta[[1]]$n)), labels=pta[[1]]$n)
axis(side = 2, at=seq(0,1, length.out=length(pta[[3]]$n)), labels=c(pta[[3]]$n),
     las=1, cex.axis=0.6)
box(which="plot")
mtext(text = PT_axis[3], font=2, line = 0.4)
mtext(text = "Time Step", side = 1, line = 2.25, cex=0.85)
box(which="plot")

plot(pta[[1]]$v[keep[3],] ~ as.numeric(pta[[1]]$n),
     pch=19, cex=2, type='l', las=1, xaxt='n',
     ylab="Relative Anomaly",
     xlab="")
axis(side = 1, at=as.numeric(pta[[1]]$n), las=1)
axis(side = 3, at = seq(from=1, to=152, by=4), labels=1980:2017)
abline(h=0, col=gray(0.5,0.3))
points(pta[[1]]$v[keep[3],] ~ as.numeric(pta[[1]]$n), 
       col=hex_col[demersal_summary[,1]], pch=19)
legend("topleft", legend = c("Spring", "Summer", "Fall", "Winter"),
       pch=19, col = hex_col, bty='n', inset=0.025)
mtext(text = "Time Step", side = 1, line = 2.25, cex=0.85)
mtext(text = "Year", side = 3, line = 2.25, cex=0.85)
# dev.off()

##################################################################
#~~~~~~~~~~~~~~~~~~
# Spatial cluster
#~~~~~~~~~~~~~~~~~~
##################################################################
#Create the matrix with the projection of species on the 4 PT
#Create the matrix with the projection of species on the 4 PT

spat_coo<-t(pta[[3]]$v[c(keep),])
rownames(spat_coo) = pta[[3]]$n
spat_labkeep <- paste0(pta[[3]]$vsnam[keep], " - ", round((100 * (pta[[3]]$d[keep])^2)/pta[[3]]$ssX[1],1), "%")
#1. Compute the distance between species
dist_spatial=dist(spat_coo, method = "euclidean")
#2. Build a tree with Ward linkage
spatial_den=hclust(dist_spatial, method='average')
#3. Plot the dendogram
spatial_dendro = as.dendrogram(spatial_den)
spatial_clusters <- cutree(spatial_dendro, k=10)[order.dendrogram(spatial_dendro)]
spatial_col_id = c(hcl.colors(n=9, palette = "Viridis")[c(1,3,7)],
                   hcl.colors(n=9, palette = "Heat")[c(5)],
                   hcl.colors(n=10, palette = "RdBu")[c(10,2)],
                   hcl.colors(n=7, palette = "Dark 2")[c(3,5)],
                   hcl.colors(n=10, palette = "BrBG")[c(2,10)],
                   hcl.colors(n=10, palette = "PiYG")[c(2,10)])


spatial_dendro = spatial_dendro %>%
  branches_attr_by_clusters(spatial_clusters, values = spatial_col_id) %>%
  set("branches_lwd", 2) %>%
  set("labels_cex", 0.8)
par(mfrow=c(1,1), mar=c(5,1.5,1.5,1.5)+0.1)
plot(spatial_dendro, xaxt='n', axes=FALSE)
rect.hclust(spatial_den, k=7, border=gray(0.2,0.2))

####### Identify optimal clusters
spat_asw = numeric(nrow(spat_coo))

# write values
for(k in 2:(nrow(spat_coo)-1)) {
  sil = silhouette(cutree(spatial_den, k=k), dist_spatial)
  spat_asw[k] = summary(sil)$avg.width
}

#best (larggest Silhouette width)
spat_k.best.silhouette = which.max(spat_asw)

spat_kt = data.frame(k=1:nrow(spat_coo), r=0) 

for(i in 2:nrow((spat_coo)-1)) {
  gr = cutree(spatial_den, i)
  distgr = grpdist(gr)
  mt = cor(dist_spatial, distgr, method = "pearson")
  spat_kt[i,2] = mt
}

spat_kt = na.omit(spat_kt)
spat_k.best.mantel = which.max(spat_kt$r)

# pdf(file = "OpenWater_demersalson_spatial_mantel_eval.pdf", width = 11, height=5, paper = "special")
par(mfrow=c(1,2), mar=c(4.5,4.5,1.5,1.5)+0.1)
plot(1:nrow(spat_coo), spat_asw, 
     main="Silhouette - Complete", xlab="Number of Clusters",
     ylab="Average Silhouette Width", type='h', lend=1,
     las=1, lwd=2, ylim=c(0, max(spat_asw)*1.05), yaxs='i')
lines(y = c(0,max(spat_asw)), x=c(spat_k.best.silhouette, spat_k.best.silhouette),
      col="red3", lwd=4, lend=1)
axis(1, spat_k.best.silhouette,  col.axis="red3", col.ticks = "red3", font=2)
axis(1, spat_k.best.silhouette, "(optimum k)", col="red3",
     font=2, col.axis="red3", line = 0.75, tick = FALSE)

plot(spat_kt$k, spat_kt$r, 
     main="Mantel - Complete", xlab="Number of Clusters",
     ylab=bquote("Pearson's"~italic(rho)), type='h', lend=1,
     las=1, lwd=2, ylim=c(0, max(spat_kt$r)*1.05), yaxs='i')
lines(y = c(0,max(spat_kt$r)), x=c(spat_k.best.mantel, spat_k.best.mantel),
      col="red3", lwd=4, lend=1)
axis(1, spat_k.best.mantel,  col.axis="red3", col.ticks = "red3", font=2)
axis(1, spat_k.best.mantel, "(optimum k)", col="red3",
     font=2, col.axis="red3", line = 0.75, tick = FALSE)
# dev.off()



spat_nclust = 9

spat_cutg = cutree(spatial_dendro, k = spat_nclust)
spat_sil = silhouette(spat_cutg, dist_spatial)
spat_silo = sortSilhouette(spat_sil)
silo_names = data.frame(code=row.names(spat_coo)[attr(spat_silo, "iOrd")])
row.names(spat_silo)=silo_names$code


# pdf(file = "PTA_spatial_sil_demersal_year_round.pdf", width = 6, height=7, paper = "special")
#quartz(width = 6, height=7)
par(mar=c(5.5, 5, 3,2))
plot(spat_silo, border=1, 
     main="Demersal, year-round",
     cex.names=0.5, nmax.lab = 100, max.strlen = 10)
# dev.off()

spat_clust_3D <- as.factor(cutree(spatial_dendro, k=spat_nclust))
dendro_spat_df = data.frame(site = unlist(partition_leaves(spatial_dendro)[1]),
                            dendr_order=1:(length(unlist(partition_leaves(spatial_dendro)[1]))))
spat_cluster_df=data.frame( site = row.names(spat_coo), cluster = spat_clust_3D)
spat_cluster_df=merge(spat_cluster_df ,dendro_spat_df)
spat_cluster_df=spat_cluster_df[order(spat_cluster_df$dendr_order),]
clust_col_order=spat_cluster_df$cluster[!duplicated(spat_cluster_df$cluster)]

spat_cluster_df = spat_cluster_df[order(spat_cluster_df$cluster),]
spat_dendro_black=as.dendrogram(spatial_den)
spat_dendro_black = spat_dendro_black %>%
  set("branches_lwd", 1.5) %>%
  set("labels_cex", 0.85) 


height_df=data.frame(object=as.character(get_nodes_attr(spat_dendro_black, "label")),
                     height = get_nodes_attr(spat_dendro_black, "height"))
par(mar=c(10,6,0,0))
plot(spat_dendro_black)

object=c()
height=c()
for(i in 2:dim(height_df)[1]){
  if(i==2){
    if(height_df$height[i]==0 & height_df$height[i-1]>0){
      object=c(object, as.character(height_df$object[i]))
      height=c(height, height_df$height[i-1])
    }
  }
  if(i>2){
    if(height_df$height[i]==0 & height_df$height[i-1]==0){
      object=c(object, as.character(height_df$object[i]))
      height=c(height, height_df$height[i-2])
    }
    if(height_df$height[i]==0 & height_df$height[i-1]>0){
      object=c(object, as.character(height_df$object[i]))
      height=c(height, height_df$height[i-1])
    }
  }
}


nearest_node_height = data.frame(station=object, nearest_node_height=height)

write.csv(x = nearest_node_height, file = "dyr_nearest_node_height.csv")



# pdf(file = "Demersal_PTA_rect_spatial_dendo.pdf", width = 12, height=7, paper = "special")
# windows(width = 7, height=9)
# quartz(width = 10, height=7)
par(mar=c(1,3,1,1))
par(mfrow=c(1,1), mar=c(5,1.5,1.5,1.5)+0.1, oma=c(0,0,0,0))
plot(spat_dendro_black, xaxt='n', axes=FALSE, lwd=5, cex.lab=0.5)
rect.dendrogram(spat_dendro_black, k=spat_nclust,
                border=hcl.colors(n=spat_nclust+1,
                                  palette = "Viridis")[clust_col_order],
                lwd=1.5)
# dev.off()

par(mar=c(4.5,4.5,1,1)+0.1, 
    fig = c(0.5,1, 0.55, 1), new = T)  
plot(spat_kt$k, spat_kt$r, 
     main="Mantel - Complete", xlab="Number of Clusters",
     ylab=bquote("Pearson's"~italic(rho)), type='h', lend=1,
     las=1, lwd=2, ylim=c(0, max(spat_kt$r)*1.05), yaxs='i')
lines(y = c(0,max(spat_kt$r)), x=c(spat_k.best.mantel, spat_k.best.mantel),
      col="red3", lwd=4, lend=1)
axis(1, spat_k.best.mantel,  col.axis="red3", col.ticks = "red3", font=2)
axis(1, spat_k.best.mantel, "(optimum k)", col="red3",
     font=2, col.axis="red3", line = 0.75, tick = FALSE)
# dev.off()



