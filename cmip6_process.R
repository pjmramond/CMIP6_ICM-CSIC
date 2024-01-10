############################################################################################
# Generate and unite global predictions for Oceanic variables from CMIP6
############################################################################################

# 1/ Data selection was made there: https://aims2.llnl.gov/search

# - Activity ID: ScenarioMIP (model predictions for different ICPP scenarios of climate change)
# - Resolution: 100 km
# - Variant Label: r1i1p1f1 (this is based on an email from Pablo Ortega)
# - Frequency: yr (for most variable) and mon (for temperature/thetao)
# - Variable ID: https://redmine.dkrz.de/attachments/download/5773/AR6%20WG1%20priority%20variables.xlsx
# - Result Type: Originals Only (avoid duplicates)
# => display 20 files per page (for more files the wget script won't be generated), add to cart, download wget script, copy the http links to new text file, do this for all files for each variable
# => I then subset the http links based on the presence of year 2100 in the filename (this is consistent across files)

# 2) On ada, I then run the download with command: for i in `awk '{print $2}' thetao_2100.txt`; do wget $i; done
# 3) the processing of the data is then performed with this script.

## Check other variables in TARA
# Info on TARA samples
info<-read.table("/export/lv6/user/pramond/Ocean_Microbiomics_Database/omd_env.csv",sep = ";", header = TRUE)
info$Ocean.region[is.na(info$Ocean.region)]<-"[AO] Arctic Ocean (MRGID:1906)"
info<-info[!duplicated(info$PANGAEA), ]
table(info$sf)
info.s<-info[info$Layer == "SRF",]

# env.
na_env<-sapply(info.s, function(y) sum(length(which(is.na(y)))))
env<-na.omit(info[, names(na_env[na_env<10]) ])
env<-env[!duplicated(env$PANGAEA), ]
rownames(env)<-env$metagID

# enviornmental variables in TARA:
predictor_variables=colnames(env)[14:ncol(env)]
head(env[, predictor_variables])

############################################################################################

# load libraries
library(ncdf4)
library(tidyverse)
library(raster)
library(rasterVis)
library(viridis)
library("RNetCDF")
library("processNC")
library(ggplot2)
library(ggpubr)
library(scales)

#############################################################################################
# Import and format data of all variables
#############################################################################################

## load all nectcdf of all monthly variables  
########################################################

# list all variables and scenarios
scenario<-c("ssp119","ssp126","ssp245","ssp370","ssp434","ssp460","ssp534-over","ssp585") 
var<-grep("[.]|the|so",list.files("/export/lv6/user/pramond/CMIP6/"), invert = TRUE, value=TRUE) # list all variables with monthly frequency

# Loop to import and transform all netcdfs for each variable and save it in a list
list_2100_nc<-list();list_2100_wg<-list();list_2100_df<-list();list_dl<-list() # create empty lists to collect loop results
for (v in var){ # loop for each variable
  # retrieve information on files for each variable
  chl=list.files(paste("/export/lv6/user/pramond/CMIP6/", v, '/', sep = "" ), pattern = "*.nc") # list all netcdfs for this variable
  # Collect information on netcdf files based on their names
  dl=data.frame(files = chl,
                model = stringr::str_split_fixed(chl, "_", n = 7)[,3],
                scenarios = stringr::str_split_fixed(chl, "_", n = 7)[,4],
                year = stringr::str_split_fixed(chl, "_", n = 7)[,7]
  )
  dl$start<-as.numeric(gsub("-2100.nc","", dl$year))
  dl$bands<-2100-dl$start+1 # get number of bands (years in netcdf)
  # filter out models with different grid and coordinate systems
  dl<-dl[dl$model != "CESM2-WACCM", ]
  dl<-dl[grep("NorESM2", dl$model, invert = TRUE), ]
  dl<-dl[grep("cmip5", dl$scenarios, invert = TRUE), ]
  
  # load all netcdf for a variable
  for (i in 1:nrow(dl)){
    # load raw raster
    path<-paste("/export/lv6/user/pramond/CMIP6/",v, "/", dl$files[i], sep = "") # make path
    r<-raster(x = path, band = dl$bands[i], level = 1,  stopIfNotEqualSpaced = FALSE) # read/import it (level 1, correspond to the surface layer)
    # convert to normalize coordinates and convert to WGS84
    rtp<-as.data.frame(rasterToPoints(r)) # convert raster to data.frame
    rtp$x<-scales::rescale(rtp$x, to = c(0, 360)) # rescale model grid
    rtp$y<-scales::rescale(rtp$y, to = c(-90, 90)) # rescale model grid
    colnames(rtp)[3]<-paste(v,dl$model[i], dl$scenarios[i], sep = ".") # change raster variable name
    r2<-rasterFromXYZ(rtp) # convert back to raster with new coordinates
    crs(r2) <- CRS('+init=EPSG:4326') # force WGS84 crs
    # convert to data.frame
    rdf<-as.data.frame(rasterToPoints(r2)) # turn back to data.frame
    # save outputs in lists
    list_2100_nc[[v]][[dl$scenarios[i]]][[dl$model[i]]]<-r # save original rasters
    list_2100_wg[[v]][[dl$scenarios[i]]][[dl$model[i]]]<-r2 # save rasters with modified coordinates
    list_2100_df[[v]][[dl$scenarios[i]]][[dl$model[i]]]<-rdf # save data.frames with standardized coordinates
    #remove(r)
  }
  # save infos on netcdf saved in lists
  list_dl[[v]]<-dl
}  
dl<-bind_rows(list_dl, .id = "variable")

# explore raw model output
var_check="no3";var # select variable
sc="ssp585" # select climate scenario
list_plots<-list() # empty list to collect all plots loaded
for (i in unique(dl[dl$variable %in% var_check & dl$scenarios %in% sc, "model"] ) ){
  list_plots[[i]]<-gplot( list_2100_wg[[var_check]][[sc]][[i]] ) + geom_tile(aes(fill = value))+
    ggtitle(paste(i))+
    scale_fill_viridis(option = "B")
}
do.call("ggarrange", c(list_plots)) # plot the plot-list 

## load all nectcdf of all monthly variables  
########################################################

# Loop to import and transform all netcdfs for each variable and save it in a list
list_mon_nc<-list();list_mon_wg<-list();list_mon_df<-list();list_dlt<-list() # create empty lists to collect loop results
for (v in c("thetao","so")){ # loop for each variable
  # retrieve information on files for each variable
  chl=list.files(paste("/export/lv6/user/pramond/CMIP6/", v, '/', sep = "" ), pattern = "*.nc") # list all netcdfs for this variable
  # Collect information on netcdf files based on their names
  dlt=data.frame(files = chl,
                model = stringr::str_split_fixed(chl, "_", n = 7)[,3],
                scenarios = stringr::str_split_fixed(chl, "_", n = 7)[,4],
                year = stringr::str_split_fixed(chl, "_", n = 7)[,7]
  )
  dlt$start<-as.numeric(gsub("01-210012.nc","", dlt$year))
  dlt<-na.omit(dlt)
  dlt$bands<-2100-dlt$start+1
  dlt<-unique(dlt)
  # filter out models with different grid and coordinate systems
  dlt<-dlt[!dlt$model %in% c("CESM2-WACCM", "BCC-CSM2-MR","CAS-ESM2-0",
                             "CESM2-WACCM", "E3SM-1-1" , "E3SM-1-1-ECA", 
                             "FGOALS-f3-L", "FGOALS-g3", "INM-CM4-8", "INM-CM5-0",
                             "MRI-ESM2-0", "TaiESM1", "MIROC6"), ]
  rownames(dlt)=NULL
  
  # load all netcdf for a variable
  for (i in 1:nrow(dlt)){
    # load raw raster
    path<-paste("/export/lv6/user/pramond/CMIP6/",v,"/",dlt$files[i], sep = "")
    r<-brick(x = path,stopIfNotEqualSpaced = FALSE, level =1)
    rm<-calc(r[[names(r)[(nbands(r)-11):nbands(r)]]], fun = mean)
    # convert to normalize coordinates and convert to WGS84
    rtp<-as.data.frame(rasterToPoints(rm)) # convert raster to data.frame
    rtp$x<-scales::rescale(rtp$x, to = c(0, 360)) # rescale model grid
    rtp$y<-scales::rescale(rtp$y, to = c(-90, 90)) # rescale model grid
    colnames(rtp)[3]<-paste(v,dlt$model[i], dlt$scenarios[i], sep = ".") # change raster variable name
    r2<-rasterFromXYZ(rtp) # convert back to raster with new coordinates
    crs(r2) <- CRS('+init=EPSG:4326') # force WGS84 crs
    # convert to data.frame
    rdf<-as.data.frame(rasterToPoints(r2)) # turn back to data.frame
    # save outputs in lists
    list_mon_nc[[v]][[dlt$scenarios[i]]][[dlt$model[i]]]<-r # save original rasters
    list_mon_wg[[v]][[dlt$scenarios[i]]][[dlt$model[i]]]<-r2 # save rasters with modified coordinates
    list_mon_df[[v]][[dlt$scenarios[i]]][[dlt$model[i]]]<-rdf # save data.frames with standardized coordinates
    print(paste(v, " / ",i, " / ",nrow(dlt), " models",  sep = ""))
  }
  # save infos on netcdf saved in lists
  list_dlt[[v]]<-dlt
}  
dlt<-bind_rows(list_dlt, .id = "variable")

# explore raw model output
var_check="so";var # select variable
sc="ssp585" # select climate scenario
list_plots<-list() # empty list to collect all plots loaded
for (i in unique(dlt[dlt$variable %in% var_check & dlt$scenarios %in% sc, "model"] ) ){
  list_plots[[i]]<-gplot( list_mon_wg[[var_check]][[sc]][[i]] ) + geom_tile(aes(fill = value))+
    ggtitle(paste(i))+
    scale_fill_viridis(option = "B")
}
do.call("ggarrange", c(list_plots)) # plot the plot-list 

## normalize model grids, format and summarize values per scenarios
###################################################################

citation("raster")

# explore available data for each variable scenarios
list_dl$thetao<-dlt[dlt$variable == "thetao",];list_dl$so<-dlt[dlt$variable == "so",]
netcdf<-bind_rows(list_dl, .id = "variable")
table(netcdf$variable,netcdf$scenarios)
list_df<-c(list_2100_df, list_mon_df)
my_merge <- function(df1, df2){merge(df1, df2, by = c("x", "y"), all = TRUE)}

# loop to get mean of all the models for each variable and selected scenarios
# takes about 1h with the round loop (to get similar behavior for positive or negative coordinates)
list_2100_new<-list();list_2100_av<-list();list_2100_fin<-list()
for (v in names(list_dl)){
  for (k in scenario ){
    mod=names(list_df[[v]][[k]])
    for (z in mod){
      for (n in 1:nrow(list_df[[v]][[k]][[z]])){
        if (list_df[[v]][[k]][[z]]$y[n] < 0){
          list_df[[v]][[k]][[z]]$y[n]<- -round(abs(list_df[[v]][[k]][[z]]$y[n]) ) }else{
            list_df[[v]][[k]][[z]]$y[n]<- round(list_df[[v]][[k]][[z]]$y[n])}
        if (list_df[[v]][[k]][[z]]$x[n] < 0){
          list_df[[v]][[k]][[z]]$x[n]<- -round(abs(list_df[[v]][[k]][[z]]$x[n]) ) }else{
            list_df[[v]][[k]][[z]]$x[n]<- round(list_df[[v]][[k]][[z]]$x[n])}
      #list_df[[v]][[k]][[z]]$y<-round(list_df[[v]][[k]][[z]]$y)
      #list_df[[v]][[k]][[z]]$x<-round(list_df[[v]][[k]][[z]]$x) }
      }
      list_2100_new[[v]][[k]][[z]]<-aggregate(list_df[[v]][[k]][[z]][,3] ,data = list_df[[v]][[k]][[z]], by = list(list_df[[v]][[k]][[z]]$x,list_df[[v]][[k]][[z]]$y),FUN = mean )
      colnames(list_2100_new[[v]][[k]][[z]])<-colnames(list_df[[v]][[k]][[z]])
    }
    list_2100_av[[k]][[v]]<-Reduce(my_merge, list_2100_new[[v]][[k]])
    
    if(ncol(list_2100_av[[k]][[v]]) > 3){ list_2100_fin[[k]][[v]]<-data.frame(x=list_2100_av[[k]][[v]]$x, y=list_2100_av[[k]][[v]]$y, z = apply(list_2100_av[[k]][[v]][,-c(1:2)], 1, mean) );colnames(list_2100_fin[[k]][[v]])[3]<-v  
    }else{list_2100_fin[[k]][[v]]<-list_2100_av[[k]][[v]];colnames(list_2100_fin[[k]][[v]])[3]<-v
    }
  }
  print(paste("Done with ", v))
}

# extract data per scenario
ssp119<-Reduce(my_merge, list_2100_fin$ssp119) 
ssp126<-Reduce(my_merge, list_2100_fin$ssp126) 
ssp245<-Reduce(my_merge, list_2100_fin$ssp245) 
ssp370<-Reduce(my_merge, list_2100_fin$ssp370) 
ssp434<-Reduce(my_merge, list_2100_fin$ssp434)
ssp460<-Reduce(my_merge, list_2100_fin$ssp460)
ssp585<-Reduce(my_merge, list_2100_fin$ssp585)
list_sc_oce_2100<-list(ssp119,ssp126,ssp245,ssp370,ssp434, ssp460,ssp585)
names(list_sc_oce_2100)<-c("ssp119","ssp126","ssp245","ssp370","ssp434", "ssp460","ssp585")
saveRDS(list_sc_oce_2100, "/export/lv6/user/pramond/CMIP6/list_sc_oce_2100_360.RData")

list_sc_oce_2100<-readRDS("/export/lv6/user/pramond/CMIP6/list_sc_oce_2100_360.RData")


# explore final averaged model output
list_plots<-list() # empty list to collect all plots loaded
for (i in names(list_dl) ){
  list_plots[[i]]<-ggplot( data = list_sc_oce_2100$ssp245, aes_string(x="x", y="y", fill = i) ) +
    geom_tile()+
    ggtitle(paste(i))+
    scale_fill_viridis(option = "B")
}
do.call("ggarrange", c(list_plots)) # plot the plot-list 


############################################################################################
# Extract and unite data for contemporaneous ocean
############################################################################################

## load all nectcdf of all yearly variables  
########################################################

# list all variables and scenarios
scenario<-c("ssp119","ssp126","ssp245","ssp370","ssp434","ssp460","ssp534-over","ssp585") 
var<-grep("[.]|the|so",list.files("/export/lv6/user/pramond/CMIP6/"), invert = TRUE, value=TRUE) # list all variables with monthly frequency

# Loop to import and transform all netcdfs for each variable and save it in a list
list_2015_nc<-list();list_2015_wg<-list();list_2015_df<-list();list_dl<-list() # create empty lists to collect loop results
for (v in var){ # loop for each variable
  # retrieve information on files for each variable
  chl=list.files(paste("/export/lv6/user/pramond/CMIP6/", v, '/', sep = "" ), pattern = "*.nc") # list all netcdfs for this variable
  chl<-chl[grep("2015", chl)]
  # Collect information on netcdf files based on their names
  dl=data.frame(files = chl,
                model = stringr::str_split_fixed(chl, "_", n = 7)[,3],
                scenarios = stringr::str_split_fixed(chl, "_", n = 7)[,4],
                year = stringr::str_split_fixed(chl, "_", n = 7)[,7]
  )
  dl$start<-as.numeric(gsub("-2100.nc","", dl$year))
  dl$bands<-2100-dl$start+1 # get number of bands (years in netcdf)
  # filter out models with different grid and coordinate systems
  dl<-dl[dl$model != "CESM2-WACCM", ]
  dl<-dl[grep("NorESM2", dl$model, invert = TRUE), ]
  dl<-dl[grep("cmip5", dl$scenarios, invert = TRUE), ]
  
  # load all netcdf for a variable
  for (i in 1:nrow(dl)){
    # load raw raster
    path<-paste("/export/lv6/user/pramond/CMIP6/",v, "/", dl$files[i], sep = "") # make path
    r<-raster(x = path, band = 1, level = 1,  stopIfNotEqualSpaced = FALSE) # read/import it (level 1, correspond to the surface layer)
    # convert to normalize coordinates and convert to WGS84
    rtp<-as.data.frame(rasterToPoints(r)) # convert raster to data.frame
    rtp$x<-scales::rescale(rtp$x, to = c(0, 360)) # rescale model grid
    rtp$y<-scales::rescale(rtp$y, to = c(-90, 90)) # rescale model grid
    colnames(rtp)[3]<-paste(v,dl$model[i], dl$scenarios[i], sep = ".") # change raster variable name
    r2<-rasterFromXYZ(rtp) # convert back to raster with new coordinates
    crs(r2) <- CRS('+init=EPSG:4326') # force WGS84 crs
    # convert to data.frame
    rdf<-as.data.frame(rasterToPoints(r2)) # turn back to data.frame
    # save outputs in lists
    list_2015_nc[[v]][[dl$scenarios[i]]][[dl$model[i]]]<-r # save original rasters
    list_2015_wg[[v]][[dl$scenarios[i]]][[dl$model[i]]]<-r2 # save rasters with modified coordinates
    list_2015_df[[v]][[dl$scenarios[i]]][[dl$model[i]]]<-rdf # save data.frames with standardized coordinates
    #remove(r)
  }
  # save infos on netcdf saved in lists
  list_dl[[v]]<-dl
}  
dl<-bind_rows(list_dl, .id = "variable")

# explore raw model output
var_check="no3";var # select variable
sc="ssp585" # select climate scenario
list_plots<-list() # empty list to collect all plots loaded
for (i in unique(dl[dl$variable %in% var_check & dl$scenarios %in% sc, "model"] ) ){
  list_plots[[i]]<-gplot( list_2015_wg[[var_check]][[sc]][[i]] ) + geom_tile(aes(fill = value))+
    ggtitle(paste(i))+
    scale_fill_viridis(option = "B")
}
do.call("ggarrange", c(list_plots)) # plot the plot-list 

## load all nectcdf of all yearly variables  
########################################################

# Loop to import and transform all netcdfs for each variable and save it in a list
list_mon15_nc<-list();list_mon15_wg<-list();list_mon15_df<-list();list_dlt<-list() # create empty lists to collect loop results
for (v in c("thetao","so")){ # loop for each variable
  # retrieve information on files for each variable
  chl=list.files(paste("/export/lv6/user/pramond/CMIP6/", v, '/', sep = "" ), pattern = "*.nc") # list all netcdfs for this variable
  chl<-chl[grep("2015",chl)]
  # Collect information on netcdf files based on their names
  dlt=data.frame(files = chl,
                 model = stringr::str_split_fixed(chl, "_", n = 7)[,3],
                 scenarios = stringr::str_split_fixed(chl, "_", n = 7)[,4],
                 year = stringr::str_split_fixed(chl, "_", n = 7)[,7]
  )
  dlt$start<-as.numeric(gsub("01-.*","", dlt$year))
  dlt$bands<-2100-dlt$start+1 # get number of bands (years in netcdf)
  dlt<-unique(dlt)
  # filter out models with different grid and coordinate systems
  dlt<-dlt[!dlt$model %in% c("CAMS-CSM1-0","CESM2-WACCM", "BCC-CSM2-MR","CAS-ESM2-0",
                             "CESM2-WACCM", "E3SM-1-1", "E3SM-1-0" , "E3SM-1-1-ECA", 
                             "FGOALS-f3-L", "FGOALS-g3", "INM-CM4-8", "INM-CM5-0",
                             "MRI-ESM2-0", "TaiESM1", "MIROC6"), ]
  rownames(dlt)=NULL
  
  # load all netcdf for a variable
  for (i in 1:nrow(dlt)){
    # load raw raster
    path<-paste("/export/lv6/user/pramond/CMIP6/",v,"/",dlt$files[i], sep = "")
    r<-brick(x = path,stopIfNotEqualSpaced = FALSE, level =1)
    rm<-calc(r[[names(r)[1:12]]], fun = mean)
    # convert to normalize coordinates and convert to WGS84
    rtp<-as.data.frame(rasterToPoints(rm)) # convert raster to data.frame
    rtp$x<-scales::rescale(rtp$x, to = c(0, 360)) # rescale model grid
    rtp$y<-scales::rescale(rtp$y, to = c(-90, 90)) # rescale model grid
    colnames(rtp)[3]<-paste(v,dlt$model[i], dlt$scenarios[i], sep = ".") # change raster variable name
    r2<-rasterFromXYZ(rtp) # convert back to raster with new coordinates
    crs(r2) <- CRS('+init=EPSG:4326') # force WGS84 crs
    # convert to data.frame
    rdf<-as.data.frame(rasterToPoints(r2)) # turn back to data.frame
    # save outputs in lists
    list_mon15_nc[[v]][[dlt$scenarios[i]]][[dlt$model[i]]]<-r # save original rasters
    list_mon15_wg[[v]][[dlt$scenarios[i]]][[dlt$model[i]]]<-r2 # save rasters with modified coordinates
    list_mon15_df[[v]][[dlt$scenarios[i]]][[dlt$model[i]]]<-rdf # save data.frames with standardized coordinates
    print(paste(v, " / ",i, " / ",nrow(dlt), " models",  sep = ""))
  }
  # save infos on netcdf saved in lists
  list_dlt[[v]]<-dlt
}  
dlt<-bind_rows(list_dlt, .id = "variable")

# explore raw model output
var_check="thetao";var # select variable
sc="ssp126";scenario # select climate scenario
list_plots<-list() # empty list to collect all plots loaded
for (i in unique(dlt[dlt$variable %in% var_check & dlt$scenarios %in% sc, "model"] ) ){
  list_plots[[i]]<-gplot( list_mon15_wg[[var_check]][[sc]][[i]] ) + geom_tile(aes(fill = value))+
    ggtitle(paste(i))+
    scale_fill_viridis(option = "B")
}
do.call("ggarrange", c(list_plots)) # plot the plot-list 

## normalize model grids, format and summarize values per scenarios
###################################################################

# explore available data for each variable scenarios
list_dl$thetao<-dlt[dlt$variable == "thetao",];list_dl$so<-dlt[dlt$variable == "so",]
netcdf<-bind_rows(list_dl, .id = "variable")
table(netcdf$variable,netcdf$scenarios)
list_df<-c(list_2015_df, list_mon15_df)
my_merge <- function(df1, df2){merge(df1, df2, by = c("x", "y"), all = TRUE)}

# loop to get mean of all the models for each variable and selected scenarios
# takes about 1h with the round loop (to get similar behavior for positive or negative coordinates)
list_2015_new<-list();list_2015_av<-list();list_2015_fin<-list()
for (v in names(list_dl)){
  for (k in c("ssp119", "ssp126","ssp245", "ssp370", "ssp585") ){
    mod=names(list_df[[v]][[k]])
    for (z in mod){
      for (n in 1:nrow(list_df[[v]][[k]][[z]])){
        if (list_df[[v]][[k]][[z]]$y[n] < 0){
          list_df[[v]][[k]][[z]]$y[n]<- -round(abs(list_df[[v]][[k]][[z]]$y[n]) ) }else{
            list_df[[v]][[k]][[z]]$y[n]<- round(list_df[[v]][[k]][[z]]$y[n])}
        if (list_df[[v]][[k]][[z]]$x[n] < 0){
          list_df[[v]][[k]][[z]]$x[n]<- -round(abs(list_df[[v]][[k]][[z]]$x[n]) ) }else{
            list_df[[v]][[k]][[z]]$x[n]<- round(list_df[[v]][[k]][[z]]$x[n])}
        #list_df[[v]][[k]][[z]]$y<-round(list_df[[v]][[k]][[z]]$y)
        #list_df[[v]][[k]][[z]]$x<-round(list_df[[v]][[k]][[z]]$x) }
      }
      list_2015_new[[v]][[k]][[z]]<-aggregate(list_df[[v]][[k]][[z]][,3] ,data = list_df[[v]][[k]][[z]], by = list(list_df[[v]][[k]][[z]]$x,list_df[[v]][[k]][[z]]$y),FUN = mean )
      colnames(list_2015_new[[v]][[k]][[z]])<-colnames(list_df[[v]][[k]][[z]])
    }
    list_2015_av[[k]][[v]]<-Reduce(my_merge, list_2015_new[[v]][[k]])
    
    if(ncol(list_2015_av[[k]][[v]]) > 3){ list_2015_fin[[k]][[v]]<-data.frame(x=list_2015_av[[k]][[v]]$x, y=list_2015_av[[k]][[v]]$y, z = apply(list_2015_av[[k]][[v]][,-c(1:2)], 1, mean) );colnames(list_2015_fin[[k]][[v]])[3]<-v  
    }else{list_2015_fin[[k]][[v]]<-list_2015_av[[k]][[v]];colnames(list_2015_fin[[k]][[v]])[3]<-v
    }
  }
  print(paste("Done with ", v))
}

# extract data per scenario
ssp119<-Reduce(my_merge, list_2015_fin$ssp119) 
ssp126<-Reduce(my_merge, list_2015_fin$ssp126) 
ssp245<-Reduce(my_merge, list_2015_fin$ssp245) 
ssp370<-Reduce(my_merge, list_2015_fin$ssp370) 
ssp585<-Reduce(my_merge, list_2015_fin$ssp585)
list_sc_oce_2015<-list(ssp119,ssp126,ssp245,ssp370,ssp585)
names(list_sc_oce_2015)<-c("ssp119","ssp126","ssp245","ssp370","ssp585")
saveRDS(list_sc_oce_2015, "/export/lv6/user/pramond/CMIP6/list_sc_oce_2015_360.RData")

# merge data for scenarios in 2015 (is quite similar calibration)
oce_2015<-list_sc_oce_2015
for (i in names(oce_2015)){colnames(oce_2015[[i]])[-c(1:2)]<-paste(i, colnames(oce_2015[[i]])[-c(1:2)], sep = "_")}
my_merge <- function(df1, df2){merge(df1, df2, by = c("x", "y"), all = TRUE)}
oce_2015<-Reduce(my_merge, oce_2015)
coord<-oce_2015[,1:2]
for (v in names(list_dl)){coord<-cbind(coord, apply(oce_2015[, grep(v, colnames(oce_2015)) ] ,MARGIN = 1, FUN = mean) )}
colnames(coord)<-colnames(list_sc_oce_2015$ssp119)
oce_2015<-coord
saveRDS(oce_2015, "/export/lv6/user/pramond/CMIP6/oce_2015_360.RData")

# explore final averaged model output
list_plots<-list() # empty list to collect all plots loaded
for (i in names(list_dl) ){
  list_plots[[i]]<-ggplot( data = oce_2015, aes_string(x="x", y="y", fill = i) ) +
    geom_tile()+
    ggtitle(paste(i))+
    scale_fill_viridis(option = "B")
}
do.call("ggarrange", c(list_plots)) # plot the plot-list 

