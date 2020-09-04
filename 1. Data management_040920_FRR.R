#- DATA MANAGEMENT FILE

#---------------------------------------------------------------
#- Remove everything in memory
rm(list=ls())

#---------------------------------------------------------------
#- set the working directory
setwd("H:\\My Documents\\MSc Statistics\\Msc Dissertation\\Data")

#load required packages
library(tidyverse)
library(GSODR)
library(spatstat)

#---------------------------------------------------------------------------------
###READ DATA
data <- read.csv("data_10062020.csv")
load("H:/My Documents/MSc Statistics/Msc Dissertation/Data/chileka.RData")

library(lubridate)
data$date <- dmy(data$Date.of.BC.collection) #1 failed to parse/NA date value

#US National Centres for Environmental Information 
chileka$date <- chileka$YEARMODA

#merge data with climate data
library(dplyr)
cdata <- left_join(data, chileka, by=c("date") )

# -------------------------------------------------------------------------------
#Keep S Typhi cases in STRATAA
newd <- cdata[ which(cdata$S..Typhi.NTS == "S.TYPHI" & 
                       cdata$TyVAC.STRATAA == "STRATAA" ),  ]
save(newd, file = "newd.RData")

# --------------------------------------------------------------------------------
#Data cleaning

#extract coordinates for typhoid
data_coords <- dplyr::select(newd, "Longitude", "Latitude")

range(data_coords$Longitude)
range(data_coords$Latitude)#some of the values are outside Ndirande range

#- remove coordinates outside of Ndirande
data_coords2 <- subset(data_coords, Longitude > 35.03)
newd2 <- subset(newd, Longitude >35.03)

#re-project the coordinates data
#convert coordinates to spatial points
library(sp)
data_coords2 <- SpatialPoints(data_coords2) #points converted into spatial object
summary(data_coords2)$is.projected
coordinates(data_coords2)
proj4string(data_coords2) <- CRS("+proj=longlat") #reproject the points
#https://epsg.io/?q=Malawi.+Zambia+and+Zimbabwe
data.UTM <- spTransform(data_coords2, CRS("+init=epsg:32736"))
data.UTM
plot(data.UTM); axis(1); axis(2)
data.UTM.df <- as.data.frame(data.UTM)

par(mfrow = c(1, 2))
plot(data_coords2); axis(1); axis(2)
plot(data.UTM); axis(1); axis(2)

#--------------------------------------------------------------------------------
#read census (HH & individual) data
census_hh <- read.csv("census_household_level_raw.csv")
census_m <- read.csv("census_member_level_raw.csv")

#-------------------------------------------------------------------------------
#merge census_mm to census_h which has geolocations
library(tidyr)
names(census_m)
#columns m02_hhidm, m01_hhidb & pid are incomplete
#Split mo3_hhmid so that the hhid can be merged with census_hh
census_mm <- separate(data = census_m, col = m03_hhmid, 
                      into = c("hhid", "indid"), 
                      sep = "-")

#columns h02_hhidm, h01_hhidb & pid are incomplete
census_hh$hhid <- paste0(census_hh$h02_hhidb, census_hh$h03_hhidm)

#- merge census member to census household
census_indv <- merge(census_hh, census_mm, by.x='hhid', by.y='hhid',  
                     sort=TRUE)

#extract coordinates from census data
census_coords <- dplyr::select(census_indv, "h07_gpslon", "h07_gpslat")

# remove rows with NA  
census_coords2 <- census_coords[complete.cases(census_coords), ]
census_coords3 <- subset(census_coords2, h07_gpslon > 35)

#convert coordinates to spatial points
#http://rstudio-pubs-static.s3.amazonaws.com/19879_7e13ab80d5ed416c8e235bd6bb93cf3e.html
census.dec <- SpatialPoints(census_coords3, proj4string = CRS("+proj=longlat"))
census.UTM <- spTransform(census.dec, CRS("+init=epsg:32736 "))
census.dec.df <- as.data.frame(census.dec)
census.UTM.df <- as.data.frame(census.UTM)

#Compare the two:
par(mfrow=c(1,2))
plot(census.dec); axis(1); axis(2)
plot(census.UTM); axis(1); axis(2)

#compare whether shape of data for dec and UTM is the same
plot(census.dec.df$h07_gpslon, census.dec.df$h07_gpslat); axis(1); axis(2)
plot(census.UTM.df$h07_gpslon, census.UTM.df$h07_gpslat); axis(1); axis(2)

#---------------------------------------------------------------
#- STUDY HAD NO DEFINED POLYGON FOR STUDY AREA
#- CREATE A POLYGON FROM THE CENSUS DATA

#Create a polygon from the census data
#create a polygon
census_poly <- clickpoly(add = TRUE)
class(census_poly)
plot(census_poly); axis(1); axis(2)
points(census.UTM)

plot(census.UTM.df$h07_gpslon, census.UTM.df$h07_gpslat); axis(1); axis(2)
census_poly2 <- clickpoly(add = TRUE)
class(census_poly2)
plot(census_poly2); axis(1); axis(2)
points(census.UTM)
plot(census_poly); axis(1); axis(2)
plot(census_poly2); axis(1); axis(2)

#convert the owin to spatialpolygons object and save the file
# convert spatstat objects to sp classes
poly <- as(census_poly, "SpatialPolygons")
poly
plot(poly)

IDs <- sapply(slot(poly, "polygons"), function(x) slot(x, "ID"))
dfpoly <- data.frame(rep(0, length(IDs)), row.names=IDs)
polyf <- SpatialPolygonsDataFrame(poly, dfpoly)

library(rgdal)
writeOGR(polyf,  layer = "npoly", dsn = '.',  driver = "ESRI Shapefile") #save to disk

#save the second polygon
poly2 <- as(census_poly2, "SpatialPolygons")
plot(poly2); axis(1); axis(2)

IDs <- sapply(slot(poly, "polygons"), function(x) slot(x, "ID"))
dfpoly <- data.frame(rep(0, length(IDs)), row.names=IDs)
polyf2 <- SpatialPolygonsDataFrame(poly2, dfpoly)

library(rgdal)
writeOGR(polyf2,  layer = "ndirandepoly", dsn = '.',  driver = "ESRI Shapefile") #save to disk

#-----------------------------
#- Read shapefile
poly <- readOGR(dsn=".", "ndirandepoly")
plot(poly)

#--------------------------------------------------------------------------------
#See if household WASH indicators can be included as spatial covariates
#this is because they have been observed at other locations in addition
#to the locations of the cases
wash <- read.csv("wash_04022020.csv")
names(wash)

#merge census census hh to wash data
#Try merging the wash data to the census data
#merge census census hh to wash data
newd.wash <- merge(census_hh,wash, by.x='hhid',
                   by.y='household_id',sort=TRUE)
#reduces number of points to 14,180
#attaching wash data to indvs reduces # of cases from 165 to 65

newd.wash2 <- newd.wash[!is.na(newd.wash$h07_gpslon),]
newd.wash3 <- newd.wash[!is.na(newd.wash$h07_gpslat),]
plot(newd.wash3$h07_gpslon, newd.wash3$h07_gpslat)
newd.wash4 <- subset(newd.wash3, h07_gpslon > 35)

#convert coordinates to spatial points
#http://rstudio-pubs-static.s3.amazonaws.com/19879_7e13ab80d5ed416c8e235bd6bb93cf3e.html
coordinates(newd.wash4) <- c("h07_gpslon", "h07_gpslat")
proj4string(newd.wash4) <- CRS("+proj=longlat") #decimal/latlon crs
wash.UTM <- spTransform(newd.wash4, CRS("+init=epsg:32736"))
wash.UTM.df <- as.data.frame(wash.UTM)

################################################################
#- GET/GENERATE ENVIRONMENTAL COVARIATES FOR STAT INFERENCE

#- WASH score
# Derived as per the PCA and linear geostatistical file
#- As generated above using spatial prediction from a linear  geostatistical model

#---------------------------------------------------------------------------
#- Distance to the health facility (HF) in meters
#- https://rpubs.com/ricardo_ochoa/415839

#- In order to have our distance measured in meters I converted
#- the coordinates from longitude/latitude to easting/northing

names(data.UTM.df)

#create a dummy/empty raster
r <- raster(ncols=100, nrows=100,
            xmn = 716000, xmx = 721000, ymn = 8253000, ymx = 8258000)
r <- rasterize(data.UTM, r, fun=function(x,...)length(x))
r
plot(r)

#set coordinates for Ndirande health center
hclon <- 35.0395908
hclat <- -15.7780157
ndirande_hc <- data.frame(hclon, hclat)

#assign projection to HF coords and re-project
cord.dec <- SpatialPoints(ndirande_hc, proj4string = CRS("+proj=longlat"))
cord.UTM <- spTransform(cord.dec, CRS("+init=epsg:32736"))
cord.UTM

# Distance from points to each raster cell
distances <- distanceFromPoints(object = r, xy = cord.UTM)
plot(distances)
distances 
distances2 <- crop(distances, extent(poly))
plot(distances2)
distancesf <- mask(distances2, poly)
plot(distancesf)

#-------------------------------------------------------------
#- Elevation

#dDownload elevation data from the rasterpackage
#- Raster data acquisition
#- https://www.gis-blog.com/r-raster-data-acquisition/
#get raster data for Ndirande
m <- readOGR(dsn=".", "ndirandeShape")
srtm <- getData('SRTM', lon=35.005833, lat=-15.786111)
plot(srtm, xlim = c(35,36))
plot(srtm, xlim = c(35,36), ylim=c(-16,-15))
plot(m, add=TRUE)
plot(srtm, xlim = c(35,35.3), ylim=c(-16,-15.6))
plot(m, add=TRUE)
plot(srtm, xlim = c(35,35.1), ylim=c(-15.85,-15.7))
plot(m, add=TRUE)
plot(srtm, xlim = c(35.02,35.07), ylim=c(-15.80,-15.73))
plot(m, add=TRUE)

sr <- "+init=epsg:32736"
elevation <- projectRaster(srtm, crs=sr)
elevation2 <- crop(elevation, poly)
plot(elevation2)
elevation.f <- rasterize(poly, elevation2, mask=TRUE) # plot the clipped raster plot(r.lcmask)
plot(elevation.f)

#-------------------------------------------------------------------------------
#- Create a PPP object for the Typhoid cases
#-----------------------------------------------------------------------

##create pixel image for the population density
library(maptools)
poly.owin <- as.owin(poly)
x_typh <- ppp(x = data.UTM.df$Longitude, y = data.UTM.df$Latitude, window = poly.owin)
x_typh <- as.ppp(x_typh) #to remove points falling outside window

#--------------------------------------------------------------------
## GO TO EXPLORATORY DATA ANALYSIS FILE
#----------------------------------------
#- 
x_census <- ppp(census.UTM.df$h07_gpslon, y=census.UTM.df$h07_gpslat,
                window = poly.owin) #452 points rejected
x_census <- as.ppp(x_typh) #to remove points falling outside window
plot(x_census)

#convert the ppp object to a pixel
census_pop <- as.im(x_census)

#work with as.im image and replace zero values with a very small number
# This is for the sake of the offset component which would not allow log of zeroes
census_im <- census_pop
census_im[census_im==0] <- 0.5*min(census_im[census_im!=0])
image(census_im)
image(log(census_im))
ppm(x_typh ~ offset(log(census_im)))
