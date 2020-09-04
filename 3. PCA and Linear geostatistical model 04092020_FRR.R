#-CREATING WASH SCORE USING PCA & LINEAR GEOSTATISTICAL MODEL

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
newd.wash4 <- subset(newd.wash3, h07_gpslon > 35) #remove obs outside Ndirande

#convert coordinates to spatial points
#http://rstudio-pubs-static.s3.amazonaws.com/19879_7e13ab80d5ed416c8e235bd6bb93cf3e.html
coordinates(newd.wash4) <- c("h07_gpslon", "h07_gpslat")
proj4string(newd.wash4) <- CRS("+proj=longlat") #decimal/latlon crs
wash.UTM <- spTransform(newd.wash4, CRS("+init=epsg:32736"))
wash.UTM.df <- as.data.frame(wash.UTM)

#carry out PCA of wash data
#-----------------------------
library(FactoMineR)
names(wash.UTM.df)
wash.active <- dplyr::select(wash.UTM.df, 
                             "h07_gpslon", "h07_gpslat", "hrooms" ,
                             "toilet_kind", "toilet_type", "wsource" )
#reverse ordering of categorical variables so that the highest level is the "best"
wash.active$rooms <- factor(recode(wash.active$hrooms,
                                   "1"=1, "2"=2, "3"=3, "4"=4, "5+"=5) )
wash.active$toi_kind <- factor(recode(wash.active$toilet_kind,
                                      "1"=4, "2"=3, "3"=2, "4"=1) )
wash.active$toi_type <- factor(recode(wash.active$toilet_type,
                                      "1"=5, "2"=4, "3"=3, "4"=2, "5"=1))
wash.active$waters <- factor(recode(wash.active$wsource,
                                    "1"=6, "2"=4, "3"=5, "4"=3, "5"=2, "6"=1))
#2 Public tap outside house, 3=Private tap outside house
wash.active2 <- dplyr::select(wash.active, "rooms", "toi_kind", "toi_type",
                              "waters")

#convert the variables to numeric
rooms <- as.numeric(levels(wash.active2$rooms))[wash.active2$rooms]
toi_kind <- as.numeric(levels(wash.active2$toi_kind))[wash.active2$toi_kind]
toi_type <- as.numeric(levels(wash.active2$toi_type))[wash.active2$toi_type]
waters <- as.numeric(levels(wash.active2$waters))[wash.active2$waters]

#- CREATE DUMMY VARIABLES FOR THE PCA

#------------------------------------------------------------
#- VARIABLE DEFINITIONS IN QUESTIONNAIRE

#1. Number of rooms in the house: 1=1, 2=2, 3=3, 4=4, 5=5+)
#2 - Kind of toilet being used by the household: 1 = They do not have one,
#-  2 = Public,, 3 = Shared with neighbours, 4 = Household use only
#3- Type of toilet used: (1 = Open defecation, 2 = Pit latrine with wood/soil 
 #- floor, 3 = Pit latrine with slab, 4 = Flush/pour flush
#4- Primary source of drinking water: (1 = Unprotected, 2 = Protected borehole/well, 3 = Public tap/standpipe,
  #4 = Public tap outside of house, 5 = Private tap outside of house, 6 = Piped to house

#------------------------------------------------------------
#- Number of rooms in house: left as continuous variable: 1 to 5
table(rooms)

#---------------------
#- Kind of toilet used in the household

#- Household use only
toilet_hh_useonly <- ifelse(toi_kind==4, 1, 0)

#- Shared with neighbours
toilet_shared_neighb <- ifelse(toi_kind==3, 1, 0)

#- Public toilet
toilet_public <- ifelse(toi_kind==2, 1, 0)

#-----------------------
#- Type of toilet being used
#- NB, There was no person that had a level 1 answer
#- Level 1 = Refused to answer, 2 = Open defaecation

#- Flush or pour toilet
toilet_flush_pour <- ifelse(toi_type==5, 1, 0)

#-Pit latrine with slab
toilet_pitwithslab <- ifelse(toi_type==4, 1, 0)

#- Pit latrine with wood/soil floor
toilet_pitwithwoodsoilfloor <- ifelse(toi_type==3, 1, 0)

#-------------------------
#- Primary source of drinking water

#- Piped tp the household
water_pipedtohh <- ifelse(waters==6, 1, 0)

#- Private tap outside of house
water_privatetap <- ifelse(waters==5, 1, 0)

#- Public tap outside of house
water_publictap <- ifelse(waters==4, 1, 0)

#- Public tap/standpipe
water_publicstandpipe <- ifelse(waters==3, 1, 0)

#- Protected well/borehole
water_protectwellborehole <- ifelse(waters==2, 1, 0)

#-------------------------------------------------------------

#- Create a data frame from the continous variable and dummy variables
library(FactoMineR)
library(factoextra)
wash.active3 <- data.frame(rooms, 
                           toilet_hh_useonly, toilet_shared_neighb, toilet_public,
                           toilet_flush_pour, toilet_pitwithslab, toilet_pitwithwoodsoilfloor,
        water_pipedtohh, water_privatetap, water_publictap, water_publicstandpipe, water_protectwellborehole)
head(wash.active3)
wash.pca <- prcomp(wash.active3, scale = T)
#wash.pca2 <- PCA(wash.active3, scale.unit = TRUE, ncp = 5, graph = TRUE)

#Scree plot
fviz_eig(wash.pca)

wash.score <- wash.pca$x[,1]
range(wash.score)
eig.val <- get_eigenvalue(wash.pca)
eig.val

#-------------------------------------------------------
#Run a linear geostatistical model in PrevMap
library(PrevMap)

#Running the linear geostatistical model in PrevMap was 

#computantionally intensive (Takes more than 2 days to converge)
#- because WASH data has 14,136 data points

#Considered using a stochastic partial differential equation (SPDE)
#approach to reduce computaional time

#Creating a mesh using R-INLA for the linear geostatistical model

#function to extract coordinates from polygon
#extract coordinates from the polygon
library(rgdal)
poly <- readOGR(dsn = ".", layer = "ndirandepoly")

getcoords <- function(polygon)
{
  xycoords <- list()
  for(i in 1:length(polygon@polygons[[1]]@Polygons))
  {
    xycoords[[i]] <- polygon@polygons[[1]]@Polygons[[i]]@coords
  }
  xycoords<- Reduce(rbind, xycoords)
  xycoords
}

ndi.coords <- getcoords(poly)
plot(ndi.coords)
lines(ndi.coords)

##Creating a mesh using R-INLA for the linear geostatistical model
loc <- ndi.coords
library(INLA)
mesh <- inla.mesh.2d(loc, offset=c(500,1000),
                     cutoff = 1,  max.edge=c(300, 600))
plot(mesh)
points(loc, col= "red")
class(mesh)

#check which transformation to use for the data
par(mfrow=c(1,2))
hist(wash.score, xlab="Wash score", main = "", 
     xlim = c(-3,4), breaks=20) #identity transformation is the best
box()

#Preliminary diagnostic for the wash spatial data
#diag.cor <- spat.corr.diagnostic(wash.score~1,coords=~h07_gpslon+h07_gpslat,
#  data=wash.UTM.df,likelihood = "Gaussian")
#Takes a long time to converge - too many data points

#Fit the linear geostatistical model using the SPDE approach
geo.fit.wash <- linear.model.MLE(wash.score ~1,
                                 coords=~h07_gpslon+h07_gpslat,
                                 start.cov.pars=c(1,2),
                                 kappa=1, #Matern function=1
                                 data=wash.UTM.df,
                                 # fixed.rel.nugget = 0,
                                 SPDE = TRUE, mesh=mesh,
                                 method="nlminb")

#variogram to assess model
#diag.fit <- variog.diagnostic.lm(geo.fit.wash, n.sim=50)
#Takes a long time to converge - too many data points

#-----------------------------------------------
#- Goodness of fit of model
#--------------------------------------------------------------------
#Variograms for linear geostatistical model

#- Use geoR to compute variograms - quicker than PrevMap for 
#scenario with many data points
library(geoR)

#create geoR data
wash.geor <- data.frame(wash.UTM.df$h07_gpslon, wash.UTM.df$h07_gpslat, wash.score)

#convert data to geoR format
wash.geor2 <- as.geodata(wash.geor)

#expoloratory analysis for geoR data
plot(wash.geor2, lowess=TRUE)

#variogram(wash.geor, wash.score, 
#coords = ~(wash.UTM.df.h07_gpslon + wash.UTM.df.h07_gpslon.1 ) )
#PrevMap variogram takes a while to fit

#Compute empirical variogram
wash.empvar <- variog(wash.geor2)
plot(wash.empvar, type="b")
lines(wash.empvar, type = "b", lty = 2)
eyefit(wash.empvar)

#to check for spatial correlation
wash.empvar <- variog(wash.geor2)
plot(wash.empvar, col="red")
eyefit(wash.empvar)

plot(wash.empvar, col = "red")
lines.variomodel(wash.empvar, cov.model = "exponential", 
                 cov.pars = c(4.32, 1474.36),
                 nugget = 0.54)
wash.empvar.env <- variog.mc.env(wash.geor2, obj.variog = wash.empvar, nsim = 999)

par(mfrow=c(1,1))
plot(wash.empvar, envelope = wash.empvar.env, col = "red", lwd=2)
#lines.variomodel(wash.empvar, cov.model = "exponential", 
#                cov.pars = c(4.32, 1474.36),
#               nugget = 0.54, col = "green", lwd=2)

#diag.cor <- spat.corr.diagnostic(wash.score~1,coords=~h07_gpslon+h07_gpslat,
#data=wash.UTM.df,likelihood = "Gaussian")

#INTERPRETATION
#- Diagnostics show presence of spatial autocorrelation

#-----------------------------------------------------------------------

#--------------------------------------------
#- interpolate/ make prediction
#Failed to use the polygon as grid so decided to extract the coordinates from polygon

library(splancs)
ndirandegrid <- gridpts(as.matrix(ndi.coords),xs=100,ys=100)
plot(ndirandegrid)

#use PrevMap tp make prediction
#Use logit as scale of prediction because it is continuous
pred.wash <- spatial.pred.linear.MLE(geo.fit.wash, grid.pred = ndirandegrid,
                                     scale.predictions="logit",
                                     standard.errors = TRUE)
plot(pred.wash, type="logit", xlab = "Easting", 
     ylab = "Northing")
class(pred.wash)

#Extract the logit values for the model
ndirandegrid2 <- as.data.frame(ndirandegrid)
ndirandegrid2$pred.wash.logit <- as.numeric(unlist(pred.wash$logit[1]))

#convert the result to a raster
#- https://gis.stackexchange.com/questions/188070/converting-data-from-data-frame-into-raster-file-using-r
class(ndirandegrid2)
library(sp)
library(raster)
coordinates(ndirandegrid2) <- ~V1+V2
gridded(ndirandegrid2) <- TRUE
str(ndirandegrid2@data)
wash.score.r <- raster(ndirandegrid2, "pred.wash.logit") 
class(wash.score.r)
plot(wash.score.r, xlab = "Easting", ylab="Northing")

#------------------------------------------------------------------
