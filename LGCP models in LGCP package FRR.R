#- This script file presents a summary of spatial and spatio-temporal
#- LGCP models fitted in the LGCP package

#- Key messages
#- Models fitted in this package use a Bayesian estimation
#- Models are computationally demanding (takes approximately 48 hours
# for full model to converge)
#- Spatial LGCP model converged. Results were not reported in the dissertation

#- Temporal data was very sparse (only a few cases were observed every few months)
#- The package, according to the published paper, seems to work well for
#- uninterrupted times (i.e. uninterrupted sequence e.g. 1,2,3,4 not 1,4)

#- Spatio-temporal model only converged for aggregated times (e.g.
# years, months or quarters). 
#- All the converged spatio-temporal models showed presence of both spatial
# and temporal correlation

#-However, the main interest in modelling typhoid fever is in looking at
# temporal trends for time in days
#- My future work will look into developing methods for very sparse
# spatio-temporal data

#-----------------------------------------------------------
#LGCP spatial and spatio-temporal models in the LGCP package

#------------------
#- Spatial Model
#- NOTE: Data used similar to the spatstat model data

##################################################################################
#- LGCP Models in LGCP package
#--------------------------------

#- LGCP package mostly uses grouped/regional covariates
#- STRATAA/This study did not have grouped covariates
#- I therefore converted the covariates to rasters

#- Environmental covariates
rasterdata <- stack("spdff.tif")
names(rasterdata) <- c("pop", "dist", "elev", "amt", "amp","wash")

#unlist the stack for spatstat operations
list2env(setNames(unstack(rasterdata), names(rasterdata)), .GlobalEnv)

washnn <- resample(wash.score.r, pop, method="bilinear")

spdf <- stack(pop, dist, elev, washnn)
names(spdf) <- c("pop", "distance", "elev", "wash")
plot(spdf)

spdff <- as(spdf, "SpatialPixelsDataFrame")

library(lgcp)

# minimum contrast estimation
#- To help with selection of cell width
minimum.contrast(x_typh, model = "exponential", method = "g",
                 intens = density(x_typh), transform = log)

#choose cell width
chooseCellwidth(x_typh, cwinit = 100) #choose 100
chooseCellwidth(x_typh, cwinit = 60) #go for 60
chooseCellwidth(x_typh, cwinit = 32) #too small - choose 60

CELLWIDTH <- 60

# select extension
EXT <- 2

# perform polygon overlay operations and compute computational grid
spdff
polyolay <- getpolyol(data = x_typh, pixelcovariates = spdff,
                      cellwidth = CELLWIDTH, ext = EXT)

# set the interpolation type for each variable
spdff@data <- guessinterp(spdff@data)
spdff@data <- assigninterp(df = spdff@data, vars = c("pop"), 
                           value = "ArealWeightedSum")
class(spdff@data$pop)

#- Model with population as an offset
FORM_pop <- X ~ pop - 1
Zmat_pop <- getZmat(formula = FORM_pop, data = x_typh,
                    pixelcovariates = spdff, cellwidth = CELLWIDTH,
                    ext = EXT, overl = polyolay)
mm <- length(attr(Zmat_pop, "mcens"))
nn <- length(attr(Zmat_pop, "ncens"))
poisson.offset <- spatialAtRisk(list(X = attr(Zmat_pop, "mcens"),
                   Y = attr(Zmat_pop, "ncens"), Zm = matrix(Zmat_pop, mm, nn)))

#define the model
names(spdff)
FORM <- X ~ distance + elev + wash

#interpolate the covariate data onto the computational grid:
Zmat <- getZmat(formula = FORM, data = x_typh, pixelcovariates = spdff,
                cellwidth = CELLWIDTH, ext = EXT, overl = polyolay)

#under Poisson model, expect # of cases to be proportional to popln at risk
Zmat[, "pop"] <- log(Zmat[, "pop"])
Zmat[, "pop"][is.infinite(Zmat[, "pop"])] <- min(Zmat[, "pop"][!is.infinite(Zmat[, 
                 "pop"])])

# choose the covariance function
CF <- CovFunction(exponentialCovFct)

# plot the interpolated covariates
plot(Zmat, ask = FALSE)

# define the priors
priors <- lgcpPrior(etaprior = PriorSpec(
  LogGaussianPrior(mean = log(c(1, 50)), variance = diag(0.015, 2))),
  betaprior = PriorSpec(
    GaussianPrior(mean = rep(0, 4), variance = diag(10^6, 4))))

#- Define an offset for the model
BASEDR <- getwd()
lgcp_typh <- lgcpPredictSpatialPlusPars(formula = FORM, sd = x_typh, 
                                        Zmat = Zmat,
                                        model.priors = priors, 
                                        #model.inits = INITS, 
                                        spatial.covmodel = CF,
                                        cellwidth = CELLWIDTH, poisson.offset = poisson.offset, 
                                        mcmc.control = mcmcpars(
                                          mala.length = 1000000, burnin = 100000, retain = 900,
                                          adaptivescheme = andrieuthomsh(inith = 1, alpha = 0.5, C = 1,
                                       targetacceptance = 0.574)),
                                        output.control = setoutput(gridfunction = dump2dir(
                                          dirname = paste(BASEDR, "/lgcp.typh/", sep = ""), forceSave = TRUE)),
                                        ext = EXT)

#Model diagnostics
# plot the log target
par(mfrow=c(1,1))
plot(ltar(lgcp_typh), type = "s", xlab = "Iteration/900", ylab = "log target", ask = FALSE)

# compute and plot autocorrelations in the latent field
lagch <- c(1, 5, 15)
Sacf <- autocorr(lgcp_typh, lagch, inWindow = NULL)
library("fields")
for (i in 1:3) {
  image.plot(xvals(lgcp_typh), yvals(lgcp_typh), Sacf[, , i], zlim = c(-1, 1), axes = FALSE, 
             xlab = "", ylab = "", asp = 1, sub = paste("Lag:", lagch[i]), ask = FALSE)
  plot(spdf, add = TRUE, ask = FALSE) }

# produce traceplots of beta and eta
traceplots(lgcp_typh, ask = FALSE)

# produce autocorrelation plots of beta and eta
par(mfrow=c(2,3))
parautocorr(lgcp_typh, ask = FALSE)

# a summary table of beta and eta
parsum_lgcp_typh <- parsummary(lgcp_typh)
parsum_lgcp_typh
# the above converted to LaTeX
library("miscFuncs")
latextable(parsum_lgcp_typh, rownames = rownames(parsum_lgcp_typh), 
           colnames = c("Parameter", colnames(parsum_lgcp_typh)), 
           digits = 4)

# a text summary of the model parameters
textsummary(lgcp_typh, digits = 4)

# a plot of the prior and posterior densities
priorpost(lgcp_typh, ask = FALSE)

# the posterior covariance function
postcov(lgcp_typh, ask = FALSE)


#--------------------------------------------------------------------------
#- TEMPORAL MODELS
#- All temporal models were not run on 1,000,000 iterations as they were
# just trial runs to determine an acceptable time format that would
# help the model to converge

#-------------------
#- Time in months (months adjusted so as not to have interrupted time points)

#remove points not falling in polygon for the tempoeral covariates
#- https://stackoverflow.com/questions/47696382/removing-data-outside-country-map-boundary-in-r
newd2.UTM1 <- newd2.UTM[!is.na(over(newd2.UTM, as(poly, "SpatialPolygons"))),]
newd2.UTM2 <- as.data.frame(newd2.UTM1)

#sort the data frame according to date 
newd2.UTM2 <- dplyr::arrange(newd2.UTM2, date)
class(newd2.UTM2$date)
newd2.UTM2$YEAR <- year(newd2.UTM2$date)
newd2.UTM2$MONTH <- month(newd2.UTM2$date)

#Join this data to the climate data
newd2.UTM2 <- left_join(newd2.UTM2, chileka_clim, by = c("YEAR", "MONTH"))

#generate start date for study
startdate = '2016/10/05'
newd2.UTM2$startdate = as.Date(startdate,'%Y/%m/%d')

#calculate the difference in days
newd2.UTM2$t <- difftime(newd2.UTM2$date, newd2.UTM2$startdate, units = c("days"))
newd2.UTM3 <- newd2.UTM2[complete.cases(newd2.UTM2), ]

#------------------------------------------------------------
#- ST MODEL WITH TIME IN QUARTERS 

#---------------------------------------------------------------------------------------------------------------
#- SPATIO-TEMPORAL LGCP MODEL IN LGCP PACKAGE
#--------------------------------------------------------------------------------
#- Aggregate at quartery level
#calculate the difference in days
AB$t <- difftime(AB$date, AB$startdate, units = c("days"))
AB2Q<- AB[complete.cases(AB), ]
AB2Q$tQ <- ifelse(AB2$MONTH.x <= 3, 1,
                  ifelse(AB2$MONTH.x <= 6, 2,
                         ifelse(AB2$MONTH.x <= 9, 3,
                                ifelse(AB2$MONTH.x <= 12, 4, NA))))

xQ <- AB2$Longitude; yQ <- AB2$Latitude; tQ <- AB2$tQ
stdataQ <- cbind(xQ,yQ,tQ)
stdata1Q <- stdataQ[order(stdataQ[,ncol(stdataQ)]),] #sort the data
tlim1Q <- as.integer(c(1,4))

library(lgcp)
xytfQ <- stppp(list(data = stdata1Q, tlim = tlim1Q, window = poly.owin))
xytfQ

# perform polygon overlay operations and compute computational grid
polyolay_xytQ <- getpolyol(data = xytfQ, pixelcovariates = spdff, 
                           cellwidth = CELLWIDTH_xyt, 
                           ext = EXT_xyt)

# the statistical model for the main effets, $\beta$
FORM_xytQ <- X ~ dist + elev + wash + amt + amp + mean_temp + mean_prec + quarter
FORM.spatialQ <- X ~ dist + elev  + wash + amt + amp 
TEMPORAL.FORMULAQ <- t ~ quarter + mean_temp + mean_prec - 1

# perform interpolation of the covariates onto the
Zmat_xytQ <- getZmat(formula = FORM.spatialQ, data = xytfQ, 
                     pixelcovariates = spdff, 
                     cellwidth = CELLWIDTH_xyt, 
                     ext = 2, overl = polyolay_xytQ)

#create a dataframe for time
tdataQ <- dplyr::select(AB2Q, "tQ", "mean_temp", "mean_prec")
names(tdataQ)[names(tdataQ) == "tQ"] <- "t"

tdata2Q <- tdataQ %>% group_by(t) %>% summarise_all(funs(mean))

tdata2Q$quarter <- tdata2Q$t

tdata3Q <- dplyr::select(tdata2Q, "t", "quarter", "mean_temp", "mean_prec")

# choose last time point and number of proceeding time-points to include
TQ <- 4
LAGLENGTHQ <- 3

# bolt on the temporal covariates
ZmatList_xytQ <- addTemporalCovariates(temporal.formula = TEMPORAL.FORMULAQ, 
                                       T = TQ, laglength = LAGLENGTHQ, 
                                       tdata = tdata3Q, Zmat = Zmat_xytQ)
#define the priors
lgprior_xytQ <- PriorSpec(LogGaussianPrior(mean = log(c(1, 1, 1)),
                                           variance = diag(1, 3)))
gprior_xytQ <- PriorSpec(GaussianPrior(mean = rep(1*10^-15, 9), 
                                       variance = diag(0.01, 9)))
priors_xytQ <- lgcpPrior(etaprior = lgprior_xytQ, betaprior = gprior_xytQ)

#Select initial values for algorithm
#INITSQ<- lgcpInits(etainit = log(c(sqrt(2.8), 1286, 0.5)), betainit = NULL)

# choose the covariance function
CF_xytQ <- CovFunction(exponentialCovFct)

# run the MCMC algorithm, note that as this takes a long time, we present a
DIRNAME <- getwd()

lg_xytQ <- lgcpPredictSpatioTemporalPlusPars(formula = FORM_xytQ, xyt = xytfQ, T = TQ, 
                                             laglength = LAGLENGTHQ, 
                                             ZmatList = ZmatList_xytQ, model.priors = priors_xytQ, 
                                             spatial.covmodel = CF_xytQ, 
                                             cellwidth = CELLWIDTH_xyt, 
                                             mcmc.control = mcmcpars(mala.length = 1000, burnin = 100, 
                                retain = 9, adaptivescheme = andrieuthomsh(inith = 1, alpha = 0.5, C = 1, 
                                     targetacceptance = 0.574)), 
                                             output.control = setoutput(gridfunction = dump2dir(dirname = file.path(DIRNAME, 
                                                                                                                    "strataa_xyt22Q"), forceSave = TRUE)), ext = EXT_xyt)

parsum <- parsummary(lg_xytQ)
textsummary(lg_xytQ, digits = 4)
##CODE ABOVE WORKS
#------------------------------------------------------------------------

#Construct a sinusoidal (temporal) covariate
# https://stackoverflow.com/questions/20104680/sine-curve-fit-using-lm-and-nls-in-r
newd2.UTM3$sin1 <-cos(2*pi*t3/366)
newd2.UTM3$sin2 <-sin(2*pi*t3/366)

#- Create an STPPP object
x <- as.integer(newd2.UTM3$Longitude)
y <- as.integer(newd2.UTM3$Latitude)
t <- as.integer(newd2.UTM3$t)

stdata <- cbind(x,y,t)
stdata2 <- stdata[complete.cases(stdata), ] #just making sure that they are no missing data
tlim <- as.integer(c(1,1230))

x_typh

library(lgcp)
xyt <- stppp(list(data = stdata2, tlim = tlim, window = poly.owin))
#this step requires that there not be missing data
xyt #no of points reduced from 165 to 162. 1 missing location, 2 missing date

#------ Computations
minc <- minimum.contrast.spatiotemporal(xyt, model = "exponential",
                                        method = "g", transform = log, 
                                        spatial.dens = density.ppp(xyt),
                                        temporal.intens = muEst(xyt))

#- Create a new time variable/group observations by number of months since study commenced
#- Some time points deliberately grouped into one month to avoid having interrupted points
t3 <- ifelse(t <= 30, 1, ifelse (t <= 60, 2, ifelse(t <=91, 3, ifelse( t<= 120, 4,
           ifelse(t <= 150, 5, ifelse (t <= 180, 6, ifelse (t <= 210, 8, ifelse(t <= 240, 8,
           ifelse(t <= 273, 9, ifelse(t <= 300, 10, ifelse(t <= 330, 11, ifelse(t <= 360, 12, 
          ifelse(t <= 390, 13, ifelse(t <= 420, 14, ifelse(t <= 450, 15, ifelse(t <= 480, 15, 
          ifelse(t <= 510, 17,  ifelse(t <= 540, 18,  ifelse(t <= 570, 19,  ifelse(t <= 600, 20, 
          ifelse(t <= 630, 21,  ifelse(t <= 660, 22,  ifelse(t <= 690, 23,  ifelse(t <= 720, 24, 
          ifelse(t <= 750, 25,  ifelse(t <= 780, 26, ifelse(t <= 810, 27, ifelse(t <= 840, 28, 
          ifelse(t <= 870, 29,  ifelse(t <= 900, 30,  ifelse(t <= 930, 31,  ifelse(t <= 960, 32, 
          ifelse(t <= 990, 33, ifelse(t <= 1020, 34, ifelse(t <= 1050, 35, ifelse(t <= 1080, 36, 
          ifelse(t <= 1110, 37, ifelse(t <= 1140, 38, ifelse(t <= 1170, 39, ifelse(t <= 1120, 40, 
          ifelse(t <= 1150, 41,  ifelse(t <= 1180, 42,  ifelse(t <= 1210, 41,  ifelse(t <= 1240, 42,"NA"
                      )))) )))) )))) ) ))))))))))))))))))) ))))))))) )))

#remove NA rows
stdata3 <- cbind(x,y,t3)
stdata4 <- stdata3[complete.cases(stdata3), ] #just making sure that they are no missing data
tlim <- as.integer(c(1,1230))

x_typh

library(lgcp)
xyt2 <- stppp(list(data = stdata4, tlim = tlim, window = poly.owin))

# minimum contrast estimation
minc2 <- minimum.contrast.spatiotemporal(xyt2, model = "exponential", 
                                         method = "g", 
                                         transform = log, spatial.dens = density.ppp(xyt2), 
                                         temporal.intens = muEst(xyt2))
minc2

# select cell width
chooseCellwidth(xyt2, cwinit = 1400)
chooseCellwidth(xyt2, cwinit = 410)
chooseCellwidth(xyt2, cwinit = 50)
chooseCellwidth(xyt2, cwinit = 80)
chooseCellwidth(xyt2, cwinit= 120)
CELLWIDTH_xyt <- 80

# select extension
EXT_xyt <- 5

# perform polygon overlay operations and compute computational grid
polyolay_xyt <- getpolyol(data = xyt2, pixelcovariates = spdff, cellwidth = CELLWIDTH_xyt, 
                          ext = EXT_xyt)

# the statistical model for the main effets, $\beta$
names(spdff)
names(newd2.UTM3)
FORM_xyt <- X ~ dist + elev + amt + amp + wash + sin1 + sin2 + mean_temp + mean_prec + mean_dew + mean_winds 
FORM.spatial <- X ~ pop + dist + elev + amt + amp + wash
TEMPORAL.FORMULA <- t ~ sin1 + sin2 + mean_temp + mean_prec + mean_dew + mean_winds - 1

# set the interpolation type for each variable
spdff@data <- guessinterp(spdff@data)
spdff@data <- assigninterp(df = spdff@data, vars = c("pop"), 
                           value = "ArealWeightedSum")
class(spdff@data$pop)

# perform interpolation of the covariates onto the
Zmat_xyt <- getZmat(formula = FORM.spatial, data = xyt2, 
                    pixelcovariates = spdff, 
                    cellwidth = CELLWIDTH_xyt, 
                    ext = EXT_xyt, overl = polyolay_xyt)

#Create the offset variable
FORM_pop_xyt <- X ~ pop - 1
Zmat_pop_xyt <- getZmat(formula = FORM_pop_xyt, data = x_typh,
                        pixelcovariates = spdff, cellwidth = CELLWIDTH_xyt,
                        ext = EXT_xyt, overl = polyolay_xyt)
mm_xyt <- length(attr(Zmat_pop_xyt, "mcens"))
nn_xyt <- length(attr(Zmat_pop_xyt, "ncens"))
poisson.offset_xyt <- spatialAtRisk(list(X = attr(Zmat_pop_xyt, "mcens"),
                                         Y = attr(Zmat_pop_xyt, "ncens"), Zm = matrix(Zmat_pop_xyt, 
                                       mm_xyt, nn_xyt)))

# plot the spatial interpolated covariates
plot(Zmat_xyt, ask = FALSE)
plot(Zmat_pop_xyt)

#create a dataframe for time
tdata <- dplyr::select(newd2.UTM3, "mean_temp", "mean_prec", 
                       "mean_dew", "mean_winds", "sin1", "sin2")
#tdata$t <- t3
tvec <- xyt2$tlim[1]:xyt2$tlim[2]
tdata$t <- 1:162

# choose last time point and number of proceeding time-points to include
#T <- 595
#LAGLENGTH <- 13
T <- 43 #Last time point is 42 months
LAGLENGTH <- 1 #Lag of 2 months/50 days

# bolt on the temporal covariates
ZmatList_xyt <- addTemporalCovariates(temporal.formula = TEMPORAL.FORMULA, T = 43, laglength = 6, 
                                      tdata = tdata, Zmat = Zmat_xyt)

# define the priors
lgprior_xyt <- PriorSpec(LogGaussianPrior(mean = log(c(1, 1000, 1)),
                                          variance = diag(1, 3)))
gprior_xyt <- PriorSpec(GaussianPrior(mean = rep(0.004, 13), variance = diag(1e+06, 13)))
priors_xyt <- lgcpPrior(etaprior = lgprior_xyt, betaprior = gprior_xyt)

# set initial values for the algorithm
INITS <- lgcpInits(etainit = log(c(sqrt(2.8), 1286, 0.5)), 
                   #betainit = NULL
                   betainit = rep(0.001, 13))

# choose the covariance function
CF_xyt <- CovFunction(exponentialCovFct)

# run the MCMC algorithm, note that as this takes a long time, I present a
# shorter run here 
DIRNAME <- getwd()
lg_xyt <- lgcpPredictSpatioTemporalPlusPars(formula = FORM_xyt, xyt = xyt2, T = 43, laglength = 6, 
                                            ZmatList = ZmatList_xyt, model.priors = priors_xyt, 
                                            model.inits = INITS, 
                                            spatial.covmodel = CF_xyt, 
                                            cellwidth = CELLWIDTH_xyt, mcmc.control = mcmcpars(mala.length = 1000, burnin = 100, 
                                      retain = 9, adaptivescheme = andrieuthomsh(inith = 2, alpha = 0.7, C = 2, 
                                        targetacceptance = 0.574)), 
                                            output.control = setoutput(gridfunction = dump2dir(dirname = file.path(DIRNAME, 
                                         "strataa_xyt22"), forceSave = TRUE)), ext = EXT_xyt)

# plot the log target
plot(ltar(lg_xyt), type = "s", xlab = "Iteration/900", ylab = "log target", ask = FALSE)

# produce traceplots of beta and eta
traceplots(lg_xyt, ask = FALSE)

# produce autocorrelation plots of beta and eta
parautocorr(lg_xyt, ask = FALSE)

# a summary table of beta and eta
parsum <- parsummary(lg_xyt)
# the above converted to LaTeX
library("miscFuncs")
latextable(parsum, rownames = rownames(parsum), colnames = c("Parameter", colnames(parsum)), 
           digits = 4)

# a text summary of the model parameters
textsummary(lg_xyt, digits = 4)

# a plot of the prior and posterior densities
priorpost(lg, ask = FALSE)

# the posterior covariance function
postcov(lg, ask = FALSE)

#############################################################################
#- Aggregate at quartery level
#calculate the difference in days
AB$t <- difftime(AB$date, AB$startdate, units = c("days"))
AB2Q<- AB[complete.cases(AB), ]
AB2Q$tQ <- ifelse(AB2$MONTH.x <= 3, 1,
                  ifelse(AB2$MONTH.x <= 6, 2,
                         ifelse(AB2$MONTH.x <= 9, 3,
                                ifelse(AB2$MONTH.x <= 12, 4, NA))))

xQ <- AB2$Longitude; yQ <- AB2$Latitude; tQ <- AB2$tQ
stdataQ <- cbind(xQ,yQ,tQ)
stdata1Q <- stdataQ[order(stdataQ[,ncol(stdataQ)]),] #sort the data
tlim1Q <- as.integer(c(1,4))

library(lgcp)
xytfQ <- stppp(list(data = stdata1Q, tlim = tlim1Q, window = poly.owin))
xytfQ

# perform polygon overlay operations and compute computational grid
polyolay_xytQ <- getpolyol(data = xytfQ, pixelcovariates = spdff, 
                           cellwidth = CELLWIDTH_xyt, 
                           ext = EXT_xyt)

# the statistical model for the main effets, $\beta$
FORM_xytQ <- X ~ dist + elev + wash + amt + amp + mean_temp + mean_prec + quarter
FORM.spatialQ <- X ~ dist + elev  + wash + amt + amp 
TEMPORAL.FORMULAQ <- t ~ quarter + mean_temp + mean_prec - 1

# perform interpolation of the covariates onto the
Zmat_xytQ <- getZmat(formula = FORM.spatialQ, data = xytfQ, 
                     pixelcovariates = spdff, 
                     cellwidth = CELLWIDTH_xyt, 
                     ext = 2, overl = polyolay_xyt)

#create a dataframe for time
tdataQ <- dplyr::select(AB2Q, "tQ", "mean_temp", "mean_prec")
names(tdataQ)[names(tdataQ) == "tQ"] <- "t"

tdata2Q <- tdataQ %>% group_by(t) %>% summarise_all(funs(mean))

tdata2Q$quarter <- tdata2Q$t

tdata3Q <- dplyr::select(tdata2Q, "t", "quarter", "mean_temp", "mean_prec")

# choose last time point and number of proceeding time-points to include
TQ <- 4
LAGLENGTHQ <- 3

# bolt on the temporal covariates
ZmatList_xytQ <- addTemporalCovariates(temporal.formula = TEMPORAL.FORMULAQ, 
                                       T = TQ, laglength = LAGLENGTHQ, 
                                       tdata = tdata3Q, Zmat = Zmat_xytQ)
#define the priors
lgprior_xytQ <- PriorSpec(LogGaussianPrior(mean = log(c(1, 1, 1)),
                                           variance = diag(1, 3)))
gprior_xytQ <- PriorSpec(GaussianPrior(mean = rep(1*10^-15, 9), 
                                       variance = diag(0.01, 9)))
priors_xytQ <- lgcpPrior(etaprior = lgprior_xytQ, betaprior = gprior_xytQ)

#Select initial values for algorithm
#INITSQ<- lgcpInits(etainit = log(c(sqrt(2.8), 1286, 0.5)), betainit = NULL)

# choose the covariance function
CF_xytQ <- CovFunction(exponentialCovFct)

# run the MCMC algorithm, note that as this takes a long time, we present a
DIRNAME <- getwd()

lg_xytQ <- lgcpPredictSpatioTemporalPlusPars(formula = FORM_xytQ, xyt = xytfQ, T = TQ, 
                                             laglength = LAGLENGTHQ, 
                                             ZmatList = ZmatList_xytQ, model.priors = priors_xytQ, 
                                             spatial.covmodel = CF_xytQ, 
                                             cellwidth = CELLWIDTH_xyt, 
                                             mcmc.control = mcmcpars(mala.length = 1000, burnin = 100, 
                                             retain = 9, adaptivescheme = andrieuthomsh(inith = 1, alpha = 0.5, C = 1, 
                                              targetacceptance = 0.574)), 
                                             output.control = setoutput(gridfunction = dump2dir(dirname = file.path(DIRNAME, 
                                        "strataa_xyt22Q"), forceSave = TRUE)), ext = EXT_xyt)

parsum <- parsummary(lg_xytQ)
textsummary(lg_xytQ, digits = 4)
##CODE ABOVE WORKS

############################################################################################################
#----------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#- Aggregate at quartery level
#calculate the difference in days
AB$t <- difftime(AB$date, AB$startdate, units = c("days"))
AB3<- AB[complete.cases(AB), ]
AB3$tQ <- ifelse(AB2$MONTH.x <= 3, 1,
                 ifelse(AB2$MONTH.x <= 6, 2,
                        ifelse(AB2$MONTH.x <= 9, 3,
                               ifelse(AB2$MONTH.x <= 12, 4, NA))))

xQ3 <- AB2$Longitude; yQ3 <- AB2$Latitude; tQ3 <- AB2$tQ
stdataQ3 <- cbind(xQ,yQ,tQ)
stdata1Q3 <- stdataQ[order(stdataQ[,ncol(stdataQ)]),] #sort the data
tlim1Q3 <- as.integer(c(1,4))

library(lgcp)
xytQ3 <- stppp(list(data = stdata1Q3, tlim = tlim1Q3, window = poly.owin))
xytQ3

# perform polygon overlay operations and compute computational grid
polyolay_xytQ3 <- getpolyol(data = xytQ3, pixelcovariates = spdff, 
                            cellwidth = CELLWIDTH_xyt, 
                            ext = EXT_xyt)

# the statistical model for the main effets
FORM_xytQ3 <- X ~ dist + elev + wash + amt + amp + 
  mean_temp + mean_prec + mean_dew + mean_winds +
  cos + sin + quarter
FORM.spatialQ3 <- X ~ dist + elev  + wash + amt + amp 
TEMPORAL.FORMULAQ3 <- t ~ mean_temp + mean_prec + mean_dew + mean_winds +
  cos + sin + quarter - 1

# perform interpolation of the covariates onto the
Zmat_xytQ3 <- getZmat(formula = FORM.spatialQ, data = xytQ3, 
                      pixelcovariates = spdff, 
                      cellwidth = CELLWIDTH_xyt, 
                      ext = 2, overl = polyolay_xytQ3)

#create a dataframe for time
tdataQ3 <- dplyr::select(AB2Q, "tQ", "mean_temp", "mean_prec", "mean_dew", "mean_winds")
names(tdataQ3)[names(tdataQ3) == "tQ"] <- "t"

tdata2Q3 <- tdataQ3 %>% group_by(t) %>% summarise_all(funs(mean))

tdata2Q3$quarter <- tdata2Q3$t
tdata2Q3$cos <- cos((2*pi*tdata2Q3$t)/4)
tdata2Q3$sin <- sin((2*pi*tdata2Q3$t)/4)

# choose last time point and number of proceeding time-points to include
TQ3 <- 4
LAGLENGTHQ3 <- 3

# bolt on the temporal covariates
ZmatList_xytQ3 <- addTemporalCovariates(temporal.formula = TEMPORAL.FORMULAQ3, 
                                        T = TQ3, laglength = LAGLENGTHQ3, 
                                        tdata = tdata2Q3, Zmat = Zmat_xytQ3)
#define the priors
lgprior_xytQ3 <- PriorSpec(LogGaussianPrior(mean = log(c(1, 140, 1)),
                                            variance = diag(0.1, 3)))
gprior_xytQ3 <- PriorSpec(GaussianPrior(mean = rep(1*10^-150, 13), 
                                        variance = diag(0.5, 13)))
priors_xytQ3 <- lgcpPrior(etaprior = lgprior_xytQ3, betaprior = gprior_xytQ3)

#betainit <- c("intercept" = 0, "dist" = 0, "elev" = 0, "wash" = 0, "amt" = 0, "amp" = 0,
#"mean_temp" = 0, "mean_prec" = 0, "mean_dew"= 0, "mean_winds" = 0, 
#"cos" = 0, "sin" = 0, "quarter" = 0)
#INITSQ3 <- lgcpInits(etainit = log(c(sqrt(2.8), 1286, 0.5)), betainit = rep(-0.000001,13))

# choose the covariance function
CF_xytQ3 <- CovFunction(exponentialCovFct)

# run the MCMC algorithm
DIRNAME <- getwd()
lg_xytQ3 <- lgcpPredictSpatioTemporalPlusPars(formula = FORM_xytQ3, xyt = xytQ3, T = TQ3, 
                                              laglength = LAGLENGTHQ3,
                                              #model.inits = INITSQ3,
                                              ZmatList = ZmatList_xytQ3, model.priors = priors_xytQ3, 
                                              spatial.covmodel = CF_xytQ3, cellwidth = CELLWIDTH_xyt, 
                                              mcmc.control = mcmcpars(mala.length = 1000, burnin = 100, 
                                               retain = 9, adaptivescheme = andrieuthomsh(inith = 1, alpha = 0.5, C = 1, 
                                               targetacceptance = 0.574)), 
                                              output.control = setoutput(gridfunction = dump2dir(dirname = file.path(DIRNAME, 
                                             "strataa_xytQ3"), forceSave = TRUE)), ext = EXT_xyt)

parsum <- parsummary(lg_xytQ3)
textsummary(lg_xytQ3, digits = 4)
#CODE ABOVE DOES NOT WORK

#--------------------------------------------------------------------------------
#- Aggregate at quartery level
#calculate the difference in days
AB$t <- difftime(AB$date, AB$startdate, units = c("days"))
AB4Q<- AB[complete.cases(AB), ]
AB4Q$tQ <- ifelse(AB2$MONTH.x <= 3, 1,
                  ifelse(AB2$MONTH.x <= 6, 2,
                         ifelse(AB2$MONTH.x <= 9, 3,
                                ifelse(AB2$MONTH.x <= 12, 4, NA))))

x4Q <- AB2$Longitude; y4Q <- AB2$Latitude; t4Q <- AB2$tQ
stdata4Q <- cbind(x4Q,y4Q,t4Q)
stdata14Q <- stdata4Q[order(stdata4Q[,ncol(stdata4Q)]),] #sort the data
tlim14Q <- as.integer(c(1,4))

library(lgcp)
xytf4Q <- stppp(list(data = stdata14Q, tlim = tlim14Q, window = poly.owin))
xytf4Q

# perform polygon overlay operations and compute computational grid
polyolay_xyt4Q <- getpolyol(data = xytf4Q, pixelcovariates = spdff, 
                            cellwidth = CELLWIDTH_xyt, 
                            ext = EXT_xyt)

# the statistical model for the main effets, $\beta$
FORM_xyt4Q <- X ~ dist + elev + wash +
  mean_temp + mean_prec + quarter 
FORM.spatial4Q <- X ~ dist + elev  + wash 
TEMPORAL.FORMULA4Q <- t ~ quarter + mean_temp + mean_prec - 1

# perform interpolation of the covariates onto the
Zmat_xyt4Q <- getZmat(formula = FORM.spatial4Q, data = xytf4Q, 
                      pixelcovariates = spdff, 
                      cellwidth = CELLWIDTH_xyt, 
                      ext = EXT_xyt, overl = polyolay_xyt4Q)

#create a dataframe for time
tdata4Q <- dplyr::select(AB2Q, "tQ", "mean_temp", "mean_prec")
names(tdata4Q)[names(tdata4Q) == "tQ"] <- "t"

tdata24Q <- tdata4Q %>% group_by(t) %>% summarise_all(funs(mean))

tdata24Q$quarter <- tdata2Q$t

tdata34Q <- dplyr::select(tdata24Q, "t", "quarter", "mean_temp", "mean_prec")
tdata34Q$cos <- cos((2*pi*tdata34Q$t) /4)
tdata34Q$sin <- sin((2*pi*tdata34Q$t) /4)

# choose last time point and number of proceeding time-points to include
T4Q <- 4
LAGLENGTH4Q <- 3

# bolt on the temporal covariates
ZmatList_xyt4Q <- addTemporalCovariates(temporal.formula = TEMPORAL.FORMULA4Q, 
                                        T = T4Q, laglength = LAGLENGTH4Q, 
                                        tdata = tdata34Q, Zmat = Zmat_xyt4Q)
#define the priors
lgprior_xyt4Q <- PriorSpec(LogGaussianPrior(mean = log(c(1, 1, 1)),
                                            variance = diag(1, 3)))
gprior_xyt4Q <- PriorSpec(GaussianPrior(mean = rep(3.0*10 *-17, 7), 
                                        variance = diag(0.01, 7)))
priors_xyt4Q <- lgcpPrior(etaprior = lgprior_xyt4Q, betaprior = gprior_xyt4Q)

#Select initial values for algorithm
#INITSQ<- lgcpInits(etainit = log(c(sqrt(2.8), 1286, 0.5)), betainit = NULL)

# choose the covariance function
CF_xyt4Q <- CovFunction(exponentialCovFct)

# run the MCMC algorithm, note that as this takes a long time, we present a
DIRNAME <- getwd()

lg_xyt4Q <- lgcpPredictSpatioTemporalPlusPars(formula = FORM_xyt4Q, xyt = xytf4Q, T = T4Q, 
                                              laglength = LAGLENGTH4Q, 
                                              ZmatList = ZmatList_xyt4Q, model.priors = priors_xyt4Q, 
                                              spatial.covmodel = CF_xyt4Q, 
                                              cellwidth = CELLWIDTH_xyt, 
                                              mcmc.control = mcmcpars(mala.length = 1000, burnin = 100, 
                                               retain = 9, adaptivescheme = andrieuthomsh(inith = 1, alpha = 0.5, C = 1, 
                                                targetacceptance = 0.574)), 
                                              output.control = setoutput(gridfunction = dump2dir(dirname = file.path(DIRNAME, 
                                            "strataa_xyt22Q"), forceSave = TRUE)), ext = EXT_xyt)

parsum <- parsummary(lg_xyt4Q)
textsummary(lg_xyt4Q, digits = 4)
##CODE ABOVE DOES NOT WORK

#########################################################################################
################################################################################
#------------------------------------------------------------------------------
#- Try removing other variables (cos, sin ) CODE WORKS!!
#--------------------------------------------------------------------------------
#- Aggregate at quartery level
#calculate the difference in days
AB$t <- difftime(AB$date, AB$startdate, units = c("days"))
AB4Q<- AB[complete.cases(AB), ]
AB4Q$tQ <- ifelse(AB2$MONTH.x <= 3, 1,
                  ifelse(AB2$MONTH.x <= 6, 2,
                         ifelse(AB2$MONTH.x <= 9, 3,
                                ifelse(AB2$MONTH.x <= 12, 4, NA))))

x4Q <- AB2$Longitude; y4Q <- AB2$Latitude; t4Q <- AB2$tQ
stdata4Q <- cbind(x4Q,y4Q,t4Q)
stdata14Q <- stdata4Q[order(stdata4Q[,ncol(stdata4Q)]),] #sort the data
tlim14Q <- as.integer(c(1,4))

library(lgcp)
xytf4Q <- stppp(list(data = stdata14Q, tlim = tlim14Q, window = poly.owin))
xytf4Q

# perform polygon overlay operations and compute computational grid
polyolay_xyt4Q <- getpolyol(data = xytf4Q, pixelcovariates = spdff, 
                            cellwidth = CELLWIDTH_xyt, 
                            ext = EXT_xyt)

# the statistical model for the main effets, $\beta$
FORM_xyt4Q <- X ~ dist + elev + wash +
  mean_temp + mean_prec + quarter +cos + sin
FORM.spatial4Q <- X ~ dist + elev  + wash 
TEMPORAL.FORMULA4Q <- t ~ quarter + mean_temp + mean_prec + cos + sin- 1

# perform interpolation of the covariates onto the
Zmat_xyt4Q <- getZmat(formula = FORM.spatial4Q, data = xytf4Q, 
                      pixelcovariates = spdff, 
                      cellwidth = CELLWIDTH_xyt, 
                      ext = EXT_xyt, overl = polyolay_xyt4Q)

#create a dataframe for time
tdata4Q <- dplyr::select(AB2Q, "tQ", "mean_temp", "mean_prec")
names(tdata4Q)[names(tdata4Q) == "tQ"] <- "t"

tdata24Q <- tdata4Q %>% group_by(t) %>% summarise_all(funs(mean))

tdata24Q$quarter <- tdata2Q$t

tdata34Q <- dplyr::select(tdata24Q, "t", "quarter", "mean_temp", "mean_prec")
tdata34Q$cos <- cos((2*pi*tdata34Q$t) /4)
tdata34Q$sin <- sin((2*pi*tdata34Q$t) /4)

# choose last time point and number of proceeding time-points to include
T4Q <- 4
LAGLENGTH4Q <- 3

# bolt on the temporal covariates
ZmatList_xyt4Q <- addTemporalCovariates(temporal.formula = TEMPORAL.FORMULA4Q, 
                                        T = T4Q, laglength = LAGLENGTH4Q, 
                                        tdata = tdata34Q, Zmat = Zmat_xyt4Q)
#define the priors
lgprior4Q <- PriorSpec(LogGaussianPrior(mean = log(c(1, 20, 1) ) , 
                                        variance = diag(0.2, 3) ) )
gprior4Q <- PriorSpec(GaussianPrior(mean = rep(0, 9), variance = diag(0.1, 9)))
priors4Q<- lgcpPrior(etaprior = lgprior4Q, betaprior = gprior4Q)

#Select initial values for algorithm
test <- setNames(rep(-0.0001,9), 
                 c("(Intercept)", 
                   "dist", "elev", "wash", 
                   "quarter", "mean_temp", "mean_prec", "cos", "sin"))
INITS4Q <- lgcpInits(etainit = log(c(sqrt(0.1), 20, 1)), betainit = test)

# choose the covariance function
CF_xyt4Q <- CovFunction(exponentialCovFct)

# run the MCMC algorithm, note that as this takes a long time, we present a
DIRNAME <- getwd()

lg_xyt4Q <- lgcpPredictSpatioTemporalPlusPars(formula = FORM_xyt4Q, xyt = xytf4Q, T = T4Q, 
                                              laglength = LAGLENGTH4Q, 
                                              model.inits = INITS4Q,
                                              ZmatList = ZmatList_xyt4Q, model.priors = priors_xyt4Q, 
                                              spatial.covmodel = CF_xyt4Q, 
                                              cellwidth = CELLWIDTH_xyt, 
                                              mcmc.control = mcmcpars(mala.length = 1000, burnin = 100, 
                                              retain = 9, adaptivescheme = andrieuthomsh(inith = 1, alpha = 0.5, C = 1, 
                                              targetacceptance = 0.574)), 
                                              output.control = setoutput(gridfunction = dump2dir(dirname = file.path(DIRNAME, 
                                             "strataa_xyt22Q"), forceSave = TRUE)), ext = EXT_xyt)

parsum <- parsummary(lg_xyt4Q)
textsummary(lg_xyt4Q, digits = 4)
##CODE ABOVE WORKS


#########################################################################################
################################################################################
#------------------------------------------------------------------------------
#- Try including an offset - CODE WORKS WITH OFFSET FROM IMAGE
#--------------------------------------------------------------------------------
#- Aggregate at quartery level
#calculate the difference in days
AB$t <- difftime(AB$date, AB$startdate, units = c("days"))
AB5Q<- AB[complete.cases(AB), ]
AB5Q$tQ <- ifelse(AB2$MONTH.x <= 3, 1,
                  ifelse(AB2$MONTH.x <= 6, 2,
                         ifelse(AB2$MONTH.x <= 9, 3,
                                ifelse(AB2$MONTH.x <= 12, 4, NA))))

x5Q <- AB5Q$Longitude; y5Q <- AB5Q$Latitude; t5Q <- AB5Q$tQ
stdata5Q <- cbind(x5Q,y5Q,t5Q)
stdata15Q <- stdata4Q[order(stdata5Q[,ncol(stdata5Q)]),] #sort the data
tlim15Q <- as.integer(c(1,4))

library(lgcp)
xytf5Q <- stppp(list(data = stdata15Q, tlim = tlim15Q, window = poly.owin))
xytf5Q

# perform polygon overlay operations and compute computational grid
polyolay_xyt5Q <- getpolyol(data = xytf5Q, pixelcovariates = spdff, 
                            cellwidth = CELLWIDTH_xyt, 
                            ext = EXT_xyt)

# the statistical model for the main effets, $\beta$
FORM_xyt5Q <- X ~ dist + elev + wash +
  mean_temp + mean_prec + quarter 
FORM.spatial5Q <- X ~ dist + elev  + wash 
TEMPORAL.FORMULA5Q <- t ~ quarter + mean_temp + mean_prec - 1

# perform interpolation of the covariates onto the
Zmat_xyt5Q <- getZmat(formula = FORM.spatial5Q, data = xytf5Q, 
                      pixelcovariates = spdff, 
                      cellwidth = CELLWIDTH_xyt, 
                      ext = EXT_xyt, overl = polyolay_xyt5Q)

#create a dataframe for time
tdata5Q <- dplyr::select(AB5Q, "tQ", "mean_temp", "mean_prec")
names(tdata5Q)[names(tdata5Q) == "tQ"] <- "t"

tdata25Q <- tdata5Q %>% group_by(t) %>% summarise_all(funs(mean))

tdata25Q$quarter <- tdata25Q$t

tdata35Q <- dplyr::select(tdata25Q, "t", "quarter", "mean_temp", "mean_prec")

# choose last time point and number of proceeding time-points to include
T5Q <- 4
LAGLENGTH5Q <- 3

# bolt on the temporal covariates
ZmatList_xyt5Q <- addTemporalCovariates(temporal.formula = TEMPORAL.FORMULA5Q, 
                                        T = T5Q, laglength = LAGLENGTH5Q, 
                                        tdata = tdata35Q, Zmat = Zmat_xyt5Q)
#define the priors
lgprior5Q <- PriorSpec(LogGaussianPrior(mean = log(c(1, 20, 1) ) , 
                                        variance = diag(0.2, 3) ) )
gprior5Q <- PriorSpec(GaussianPrior(mean = rep(0, 7), variance = diag(0.1, 7)))
priors5Q<- lgcpPrior(etaprior = lgprior5Q, betaprior = gprior5Q)

#Select initial values for algorithm
test <- setNames(rep(-0.0001,9), 
                 c("(Intercept)", 
                   "dist", "elev", "wash", 
                   "quarter", "mean_temp", "mean_prec", "cos", "sin"))
INITS5Q <- lgcpInits(etainit = log(c(sqrt(0.1), 20, 1)), betainit = test)

# choose the covariance function
CF_xyt5Q <- CovFunction(exponentialCovFct)

# run the MCMC algorithm, note that as this takes a long time, we present a
DIRNAME <- getwd()

#Create the offset variable
#FORM_pop_5Q <- X ~ pop - 1
#Zmat_pop_5Q <- getZmat(formula = FORM_pop_5Q, data = xytf5Q,
#                        pixelcovariates = spdff, cellwidth = CELLWIDTH_xyt,
#                         ext = EXT_xyt, overl = polyolay_xyt5Q)
#mm_xyt5Q <- length(attr(Zmat_pop_5Q, "mcens"))
#nn_xyt5Q<- length(attr(Zmat_pop_5Q, "ncens"))
#poisson.offset_5Q <- spatialAtRisk(list(X = attr(Zmat_pop_5Q, "mcens"),
#                                         Y = attr(Zmat_pop_5Q, "ncens"), Zm = matrix(Zmat_pop_5Q, 
#                                         mm_xyt5Q, nn_xyt5Q)))
#range(poisson.offset_5Q$Zm)


#create list of ofsets
#length(poisson.offset_5Q)
#offset.list5Q <- list()
#for (i in 1:4) {
# offset.list5Q[[i]] <- poisson.offset_5Q }

#length(offset.list5Q)
#Using the offset above gave the error below
#Error in glm.fit(x = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  : 
#NA/NaN/Inf in 'y'

#Opted to create an offset from the image of the census
offs <- spatialAtRisk(census_im)
offs2 <- list()
for(i in 1:4) {offs2[[i]] <- offs}

lg_xyt5Q <- lgcpPredictSpatioTemporalPlusPars(formula = FORM_xyt5Q, xyt = xytf5Q, T = T5Q, 
                                              laglength = LAGLENGTH5Q, 
                                              # model.inits = INITS5Q,
                                              ZmatList = ZmatList_xyt5Q, model.priors = priors5Q, 
                                              spatial.covmodel = CF_xyt5Q, 
                                              cellwidth = CELLWIDTH_xyt, poisson.offset = offs2,
                                              mcmc.control = mcmcpars(mala.length = 1000, burnin = 100, 
                                              retain = 9, adaptivescheme = andrieuthomsh(inith = 1, alpha = 0.5, C = 1, 
                                               targetacceptance = 0.574)), 
                                              output.control = setoutput(gridfunction = dump2dir(dirname = file.path(DIRNAME, 
                                             "strataa_xyt22Q"), forceSave = TRUE)), ext = EXT_xyt)

parsum <- parsummary(lg_xyt5Q)
textsummary(lg_xyt5Q, digits = 4)
##CODE ABOVE WORKS

##################################################################################################################
###Try unique time points
#- Aggregate at quartery level
#calculate the difference in days
AB$t <- difftime(AB$date, AB$startdate, units = c("days"))
AB6D<- AB[complete.cases(AB), ]
AB6D$td <- AB6D$t

x6D <- AB6D$Longitude; y6D <- AB6D$Latitude; t6D <- AB6D$td
stdata6D <- cbind(x6D,y6D,t6D)
stdata16D <- stdata6D[order(stdata6D[,ncol(stdata6D)]),] #sort the data
tlim16D <- as.integer(c(1,1240))

library(lgcp)
xytf6D <- stppp(list(data = stdata16D, tlim = tlim16D, window = poly.owin))
xytf6D

# perform polygon overlay operations and compute computational grid
polyolay_xyt6D <- getpolyol(data = xytf6D, pixelcovariates = spdff, 
                            cellwidth = CELLWIDTH_xyt, 
                            ext = EXT_xyt)

# the statistical model for the main effets, $\beta$
FORM_xyt6D <- X ~ dist + elev + wash +
  mean_temp + mean_prec  
FORM.spatial6D <- X ~ dist + elev  + wash 
TEMPORAL.FORMULA6D <- t ~  mean_temp + mean_prec - 1

# perform interpolation of the covariates onto the
Zmat_xyt6D <- getZmat(formula = FORM.spatial6D, data = xytf6D, 
                      pixelcovariates = spdff, 
                      cellwidth = CELLWIDTH_xyt, 
                      ext = EXT_xyt, overl = polyolay_xyt6D)

#create a dataframe for time
tdata6D <- dplyr::select(AB6D, "td", "mean_temp", "mean_prec")
names(tdata6D)[names(tdata6D) == "td"] <- "t"

# choose last time point and number of proceeding time-points to include
T6D <- 1230
LAGLENGTH6D <- 0

# bolt on the temporal covariates
#class(tdata26D$t)
#tD <- as.vector(tdata26D$t)
#class(tD) ; tdata26D$t <- tD; class(tdata26D$t)

ZmatList_xyt6D <- addTemporalCovariates(temporal.formula = TEMPORAL.FORMULA6D, 
                                        T = T6D, laglength = LAGLENGTH6D, 
                                        tdata = tdata6D, Zmat = Zmat_xyt6D)
#ERROR
#Error in modmat[i, ] : subscript out of bounds
#Only allowing laglength of 0

#define the priors
lgprior6D <- PriorSpec(LogGaussianPrior(mean = log(c(1, 20, 1) ) , 
                                        variance = diag(0.2, 3) ) )
gprior6D <- PriorSpec(GaussianPrior(mean = rep(0, 6), variance = diag(0.1, 6)))
priors6D<- lgcpPrior(etaprior = lgprior6D, betaprior = gprior6D)

# choose the covariance function
CF_xyt6D <- CovFunction(exponentialCovFct)

# run the MCMC algorithm, 
DIRNAME <- getwd()

#Opted to create an offset from the image of the census
offs6D <- spatialAtRisk(census_im)
offs26D <- list()
for(i in 1:2) {
  offs26D[[i]] <- offs6D}

betainit6D <- c("intercept" = 0.1, "dist" = 0.1, "elev" = 0.1, "wash" = 0.1,
                "mean_temp" = 0.1, "mean_prec" = 0.1)
INITS6D <- lgcpInits(etainit = log(c(sqrt(2.8), 1286, 0.5)), 
                     betainit = betainit6D )

lg_xyt6D <- lgcpPredictSpatioTemporalPlusPars(formula = FORM_xyt6D, xyt = xytf6D, 
                                              T = T6D, laglength = LAGLENGTH6D, 
                                              model.inits = INITS6D,
                                              ZmatList = ZmatList_xyt6D, model.priors = priors6D, 
                                              spatial.covmodel = CF_xyt6D, 
                                              cellwidth = CELLWIDTH_xyt, poisson.offset = offs26D,
                                              mcmc.control = mcmcpars(mala.length = 1000, burnin = 100, 
                            retain = 9, adaptivescheme = andrieuthomsh(inith = 1, alpha = 0.5, C = 1, 
                                  targetacceptance = 0.574)), 
                             output.control = setoutput(gridfunction = dump2dir(dirname = file.path(DIRNAME, 
                        "strataa_xyt6D"), forceSave = TRUE)), ext = EXT_xyt)

parsum <- parsummary(lg_xyt6D)
textsummary(lg_xyt6D, digits = 4)
##CODE ABOVE DOES NOT WORK

#-------------------------------------------------------------
