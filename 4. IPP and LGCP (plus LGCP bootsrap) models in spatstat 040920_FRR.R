#-------------------------------------------------------------------------------
#- Inhomogeneous Poisson process and LGCP models
#-----------------------------------------------------------------------

#- Create planar point pattern (PPP) objects for the typhoid cases
# and census data

library(maptools)
poly.owin <- as.owin(poly)
x_typh <- ppp(x = data.UTM.df$Longitude, y = data.UTM.df$Latitude, window = poly.owin)
x_typh <- as.ppp(x_typh) #to remove points falling outside window

#- ##create pixel image for the population density
#- Spatstat does not allow a planar point pattern to be used as covariate
x_census <- ppp(census.UTM.df$h07_gpslon, y=census.UTM.df$h07_gpslat,
                window = poly.owin) #452 points rejected
x_census <- as.ppp(x_typh) #to remove points falling outside window
plot(x_census)

#convert the ppp object to a pixel
census_pop <- as.im(x_census)

#Some areas had zero population
#This would be problematic when specifying the offset
#work with as.im image and replace low values
census_im <- census_pop
census_im[census_im==0] <- 0.5*min(census_im[census_im!=0])
image(census_im)
image(log(census_im))
ppm(x_typh ~ offset(log(census_im)))

#------------------------------------------------------------------------------
#- Investigate effect of marks on occurrence of typhoid
#- AIC not used for marks because it does not apply to marks (Baddeley 2015)
#- ppm can handle categorical marks
#Gender
newd2$Gender <- as.factor(newd2$Gender)
x_gender <- ppp(x=data.UTM.df$Longitude, y=data.UTM.df$Latitude, marks=newd2$Gender,
                window = poly.owin)

m1 <- ppm(x_gender ~ offset(log(census_im)) +marks)
auc(m1) #87% being explained by coordinates and sex
roc(m1)

#Categorize age (ppm is not yet implemented for marked point patterns)
newd2$age_cat <- cut(newd2$Age, 
                     breaks=c(-Inf, 5, 17, Inf), 
                     labels=c("0-5","6-17","18+"))
x_age2 <- ppp(x=data.UTM.df$Longitude, y=data.UTM.df$Latitude,  marks=newd2$age,
              window = poly.owin)
m2 <- ppm(x_age2 ~ offset(log(census_im)) + marks )
auc(m2) #91% of variation explained by age and coordinates
roc(m2)

#Combine age and gender into one variable/mark
#Gender age
newd2$gender.age <- paste(newd2$Gender, newd2$age_cat, sep = "")
class(newd2$gender.age)
newd2$gender.age <- as.factor(newd2$gender.age)
class(newd2$gender.age)
x_g.age <- ppp(x=data.UTM.df$Longitude, y=data.UTM.df$Latitude, 
               marks=newd2$gender.age, window = poly.owin)
ppm(x_g.age ~ offset(log(census_im)) + marks )

#----------------------------------------------------
#- Environmental covariates

#convert to images for ppm analyses
#- Spatstat can not handle raster covariates
#- Convert the rasters to images 
elev.im <- as.im(elevation.f); plot(elev.im)
dist.im <- as.im(distancesf); plot(dist.im)
wash.im <- as.im(wash.score.r); plot(wash.im)

#-----------------------------------------------------------------------
#- Fit the IPP model with marks
ppm.typh.marks_null <- ppm(x_g.age ~ 1 )
AIC(ppm.typh.marks_null) #4373.956

ppm.typh.marks_null_offset <- ppm(x_g.age ~offset(log(census_im)) )
AIC(ppm.typh.marks_null_offset) #4158.503
#- AIC lower than model without offset
#- model with offset better

ppm.typh.marks_a <- ppm(x_g.age, ~ offset(log(census_im)) + marks +
                          elev.im  ,
                        covariates = list(dist = elev.im))
AIC(ppm.typh.marks_a) #4103.448
#- Although not statistically significant, elevation
#- contributes affects the outcome (AIC lower than null model)

ppm.typh.marks_b <- ppm(x_g.age, ~ offset(log(census_im)) + marks +
                          dist.im +elev.im ,
                        covariates = list(dist = dist.im, elev=elev.im))
AIC(ppm.typh.marks_b) #4083.304
#- AIC lower than that for model with elevation alone
#- Distance contributes affects the outcome

ppm.typh.marks_c <- ppm(x_g.age, ~ offset(log(census_im)) + marks +
                        dist.im + elev.im + wash.im ,
                      covariates = list(dist = dist.im, elev=elev.im,
                                        wash = wash.im))
AIC(ppm.typh.marks_c) #4018.662
#- AIC is the lowest amongst all the models
#- Include all three covariates and the offset in the models

ppm.typh.marks_c
summary(ppm.typh.marks_c)

#---------------------------------------------------
#- Fit IPP model with no marks to compare the results to the LGCP
# estimates. This is because (both spatstat and lgcp) packages
# can not handle marked point patterns at the moment
#- Fit the IPP model with marks
ppm.typh.nomarks_null <- ppm(x_typh ~ 1 )
AIC(ppm.typh.nomarks_null) #3800.593

ppm.typh.nomarks_null_offset <- ppm(x_typh ~offset(log(census_im)) )
AIC(ppm.typh.nomarks_null_offset) #3585.14
#- AIC lower than model without offset
#- model with offset better

ppm.typh.nomarks_a <- ppm(x_typh, ~ offset(log(census_im)) +
                          elev.im  ,
                        covariates = list(dist = elev.im))
AIC(ppm.typh.nomarks_a) #43574.446
#- Although not statistically significant, elevation
#- contributes affects the outcome (AIC lower than null model)
#- Elevation slightly improves the model

ppm.typh.nomarks_b <- ppm(x_typh, ~ offset(log(census_im)) +
                          dist.im +elev.im ,
                        covariates = list(dist = dist.im, elev=elev.im))
AIC(ppm.typh.nomarks_b) #43554.302
#- AIC lower than that for model with elevation alone
#- Distance contributes affects the outcome
#- Distance to HF improves the model

ppm.typh.nomarks_c <- ppm(x_typh, ~ offset(log(census_im))+
                          dist.im + elev.im + wash.im ,
                        covariates = list(dist = dist.im, elev=elev.im,
                                          wash = wash.im))
AIC(ppm.typh.nomarks_c) #3489.66
#- WASH score improves the model
#- AIC is the lowest amongst all the models
#- Include all three covariates and the offset in the models

ppm.typh.nomarks_c
summary(ppm.typh.nomarks_c)

valid.ppm(ppm.typh.nomarks_c) #model is well defined
valid.ppm(ppm.typh.nomarks_c)

ppm.typh.marks <- ppm(x_g.age, ~ offset(log(census_im)) + marks +
                          dist.im + elev.im + wash.im ,
                        covariates = list(dist = dist.im, elev=elev.im,
                                          wash = wash.im))

ppm.typh.nomarks <- ppm(x_typh, ~ offset(log(census_im))+
                            dist.im + elev.im + wash.im ,
                          covariates = list(dist = dist.im, elev=elev.im,
                                            wash = wash.im))

#-----------Model Diagnostics

#- Quantile-Quantile plot to check independence assumption
qqplot.ppm(ppm.typh.nomarks, nsim=99)
#independence assumption not met

#partial residuals
#Effect of all variables was linear on log scale
#Functional form of covariates
pres_elev <- parres(ppm.typh.nomarks, "elev.im")
plot(pres_elev)

pres_dist <- parres(ppm.typh.nomarks, "dist.im")
plot(pres_dist)

pres_wash <- parres(ppm.typh.nomarks, "wash.im")
plot(pres_wash)

#Diagnostics for the IPP model
res <- residuals(ppm.typh.marks)
res
plot(res)

res <- residuals(ppm.typh.nomarks)
res
plot(res)
#Total residual is zero up to a numerical error

#Four-panel diagnosis plot (recommended by Baddeley)
diagnose.ppm(ppm.typh.nomarks)

#Zoom in on 4th graph from diagnose.ppm
plot((residuals(ppm.typh.nomarks, type="pearson")))
plot(Smooth(residuals(ppm.typh.nomarks)))
plot(Smooth(residuals(ppm.typh.nomarks, type="pearson")),
     main="Smooth Pearson residuals for IPP model")#clear trend
axis(1); axis(2); box()

#-
#Smooth raw residuals for IPP
diagnose.ppm(ppm.typh.nomarks, which="smooth",
             main="Smoothed raw residuals for IPP model")
axis(1); axis(2); box()

pres2e_ppm <- residuals(ppm.typh.nomarks, type="pearson")
psmo2e_ppm <- Smooth(pres2e_ppm)
1/(2 * sqrt(pi) * attr(psmo2e_ppm, "sigma")) #0.0006661422
2 * 1/(2 * sqrt(pi) * attr(psmo2e_ppm, "sigma")) #0.0006661422

#------------------------------------------------------------
#- Log-Gaussian Cox Process Model
#- LGCP model can not handle marks
#- I fitted the LGCP model without marks

#- LGCP model fitted using composite likelihood methods

#- Methods for selecting and diagnosing LGCP models in spatat
# are limited
#- I therefore did not use AIC here

#- Selected models based on coefficients and 95% CIs

kppm.typh_null <- kppm(x_typh ~ 1,
                  model = "exp", "LGCP", method="clik2")
kppm.typh_null

kppm.typh_null_offset <- kppm(x_typh ~ offset(log(census_im)),
                       model = "exp", "LGCP", method="clik2")
kppm.typh_null_offset

kppm.typh_a <- kppm(x_typh ~ offset(log(census_im))+
                           elev, covariates=list(elev=elev.im),
                              model = "exp", "LGCP", method="clik2")
kppm.typh_a
summary(kppm.typh_a) 
#elevation statistically significant in unadjusted model

kppm.typh_b <- kppm(x_typh ~ offset(log(census_im))+
                      elev + dist, 
                    covariates=list(elev=elev.im, dist=dist.im),
                    model = "exp", "LGCP", method="clik2")
kppm.typh_b
summary(kppm.typh_b) 
#distance statistically significant
#elevation no longer significant
#I still chose to maintain the elevation for comparability with
#IPP model

kppm.typh_c <- kppm(x_typh ~ offset(log(census_im))+
                      elev + dist + wash, 
                    covariates=list(elev=elev.im, dist=dist.im,
                                    wash=wash.im),
                    model = "exp", "LGCP", method="clik2")
kppm.typh_c
summary(kppm.typh_c)
#- wash statistically significant
#- distance statistically significant
#- elevation no longer significant
#- Final model used all three covariates

kppm.typh <- kppm(x_typh, ~ offset(log(census_im)) +
                    dist + elev + wash ,
                  covariates = list(dist = dist.im, elev=elev.im,
                                    wash = wash.im),
                  model = "exp", "LGCP", method="clik2")
kppm.typh
summary(kppm.typh)
#unlist(parameters(kppm.typh))

#------------------------
#- Carry out parametric bootstrap to get confidence intervals
# for the cluster parameters
#----------------------------------------------------
#- PARAMETRIC BOOTSTRAP

#- Bootstrapping of cluster parameters for LGCP Model

#- Fix the seed for reproducibility
set.seed(500)

#- Simulate data 1,000 models from the fitted LGCP model
model_sim <- simulate.kppm(kppm.typh, nsim=1000)

#- Initialize lists to store simulation results
kppm_sim <- list()
para1 <- list()
para2 <- list()

for (i in 1:length(model_sim)) {
  kppm_sim[[i]] <- kppm(model_sim[[i]] ~ offset(log(census_im)) +
                          dist + elev + wash ,
                        covariates = list(dist = dist.im, elev=elev.im,
                                          wash = wash.im),
                        model = "exp", "LGCP",
                        method="clik2")  
  para1[i] <- kppm_sim[[i]]$clustpar[1]
  para2[i] <- kppm_sim[[i]]$clustpar[2]
}

para1f = as.numeric(as.character(unlist(para1)))
para2f = as.numeric(as.character(unlist(para2)))

#- Calculate confidence interval using percentile method
#calculate CI
quantile(para1f, prob=0.025)
quantile(para1f, prob=0.975)
quantile(para2f, prob=0.025)
quantile(para2f, prob=0.975)

#---------------------------------------------------------

#- model diagnostics for LGCP model
#- As previously mentioned, the methods are very limited
envLT <- envelope(kppm.typh, Lest, nsim=39)
plot(envLT)

summary(kppm.typh)

#Perform a prediction using the LGCP model
kppm.pred <- predict(kppm.typh)

#Find the maximum predicted intensity per 100,000 population
min(kppm.pred)
max(kppm.pred) * 100000

#Plot the predicted intensity
plot(kppm.pred, main = "Predicted intensity of typhoid fever"); axis(1);axis(2); box()
plot(kppm.pred)

image(kppm.pred)
points(x_typh)

#Model diagnostics for the LGCP model
#Diagnostics
res <- residuals() #close to zero which is good
plot(Smooth(res))

#-----MODEL DIAGNOSTICS FOR LGCP MODELS
#- The diagnose.ppm command does not work for the LGCP model
#- There is no similar function for LGCP model
#- That is why the design of the smoothed residuals for the IPP
# and LGCP models is different

#Residuals for LGCP model
res <- residuals(kppm.typh)
res
plot(res)

#diagnose.ppm not available for kppm models
#Plot similar plot using raw residuals
plot(Smooth(residuals(kppm.typh, type="raw")),
     main = "Smoothed raw residuals for LGCP model")
axis(1); axis(2); box()

pres2e_kppm <- residuals(kppm.typh, type="pearson")
psmo2e_kppm <- Smooth(pres2e_kppm)
1/(2 * sqrt(pi) * attr(psmo2e_kppm, "sigma")) 0.0006661422

plot(Smooth(residuals(kppm.typh, type="Pearson")),
     main = "Smoothed Pearson residuals for LGCP model")
axis(1); axis(2); box()

#---------------------------------------------------------------------
