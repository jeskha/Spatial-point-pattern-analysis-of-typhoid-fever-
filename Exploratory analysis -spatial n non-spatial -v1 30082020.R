#-------------------------------------------------------------------
#- DESCRIPTIVE ANALYSIS

#- Non-Spatial expoloratory analysis
table(data$S..Typhi.NTS) #-get numbers for data

table(newd2$Gender)

class(newd2$Age)
summary(newd2$Age)

#- Summarize age by gender
ggplot(newd2, aes(x=Gender, y=Age)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=1,
               outlier.size=4, 
               fill = "white", colour = "#3366FF") + 
  # coord_flip() + 
  theme_bw() +
  labs(y = "Age (Years)", x = "") +
  theme(axis.text = element_text(size = 12, color= "black"),
        axis.title = element_text(size=14, colour = "black", face = "plain"),
        legend.text = element_text(size = 12),
        legend.title=element_blank()) +
   geom_jitter(width = 0.15) + 
  stat_summary(geom = "text", fun = quantile,
               aes(label=sprintf("%1.1f", ..y..) ),
               position = position_nudge((x=0.43)), size = 4.5)

#----------------------------------------------------------------------
#- Draw a map of the vases and the population census

library(ggmap)
library(tidyverse)

#plot Ndirande with points, census and health facility
head(newd2)
head(census_coords3)
hclon <- 35.0395908
hclat <- -15.7780157
ndirande_hc <- data.frame(hclat, hclon)

#Plot Ndirande typhoid fever cases
register_google(key = "AIzaSyC1edoiw3S_Gw0TI1iQ9VCyll91xJ8veMI")

ndirande_loc <- c(left = 35.030, bottom = -15.775, right = 35.050, top = -15.755)
ndirande_roadmap <- get_map(location = ndirande_loc,
                           zoom = 13, 
                           maptype = "roadmap",
                           source = "google")
ggmap(ndirande_roadmap)

ndirande_satellite <- get_map(location = ndirande_loc,
                            zoom = 12, 
                            maptype = "satellite",
                            source = "google")
ggmap(ndirande_satellite)

ndirande_loc2 <- c(left = 35.02, bottom = -15.80, right = 35.07, top = -15.74)
ndirande_terrain <- get_map(location = ndirande_loc2,
                              zoom = 13, 
                              maptype = "satellite", source = "google")
ggmap(ndirande_terrain)

x <- c(cases$x, contr$x)
y <- c(cases$y, contr$y)
cc <- c(rep("case", cases$n), rep("control", contr$n))

cc <- c(rep("Typhoid fever case", count(newd2)),
        rep("Population census", count(census_coords3)),
        rep("Ndirande health center", 1))

census_mapd <- data.frame(Longitude = census_coords3$h07_gpslon, Latitude = census_coords3$h07_gpslat,
                          Key = "Population census")
cases_mapd <- data.frame(Longitude = newd2$Longitude, Latitude = newd2$Latitude, Key = "Typhoid fever case")
hc_mapd <- data.frame(Longitude = 35.0395908, Latitude = -15.7780157, Key = "Ndirande health center" )

map_d <- rbind(census_mapd, cases_mapd, hc_mapd)

ggmap(ndirande_terrain) + geom_point(data = map_d, 
                                     aes(x = Longitude, y = Latitude, col = Key), 
                                     size = 2) +
  xlab("Longitude") + ylab("Latitude")

  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        rect = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank())
  
#-------------------------------------------------------------------
#- Two different maps for the cases and for the population census
#- TYPHOID CASES
cc <- c(rep("Typhoid fever case", count(newd2)),
          rep("Ndirande health center", 1))
  
census_mapd <- data.frame(Longitude = census_coords3$h07_gpslon, Latitude = census_coords3$h07_gpslat,
                            Key = "Population census")
cases_mapd <- data.frame(Longitude = newd2$Longitude, Latitude = newd2$Latitude, Key = "Typhoid fever case")
hc_mapd <- data.frame(Longitude = 35.0395908, Latitude = -15.7780157, Key = "Ndirande health center" )
  
map_d <- rbind(census_mapd, cases_mapd, hc_mapd)


map_d <- rbind(cases_mapd, hc_mapd)
map_d2 <- rbind(census_mapd, hc_mapd)

#par(mfrow=c(1,2))
ggmap(ndirande_terrain) + geom_point(data = map_d, 
                                       aes(x = Longitude, y = Latitude, col = Key), 
                                       size = 1) +
  scale_color_manual(values = c("Ndirande health center" = "black", "Typhoid fever case" = "red"))+
  labs(x = "Longitude", y = "Latitude") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_point(data = hc_mapd, aes(x=Longitude, y = Latitude), size = 5) +
  guides(colour = guide_legend(override.aes = list(size=5)))

#- census and health center
ggmap(ndirande_terrain) + geom_point(data = map_d2, 
                                     aes(x = Longitude, y = Latitude, col = Key), 
                                     size = 0.5) +
  scale_color_manual(values = c("Ndirande health center" = "black", "Population census" = "red"))+
  labs(x = "Longitude", y = "Latitude") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_point(data = hc_mapd, aes(x=Longitude, y = Latitude), size = 5) +
  guides(colour = guide_legend(override.aes = list(size=5)))

####--------------------------------------------------------
#- Homogeneous K-function
#- Read shapefile
poly <- rgdal::readOGR(dsn=".", "ndirandepoly")
plot(poly)

library(maptools)
poly.owin <- as.owin(poly)
x_typh <- ppp(x = data.UTM.df$Longitude, y = data.UTM.df$Latitude, window = poly.owin)
x_typh <- as.ppp(x_typh) 

typh_kest <- Kest(x_typh)
plot(typh_kest)
plot(typh_kest, iso~r)
plot(typh_kest, cbind(iso, theo) ~ r, xlab= "r (metres)")
plot(typh_kest, cbind(iso, trans, theo) ~ r)

#plot the K-function with an envelope
E <- envelope(x_typh, Kest, nsim=162, fix.n=TRUE)
plot(E, xlab="r (metres)", main= "K-function for typhoid fever")

#use more formal test
mad.test(x_typh, Lest, nsim=162, rmax=2, use.theo=TRUE)
dclf.test(x_typh, Lest, nsim=162, rmax=2, use.theo=TRUE)$p.value

#--------------------------------------------------------------------
#- Inhomogeneous K function
  
# (1) intensity function estimated by model-fitting
fit <- ppm(x_typh, ~ offset(log(census_im)), Poisson())
# (a) predict intensity values at points themselves,
#     obtaining a vector of lambda values
lambda <- predict(fit, locations=x_typh, type="trend")
# inhomogeneous K function
Ki <- Kinhom(x_typh, lambda)
plot(Ki, main = "Inhomogeneous K-function",
     xlab = "r (meters)")
Kenv <- envelope(x_typh, Kest, nrank=2, nsim=99)

#Inhomogeneous K function
mydense<-density(x_typh); plot(mydense)  
Ki <- Kinhom(x_typh, mydense)
plot(Ki, main = "Inhomogeneous K function")
Kenv2<-envelope(x_typh, Kinhom, simulate = expression(rpoispp(mydense)),nrank=2, nsim=99) # 95 C.I
plot(Kenv2)

fit <- ppm(x_typh ~ offset(log(census_im)))
lam <- predict(fit, locations = x_typh)
Ki <- Kinhom(x_typh, lam, correction = c("border", "bord.modif"))
plot(Ki, main = "Inhomogeneous K function", xlab = "r (meters)")

#-----------------------------------------------------------------
#- Inhomogeneous K-function with confidence intervals
#Inhomogeneous K function
mydense <- density(x_census); plot(mydense)
Ki <- Kinhom(x_typh, mydense)
plot(Ki, main = "Inhomogeneous K function")
Kenv2 <- envelope(x_typh, Kinhom, 
                  simulate = expression(rpoispp(mydense)),
                  nsim=999, correction="none") # 95 C.I
par(mfrow=c(1,2))
par(mar=c(5.1, 4.1, 4.1, 2.1))
plot(Kenv2, xlab = "r (one unit = 1 metre)", ylab = "K(r)",
     main = "Inhomogeneous K-function")

plot(Kenv2, xlab = "r (one unit = 1 metre)", ylab = "K(r)",
     main = "b",
     xlim=c(0,400), xaxt = 'n', ylim=c(0,600000))
axis(side=1, at = seq(100, 400, by = 100))

#zoom in on the smaller distances
plot(Kenv2, xlab = "r (one unit = 1 metre)", ylab = "K(r)",
     main = "c",
     xlim=c(0,200))


plot(Kenv2, xlab = "r (one unit = 1 metre)", ylab = "K(r)",
     main = "Inhomogeneous K function",
     xlim=c(200,400))


#--------------------------------------
plot(envelope(x_typh, Kest, nsim=39))


