######################################
##                                  ##
##  Spatial Econometrics - Project  ##
##                                  ##
######################################


#rm(list=ls())

# Install (if necessary) and import packages
library(spdep)
library(rgdal)
library(maptools)
library(sp)
library(RColorBrewer)
library(classInt)
library(GISTools)
library(maps)
library(geosphere)
library(moments)
library(dplyr)
library(tseries)
library(epiDisplay)
library(tidyverse)
library(car)
library(lmtest)
library(MASS)
library(olsrr)
library(haven)
library(ggplot2)
library(ggpubr)

# Reading maps with rgdal
ward <- readOGR(".", "London_Ward")
ward_merged <- readOGR(".", "London_Ward_CityMerged")
borough <- readOGR(".", "London_Borough_Excluding_MHW")

# Changing projections
ward <- spTransform(ward, CRS("+proj=longlat +datum=NAD83"))
#degAxis(1)  
#degAxis(2)
ward_merged <- spTransform(ward_merged, CRS("+proj=longlat +datum=NAD83"))
#degAxis(1)
#degAxis(2)
borough <- spTransform(borough, CRS("+proj=longlat +datum=NAD83"))
#degAxis(1)
#degAxis(2)

dim(ward)
head(ward)
tail(ward)
dim(ward_merged)
head(ward_merged)
tail(ward_merged)
dim(borough)
head(borough)
tail(borough)

# Now let's plot the data
par(mar=c(1,1,1,1))  # to enlarge the plot
#par(mar=c(5,4,4,2))
plot(ward)
plot(ward_merged)
plot(borough)

# Save plot
#pdf('London wards vs London boroughs.pdf')
#par(mfrow=c(1,2))
#plot(ward_merged)
#plot(borough)
#dev.off()
#par(mfrow=c(1,1))


## Let's consider the single variables and plot a better visualisation

## WARD

# Let's check the type of our data and transform them into data.frame type
class(ward)
ward.df <- as.data.frame(ward)
head(ward.df)
#ward@bbox  # Min and Max values

# Let's now refine our visualization using "choropleth()"
var_ward <- ward.df$HECTARES
summary(var_ward)
head(var_ward)
dim(var_ward)
class(var_ward)

choropleth(ward, var_ward)
shades <- auto.shading(var_ward, n=6,
                       cols=add.alpha(brewer.pal(6,"Blues"),0.5))
choropleth(ward, var_ward, shading=shades)
#degAxis(1)
#degAxis(2)
choro.legend(0.2, 51.45, shades, between="to", fmt="%g", cex=0.65, bty="n")
#?choro.legend
#brks_ward <- c(90,130,180,240,350)
#size <- (brks_ward/800)*1.2
#legend(0.25, 51.45, legend=brks_ward, bg=brewer.pal(6,"Blues"),
#       cex=0.5, bty="n", pch=21)
#?legend
# It works but the legend's levels are too far from the plot, so theyso they do
# not appear in the plot (enlarged or downloaded)
#title(main="Hectares of London wards (not merged)")

# Save plot
#pdf('Hectares of London wards (not merged) + legend.pdf')
#shades <- auto.shading(var_ward, n=6,
#                       cols=add.alpha(brewer.pal(6,"Blues"),0.5))
#choropleth(ward, var_ward, shading=shades)
#choro.legend(0.19, 51.38, shades, between="to", fmt="%g", cex=0.65, bty="n")
#dev.off()


## WARD_MERGED

# Let's check the type of our data and transform them into data.frame type
class(ward_merged)
ward_merged.df <- as.data.frame(ward_merged)
head(ward_merged.df)
#ward_merged@bbox  # Min and Max values

# Let's now refine our visualization using "choropleth()"
var_ward_merged <- ward_merged.df$HECTARES
summary(var_ward_merged)
head(var_ward_merged)
dim(var_ward_merged)
class(var_ward_merged)

choropleth(ward_merged, var_ward_merged)
shades <- auto.shading(var_ward_merged, n=6,
                       cols=add.alpha(brewer.pal(6,"Blues"),0.5))
choropleth(ward_merged, var_ward_merged, shading=shades)
choro.legend(0.2, 51.45, shades, between="to", fmt="%g", cex=0.65, bty="n")
#title(main="Hectares of London wards (merged)")

# Save plot
#pdf('Hectares of London wards (merged) + legend.pdf')
#shades <- auto.shading(var_ward_merged, n=6,
#                       cols=add.alpha(brewer.pal(6,"Blues"),0.5))
#choropleth(ward_merged, var_ward_merged, shading=shades)
#choro.legend(0.19, 51.38, shades, between="to", fmt="%g", cex=0.65, bty="n")
#dev.off()


## Let's do the same analyses but considering only the boroughs (without wards)

## BOROUGH

# Let's check the type of our data and transform them into data.frame type
class(borough)
borough.df <- as.data.frame(borough)
head(borough.df)

# Let's now refine our visualization using "choropleth()"
var_borough <- borough.df$HECTARES
choropleth(borough, var_borough)
shades <- auto.shading(var_borough, n=6,
                       cols=add.alpha(brewer.pal(6,"Blues"),0.5))
choropleth(borough, var_borough, shading=shades)
choro.legend(0.2, 51.45, shades, between="to", fmt="%g", cex=0.65, bty="n")
#title(main="Hectares of London wards (excluding boroughs)")
#degAxis(1)
#degAxis(2)

# Save plot
#pdf('Hectares of London wards (excluding boroughs) + legend.pdf')
#shades <- auto.shading(var_borough, n=6,
#                       cols=add.alpha(brewer.pal(6,"Blues"),0.5))
#choropleth(borough, var_borough, shading=shades)
#choro.legend(0.16, 51.37, shades, between="to", fmt="%g", cex=0.65, bty="n")
#dev.off()

# According to all these plots it looks like the area (hectares) of boroughs
# increase with the increasing of the distance to the city centre.
# Let's check this hypothesis calculating, as first, the centroids of each
# borough

crds_borough <- coordinates(borough)
crds_borough

plot(borough, lwd=1)
points(crds_borough, pch=21, bg="red", cex=0.8)
#identify(crds_borough, labels=as.character(borough.df$NAME))
# Does not work properly
pointLabel(crds_borough, borough.df$NAME, cex=0.6)

# Ok, let's take "City of London" as a city centre reference

# Now let's calculate the distances between two coordinates. I will use a for
# loop to make it quicker
#install.packages('geosphere')
#library(geosphere)
crds_borough

dist <- rep(NA,33)
dist
for (i in 1:33){
    dist[i] <- distm(c(crds_borough[i,1], crds_borough[i,2]),
          c(crds_borough[33,1], crds_borough[33,2]), fun=distHaversine)
}
dist
# Distances are in meters

dist.df <- as.data.frame(dist)
borough.df.ext <- cbind(borough.df, dist.df)
#View(borough.df.ext)

# I would like to remove NA values from my dataframe for a better
# visualisation. It seems that the columns "SUB_2009" and "SUB_2006" have
# only NAs. Let's check
all(is.na(borough.df.ext$SUB_2009))
all(is.na(borough.df.ext$SUB_2006))
# Yes. Let's remove them
#data <- data[-c(4, 6:62)]
borough.df.ext <- borough.df.ext[-c(6:7)]
#View(borough.df.ext)
dim(borough.df.ext)
# Ok

# Finally I plot the data regarding the distances and the hectares, with the
# aim of verifying the initial intuition
par(mar=c(5,4,4,2))
plot(borough.df.ext$dist, borough.df.ext$HECTARES, xlab="Distance",
     ylab="Hectares", pch=16, col=4, cex=1)
#?plot
# As I thought, in general, the further away is the borough, the higher is
# its area (hectares)

# Save the plot
#pdf('Distance to the city centre vs Hectares of London boroughs.pdf')
#par(mar=c(5,4,4,2))
#plot(borough.df.ext$dist, borough.df.ext$HECTARES, xlab="Distance",
#     ylab="Hectares", pch=16, col=4, cex=1)
#dev.off()


## MODEL

# Mention about direct effect --> OLS

model1 <- lm(HECTARES ~ dist, data=borough.df.ext)
summary(model1)

# The intercept is not significant --> I remove it (with "- 1")
model2 <- lm(HECTARES ~ dist - 1, data=borough.df.ext)
summary(model2)
# Better: more significance against the null hypothesis and an higher R^2
# Comments about summary


## DIAGNOSTICS

# Reset test
resettest(model2, power=2:3, type=c("fitted"))
# H0 not rejected --> the model is correctly represented

# Normality of residuals
summary(model2$residuals)
skewness(model2$residuals)
kurtosis(model2$residuals)

# Histogram + normal curve
par(mar=c(5,4,4,2))
h <- hist(model2$residuals, col="blue", xlab="Residuals", freq=F, nclass=50)
xfit <- seq(min(model2$residuals), max(model2$residuals), length=50)
yfit <- dnorm(xfit, mean=mean(model2$residuals), sd=sd(model2$residuals))
lines(xfit, yfit, col="red", lwd=2)

# Q-Q Plot for standardized residuals
plot(model2, which=2)

res.std <- rstandard(model2)
boxplot(res.std)

jarque.bera.test(model2$residuals)
# H0 rejected --> residuals are not normally distributed

# Homoscedasticity
# Residuals vs fitted - we should have randomly located points
plot(model2, which=1)

# Breusch-Pagan test
bptest(model2, studentize=FALSE)
# Error: it needs intercept --> if you want try model1

# Akaike Information Criterion
#AIC(model1)
AIC(model2)


## Trying again using the robust variance-covariance matrix

# Robust regression
reg_robust <- rlm(HECTARES ~ dist - 1, data=borough.df.ext)
summary(reg_robust)

# Reset test
resettest(reg_robust, power=2:3, type=c("fitted"))

# Homoscedasticity
plot(reg_robust, which=1)
bptest(reg_robust, studentize=FALSE)

# Residuals analysis
summary(reg_robust$residuals)
skewness(reg_robust$residuals)
kurtosis(reg_robust$residuals)

h <- hist(reg_robust$residuals, col="blue", xlab="Residuals", freq=F,
              nclass=50)
xfit <- seq(min(reg_robust$residuals), max(reg_robust$residuals),
                length=50)
yfit <- dnorm(xfit, mean=mean(reg_robust$residuals),
              sd=sd(reg_robust$residuals))
lines(xfit, yfit, col="red", lwd=2)

plot(reg_robust, which=2)

res.std <- rstandard(reg_robust)
boxplot(res.std)

jarque.bera.test(reg_robust$residuals)
# H0 rejected again

# --> NO robust regression


## Try with logarithm

model3 <- lm(log(HECTARES) ~ dist, data=borough.df.ext)
summary(model3)

# Reset test
resettest(model3, power=2:3, type=c("fitted"))
# H0 rejected --> the model is not correctly represented

# Normality of residuals
summary(model3$residuals)
skewness(model3$residuals)
kurtosis(model3$residuals)

# Histogram + normal curve
par(mar=c(5,4,4,2))
h <- hist(model3$residuals, col="blue", xlab="Residuals", freq=F, nclass=50)
xfit <- seq(min(model3$residuals), max(model3$residuals), length=50)
yfit <- dnorm(xfit, mean=mean(model3$residuals), sd=sd(model3$residuals))
lines(xfit, yfit, col="red", lwd=2)

# Q-Q Plot for standardized residuals
plot(model3, which=2)

res.std <- rstandard(model3)
boxplot(res.std)

jarque.bera.test(model3$residuals)
# H0 not rejected --> Normality of the residuals

# Homoscedasticity
# Residuals vs fitted - we should have randomly located points
plot(model3, which=1)

# Breusch-Pagan test
bptest(model3, studentize=FALSE)
# H0 not rejected --> Homoscedasticity of the residuals

# Akaike Information Criterion
AIC(model3)
# Very low
