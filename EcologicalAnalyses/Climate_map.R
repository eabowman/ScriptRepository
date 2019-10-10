## Script adapted from script created by Yu-Ling Huang
## modified by Liz Bowman (eabowman on github) Oct. 1, 2019

#install.packages(c('raster','rgdal','scales','maptools','dplyr','tidyr'))
library(raster)
library(rgdal)
library(scales)
library(maptools)
library(dplyr)
library(tidyr)

sites <- read.csv('data/20190923_sampling_data.csv', as.is = T)

# isolate site coordinates
sites %>%
  select(site, range, lat, long) %>%
  distinct(site, .keep_all = T) -> sites
sites$range <- factor(sites$range)

#get USA map from GADM, level=0 no state border, level=1 state border, level=2 conty border...
#US_0<-getData('GADM', country="USA", level=0)
US_1<-getData('GADM', country="USA", level=1)
#US_2<-getData('GADM', country="USA", level=2)
#US_3<-getData('GADM', country="USA", level=3)
#US_4<-getData('GADM', country="USA", level=4)

#check the data class
class(US_1)
#check data fram names
names(US_1)
#list the all Name of states
unique(US_1$NAME_1)
#get the Arizona state border map from the whole USA map
Arizona<-subset(US_1, NAME_1=="Arizona")
plot(Arizona)

#read the downloaded WorldClm data (tile12) by raster
AnnPrec<-raster('~/../../Volumes/Cenococcum/PhD/worldclim/prec_12/wc2.0_30s_prec_12.tif')
AnnMeanTemp<-raster('~/../../Volumes/Cenococcum/PhD/worldclim/tavg_12/wc2.0_30s_tavg_12.tif')
plot(AnnPrec)

# or can directly download the whole bio data using the function below:
#var can be bio1 to bio12, bio is the whole climate variables data set
#res = resolutions of the data
#lon, lat = specify a point then it will download the data tile in that area
Bio<-getData("worldclim", var="bio", res=0.5, lon=-110, lat=33)
#crop the cilmate data with specific region, in this case it's Arizona map we get from US_1
Bio.sub<-crop(Bio, extent(Arizona))
#mask function make sure the climate data only within the Arizona border
Bio.sub.mask<-mask(Bio.sub, Arizona)


# << Annual temperature map >> -----------
pdf('./figures/TemperatureMap.pdf', height = 7, width = 10)

#add all points
color.vec<-c("blue", "red","gray","yellow")
Temp.col<-colorRampPalette(c("blue","deepskyblue","cyan", "yellow","orange", "darkorange", "red"))
plot(Bio.sub.mask$bio1_12/10, main="",legend=FALSE,axes=FALSE, box=FALSE, col=Temp.col(9),
     xlab= "", ylab="")
plot(Arizona, add = T)
axis(1, pos=31, cex.axis=2, tck=-0.02)
axis(2, las=2, cex.axis=2, tck=-0.02)
# mtext(text = 'Longitude', side = 1, line = 4.5, cex = 3)
# mtext(text = 'Latitude', side = 2, line = 3, cex = 3)
points(sites$long, sites$lat, pch=19, cex=2)
# temperature legend
par(xpd=TRUE)
Bio.sub.mask.range<-c(0,20)
plot(Bio.sub.mask$bio1_12/10, legend.only=TRUE,col=Temp.col(9),
     legend.width=1, legend.shrink=1,
     axis.args=list(at=seq(Bio.sub.mask.range[1], Bio.sub.mask.range[2], 10),
                    labels=seq(Bio.sub.mask.range[1], Bio.sub.mask.range[2],10),
                    cex.axis=1.25),
     legend.args=list(text=' Temperature (°C)', side=4, font=2, line=2.5, cex=1.25))

dev.off()

# << Annual Precipitation map >> -------------
pdf('./figures/PrecipitationMap.pdf', height = 7, width = 10)

# add sample points to the ann. prec. map, group them by colors based on biotic communities
color.vec<-c("green", "gray","blue")
Ann.prec.col<-colorRampPalette(c("gold","yellow","cyan","deepskyblue","lightblue","cornflowerblue", "blue"))
plot(Bio.sub.mask$bio12_12, legend=FALSE,axes=FALSE, box=FALSE, col=Ann.prec.col(7),
     xlab= "", ylab="")
plot(Arizona,add=TRUE)
axis(1, pos=31, cex.axis=2, tck=-0.02)
axis(2, las=2, cex.axis=2, tck=-0.02)
# mtext(text = 'Longitude', side = 1, line = 4)
# mtext(text = 'Latitude', side = 2, line = 3)
points(sites$long, sites$lat, pch=19, cex=2)

# precipitation legend
par(xpd=TRUE)
Bio.sub.mask.range<-c(100,1000)
plot(Bio.sub.mask$bio12_12, legend.only=TRUE,col=Ann.prec.col(7),
     legend.width=2, legend.shrink=1,
     axis.args=list(at=seq(Bio.sub.mask.range[1], Bio.sub.mask.range[2], 200),
                    labels=seq(Bio.sub.mask.range[1], Bio.sub.mask.range[2],200), 
                    cex.axis=1.25),
     legend.args=list(text=' Precipitation (mm)', side=4, font=2, line=3, cex=1.25))

dev.off()

# << make elevation map with sample site points >> ------------
pdf('./figures/ElevationMap.pdf', height = 8, width = 10)

Alt<-getData("worldclim", var="alt", res=0.5, lon=-110, lat=33)
Alt.sub<-crop(Alt, extent(Arizona))
Alt.sub.mask<-mask(Alt.sub, Arizona)
Alt.col<-colorRampPalette(c("white", "gray90","gray80","gray70","gray60","gray50","gray30","black"))

# point color indicate geographic area
color.vec<-c("red", "yellow","blue","green")
plot(Alt.sub.mask$alt_12,legend=FALSE,axes=FALSE, box=FALSE, col=Alt.col(8),
     xlab= "", ylab="")
plot(Arizona,add=TRUE)
axis(1, pos=31, cex.axis=2, tck=-0.02)
axis(2, las=2, cex.axis=2, tck=-0.02)
# mtext(text = 'Longitude', side = 1, line = 4)
# mtext(text = 'Latitude', side = 2, line = 3)
points(sites$long, sites$lat, pch=19, cex=2)

# Legend
par(xpd=TRUE, mar = c(5,4,4,0))
Alt.sub.mask.range<-c(500,maxValue(Alt.sub.mask))
plot(Alt.sub.mask$alt_12, legend.only=TRUE,col=Alt.col(8),
     legend.width=2, legend.shrink=1,
     axis.args=list(at=seq(Alt.sub.mask.range[1], Alt.sub.mask.range[2], 1000),
                    labels=seq(Alt.sub.mask.range[1], Alt.sub.mask.range[2],1000), 
                    cex.axis=1.25),
     legend.args=list(text='Elevation (m)', side=4, font=2, line=3.5, cex=1.25))

dev.off()

