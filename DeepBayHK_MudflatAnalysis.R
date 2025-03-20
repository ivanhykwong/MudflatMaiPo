library(terra)
library(dplyr)

setwd("") # Result directory
mud1991 <- "MudflatElevation_1991-2000.tif"
mud2001 <- "MudflatElevation_2001-2010.tif"
mud2011 <- "MudflatElevation_2011-2020.tif"
mud2020 <- "MudflatElevation_2020-2021.tif"
mud2020 <- resample(mud2020, mud2011)
aoi <- vect("") # shapefile for cropping all raster to same extent
suppvectdir <- "" # supplementary vector directory


#### Figure 9 area coverage ####

pix.area.ha <- function(rasname, height){
  ras <- rast(rasname)
  ras <- crop(ras, aoi)
  if (res(ras)[1]==3) {
    ras <- aggregate(ras, fact=10)  # pixel size 3 to 30
  }
  ras <- ifel(ras<=height, NA, (res(ras)[1]^2)/1000) # pixel area in ha
  return(global(ras,"sum",na.rm=TRUE)[[1]])
}

df <- data.frame()
yearlist <- c(1991,2001,2011,2020)
raslist <- c(mud1991,mud2001,mud2011,mud2020)
for (i in 1:4){
  rasname <- raslist[i]
  df1 <- data.frame(year=yearlist[i], area60=pix.area.ha(rasname,60),
                    area90=pix.area.ha(rasname,90), 
                    area120=pix.area.ha(rasname,120), 
                    area145=pix.area.ha(rasname,145),
                    area170=pix.area.ha(rasname,170),
                    area200=pix.area.ha(rasname,200))
  df <- rbind(df,df1)
}

veg.area.ha <- function(rasname){
  ras <- rast(rasname)
  ras <- crop(ras, aoi)
  if (res(ras)[1]==3) {
    ras <- aggregate(ras, fact=10)  # pixel area in ha
  }
  ras <- (ras*0+1)*(res(ras)[1]^2)/1000 # pixel area in ha
  return(global(ras,"sum",na.rm=TRUE)[[1]])
}
for (r in file.path(suppvectdir, paste0("veg_", c("1991-2000.tif","2001-2010.tif","2011-2020.tif","2020-2021.tif")))){
  print(veg.area.ha(r))
}


#### Figure 10 Transect ####

# Draw transect line in ArcGIS 
# -> "Generate Points Along Lines" tool to create points at 3-m interval
transect_pt <- vect("transect_pt.shp")
df <- as.data.frame(transect_pt)
yearlist <- c(1991,2001,2011,2020)
raslist <- c(mud1991,mud2001,mud2011,mud2020)
for (i in 1:4){
  rasname <- raslist[i]
  ras <- rast(rasname)
  ras <- crop(ras, aoi)
  df1 <- terra::extract(ras, transect_pt, method="bilinear", ID=FALSE)
  colnames(df1) <- yearlist[i]
  df <- cbind(df, df1)
}
head(df)
write.csv(df, "transect_result.csv")


#### Figure 8 rate of change ####

mud2020_30m <- resample(mud2020, mud2011)
change1991 <- (mud2001-mud1991)/10 # 10 years
change2001 <- (mud2011-mud2001)/10
change2011 <- (mud2020_30m-mud2011)/5 # 5 years
changeoverall <- c(change1991,change1991,change2001,change2001,change2011)
changeoverall <- mean(changeoverall)
writeRaster(changeoverall, "changeoverall.tif")


