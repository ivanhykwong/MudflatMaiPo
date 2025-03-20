library(terra)
library(sf)
library(smoothr)
library(gstat)
library(stars)
library(dplyr)

#### input data ####
indir <- "" # Landsat input image directory (pre-processed)
imglist <- list.files(indir, pattern = ".tif$")
aoi <- vect("") # area of interest (shapefile)
outdir <- "" # output directory
OHTBTpath <- "" # Tide height data requested from HKO
imgmeta <- read.csv("") # a csv file storing metadata of the Landsat imagery, with at least 3 columns (Landsat image index, Date, Time) 
land <- vect("") # a land area shapefile for masking

#### crop AOI ####
for (imgname in imglist){
  img <- rast(file.path(indir, imgname))
  img <- crop(img, aoi)
  if (is.na(global(img[[1]],max,na.rm=TRUE)[[1]])==FALSE){ # not empty raster
    writeRaster(img, file.path(outdir, imgname))
  }
}

#### get height from time ####

get.sea.level <- function(imgtime){
  year <- format(imgtime,"%Y")
  tbt <- read.fwf(paste0(OHTBTpath,year), widths=c(2,2,rep(4,24)))
  colnames(tbt) <- c("month","day", seq(1,24))
  month <- as.numeric(format(imgtime,"%m"))
  day <- as.numeric(format(imgtime,"%d"))
  tbt_row <- tbt[tbt$month==month & tbt$day==day,]
  tbt_row <- tbt_row[,3:26]
  hour <- as.numeric(format(imgtime,"%H"))+8
  min1 <- as.numeric(format(imgtime,"%M"))
  min2 <- 60-min1
  height1 <- tbt_row[,hour]
  height2 <- tbt_row[,hour+1]
  if (height1==9999 | height2==9999){ return(9999) }
  height <- height1*min2/60 + height2*min1/60
  return(height)
}

#### create.waterline ####

imgdir <- outdir
setwd(imgdir)
imglist <- list.files(imgdir, pattern="tif$")
linedir <- file.path(imgdir, "waterline")
dir.create(linedir, showWarnings = F)

create.waterline <- function(imgname){
  # ignore <1991 and >2020
  imgyear <- as.integer(substr(imgname, 18, 21))
  if (imgyear<1991 | imgyear>2020) {return()}
  imgindex <- substr(imgname, 6, 25)
  imgdate <- imgmeta[imgmeta$index==imgindex,"Date"]
  imgtime <- imgmeta[imgmeta$index==imgindex,"Time"]
  imgtime <- paste(imgdate, imgtime)
  imgtime <- strptime(imgtime, tz = "UTC", format = '%d/%m/%Y %H:%M:%S')
  height <- get.sea.level(imgtime)
  if (height==9999) {return()} # end function if no valid height
  if (height > 215) {return()} # average high tide, no valid mudflat pixel if higher than this value
  
  img <- rast(imgname)
  if (substr(imgname,3,4)=="TM" | substr(imgname,3,4)=="ET"){
    names(img) <- c("B","G","R","NIR","SW1","SW2")
  } else if (substr(imgname,3,4)=="OL"){
    names(img) <- c("C","B","G","R","NIR","SW1","SW2")
  }
  img_out <- img
  brightmask <- img$R < 0.16
  img <- mask(img, brightmask, maskvalues=0)
  img <- mask(img, land, inverse=TRUE)
  ndvi <- (img$NIR-img$R)/(img$NIR+img$R)
  ndvi03 <- ndvi<0.3
  ndvi03 <- focal(ndvi03, 3, "min") # expand vegetation layer for 1 pixel
  img <- mask(img, ndvi03, maskvalues=0)
  
  ndwi <- (img$G-img$NIR)/(img$G+img$NIR)
  ndwi0 <- ndwi>0
  ndwi0_poly <- st_as_sf(as.polygons(ndwi0))
  if (nrow(ndwi0_poly)<=1) {return()} # 0 or 1 polygon, i.e. no water-land separation
  ndwi0_poly <- st_cast(ndwi0_poly,"POLYGON")
  ndwi0_poly$area <- st_area(ndwi0_poly)
  colnames(ndwi0_poly) <- c("A","geometry","area")
  mudflatline <- st_cast(ndwi0_poly[ndwi0_poly$A==0,],"LINESTRING")
  waterline <- st_cast(ndwi0_poly[ndwi0_poly$A==1,],"LINESTRING")
  interline <- st_intersection(mudflatline, waterline)
  interline <- st_collection_extract(interline, "LINESTRING")
  if (nrow(interline)==0) {return()} # 0 valid line
  
  interline$date <- as.character(imgtime)
  interline$height <- height
  interline <- interline[,c("date","height")]
  write_sf(interline,
           file.path(linedir,gsub("hk1980.tif","waterline.shp",imgname)),
           layer_options = "SHPT=ARC")
  writeRaster(img_out, file.path(linedir,imgname))
}
lapply(imglist, create.waterline)

#### contour.to.point ####

setwd(linedir)
shplist <- list.files(linedir, pattern=".shp$")
ptdir <- file.path(linedir,"pt")
dir.create(ptdir, showWarnings = F)

contour.to.point <- function(shpname){
  shp <- read_sf(shpname)
  if (st_geometry_type(shp, by_geometry = FALSE) != "MULTILINESTRING"){
    shp <- st_collection_extract(shp, "LINESTRING")
  }
  shp <- shp[,"height"]
  pt <- shp %>% st_cast("LINESTRING") %>% st_cast("MULTIPOINT") %>% st_cast("POINT")
  write_sf(pt, paste0("pt/",gsub("_waterline.shp", "_pt.shp", shpname)), 
           layer_options = "SHPT=POINT")
}
lapply(shplist, contour.to.point)

#### merge.point ####

setwd(ptdir)
shplist <- list.files(ptdir, pattern=".shp$")
resultdir <- file.path(outdir, "result")
dir.create(resultdir, showWarnings = F)

merge.shp.year <- function(shplist, y1, y2, outname){
  shpl <- substr(shplist, 18, 21)
  shpl <- as.integer(shpl)
  shpl <- shpl>=y1 & shpl<=y2
  shplist_y <- shplist[shpl]
  shplist_y <- lapply(shplist_y, read_sf)
  shp_merge <- do.call(rbind, shplist_y)
  shp_merge <- cbind(shp_merge, st_coordinates(shp_merge))
  df <- st_drop_geometry(shp_merge)
  df <- df %>% group_by(X,Y) %>% summarise(height=median(height))
  shp <- st_as_sf(df, coords=c("X","Y"), crs=st_crs(shp_merge))
  write_sf(shp, file.path(resultdir,outname))
}
merge.shp.year(shplist, 1991, 2000, "pt_1991-2000.shp")
merge.shp.year(shplist, 2001, 2010, "pt_2001-2010.shp")
merge.shp.year(shplist, 2011, 2020, "pt_2011-2020.shp")

#### kriging ####
# ArcGIS is used in this step
# Open pt_year.shp in ArcGIS Pro
# Geoprocessing tool -> Optimized Outlier Analysis -> Remove HL and LH points
# -> Ordinary Kriging in Geoprocessing
# Save result "pt_year_kriging.tif" inside resultdir (outdir/result)


#### get date from height ####

get.date <- function(y1, y2, h){
  shpdir <- linedir
  shplist <- list.files(shpdir, pattern=".shp$")
  shpl <- substr(shplist, 18, 21)
  shpl <- as.integer(shpl)
  shpl <- shpl>=y1 & shpl<=y2
  shplist_y <- shplist[shpl]
  for (shpname in shplist_y){
    shp <- st_read(file.path(shpdir,shpname), quiet=TRUE)
    if(round(shp$height[[1]],1)==h){
      print(shpname)
    }
  }
}
# Example: get.date(1991, 2000, 102.4)

#### create water layer ####

setwd(imgdir)
imglist <- list.files(imgdir, pattern="tif$")
suppvectdir <- file.path(outdir,"suppvect")
dir.create(suppvectdir, showWarnings = F)

get.wat.mask <- function(imgname){
  img <- rast(imgname)
  if (substr(imgname,3,4)=="TM" | substr(imgname,3,4)=="ET"){
    names(img) <- c("B","G","R","NIR","SW1","SW2")
  } else if (substr(imgname,3,4)=="OL"){
    names(img) <- c("C","B","G","R","NIR","SW1","SW2")
  } else if (substr(imgname,3,4)=="PL") {
    names(img) <- c("B","G","R","NIR")
    img <- img/10000
  }
  brightmask <- img$R < 0.16
  img <- mask(img, brightmask, maskvalues=0)
  ndwi <- (img$G-img$NIR)/(img$G+img$NIR)
  ndwi0 <- ndwi>=0
  vegl <- mask(ndwi0, land, inverse=TRUE)
  return(vegl)
}
watlist <- lapply(imglist, get.wat.mask)

merge.wat <- function(watlist, y1, y2, outname){
  imgl <- substr(imglist, 18, 21)
  imgl <- as.integer(imgl)
  imgl <- imgl>=y1 & imgl<=y2
  watlist_y <- watlist[imgl]
  ras <- rast(watlist_y)
  ras <- mean(ras, na.rm=TRUE) > 0.5
  ras[ras==0] <- NA
  writeRaster(ras, file.path(suppvectdir, outname))
}
merge.wat(watlist, 1991, 2000, "water_1991-2000.tif")
merge.wat(watlist, 2001, 2010, "water_2001-2010.tif")
merge.wat(watlist, 2006, 2020, "water_2011-2020.tif")

#### create intertidal vegetation layer ####

get.veg.mask <- function(imgname){
  img <- rast(imgname)
  if (substr(imgname,3,4)=="TM" | substr(imgname,3,4)=="ET"){
    names(img) <- c("B","G","R","NIR","SW1","SW2")
  } else if (substr(imgname,3,4)=="OL"){
    names(img) <- c("C","B","G","R","NIR","SW1","SW2")
  }
  brightmask <- img$R < 0.16
  img <- mask(img, brightmask, maskvalues=0)
  ndvi <- (img$NIR-img$R)/(img$NIR+img$R)
  ndvi03 <- ndvi>=0.3
  vegl <- mask(ndvi03, land, inverse=TRUE)
  return(vegl)
}
veglist <- lapply(imglist, get.veg.mask)

merge.veg <- function(veglist, y1, y2, outname){
  imgl <- substr(imglist, 18, 21)
  imgl <- as.integer(imgl)
  imgl <- imgl>=y1 & imgl<=y2
  veglist_y <- veglist[imgl]
  ras <- rast(veglist_y)
  ras <- mean(ras, na.rm=TRUE) > 0.5
  ras[ras==0] <- NA
  writeRaster(ras, file.path(suppvectdir, outname))
}
merge.veg(veglist, 1991, 2000, "veg_1991-2000.tif")
merge.veg(veglist, 2001, 2010, "veg_2001-2010.tif")
merge.veg(veglist, 2011, 2020, "veg_2011-2020.tif")

#### mask water, land and vegetation ####

setwd(resultdir)
ras1991 <- rast("pt_1991-2000_kriging.tif")
ras2001 <- rast("pt_2001-2010_kriging.tif")
ras2011 <- rast("pt_2011-2020_kriging.tif")

setwd(suppvectdir)
veg1991 <- rast("veg_1991-2000.tif")
veg2001 <- rast("veg_2001-2010.tif")
veg2011 <- rast("veg_2011-2020.tif")
wat1991 <- rast("water_1991-2000.tif")
wat2001 <- rast("water_2001-2010.tif")
wat2011 <- rast("water_2011-2020.tif")

ras1991 <- ras1991 |> mask(veg1991,inverse=T) |> mask(land,inverse=T) |> mask(wat1991,inverse=T)
ras2001 <- ras2001 |> mask(veg2001,inverse=T) |> mask(land,inverse=T) |> mask(wat2001,inverse=T)
ras2011 <- ras2011 |> mask(veg2011,inverse=T) |> mask(land,inverse=T) |> mask(wat2011,inverse=T)

#### temporal smoothing ####

tempsmooth <- function(d1,d2,d3){
  if (!is.null(d1)) {
    if (res(d1)[1] != res(d2)[1]) {d1 <- resample(d1,d2,"bilinear")}
    d1 <- mask(d1, d2)
    d <- c(d1,d2,d2)
  } else {
    d <- c(d2,d2)
  }
  if (!is.null(d3)) {
    if (res(d3)[1] != res(d2)[1]) {d3 <- resample(d3,d2,"bilinear")}
    d3 <- mask(d3, d2)
    d <- c(d,d3)
  }
  d <- mean(d, na.rm=TRUE)
  return(d)
}

ras1991 <- tempsmooth(NULL, ras1991, ras2001)
ras2001 <- tempsmooth(ras1991, ras2001, ras2011)
ras2011 <- tempsmooth(ras2001, ras2011, ras2020)
ras2020 <- tempsmooth(ras2011, ras2020, NULL)

setwd(resultdir)
writeRaster(ras1991, "MudflatElevation_1991-2000.tif")
writeRaster(ras2001, "MudflatElevation_2001-2010.tif")
writeRaster(ras2011, "MudflatElevation_2011-2020.tif")
writeRaster(ras2020, "MudflatElevation_2020-2021.tif")

