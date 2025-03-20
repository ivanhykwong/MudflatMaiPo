library(terra)
library(rjson)
library(sf)
library(smoothr)
library(gstat)
library(stars)

#### input data ####
imgdir <- "" # PlanetScope input image directory (pre-processed)
outdir <- "" # output directory
OHTBTpath <- "" # Tide height data requested from HKO
land <- vect("") # a land area shapefile for masking
aoi <- st_read("") # Area of interest shapefile
imglowesttide <- rast("") # PlanetScope input image with the lowest tide level, for creating water mask

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
  height <- height1*min2/60 + height2*min1/60
  return(height)
}

#### create.waterline ####

setwd(imgdir)
imglist <- list.files(imgdir, pattern="composite.tif$")
linedir <- file.path(imgdir, "waterline")
dir.create(linedir, showWarnings = F)

create.waterline <- function(imgname){
  imgmeta <- gsub("composite.tif", "composite_metadata.json", imgname)
  imgmeta <- fromJSON(file=imgmeta)
  imgtime <- imgmeta$properties$acquired
  imgtime <- strptime(imgtime, tz = "UTC", format = '%Y-%m-%dT%H:%M:%OS')
  height <- get.sea.level(imgtime)
  
  img <- rast(imgname)
  imgmask <- gsub("composite.tif", "composite_udm2.tif", imgname)
  imgmask <- rast(imgmask)$clear
  img <- mask(img, imgmask, maskvalues=0)
  
  names(img) <- c("B","G","R","NIR")
  brightmask <- img$R < 1600
  img <- mask(img, brightmask, maskvalues=0)

  ndvi <- (img$NIR-img$R)/(img$NIR+img$R)
  ndwi <- (img$G-img$NIR)/(img$G+img$NIR)
  ndwi0 <- ndwi>0 & ndvi<0.3
  ndwi0 <- focal(ndwi0, w=5, fun="modal", na.policy="omit", na.rm=TRUE)
  
  ndwi0_poly <- st_as_sf(as.polygons(ndwi0))
  ndwi0_poly <- st_cast(ndwi0_poly,"POLYGON")
  ndwi0_poly$area <- st_area(ndwi0_poly)
  ndwi0_poly <- ndwi0_poly[ndwi0_poly$area>units::set_units(100, m^2),]
  ndwi0_poly <- fill_holes(ndwi0_poly,units::set_units(100, m^2))
  ndwi0_poly$area <- st_area(ndwi0_poly)
  
  mudflatline <- st_cast(ndwi0_poly[ndwi0_poly$focal_modal==0,],"LINESTRING")
  waterline <- st_cast(ndwi0_poly[ndwi0_poly$focal_modal==1,],"LINESTRING")
  interline <- st_intersection(mudflatline, waterline)
  interline <- st_collection_extract(interline, "LINESTRING")
  
  interline$date <- as.character(imgtime)
  interline$height <- height
  interline <- interline[,c("date","height")]
  write_sf(interline, file.path(linedir, gsub("composite.tif","waterline.shp",imgname)))
}
lapply(imglist, create.waterline)

#### contour.to.point ####

setwd(linedir)
shplist <- list.files(linedir, pattern=".shp$")
ptdir <- file.path(linedir,"pt")
dir.create(ptdir, showWarnings = F)

contour.to.point <- function(shpname){
  shp <- st_read(shpname)
  if (shp[1,]$height > 215) {return()}
  if (st_geometry_type(shp, by_geometry = FALSE) != "MULTILINESTRING"){
    shp <- st_collection_extract(shp, "LINESTRING")
  }
  pt <- shp %>% st_cast("LINESTRING") %>% st_cast("MULTIPOINT") %>% st_cast("POINT")
  if (st_crs(pt) != st_crs(aoi)){
    pt <- st_transform(pt, st_crs(aoi))
  }
  pt <- st_intersection(pt, aoi)
  pt <- st_zm(pt) # remove Z
  set.seed(1234)
  pt <- pt[sample(nrow(pt), size=1000),] # sample.point 1000
  pt <- pt[,"height"]
  write_sf(pt, file.path(ptdir, gsub("_waterline.shp", "_pt.shp", shpname)), 
           layer_options = "SHPT=POINT")
}
lapply(shplist, contour.to.point)

#### merge.point ####

setwd(ptdir)
shplist <- list.files(ptdir, pattern=".shp$")
resultdir <- file.path(outdir, "result")
dir.create(resultdir, showWarnings = F)

shplist <- list.files(ptdir, pattern=".shp$")
shplist <- lapply(shplist, st_read)
shp_merge <- do.call(rbind, shplist)
write_sf(shp_merge, file.path(resultdir,"pt.shp"))

#### kriging ####
# ArcGIS is used in this step
# Open pt.shp in ArcGIS Pro
# Geoprocessing tool -> Optimized Outlier Analysis -> Remove HL and LH points
# -> Ordinary Kriging in Geoprocessing
# Save result "pt_2020-2021_kriging.tif" inside resultdir (outdir/result)

#### create water layer ####

suppvectdir <- file.path(outdir,"suppvect")
dir.create(suppvectdir, showWarnings = F)

img <- imglowesttide
names(img) <- c("B","G","R","NIR")
brightmask <- img$R < 1600
img <- mask(img, brightmask, maskvalues=0)
ndwi <- 1.0*(img$G-img$NIR)/(img$G+img$NIR)
ndwi0 <- ndwi>=0
wat <- mask(ndwi0, land, inverse=TRUE)
wat <- focal(wat,5,"max",na.rm=TRUE,na.policy="omit")
wat[wat==0] <- NA
writeRaster(wat, file.path(suppvectdir, "water_planet.tif"))

#### create intertidal vegetation layer ####

setwd(imgdir)
imglist <- list.files(imgdir, pattern="composite.tif$")
imglist <- lapply(imglist, rast)
rsrc <- sprc(imglist)
img <- mosaic(rsrc, fun="median")
img <- crop(img,imglist[[1]])
names(img) <- c("B","G","R","NIR")
ndvi <- (img$NIR-img$R)/(img$NIR+img$R)
ndvi03 <- ndvi>0.3
veg <- mask(ndvi03, land, inverse=TRUE)
veg[veg==0] <- NA
writeRaster(veg, file.path(suppvectdir, "veg_planet.tif"))

#### mask water, land and vegetation ####

setwd(resultdir)
ras2020 <- rast("pt_2020-2021_kriging.tif")
setwd(suppvectdir)
wat2020 <- rast("water_planet.tif")
wat2020 <- crop(wat2020, ras2020)
veg2020 <- rast("veg_planet.tif")
veg2020 <- crop(veg2020, ras2020)

ras2020 <- ras2020 |> mask(veg2020,inverse=T) |> mask(land,inverse=T) |> mask(wat2020,inverse=T)
setwd(resultdir)
writeRaster(ras2020, "MudflatElevation_2020-2021.tif")

