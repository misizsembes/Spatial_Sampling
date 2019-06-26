#THIS TOOL RANDOMLY SAMPLES SHELTERS INSIDE CAMPS (COORDINATES OF SHELTERS NEEDED)
#IT CAN BE ADJUSTED TO SIMPLY DRAW RANDOM POINTS 
#IT DRAWS RANDOM POINTS AND THEN IDENTIFIES THE NEAREST SHELTER; DUPLICATED SHELTERS ARE REMOVED

#SET WORKING DIRECTORY
setwd("~/Desktop/yemen/Jordan")
if(!require(SearchTrees)){
  install.packages("SearchTrees")
  library(SearchTrees)
}
if(!require(FNN)){
  install.packages("FNN")
  library(FNN)
}
if(!require(rgdal)){
  install.packages("rgdal")
  library(rgdal)
}
if(!require(spatstat)){
  install.packages("spatstat")
  library(spatstat)
}
library(raster)
library(tidyr)
library(dplyr)
library(sp)

#OPEN DATA AND PROJECT SPATIAL DATA
#CAMP BOUNDARIES 
camp_bound <- readOGR("Camp_Boundary/zaatari.shp")
#PROJECT CAMP BOUNDARY
 camp_bound <- spTransform(camp_bound, CRS("+proj=utm +zone=37N +datum=WGS84"))
 #IMPORT SHELTER SHAPEFILE
 households <- readOGR("Households/houses_nodupes2.shp")
 #ADD SEQUENTIAL IDs
 households$IDX <-seq.int(nrow(households))
 #PROJECT SHELTER SHAPEFILE
 households <- spTransform(households, CRS("+proj=utm +zone=37N +datum=WGS84"))
 #SET THE BUFFER AROUND THE ACTUAL SAMPLE SIZE -- TO ACCOUNT FOR ERRORS/ETC...
 sample_buffer <- .20
 
 #DISPLAY SHAPEFILES TO ENSURE THEY LAYER CORRECTLY
 plot(camp_bound)+
   points(households, col='red', cex=.1)
 
 #CREATE RASTER INSIDE CAMP BOUNDARIES--UNITS IN METERS (100 METER RESOLUTION)
 r <- raster(camp_bound)
 #SET RESOLUTION OF RASTER CELLS: SIZE OF EACH SQUARE (WHAT'S A MEANINGFUL SAMPLING UNIT)
 res(r) <- rastercellsize
 r <- rasterize(camp_bound, r)
 plot(r)+plot(camp_bound,add=TRUE)
 quads <- as(r, 'SpatialPolygons')
 plot(quads, add=TRUE)
 points(households, col='red', cex=.1)
 
 #SET THE RASTER CELL SIZE: WHAT SIZE COMPRISES MEANINGFUL A SPATIAL SAMPLING UNIT 
 #HOW BIG (units defined by the map projection) SHOULD EACH SQUARE INSIDE THE CAMP OUTLINE BE?
 #FOR EXAMPLE IF THE VALUE IS 100, THEN EACH SQUARE INSIDE THE CAMP SHAPEFILE WILL BE 100x100 units (e.g.,METERS)
 rastercellsize<-100 
 #SET SAMPLE SIZE
 sampsize <-150
 samp_w_buff<-sampsize+(round(sampsize*sample_buffer,0))
 
 #COUNT # OF HOUSEHOLDS IN EACH RASTER CELL
 nc <- rasterize(coordinates(households), r, fun='count', background=0)
 plot(nc)
 plot(camp_bound, add=TRUE)
 
 numhh <- mask(nc, r)
 plot(numhh)
 plot(camp_bound, add=TRUE)
#writeRaster(numhh, filename="rast1.tif", format="GTiff", overwrite=TRUE)

#CROP RASTER BY EXTENT OF CAMP BOUNDARIES
numhh2 <- crop(numhh, extent(246862.2, 249862.2, 3574611, 3577311))
plot(numhh2)

#CONVERT RASTER CELLS TO POLYGON FEATURES, WITH CELLS VALUES AS POLYGON ATTRIBUTES
pol <- rasterToPolygons(numhh2)
pol$ID <- seq.int(nrow(pol))
#writeOGR(pol, ".", "whatever2", driver="ESRI Shapefile")
#plot(camp_bound) ; points(spsample(camp_bound, n=567, type='random'), col='red', pch=3, cex=0.5)

#SAMPLE RASTER CELLS, PROBABILITY WEIGHTED BY THE NUMBER OF HOUSEHOLDS IN THE POLYGON
surveys<-sample(pol$ID, samp_w_buff, replace = TRUE, prob = pol$layer)
surveys <- data.frame(matrix(unlist(surveys), nrow=samp_w_buff, byrow=T),stringsAsFactors=FALSE)
colnames(surveys)[1] <-"select_loc"
surveys$ones <- 1

#COUNT THE NUMBER OF TIMES EACH POLYGON WAS DRAWN IN THE RANDOM SAMPLING 
#(THIS IS THE NUMBER OF SURVEYS NEEDS IN THAT POLYGON)
points_poly <- surveys %>%
  group_by(select_loc) %>%
  summarise_all(funs(sum), na.rm = TRUE)
colnames(points_poly)[2]<- "surveys"
colnames(points_poly)[1]<- "polyid"

#SELECT THE POLYGONS THAT WERE SELECTED AT ALL
pol<-pol[pol$ID %in% points_poly$polyid ,]
plot(pol)

#GENERATE RANDOM POINTS IN EACH POLYGON, THE NUMBER DRAWN BEING THE NUMBER OF TIMES 
#THE POLYGON WAS DRAWN IN THE RANDOM SAMPLING 
xpointz <- vector("list", nrow(points_poly))
ypointz <- vector("list", nrow(points_poly))
koordz<- vector("list",nrow(points_poly))
for (i in 1:nrow(points_poly)){
  block<-subset(pol,(pol$ID == points_poly$polyid[i]))
 point <- spsample(block, n=points_poly$surveys[i], type='random')
koordz[[i]]<-point@coords
  xpointz[[i]]<-point@coords[,1]
 ypointz[[i]]<-point@coords[,2]
}

#UNLIST X & Y COORDINATES AND JOIN IN DATAFRAME
xcoord  <-  as.data.frame(matrix(unlist(xpointz), nrow=length(unlist(xpointz[1]))))
xcoord <-t(xcoord)
xcoord<-as.numeric(xcoord)
ycoord  <-  as.data.frame(matrix(unlist(ypointz), nrow=length(unlist(ypointz[1]))))
ycoord <-t(ycoord)
ycoord<-as.numeric(ycoord)

coordz<-data.frame(xcoord,ycoord)

#CHECK IF SAMPLE POINT GENERATION LOOKS CORRECT
plot(numhh2)
points(coordz, col='red', cex=.3)

#COORDINATES FOR HOUSEHOLDS
hlong<-households@coords[,2]
hlat<-households@coords[,1]
hcoordz<-data.frame(hlat,hlong)

#NEAREST NEIGHBOR (NEAREST HOUSEHOLD TO SAMPLE POINT)
nn = get.knnx(hcoordz,coordz,k=1)
nn$ID<-seq.int(nrow(coordz))
  neighbor <- data.frame(nn$ID,nn$nn.index)
  colnames(neighbor)[1] <-"order"
  colnames(neighbor)[2] <-"nghbr_id"
  
  #REMOVE DUPLICATE HOUSEHOLDS
  deduped_HHs_toSurvey <- unique( neighbor[ , 2 ] )
  deduped_HHs_toSurvey<-as.data.frame(deduped_HHs_toSurvey)
  colnames(deduped_HHs_toSurvey)[1] <-"IDX"
  deduped_HHs_toSurvey2<-sample(deduped_HHs_toSurvey$IDX, sampsize, replace = FALSE)
  deduped_HHs_toSurvey2 <- as.data.frame(deduped_HHs_toSurvey2)
  colnames(deduped_HHs_toSurvey2)[1] <-"IDX"
  #CREATE DATAFRAME FROM THE HOUSEHOLD SHAPEFILE
  hh_justnum<-data.frame(households$house_id,households$houses_dat,households$houses_d_1,households$IDX)
  colnames(hh_justnum)[1] <-"house_id"
  colnames(hh_justnum)[2] <-"long"
  colnames(hh_justnum)[3] <-"lat"
  colnames(hh_justnum)[4] <-"IDX"
  #JOIN DEDUPED SAMPLE HOUSEHOLDS WITH MASTER HOUSEHOLD LIST TO OBTAIN HOUSE ID AND COORDINATES
    HHToSurvey<-merge(x = deduped_HHs_toSurvey2, y = hh_justnum, by = "IDX", all.x = TRUE)
    nrow(HHToSurvey)
  write.csv(HHToSurvey,paste0("HHToSurvey_", sampsize,".csv"))
  
  
  #CHECK DISTRIBUTION OF SAMPLED HOUSEHOLDS 
  #(SHOULD LOOK SIMILAR TO PROBABILITY DENSITY OF HOUSEHOLDS)
  checksampdist <- rasterize(coordinates(coordz), r, fun='count', background=0)
  plot(checksampdist)
  plot(camp_bound, add=TRUE)
