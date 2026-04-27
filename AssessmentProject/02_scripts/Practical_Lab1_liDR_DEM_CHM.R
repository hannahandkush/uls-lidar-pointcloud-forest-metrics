
#--- Practical Lab1 SCRIPT R DEM e CHM Course LiDAR ISA 2026 ############################################################
#Script to generate Digital Elevation Model (DEM) and Canopy Height Model (CHM) 
#-----------------------------Instalar e Importar packages -------------------------------

#if you do not have install Rtools in your computer you have to download from here and install the last version 
#https://cran.r-project.org/bin/windows/Rtools/
install.packages("lidR")
install.packages("terra")
install.packages ("sf")
install.packages("raster")
install.packages("sp")
install.packages("stringr")
install.packages("future")
install.packages("RMCC")
install.packages("gstat")

library(lidR)
library(terra)
library(sf)
library(raster)
library(sp)
library(stringr)
library(future)
library(RMCC)

# Setting the number of cores--------------------------

plan(multisession, workers = 2L)
set_lidr_threads(2L)

#Simple rule
#Setting	Meaning
#workers	number of parallel R processes
#set_lidr_threads()	threads used inside each process
#Total CPU usage ≈ workers × threads

#1----------------------------Working directories ------------------------------------------

#define your working directories, currently aligned with Module 2 G3 assessmnet folder
mainDir <- "C:/Users/isa129381/OneDrive - Universidade de Lisboa/LiDar/module2/AssessmentProject/01_data"
ResultadosDir <- "C:/Users/isa129381/OneDrive - Universidade de Lisboa/LiDar/module2/AssessmentProject/03_processing"
shpPath <- "C:/Users/isa129381/OneDrive - Universidade de Lisboa/LiDar/module2/AssessmentProject/01_data"
##define working directories where the unclassified point cloud is stored

#Ecolyptus G3 las file entered here#
inDir <- file.path(mainDir, "20240626_CUB_EN_TEIX_I_3_100_8_REC_25.las")

#generating results folders

noiseDir <- file.path(ResultadosDir, "01_Noise"); dir.create(noiseDir, showWarnings = F)
groundDir <- file.path(ResultadosDir, "02_GroundClassification"); dir.create(groundDir, showWarnings = F)
clipDir <- file.path(ResultadosDir, "03_Clipped"); dir.create(clipDir, showWarnings = F)
nlasDir <- file.path(ResultadosDir, "04_Normalized"); dir.create(nlasDir, showWarnings = F)
dtmDir <- file.path(ResultadosDir, "05_DTM"); dir.create(dtmDir, showWarnings = F)
chmDir <- file.path(ResultadosDir, "06_CHM"); dir.create(chmDir, showWarnings = F)

#2 ---------------------------Read .laz with readLAS and also with readLAScatalog -------------------------------

las=readLAS(file.path(inDir))

ctg = readLAScatalog(file.path(inDir))

# Optimization parameters
opt_chunk_buffer(ctg) <- 5
opt_chunk_size(ctg) <- 25
opt_laz_compression(ctg) <- TRUE

opt_chunk_alignment(ctg) <- c(25, -4)

plot(ctg, chunk = T, map = F)

#Important generate an spatial index 

lidR:::catalog_laxindex(ctg)

las_check(ctg)

#3----------------------------Visualize the clouds with different options in the RGL viewer-------------------------

plot(las, color="RGB")

plot(las, color = "Intensity")

plot(las, color = "Z")

#4----------------------------Remove duplicate point and noise-------------------------------------


# Removing the duplicated points 
las_02 <- filter_duplicates(las)


# Noise classification
las_02A <- classify_noise(las_02, algorithm = sor(k = 10, m = 3, quantile = F))
las_02B <- classify_noise(las_02, algorithm = ivf(res = 1, n = 6))

# remove the point we have classified as noise

las_03A = filter_poi(las_02A, Classification !=LASNOISE)
las_03B = filter_poi(las_02B, Classification !=LASNOISE)

# Plotting the point cloud
plot(las_03A, bg = "black", axis = T, legend = T, color = "Z")   
plot(las_03B, bg = "black", axis = T, legend = T, color = "Z")    

# Saving the LAS Data

writeLAS(las_03A, file.path(noiseDir, "LugoForestal_lote1_eucalyptus_sem_noise.laz"))

# Saving memory
remove(las, las_02, las_02A, las_02B, las_03A,las_03B)

#5 ---------------------------Classify ground--------------------------------------------- two options and within each you can do PMF or MCC depending on vegetation (MCC heavier)

#two alternative to work with readLAS or readLAScatalog. To work with many .laz files is recomended work always with the catalog



#5a)working with readLAS------------------------------------ slower

#las=readLAS(file.path(noiseDir, "LugoForestal_lote1_eucalyptus_sem_noise.laz"))

#Apply algorithm to classify the ground

#1)Progressive Morphological Filter

#nube_clasificada_pmf <-classify_ground(las,pmf(ws=3, th=0.1),last_returns = TRUE)

#2)Multiscale Curvature Classification (MCC) <--- do this for eucalyptus project
#Computationally demanding, but yields better results in forested areas

#nube_clasificada_mcc <-classify_ground(las, mcc(1.5,0.3))

#export classified point clouds

#writeLAS(nube_clasificada_pmf, file.path(groundDir, "LugoForestal_lote1_eucalyptus_classificada_pmf.laz")) #change name for output .laz file 

#writeLAS(nube_clasificada_mcc, file.path(groundDir, "LugoForestal_lote1_eucalyptus_classificada_mcc.laz"))




#5b)working with readLAScatalog ------------------------------------ faster use MCC for higher accuracy

## Progressive Morphological Filter 
# PMF classification
#myfunPMF = function(cluster, ws, th){
#  las <- lidR::readLAS(cluster)
#  if (lidR::is.empty(las)) return(NULL)
#  las_02<-filter_duplicates(las)
#  las_02A <- classify_noise(las_02, algorithm = sor(k = 10, m = 3, quantile = F))
#  las_03A = filter_poi(las_02A, Classification !=LASNOISE)
#  las_04 <- classify_ground(las_03A, pmf(ws, th))
#  return(las_04)
#}

## MCC version
myfunMCC = function(cluster, ws, th){
  las <- lidR::readLAS(cluster)
  if (lidR::is.empty(las)) return(NULL)
  las_02<-filter_duplicates(las)
  las_02A <- classify_noise(las_02, algorithm = sor(k = 10, m = 3, quantile = F))
  las_03A = filter_poi(las_02A, Classification !=LASNOISE)
  las_04 <- classify_ground(las, mcc(1.5,0.3))
  return(las_04)
}


# Define the function where the results will be written to the results directory.
opt_output_files(ctg) <- file.path(groundDir, "01_Ground_MCC_{XLEFT}_{YBOTTOM}")

#Option to export the results to a temporary directory if you want to store space in your computer

#opt_output_files(ctg) <- file.path(tempdir(), "01_Ground_PMF_{XLEFT}_{YBOTTOM}")

#opt <- list(need_buffer = TRUE, automerge = TRUE)

out1 <- catalog_apply(ctg, myfunMCC, ws=3, th=0.1, .options = opt)

plot(out1)

#We can read the catalog, which is a list of .laz files, to merge them and export the result as a single .laz file.

las_merge = readLAS(out1)

writeLAS(las_merge, file.path(groundDir, "LugoForestal_lote1_eucalyptus_classificada_pmf_catalog.laz"))


#6----------------------------Read classified .laz  with catalog ----------

ctg=readLAScatalog(file.path(groundDir, "LugoForestal_lote1_eucalyptus_classificada_pmf.laz"))

# In the case of point clouds classified with more classes, the way to select 
# the classes we do not want to work with would be, for example:
# to remove classes 6 = building, 7 = noise, 12 = overlap

#ctg=readLAScatalog(file.path(groundDir, "LugoForestal_lote1_eucalyptus_classificada_pmf.laz"),select="*",filter="-drop_class 6 7 12")

opt_chunk_buffer(ctg) <- 10
opt_chunk_size(ctg) <- 250
opt_laz_compression(ctg) <- TRUE
opt_chunk_alignment(ctg) <- c(25, -4)
plot(ctg, chunk = T, map = F)


lidR:::catalog_laxindex(ctg)

las_check(ctg)

#7 ---------------------------clip our point cloud using the boundaries of our stands.-------------------------------------


Extent0=shapefile(paste0(shpPath, "//20240626_CUB_EN_TEIX_I_3_25_final.shp"))

plot(Extent0)

Extent1=raster::aggregate(Extent0)

plot(Extent1)

# export the results

opt_output_files(ctg) <- file.path(clipDir,"//LugoForestal_lote1_eucalyptus_classificada_pmf_clip")

#Recorte de la nube con el Extent1, para tenerla de base para los c?lculos posteriores

ClipExtent1= clip_roi(ctg, Extent1)


#8 ---------------------------Generating the Digital Terrain Model (DTM)   ---------------------------

#read .laz with readLAS ou readLAScatalog. 

#In this case, we can read the classified point cloud using readLAS (5a) or work with readLAScatalog (5b)

#las= readLAS(file.path(groundDir, "LugoForestal_lote1_eucalyptus_classificada_pmf.laz"))

ctg= readLAScatalog(file.path(groundDir, "LugoForestal_lote1_eucalyptus_classificada_pmf.laz"))

opt_chunk_buffer(ctg) <- 10
opt_chunk_size(ctg) <- 250
opt_laz_compression(ctg) <- TRUE
opt_chunk_alignment(ctg) <- c(25, -4)

lidR:::catalog_laxindex(ctg)

## 1 - Generete a DTM model with the TIN algorithm (triangulation, unpopular because of smoothing effect)
#DTM_TIN <- rasterize_terrain(ctg, res = 2, algorithm = tin(),use_class=2,keep_lowest = FALSE, pkg="terra")
#plot(DTM_TIN, bg = "white")

## 2 -  Invert Distance weighting (IDW) (probability decreases with distance, also unpopular, doesn't consider the environmental conditions)
#DTM_IDW <- rasterize_terrain(ctg,res=2,algorithm = knnidw(k = 10L, p = 2),pkg="terra")
#plot(DTM_IDW, bg = "white")

## 3 - Kriging (popuar choice for research because it models adaptively rather than blindly connecting points, computationally intense)
DTM_KRI <- rasterize_terrain(ctg,res=2, algorithm = kriging(k = 40))
plot_dtm3d(DTM_KRI, bg = "white")

# clip with vector .shp our stands and export. 

Extent=vect(paste0(shpPath, "/20240626_CUB_EN_TEIX_I_3_25_final.shp"))

#DTM_TIN_c=terra::mask(DTM_TIN,Extent,filename=paste0(dtmDir, "//DTM_TIN_2m.tif"),overwrite=TRUE)
#DTM_IDW_c=terra::mask(DTM_IDW,Extent,filename=paste0(dtmDir, "//DTM_IDW_2m.tif"),overwrite=TRUE)
DTM_KRI_c=terra::mask(DTM_KRI,Extent,filename=paste0(dtmDir, "//DTM_KRI_2m.tif"),overwrite=TRUE)


#9 ---------------------------Normalized the  point clouds-Set point heights above ground level.-----------------------------------------
#read clip of the stands

las= readLAS(file.path(clipDir, "LugoForestal_lote1_eucalyptus_classificada_pmf_clip.laz"))


#We normalize the point cloud. We have 3 algorithms or we can use a raster from the digital terrain model (DEM) 
#to perform the normalization.

#using tin e k-nn

#TIN = Triangulated Irregular Network: a method that creates a mesh of triangles from ground points to model the terrain.

#K-NN = K-Nearest Neighbors: a method that estimates ground height by averaging the heights of the nearest neighboring points

#las_n=normalize_height(las, tin())

#las_n2=normalize_height(las, knnidw())

#using the raster DEM we have generated previously 

las_n3=normalize_height(las, tin(),dtm=DTM_KRI)


#Filter to remove points with negative values

#las_n_f=filter_poi(las_n, Z>= 0)
#las_n2_f=filter_poi(las_n2,Z>=0)
las_n3_f=filter_poi(las_n3,Z>=0)

#plot(las_n_f, bg = "black", axis = T, legend = T, color = "Z")
#plot(las_n2_f, bg = "black", axis = T, legend = T, color = "Z")
plot(las_n3_f, bg = "black", axis = T, legend = T, color = "Z")

#export the results

writeLAS(las_n3_f, file.path(nlasDir, "LugoForestal_lote1_eucalyptus_normalizada.laz")) #filename


#10 --------------------------Generate the Canopy Height Models (CHM)------------------------------ (optional step, 2D index showing heights at different points because ABA is more informative)

ctg_norm= readLAScatalog(file.path(nlasDir, "LugoForestal_lote1_eucalyptus_normalizada.laz")) #filename

opt_chunk_buffer(ctg_norm) <- 10
opt_chunk_size(ctg_norm) <- 250
opt_laz_compression(ctg_norm) <- TRUE

plot(ctg_norm, chunk = T, map = F)

lidR:::catalog_laxindex(ctg_norm)


#Generate the Canopy Height Models (CHM) with 1 m resolutioin

#output raster the function attributes the height of the highest point found

chm1=rasterize_canopy(ctg_norm, 1,  p2r(),pkg = "terra") 

plot(chm1, col = height.colors(50),main="CHM")

# It implements the pit-free algorithm developed by Khosravipour et al. (2014), which is based on the computation of a set of classical triangulation at different heights   
# better in r than in python, this algorithm models adaptively to interpolate using multiple neighbouring points to estimate height of null data

chm2 = rasterize_canopy(ctg_norm, 1, pitfree(c(0,2,5,10,15,20,30,40), c(0,1), subcircle = 0.2))

plot(chm2)
#clip with vector .shp our stands and export package terra.

Extent=vect(paste0(shpPath, "//20240626_CUB_EN_TEIX_I_3_25_final.shp"))

chm=terra::mask(chm2,Extent,filename=paste0(chmDir, "//CHM_pit_free_1m.tif"),overwrite=TRUE) 

plot(chm, col = height.colors(50),main="CHM")


#postprocessing to fill NA data or smooth function

fill.na <- function(x, i=5) { if (is.na(x)[i]) { return(mean(x, na.rm = TRUE)) } else { return(x[i]) }}
w <- matrix(1, 3, 3)

#We can post-process the CHM using three types: original, filled, and smoothed 

filled <- terra::focal(chm[[1]], w, fun = fill.na)
smoothed <- terra::focal(chm[[1]], w, fun = mean, na.rm = TRUE)

chms <- c(chm, filled, smoothed)
names(chms) <- c("Base", "Filled", "Smoothed")

par(mfrow = c(3,1))

terra::plot(chms[[1]])
terra::plot(chms[[2]])
terra::plot(chms[[3]])

#Export the CHMs to raster using the boundary of the study area

Extent=vect(paste0(shpPath, "//20240626_CUB_EN_TEIX_I_3_25_final.shp"))

# Writing the raster
terra::mask(filled,Extent, filename=paste0(chmDir, "//CHM_pit_free_1m_filled.tif"),overwrite=TRUE)
terra::mask(smoothed, Extent, filename=paste0(chmDir, "//CHM_pit_free_1m_smooth.tif"), overwrite=TRUE)

