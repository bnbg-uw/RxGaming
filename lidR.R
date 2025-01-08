#TODO handle units!!

# September 2024
# Bryce Bartl-Geller
# bnbg@uw.edu
# Forest Resilience Lab, University of Washington

# This  script is delivered along side  RxGaming and is intended to prepare raw
# lidar point data into a post-processed form ready to be read into an RxGaming
# Project.

#Location of raw lidar point data
pointdata = "\\\\172.25.182.82\\pfc-frl\\01_LiDAR_data\\LiDAR_data_backup\\CA_Creek_Fire_2021_ASO\\laz"

#treatment units shapefile
treatmentunits = "E:/Dropbox/RxGaming Paper/map/roush15.shp"

#unique identifier column from the shapefile
uniqueid = "Poly_Name"

#location to write output files
outdir = "F:/tmpcreek/"

# Do not modify anything below this line
#-------------------------------------------------------------------------------

#Load or install required packages
library(tools)
if(!require(lidR) | !require(sf)) {
  install.packages(c("lidR", "sf"))
  library(lidR)
  library(sf)
}

if(!require(EBImage)) {
  install.packages("BiocManager") 
  BiocManager::install("EBImage")
  library(EBImage)
}

#Create output directories
dir.create(outdir)
dir.create(file.path(outdir, "las"))
dir.create(file.path(outdir, "chm"))
dir.create(file.path(outdir, "segments"))
dir.create(file.path(outdir, "layout"))
dir.create(file.path(outdir, "mask"))

#Read the treatment unit shapefile
prj = st_read(treatmentunits)

#Read the lidar data
lasctg = readLAScatalog(pointdata)

#project the tx shapefile to the lidar data and simplify the attribute table.  Write to disk.
prj = st_transform(prj, projection(lasctg))
w = which(names(prj) == uniqueid)
prj = prj[,w]
colnames(prj)[1] = "uniqueid"
st_write(prj, paste0(outdir, "/layout/layout.shp"))

#Begin clipping the las by the tx units
opt_output_files(lasctg) = paste0(outdir, "/las/lasclip_{uniqueid}")
clip_roi(lasctg, prj)

# Make jittered versions of the lidar data to better account for innacuracy in lidar points, and to break ties.
# This also ensure methodological consistency with Lapis.
f = list.files(paste0(outdir, "/las/"), full.names=T)
for(l in f) {
  .l = readLAS(l)
  .l = normalize_height(.l, tin())
  
  #Deal with outliers and noise
  .l = classify_noise(.l, sor())
  .l = .l[!.l$Withheld_flag,]
  .l = .l[!(.l$Classification %in% c(7,9,18)),]
  .l = .l[.l$Z > -8 & .l$Z < 100, ]
  for(x in c(-1,0,1)) {
    for(y in c(-1,0,1)) {
      print(paste(x,y))
      if(!x & !y) {
        offset = 0
        eps = 0
      } else if(x & y) {
        offset = sqrt(0.2^2 + 0.2^2)
        eps = -0.000002
      } else {
        offset = 0.2
        eps = -0.000001
      }
      new_las = .l
      new_las$X = new_las$X + x*offset
      new_las$Y = new_las$Y + y*offset
      new_las$Z = new_las$Z + eps
      writeLAS(new_las, paste0(file_path_sans_ext(l), "_", x, "_", y, ".las"))
    }
  }
  unlink(l)
}

#This function id's the tao highpoints.
idHighPoint = function(x) {
  if (is.na(x[5])) {
    return(NA)
  }
  if (x[5]==max(x, na.rm=T) & x[5] > 2 & !all(na.omit(x) == max(x, na.rm=T))) {
    return(1)
  }
  return(0)
}

#Create the layers needed for rxgaming
f = list.files(file.path(outdir, "/las/"), full.names=T)
for(u in prj$uniqueid) {
  print(u)
  w = grep(u, f)
  chms = list()
  for(i in 1:length(w)) {
    print(paste(i, "of 9"))
    l = readLAS(f[w[i]])
    chms[[i]] = rasterize_canopy(l, res=0.75, p2r())
    writeRaster(chms[[i]], paste0(outdir, i, ".tif"))
  }
  chm = mosaic(sprc(chms), fun = "max")
  kernel <- matrix(1, 3, 3)
  chm <- terra::focal(x = chm, w = kernel, fun = mean, na.rm = TRUE)
  terra::writeRaster(chm, paste0(outdir, "/chm/", u, "_chm.tif"))
  mask = rasterize_canopy(l, res = 30)
  terra::writeRaster(mask, paste0(outdir, "/mask/mask_", u, ".tif"))
  tops = terra::focal(chm, 3, idHighPoint)
  terra::writeRaster(tops, paste0(outdir, "/segments/", u, "_tops.tif"))
  trees = lidR::watershed(chm, tol = 0, ext = 1)()
  terra::writeRaster(trees, paste0(outdir, "/segments/", u, "_segments.tif"))
}

r = NULL
for(.r in list.files(paste0(outdir, "/mask"), full.names=T)) {
  if(is.null(r)) {
    r = terra::rast(.r)
  } else {
    r = merge(r, terra::rast(.r))
  }
}
terra::writeRaster(r, paste0(outdir, "mask/mask.tif"))
