# Read the composition model results and mask UMW domain
# Andria Dawson

library(ncdf4)
library(ncdf4.helpers)
library(raster)
library(reshape2)
library(maptools)

us.shp <- readShapeLines('data/map_data/us_alb.shp', proj4string=CRS('+init=epsg:3175'))

#
# read data
#

# comp.nc = nc_open('data/composition/PLScomposition_western_0.2-release.nc')
# comp.nc = nc_open('data/composition/composition_midwest_v0.3.nc')
comp.nc = nc_open('data/composition/stat/composition_v0.4.nc')
mask.nc = nc_open('data/composition/paleonmask.nc')

taxa = nc.get.variable.list(comp.nc)

# XXX: this uses a lot of memory... move reading to loop below?
comp=list()
for (i in 1:length(taxa)){
  comp[[i]] = ncvar_get(comp.nc, taxa[i])
}
names(comp) = taxa

# get x and y coords
comp.x = nc.get.dim.for.axis(comp.nc, "Ash", "X")$vals
comp.y = nc.get.dim.for.axis(comp.nc, "Ash", "Y")$vals

water     = ncvar_get(mask.nc, 'water')
subregion = ncvar_get(mask.nc, 'subregion')
domain    = ncvar_get(mask.nc, 'domain')

#nc_close(comp.nc)
#nc_close(mask.nc)

#
# build grid, water, and region
#

westernDomainX <- 1:146  # thru x=1093000 out of 1:296
westernDomainY <- 1:180  # full N-S

comp.x <- comp.x[westernDomainX]
comp.y <- comp.y[westernDomainY]

# For "region", the codes are:
# 1: New England except Maine
# 2: Illinois
# 3 Indiana
# 4 Maine
# 5 southern Michigan (southern portion of lower peninsula)
# 6 Minnesota
# 7 New Jersey
# 8 New York
# 9 Ohio
# 10 Pennsylvania
# 11 Wisconsin
# 12 northern Michigan (UP and northern portion of LP)

reg_codes = c('new_england', 'illinois', 'indiana', 'maine', 'michigan_south', 
              'minnesota', 'new_jersey', 'new_york', 'ohio', 'pennsylvania', 
              'wisconsin', 'michigan_north')

water = water[westernDomainX, westernDomainY]
region = subregion[westernDomainX, westernDomainY]
region[region %in% c(2, 3, 9, 0)] <- NA # NA for states we don't want

coords  = expand.grid(x=comp.x, y=comp.y)

# compute counts
pm = list()
counts = matrix(NA, nrow=length(coords[,1]), ncol=(length(comp)))
for (i in 1:length(comp)){
  pm[[i]] = apply(comp[[i]][westernDomainX, westernDomainY,], c(1,2), 'mean')  
  # reverse the y axis
#   pm[[i]] = pm[[i]][,ncol(pm[[i]]):1]
  counts[,i] = as.vector(pm[[i]])
}
names(pm) = taxa
colnames(counts) = taxa

comp.df = data.frame(x=coords$x, y=coords$y, region=as.vector(region), water=as.vector(water))
veg = cbind(comp.df, counts)
veg = veg[!is.na(veg$region),]

# do we need to exclude interior cells?
# veg2 = veg[(veg$water<100),]

# check the orientation!
library(ggplot2)
d <- ggplot(veg) + geom_point(aes(x=veg$x, y=veg$y, colour=veg$Oak))
d

# don't want the code, want the state name
veg$region = reg_codes[veg$region]

write.table(veg, file='data/composition/composition_v0.5.csv', sep=',', row.names=F)


























