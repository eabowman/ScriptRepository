## Code to calculate distance between two GPS coordinates
## Modified by Liz Bowman Sept. 24, 2019

# install.packages(c('geosphere','rlist'))
library(geosphere); library(rlist)

# Geographic coordinates of all sites
coords <- read.csv('data_output/SiteCoordinates.csv')

# Data for tree names to create distance matrix by tree
otu.data <- read.csv(paste0(dat.dir, '20190919_AllEMSequences.csv'), as.is = T)
# Create data frame with coords by tree
tree.dist <- data.frame(tree = unique(otu.data$Tree),
                        site = NA, 
                        lat = NA,
                        lon = NA,
                        pairwise.dist = NA)
for(i in tree.dist$tree){
  site.i <- unique(otu.data[otu.data$Tree == i, 'Site'])
  lat.i <- coords[coords$site == site.i, 'lat']
  lon.i <- coords[coords$site == site.i, 'long']
  tree.dist[tree.dist$tree == i, 'site'] <- site.i
  tree.dist[tree.dist$tree == i, 'lat'] <- lat.i
  tree.dist[tree.dist$tree == i, 'lon'] <- lon.i
}

#========================================================================================#
# Method 1 ------------------
#========================================================================================#

#--Function from https://blog.exploratory.io/calculating-distances-between-two-geo-coded-locations-358e65fcafae
# written by https://twitter.com/karasu)
get_geo_distance = function(long1, lat1, long2, lat2, units = "km") {
  loadNamespace("purrr")
  loadNamespace("geosphere")
  longlat1 = purrr::map2(long1, lat1, function(x,y) c(x,y))
  longlat2 = purrr::map2(long2, lat2, function(x,y) c(x,y))
  distance_list = purrr::map2(longlat1, longlat2, function(x,y) geosphere::distHaversine(x, y))
  distance_m = list.extract(distance_list, position = 1)
  if (units == "km") {
    distance = distance_m / 1000.0;
  }
  else if (units == "miles") {
    distance = distance_m / 1609.344
  }
  else {
    distance = distance_m
    # This will return in meter as same way as distHaversine function. 
  }
  distance
}

# Site level ----
#--Create matrix allowing pair-wise distances to be calculated
pairwise.dist <- data.frame(site1 = rep(coords$site, 40),
                            site2 = rep(coords$site, each = 40),
                            dist = NA)

#--Calculate the distance between each site and each site
for(i in pairwise.dist$site1){
  for(t in pairwise.dist$site1){
    dist.i <- get_geo_distance(coords[coords$site == i, 'long'],
                               coords[coords$site == i, 'lat'],
                               coords[coords$site == t, 'long'],
                               coords[coords$site == t, 'lat'])
    pairwise.dist[pairwise.dist$site1 == i & pairwise.dist$site2 == t,
                  'dist'] <- dist.i 
  }
}


write.csv(pairwise.dist, 'data_output/20191009_PairwiseSpatialDistances_site.csv',row.names = F)

# Site level ----
#--Create matrix allowing pair-wise distances to be calculated
# pairwise.dist.tree <- data.frame(tree = rep(tree.dist$tree, each = 207),
#                                  site1 = rep(unique(tree.dist$site), each = 40),
#                                  site2 = rep(unique(tree.dist$site), 40),
#                                  dist = NA)
# 
# #--Calculate the distance between each site and each site
# for(i in pairwise.dist$site1){
#   for(t in pairwise.dist$site1){
#     dist.i <- get_geo_distance(coords[coords$site == i, 'long'],
#                                coords[coords$site == i, 'lat'],
#                                coords[coords$site == t, 'long'],
#                                coords[coords$site == t, 'lat'])
#     pairwise.dist[pairwise.dist$site1 == i & pairwise.dist$site2 == t,
#                   'dist'] <- dist.i 
#   }
# }
# 
# 
# write.csv(pairwise.dist, 'data_output/20191009_PairwiseSpatialDistances.csv',row.names = F)

#========================================================================================#
# Method 2 ------------------
#========================================================================================#

# install.packages('Imap')
library(Imap)

#--Function from https://eurekastatistics.com/calculating-a-distance-matrix-for-geographic-points-using-r/
ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace){
  # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
  # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.
  
  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}

GeoDistanceInMetresMatrix <- function(df.geopoints){
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.
  
  GeoDistanceInMetres <- function(g1, g2){
    # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
    # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
    # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
    # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
    # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="km")))
    }
    return(mapply(DistM, g1, g2))
  }
  
  n.geopoints <- nrow(df.geopoints)
  
  # The index column is used to ensure we only do calculations for the upper triangle of points
  df.geopoints$index <- 1:n.geopoints
  
  # Create a list of lists
  list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})
  
  # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")
  
  # Set the row and column names
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name
  
  return(mat.distances)
}

# Site level ----
names(coords) <- c('name','lat','lon')
dist.matrix <- GeoDistanceInMetresMatrix(coords)

write.csv(dist.matrix, 'data_output/20191009_SpatialDistancesMatrix_site.csv', row.names = T)

# Tree level ----
tree.coords <- tree.dist[c('tree','lat','lon')]
names(tree.dist.matrix)[1] <- 'name' 

tree.dist.matrix <- GeoDistanceInMetresMatrix(tree.coords)

# Add tree names to columns and rows
rownames(tree.dist.matrix) <- tree.coords$tree
names(tree.dist.matrix) <- tree.coords$tree

write.csv(tree.dist.matrix, 'data_output/20191009_SpatialDistancesMatrix_tree.csv', row.names = T)
