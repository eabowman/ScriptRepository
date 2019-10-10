## Script created by Liz Bowman Oct. 9, 2019
## for analysing creating pair-wise distance matrices

#========================================================================================#
# Load libraries and data ------------------
#========================================================================================#

#install.packages(c('dplyr','tidyr','vegan','ggplot2'))
library(dplyr);library(tidyr);library(vegan);library(ggplot2);library(nlme)

#--file paths
dat.dir <- '~/Documents/PhD/3_Sky_islands-EM/data/'
dat.out <- '~/Documents/PhD/3_Sky_islands-EM/data_output/'
fig.dir <- '~/Documents/PhD/3_Sky_islands-EM/figures/'
res.dir <- '~/Documents/PhD/3_Sky_islands-EM/results/'

# otu.data <- read.csv(paste0(dat.dir, '20190919_AllEMSequences.csv'), as.is = T)
clim.data <- read.csv(paste0(dat.out, '20190923_climate_data.csv'), as.is = T)
# stoc.data <- read.csv(paste0(dat.out,'20190924_SiteDistances.csv'), as.is = T)

# Site level environmental data
clim.data %>%
  select(site, range, lat, long, masl, prec, Tavg, forest.type) %>%
  distinct(.) -> clim.data.site

#========================================================================================#
# Environmental matrix: single factor and multifactor ------------------
#========================================================================================#

#<< Multifactor >> ----
# Site level ----
pca.env <- clim.data.site[c('prec','forest.type')]
pca.env$forest.type <- as.numeric(factor(pca.env$forest.type, labels = c(1:6))) # make numeric
rownames(pca.env) <- clim.data.site$site

pca.env.analysis <- prcomp(pca.env,
                           center = TRUE,
                           scale. = TRUE)

summary(pca.env.analysis)
env.eigenvector <- scores(pca.env.analysis, choices = c(1:4))
env.eigenvector <- data.frame(env.eigenvector)
names(env.eigenvector) <- c('PCA1','PCA2')
env.eigenvector$site <- rownames(pca.env)

clim.data.site['env.pca'] <- NA
for(i in unique(clim.data.site$site)){
  clim.data.site[clim.data.site$site == i, 'env.pca'] <- 
    env.eigenvector[env.eigenvector$site == i, 'PCA1']
}

# distance matrix
multi.env.site <- dist(clim.data.site$env.pca, method = 'euclidean')
multi.env.site <- as.matrix(multi.env.site)
rownames(multi.env.site) <- clim.data.site$site
colnames(multi.env.site) <- clim.data.site$site

write.csv(multi.env.site, 'data_output/20191009_EnvironmentDistanceMatrix_site.csv',row.names = T)

# Tree level ----
# Perform PCA of environmental factors
pca.env <- clim.data[c('prec','forest.type')]
pca.env$forest.type <- as.numeric(factor(pca.env$forest.type, labels = c(1:6))) # make numeric
rownames(pca.env) <- clim.data$tree

pca.env.analysis <- prcomp(pca.env,
                           center = TRUE,
                           scale. = TRUE)

summary(pca.env.analysis)
env.eigenvector <- scores(pca.env.analysis, choices = c(1:4))
env.eigenvector <- data.frame(env.eigenvector)
names(env.eigenvector) <- c('PCA1','PCA2')
env.eigenvector$tree <- rownames(pca.env)

clim.data['env.pca'] <- NA
for(i in unique(clim.data$tree)){
  clim.data[clim.data$tree == i, 'env.pca'] <- 
    env.eigenvector[env.eigenvector$tree == i, 'PCA1']
}

# distance matrix
multi.env.tree <- dist(clim.data$env.pca, method = 'euclidean')
multi.env.tree <- as.matrix(multi.env.tree)
rownames(multi.env.tree) <- clim.data$tree
colnames(multi.env.tree) <- clim.data$tree

write.csv(multi.env.tree, 'data_output/20191009_EnvironmentDistanceMatrix_tree.csv',row.names = T)

#<< Climate >> ----
# Site level ----
# distance matrix
clim.site <- dist(clim.data.site$prec, method = 'euclidean')
clim.site <- as.matrix(clim.site)
rownames(clim.site) <- clim.data.site$site
colnames(clim.site) <- clim.data.site$site

write.csv(clim.site, 'data_output/20191009_ClimateDistanceMatrix_site.csv',row.names = T)

# Tree level ----
# distance matrix
clim.tree <- dist(clim.data$prec, method = 'euclidean')
clim.tree <- as.matrix(clim.tree)
rownames(clim.tree) <- clim.data$tree
colnames(clim.tree) <- clim.data$tree

write.csv(clim.tree, 'data_output/20191009_ClimateDistanceMatrix_tree.csv',row.names = T)

#<< Forest type??? >> ----
# Site level ----
# distance matrix
for.site <- dist(clim.data.site$forest.type, method = 'euclidean')
for.site <- as.matrix(for.site)

write.csv(for.site, 'data_output/20191009_ForestDistanceMatrix_site.csv',row.names = F)

# Tree level ----
# distance matrix
for.tree <- dist(clim.data$forest.type, method = 'euclidean')
for.tree <- as.matrix(for.tree)

write.csv(for.tree, 'data_output/20191009_ForestDistanceMatrix_tree.csv',row.names = F)
