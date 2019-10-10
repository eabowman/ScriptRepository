## Script created by Liz Bowman Oct. 10, 2019
## for community analyses (NMDS)

#----------------------------------------------------------------------------------------#
# Load data and libraries----
#----------------------------------------------------------------------------------------#
library (vegan)

#--Load site by species matrix
otu.matrix <- read.csv()

#----------------------------------------------------------------------------------------#
# Load libraries----
#----------------------------------------------------------------------------------------#


#--isolate species data by remove columns with metadata
comm.matrix <- otu.matrix[2:length(otu.matrix)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 2]
# Keep only samples with data; remove singletons sometimes can result 
# in some sites not having any data anymore
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 1, ] # remove rows with sums of 0
simmilarity.matrix <- otu.matrix[row.names(comm.matrix),] # remove the same samples from all of your data so they match

#--distance matrix of community data
# change method to the similarity index you are using
comm.dist <- vegdist(comm.matrix, method = "jaccard", binary = TRUE)

#--NMDS analysis
nmds.comm <- metaMDS(comm.dist, dist = "bray", permutations = 999,
                       try = 100, trymax = 100)

#--Check the stress
jaccard.otu$stress

#--Plot NMDS
# Create dataframe
data.scores <- as.data.frame(scores(nmds.comm)) # similarity matrix
# Add metadata required for plot; not necessary to put in this data frame, but makes it easier
data.scores$tree <- simmilarity.matrix$Tree # metadata1 needed for plot
data.scores$site <- simmilarity.matrix$Site # metadata1 needed for plot
data.scores$range <- simmilarity.matrix$Range # metadata1 needed for plot


ggplot() + 
  geom_point(data = data.scores,aes(x = NMDS1,
                                    y = NMDS2,
                                    color=range), # color or shape equal to explanatory variable.
             size=4) + # add the point markers
  coord_equal() +
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=28,margin = margin(t = 30)), # remove x-axis labels
        axis.title.y = element_text(size=28,margin = margin(r = 30)), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())
