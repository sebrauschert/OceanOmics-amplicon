#===========================================================================
# Creating basic ecology plots with phyloseq
# Jessica Pearce
# 28/07/2022
#===========================================================================

library(tidyverse)
library(phyloseq)
library(microbiome)
library(maptools)
library(ggplot2)
library(ozmaps)
library(sf)
library(ggthemes)

oz_states <- ozmaps::ozmap_states

voyage = "RSV5"
assay = "16S"

#===========================================================================
# Read in phyloseq object from LULU output

phyloseq_object <- readRDS(paste0("06-report/",voyage,"_",assay,"_phyloseq.rds"))

# phyloseq <- subset_samples(Rowleys, !Time=="Aug-21")

#===========================================================================
# Plot Map
coords <- read_csv(paste0("06-report/",voyage,"_metadata.csv")) %>%
  select(Latitude, Longitude) %>%
  rename(lat = Latitude,
         long = Longitude)
  
ggplot(oz_states) + 
  geom_sf() + 
  coord_sf(xlim=c(min(coords$long, na.rm = T)-1, max(coords$long, na.rm = T)+1), 
           ylim=c(min(coords$lat, na.rm = T)-1, max(coords$lat, na.rm = T)+1)) + 
  geom_sf_label(aes(label = NAME), label.padding = unit(1, "mm")) +
  geom_point(data = coords, 
             mapping = aes(x = long, y = lat), 
             fill = "tomato", size = 2, alpha=0.8, shape=21) + 
  xlab('Longitude') +
  ylab('Latitude') +
  theme_minimal()

#===========================================================================
# Alpha Diversity

plot_richness(phyloseq_object, x="Atoll",
              measures=c("Shannon", "Simpson","Chao1")) +
  theme_bw()

# p + scale_x_discrete(limits = c("",""))

#===========================================================================
# Ordination plot 
# Transform data to proportions as appropriate for Bray-Curtis distances
# subset <- subset_samples(phyloseq_object, sample_names(phyloseq_object) !="")

ps.prop       <- transform_sample_counts(phyloseq_object, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, 
                          method="NMDS", 
                          distance="bray")

plot_ordination(ps.prop, 
                ord.nmds.bray, 
                color = "Atoll",
                title="Rowley Shoals and Montebelloes August 2021") +
  geom_point(size=5) + 
  theme_bw()

# Subset for only Rowley Shoals sites
subset <- subset_samples(phyloseq_object, !Atoll=="Montebellos")
ps.prop       <- transform_sample_counts(subset, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, 
                          method="NMDS", 
                          distance="bray")

plot_ordination(ps.prop, 
                ord.nmds.bray, 
                color = "Atoll",
                title="Rowley Shoals August 2021") +
  geom_point(size=5) + 
  theme_bw()

# Warning message:
  # In metaMDS(veganifyOTU(physeq), distance, ...) :
  # stress is (nearly) zero: you may have insufficient data

#===========================================================================
# Top 20 Taxa: relatice abundance

# top20 <- names(sort(taxa_sums(Rowleys), decreasing=TRUE))[1:20]
# ps.top20 <- transform_sample_counts(Rowleys, function(OTU) OTU/sum(OTU))
# ps.top20 <- prune_taxa(top20, ps.top20)

plot_bar_top20 = function (phyloseq_object, x = "Sample", y = "Abundance", fill = NULL, title = NULL, 
                           facet_grid = NULL) {
  
  
  top20    <- names(sort(taxa_sums(phyloseq_object), decreasing=TRUE))[1:20]
  ps.top20 <- transform_sample_counts(phyloseq_object, function(OTU) OTU/sum(OTU))
  ps.top20 <- prune_taxa(top20, ps.top20)
  
  mdf = psmelt(ps.top20)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "fill")
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}


plot_bar_top20(phyloseq_object, x="Atoll", fill="species") + 
  labs(title = "Taxonomic Diversity by Site",
       subtitle = "Voyage: Rowley Shoals and Montebellos August 2021",
       caption = "Data Source: OceanOmics, The Minderoo Foundation") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) 
# p + scale_x_discrete(limits = c("PC_1", "PC_2", "PC_3", "PC_4", "PC_5", "PC_6", "PC_7", "PC_8", "PC_9",
#                                 "PC_10", "PC_11", "PC_12", "PC_13", "PC_14", "PC_15", "PC_16", "PC_17", "PC_18", "PC_19", "PC_20", "D1"))
