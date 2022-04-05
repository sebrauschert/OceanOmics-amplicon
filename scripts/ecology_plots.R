# Workflow with phyloseq
#
# taken from here: https://f1000research.com/articles/5-1492/v1


library(tidyverse)
library(phyloseq)
library(microbiome)
library(maptools)
library(ggplot2)
library(ozmaps)
library(sf)

oz_states <- ozmaps::ozmap_states
oz_states

# Specify paths to data 
otu_mat    <- "03-dada2/results/RS19_16S/RS19_asv_final_table_16S.csv"
tax_mat    <- "04-taxa/LCA_RS19/lca_results_16S_95.csv"
samples_df <- "05-report/raw_data_rowleys_eDNA_2019.csv"


#===========================================================================
# CUSTOM WRANGLING
meta %>%
  select(`Sample ID`, names(meta)[names(meta) != "Sample ID"]) %>%
  write_csv(samples_df)

taxa %>%
  select(OTU, names(taxa)[names(taxa) != "OTU"]) %>%
  filter(domain %in% 'Eukaryota') %>%
  filter(!(OTU %in% "ASV_4111")) -> taxa
write_csv(taxa, paste0(tax_mat))

# Tests to check if there are duplicated ASVs
taxa %>% count(OTU) %>%
  arrange(desc(n)) %>%
  filter(n >1) %>%
  select(OTU) %>%
  unlist() %>%
  as.character() -> duplicate_ASVs

names(otu) <- str_remove(as.vector(unlist(names(otu))) , ">")

write_csv(otu, otu_mat)
#===========================================================================

file.exists(otu_mat)
file.exists(tax_mat)
file.exists(samples_df)


# Read data in as a phyloseq object
rowleys2019 <- read_csv2phyloseq(
  otu.file = otu_mat,
  taxonomy.file = tax_mat,
  metadata.file = samples_df,
  sep = ",")

#===========================================================================
# Plot Map
coords <- read_csv(samples_df) %>%
  select(Lat, Long) %>%
  rename(lat = Lat,
         long = Long)
  
ggplot(oz_states) + 
  geom_sf() + 
  coord_sf(xlim=c(min(coords$long)-5, max(coords$long+5)), 
           ylim=c(min(coords$lat)-5, max(coords$lat)+5)) + 
  geom_sf_label(aes(label = NAME), label.padding = unit(1, "mm")) +
  geom_point(data = coords, 
             mapping = aes(x = long, y = lat), 
             fill = "tomato", size = 2, alpha=0.8, shape=21) + 
  xlab('Longitude') +
  ylab('Latitude') +
  theme_map()



#===========================================================================
# Remove irrelevant taxa
rowleys2019_no_mito = subset_taxa(rowleys2019, !domain=="NA" )
rowleys2019_no_mito = subset_taxa(rowleys2019_no_mito, !domain=="Eukaryota" )
rowleys2019_no_mito = subset_taxa(rowleys2019_no_mito, !phylum=="NA" )
rowleys2019_no_mito = subset_taxa(rowleys2019_no_mito, !order=="Chloroplast" )
rowleys2019_no_mito = subset_taxa(rowleys2019_no_mito, !family=="Mitochondria" )
#rowleys2019_no_mito = subset_samples(rowleys2019_no_mito, Depth < 50) # run this line if you want to remove or keep deep samples


#===========================================================================
# Alpha Diversity

plot_richness(rowleys2019_no_mito, x="Atoll", 
              measures=c("Shannon", "Simpson")) +
  theme_bw()

#===========================================================================
# Ordination plot 
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop       <- transform_sample_counts(rowleys2019_no_mito, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, 
                          method="NMDS", 
                          distance="bray")

plot_ordination(ps.prop, 
                ord.nmds.bray, 
                color="Atoll", 
                title="Rowleys 2019") +
  geom_point(size=5) + 
  theme_bw()


#===========================================================================
# Top 20 Taxa: barplot

top20 <- names(sort(taxa_sums(rowleys2019_no_mito), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(rowleys2019_no_mito, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)

plot_bar_top20 = function (physeq, x = "Sample", y = "Abundance", fill = NULL, title = NULL, 
                        facet_grid = NULL) {
  
  
  top20    <- names(sort(taxa_sums(physeq), decreasing=TRUE))[1:20]
  ps.top20 <- transform_sample_counts(physeq, function(OTU) OTU/sum(OTU))
  ps.top20 <- prune_taxa(top20, ps.top20)
  
  
  mdf = psmelt(ps.top20)
  p = ggplot(mdf, aes_string(x = x, y = y, fill = fill))
  p = p + geom_bar(stat = "identity", position = "stack")
  p = p + theme(axis.text.x = element_text(angle = -90, hjust = 0))
  if (!is.null(facet_grid)) {
    p <- p + facet_grid(facet_grid)
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  return(p)
}


plot_bar_top20(rowleys2019_no_mito, x="Atoll", fill="family") + 
  facet_wrap(~Environment, scales="free_x") +
  labs(title = "Taxonomic Diversity by Atoll and Environment",
       subtitle = "Voyage: Rowley Shoals 2019",
       caption = "Data Source: OceanOmics, The Minderoo Foundation") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) 
