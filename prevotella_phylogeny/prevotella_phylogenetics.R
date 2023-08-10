# configure the environment -----------------------------------------------

# load packages
library(tidyverse)
library(readxl)
library(ape)
library(ggtree)
library(cluster)
library(vegan)
library(FactoMineR)

# set the working directory
setwd("working directory")

# file tracking
datestring <- format(Sys.time(), "%y%m%d_%H%M")
study <- "MDCF_POC"
experiment <- "prevotella_phylo"

# custom colors
cbbPalette <- c("#999999", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "darkred", "darkblue", "darkgrey", "white", "pink")

# load gtdb taxonomic assignments assignments for mags and references. "pm" is Primary MAM study
tax.pm <- read.table("primary mam mag taxonomic assignments", sep = "\t", header = TRUE)
tax.ref <- read_excel("reference genome/mag taxonomic assignments from Tett et al., 2019")
tax.ref2 <- read_excel("reference genome/mag taxonomic assignments from PATRIC (now BV-BRC")

# load tree (generated from fastree via accompanying sbatch scripts)
tree <- read.tree("tree.fasttree")
tree$tip.label <- gsub("^NZ_", "", tree$tip.label) # some tip formatting - a fraction of the primary MAM MAGs are NZ_ prefixed
tree$tip.label <- gsub("^Prevotella_copri_|^Prevotella_bivia_", "", tree$tip.label) # some tip formatting - some tips also had a species prefix

# generate selectors for prevotella mags from each set
prevotella_in_pm <- tax.pm %>%
  select(MAG, genus, species) %>%
  filter(genus == "g__Prevotella") %>%
  mutate(species = gsub("s__", "", species))

prevotella_isolates <- read_excel("isolate taxonomy") %>%
  filter(! isolate %in% c("Prevotella_copri_P43C7771105_G2_2_JG_bc1018", "Prevotella_copri_P43C7771105_G5_2_JG_bc1002")) %>%
  mutate(isolate = gsub("Prevotella_copri_|Prevotella_bivia_", "", isolate),
         genus = "Prevotella", 
         species = "Prevotella copri") %>%
  select(isolate, genus, species)

# Trim tree to the MAGs in the table
tree.trim <- drop.tip(tree, setdiff(tree$tip.label, c(prevotella_in_pm$MAG, 
                                                      prevotella_isolates$isolate,
                                                      tax.ref$MAG,
                                                      tax.ref2$entry,
                                                      "226186.12")))

# Root the tree on Bacteroides thetaiotaomicron VPI-5482
tree.reroot = root(tree.trim, '226186.12') 

# generate a labeling structure to tag each tree tip with its corresponding MAG or reference identity and metadata
prevotella_mag_groupings <- prevotella_in_pm %>%
  select(MAG, species) %>%
  mutate(study = "PM",
         strain = MAG) %>%
  bind_rows(prevotella_isolates %>%
              mutate(MAG = isolate,
                     study = "isolate",
                     strain = isolate) %>%
              select(-isolate, -genus)) %>%
  bind_rows(tax.ref %>%
              select(MAG, Clade) %>%
              dplyr::rename(species = Clade) %>%
              mutate(study = "ref_tett",
                     strain = MAG)) %>%
  bind_rows(tax.ref2 %>%
              dplyr::rename(species = isolate) %>%
              select(MAG, species, strain) %>%
              filter(!MAG %in% c(prevotella_isolates$isolate, "P43C7771105_G5_1_JG_v061022")) %>%
              mutate(study = "ref_patric")) %>%
  bind_rows(data.frame(MAG = "226186.12",
                       species = "Bacteroides thetaiotaomicron VPI-5482",
                       strain = "VPI-5482",
                       study = "ref_patric"))
prevotella_mag_groupings$MAG <- gsub("^Prevotella_copri_", "", prevotella_mag_groupings$MAG)
rownames(prevotella_mag_groupings) <- prevotella_mag_groupings$MAG

# initial tree view, without node labels. This tree will be too large to publish, almost to parse without zooming in illustrator
p <- ggtree(tree.reroot, layout = "rectangular") %<+% prevotella_mag_groupings.format + 
  geom_treescale() +
  geom_tiplab(aes(label = label), size = 2, offset = 0.001) +
#  geom_tiplab(aes(label = species), size = 2, offset = 0.02) +
  geom_tippoint(aes(color = study)) +
  scale_color_manual(values = cbbPalette)
print(p)

# I want to cut the tree down to something relevant. We want to look at Prevotella copri specifically, so
# I'll identify the portion of the tree relevant to this species and extract it. Re-plot the tree with node label numbers
p <- ggtree(tree.reroot, layout = "rectangular") %<+% prevotella_mag_groupings.format + 
  geom_treescale() +
  geom_tiplab(aes(label = label), size = 2, offset = 0.001) +
#  geom_tiplab(aes(label = species), size = 2, offset = 0.02, align = TRUE) +
  #geom_nodelab() +
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3) +
  geom_tippoint(aes(color = study)) +
  scale_color_manual(values = cbbPalette)
print(p)

# extract the portion of the tree that contains Prevotella copri
tree.reroot.pc <- extract.clade(tree.reroot, "node number")

# plot the subtree
p <- ggtree(tree.reroot.pc, layout = "rectangular") %<+% prevotella_mag_groupings.format + 
  geom_treescale() +
  geom_tiplab(aes(label = label), size = 2, offset = 0.001) +
  geom_tiplab(aes(label = species), size = 2, offset = 0.02) +
  #geom_nodelab() +
#  geom_text2(aes(subset=!isTip, label=node), hjust=-.3) +
  geom_tippoint(aes(color = study)) +
  scale_color_manual(values = cbbPalette)
print(p)

# I want to use the phylogenetic information from the tree to group isolates. So, I'll loop through different
# thresholds for phylogenetic distance and assign "clade" designations based on clustering.

# start with the metadata from above
prevotella_mag_groupings.clades <- prevotella_mag_groupings.format

# filter to the P. copri subtree
prevotella_mag_groupings.clades <- prevotella_mag_groupings.clades %>%
  filter(MAG_ID %in% tree.reroot$tip.label)

# loop through thresholds
for (h in c(0.2, 0.1, 0.05, 0.01, 0.005, 0.004, 0.003, 0.002, 0.001)) {
  # extract phylogenetic distances
  dd <- as.dist(cophenetic.phylo(tree.reroot))
  # consider complete vs single vs average clustering
  psclust <- cutree(as.hclust(agnes(dd, method = "single")), h = h)
  #cliques = levels(factor(psclust))[tapply(psclust, factor(psclust), function(x) {length(x) > 1})]
  cliques <- unique(psclust)
  
  tree_clusters <- data.frame()
  for (i in cliques) {
    tree_clusters <- tree_clusters %>%
      bind_rows(data.frame(cluster = i,
                           tips = names(psclust)[psclust == i]))
  }
  
  prevotella_mag_groupings.clades <- prevotella_mag_groupings.clades %>%
    left_join(tree_clusters %>% 
                #                mutate(tips = gsub(".mod2", "", tips)) %>% 
                rename(!!paste0("cluster_", gsub("\\.", "", as.character(h))) := cluster), 
              by = c("MAG_ID" = "tips"))
}

# starting to look at the output
prevotella_mag_groupings.clades <- prevotella_mag_groupings.clades %>%
  distinct() %>%
  filter(MAG_ID %in% tree.reroot.pc$tip.label)

# Which species/reference genomes are in each cluster? Here, I've decided that a PD of 0.004 
# seems to have useful discriminatory power amount clusters of P. copri. So, which species and
# clade references are common in each?
prevotella_mag_groupings.clades.species <- prevotella_mag_groupings.clades %>%
  group_by(cluster_0004) %>% 
  summarize(n_mags = length(MAG_ID),
            n_species = length(unique(species)),
            species = paste(sort(unique(species)), collapse = ";"))
prevotella_mag_groupings.clades.species

# Where do our WLZ-associated MAGs fall in the PD <= 0.004 clusters?
prevotella_mag_groupings.clades %>%
  filter(MAG_ID %in% c("MAGBg0018", "MAGBg0019"))

# What else is in this cluster?
prevotella_mag_groupings.clades %>%
  filter(cluster_0004 == 1)

# Label the clades - these are my designations based on what references are in each cluster
clades_key <- data.frame(node = c(1,2,3,4,5,6,86,87,88,89,90,91,92),
                        clade = c("A1","A2","A3","B1","B2","sp900313215","Ph","D","Pc1","C1","C2","Pc2","A4"))

# create a new set of metadata that labels the tips
clades_df <- prevotella_mag_groupings.clades %>%
  select(MAG_ID, species, study, strain, cluster_0004) %>%
  left_join(clades_key, by = c("cluster_0004" = "node")) %>%
  add_count(cluster_0004)

# preparing to relabel the tree branches/tips with clade designations
clades <- split(clades_df$MAG_ID, clades_df$clade)

# group the tree tips/branches by the clade designations from the clusters above
tree.reroot.pc <- groupOTU(tree.reroot.pc, clades)

# drop everything except the MAGs/isolates I wan't to show/compare
tree.reroot.pc.trim <- drop.tip(tree.reroot.pc, prevotella_mag_groupings.format %>%
                                  filter(study %in% c("ref_tett", "ref_patric", "technical")) %>%
                                  pull(MAG_ID))
# make sure the grouping was retained
tree.reroot.pc.trim <- groupOTU(tree.reroot.pc.trim, clades)

# view the labeled tree
ggtree(tree.reroot.pc, (aes(color = group)), layout = "rectangular") %<+% prevotella_mag_groupings.format + 
  geom_treescale() +
  geom_tiplab(aes(label = label), size = 2, offset = 0.001, color = "black") +
  geom_tippoint(aes(shape = study)) +
  scale_color_manual(values = cbbPalette)

# save the clade designations/tip labels
write.table(clades_df, paste0(datestring, "_pcopri_clades.txt"), sep = "\t", row.names = FALSE)
