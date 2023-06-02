# configure environment

library(readxl)
library(ggbeeswarm)

setwd("~/Documents/Projects/human_studies/MDCF POC/220812 synthesizing fig 3 and 4/")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#names(cbbPalette) <- c("black", "orange", "sky_blue", "bluish_green", "yellow", "blue", "vermillion", "reddish_purple")
cbbPalette_mod <- c("#E69F00", "#D55E00", "#CC79A7", "#56B4E9", "#0072B2", "#009E73")

# load fgsea enrichments for transcript DE analysis
dat <- read_excel("220812_fig_3_4_enrichments.xlsx") %>%
  arrange(category, NES) %>%
  mutate(pathway = factor(pathway, levels = unique(pathway)))

# pathway enrichment plot - all pathways
ggplot(dat, aes(x = NES, y = pathway)) +
  geom_segment(aes(x = 0, y = pathway, xend = NES, yend = pathway)) +
  geom_vline(xintercept = 0) +
  geom_point(aes(fill = category, size = size, shape = term)) +
  scale_fill_manual(values = cbbPalette) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  facet_grid(. ~ analysis) +
  theme_classic()

# pathway enrichment plot - carbohydrate utilization only
pdf("comparison_substrate.pdf", width = 12)
p <- ggplot(dat %>% filter(category == "Carbohydrate utilization (uptake and catabolism)"), aes(x = NES, y = pathway)) +
  geom_segment(aes(x = 0, y = pathway, xend = NES, yend = pathway)) +
  geom_vline(xintercept = 0) +
  geom_point(aes(fill = category, size = size, shape = term)) +
  scale_fill_manual(values = cbbPalette) +
  scale_shape_manual(values = c(21, 22, 23, 24)) +
  facet_grid(. ~ analysis) +
  theme_classic()
print(p)
dev.off()
#

# load results by diet and quartile
dat.diet <- read_excel("~/Library/CloudStorage/Box-Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/Volcano plot/TableS10_rIntake.xlsx")
dat.quartile <- read_excel("~/Library/CloudStorage/Box-Box/Hibberd_Webber_et_al_MDCF_POC_MAGs/RNAseq/summer 2022 followup analysis/Volcano plot/TableS12_rIntake.xlsx")

# specific targets
diet_path_targets <- c("arabinose utilization", "alpha-arabinooligosaccharides utilization", "fucose utilization")

# volcano plot
pdf("volcano_diet.pdf", height = 4, width = 10)
p <- ggplot(dat.diet, aes(x = log2_FC, y = -1*log10(p_value))) +
  geom_hline(yintercept = -1*log10(0.05)) +
  geom_vline(xintercept = 0) +
  geom_point(data = dat.diet %>% filter(q_value >= 0.05), size = 1, color = "grey") +
  geom_point(data = dat.diet %>% filter(q_value >= 0.05 & pathway %in% diet_path_targets), size = 2, aes(fill = pathway), pch = 21) +
  geom_point(data = dat.diet %>% filter(q_value < 0.05), size = 1, color = "darkgrey") +
  geom_point(data = dat.diet %>% filter(q_value < 0.05 & pathway %in% diet_path_targets), size = 4, aes(fill = pathway), pch = 21) +
  scale_fill_manual(values = cbbPalette) +
  facet_grid(. ~ contrast, scales = "free") +
  theme_classic()
print(p)
dev.off()

# volcano plot for MAGs Bg0018/Bg0019
ggplot(dat.diet, aes(x = log2_FC, y = -1*log10(q_value))) +
  geom_hline(yintercept = -1*log10(0.05)) +
  geom_vline(xintercept = 0) +
  geom_point(data = dat.diet, size = 1, color = "grey") +
  geom_point(data = dat.diet %>% filter(MAG_ID %in% c("MAGBg0019", "MAGBg0018")), aes(fill = MAG_ID), size = 2, color = "black", pch = 21) +
  facet_grid(. ~ contrast, scales = "free") +
  theme_classic()

tally <- dat.diet %>% 
  group_by(MAG_ID) %>%
  count() %>%
  arrange(desc(n))

tally %>% 
  ungroup() %>%
  filter(n > 500) %>%
  summarize(sum(n), length(n))

nrow(dat.diet)

#15246/17088

ggplot(dat.diet, aes(x = log2_FC, y = -1*log10(q_value))) +
  geom_hline(yintercept = -1*log10(0.05)) +
  geom_vline(xintercept = 0) +
  geom_point(data = dat.diet, size = 1, color = "grey") +
  geom_point(data = dat.diet %>% 
               group_by(contrast) %>%
               slice_min(Sn_log10_p_value, n = 50), aes(fill = MAG_ID), size = 2, color = "black", pch = 21) +
  facet_grid(. ~ contrast, scales = "free") +
  theme_classic()

dat.diet %>%
  filter(! MAG_ID %in% c("MAGBg0018", "MAGBg0019")) %>%
  nrow()

dat.diet %>%
  filter(q_value < 0.05 & ! MAG_ID %in% c("MAGBg0018", "MAGBg0019")) %>%
  nrow()

# quartile plotting

# specific targets from GSEA
quartile_path_targets <- c("xylooligosaccharides utilization", 
                           "fructooligosaccharides utilization", 
                           "oligogalacturonate utilization",
                           "galactooligosaccharides utilization",
                           "galactose utilization",
                           "glucuronate utilization",
                           "galacturonate utilization",
                           "alpha-arabinooligosaccharides utilization",
                           "beta-glucosides utilization",
                           "arabinose utilization",
                           "fucose utilization",
                           "melibiose utilization")

# volcano plots
pdf("volcano_quartile.pdf", height = 4, width = 10)
p <- ggplot(dat.quartile, aes(x = log2_FC, y = -1*log10(q_value))) +
  geom_hline(yintercept = -1*log10(0.05)) +
  geom_vline(xintercept = 0) +
  geom_point(data = dat.quartile %>% filter(q_value >= 0.05), size = 1, color = "grey") +
  geom_point(data = dat.quartile %>% filter(q_value >= 0.05 & pathway %in% quartile_path_targets), size = 2, aes(fill = pathway), pch = 21) +
  geom_point(data = dat.quartile %>% filter(q_value < 0.05), size = 1, color = "darkgrey") +
  geom_point(data = dat.quartile %>% filter(q_value < 0.05 & pathway %in% quartile_path_targets), size = 4, aes(fill = pathway), pch = 21) +
#  scale_fill_manual(values = cbbPalette) +
  facet_grid(. ~ coefficent, scales = "free") +
  theme_classic()
print(p)
dev.off()

ggplot(dat.quartile, aes(x = log2_FC, y = -1*log10(q_value), shape = coefficent)) +
  geom_hline(yintercept = -1*log10(0.05)) +
  geom_vline(xintercept = 0) +
  geom_point(data = dat.quartile %>% filter(q_value >= 0.05), size = 1, color = "grey") +
  geom_point(data = dat.quartile %>% filter(q_value >= 0.05 & pathway %in% quartile_path_targets), size = 1, aes(fill = pathway), pch = 21) +
  geom_point(data = dat.quartile %>% filter(q_value < 0.05), size = 2, color = "darkgrey") +
  geom_point(data = dat.quartile %>% filter(q_value < 0.05 & pathway %in% quartile_path_targets), size = 2, aes(fill = pathway), pch = 21) +
  theme_classic()

# volcano, reduced to MAG Bg0019/Bg0018
ggplot(dat.quartile, aes(x = log2_FC, y = -1*log10(q_value))) +
  geom_hline(yintercept = -1*log10(0.05)) +
  geom_vline(xintercept = 0) +
  geom_point(data = dat.quartile, size = 1, color = "grey") +
  geom_point(data = dat.quartile %>% filter(MAG_ID %in% c("MAGBg0019", "MAGBg0018")), aes(fill = MAG_ID), size = 4, color = "black", pch = 21) +
  facet_grid(. ~ coefficent, scales = "free") +
  theme_classic()

# checking enrichment results for both analyses
df <- dat.quartile %>%
#df <- dat.diet %>%
#  filter(contrast == "arm_time") %>%
  filter(coefficent == "quartile_week") %>%
  arrange(sn_log10_p_value)

# setting up gene sets - transcripts assigned to each pathway
sets <- split(df$locus_tag, df$pathway)

# ranks determined by signed significance value
ranks <- df$sn_log10_p_value
names(ranks) <- df$locus_tag

# run GSEA
fgsea.res <- fgsea(sets,
                   ranks,
                   minSize = 10,
                   nPermSimple = 10000) %>%
  arrange(padj)

fgsea.res %>%
  filter(padj < 0.1)

sig_pathways <- fgsea.res %>%
  filter(padj < 0.1) %>%
  pull(pathway)

stats <- ranks
gseaParam = 1
rnk <- rank(-stats)
ord <- order(rnk)
statsAdj <- stats[ord]
statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
statsAdj <- statsAdj/max(abs(statsAdj))

# processing leading edges
pathway_df <- data.frame()

for (i in 1:length(sig_pathways)) {
  target_pathway <- sig_pathways[i]
  pathway <- sets[[target_pathway]]
  leading_edge <- fgsea.res %>%
    filter(pathway == target_pathway) %>%
    pull(leadingEdge) %>%
    unlist()
  
  pathway_df <- pathway_df %>%
    bind_rows(data.frame(pathway_name = target_pathway,
                           locus_tag = pathway,
                           rank = unname(as.vector(na.omit(match(pathway, names(statsAdj)))))) %>%
    mutate(leadingEdge = ifelse(locus_tag %in% leading_edge, TRUE, FALSE))) %>%
    arrange(rank)
  
}

pathway_df.anno <- pathway_df %>%
  left_join(dat.diet %>% 
              select(locus_tag, MAG_ID, genus, species) %>%
              distinct(),
            by = "locus_tag")

pathway_df.anno.filt <- pathway_df.anno %>%
#  filter(pathway_name %in% diet_path_targets)
  filter(pathway_name %in% quartile_path_targets)

pathway_df.anno.filt.bulk <- pathway_df.anno.filt.bulk %>%
  bind_rows(pathway_df.anno.filt %>%
              mutate(contrast = "quartile_week"))

target_mags <- c("MAGBg0001", "MAGBg0006", "MAGBg0010", "MAGBg0018", "MAGBg0019", "MAGBg0023")

# plotting driver transcripts for each enrichment result
pdf("221114_bulk.pdf", height = 14, width = 28)
p <- ggplot(pathway_df.anno.filt.bulk, aes(x = rank, y = pathway_name)) +
  geom_hline(aes(yintercept = pathway_name)) +
  geom_quasirandom(data = pathway_df.anno.filt.bulk %>% 
                     filter(leadingEdge == TRUE & !MAG_ID %in% target_mags), fill = "grey", pch = 21, size = 5, groupOnX = FALSE, varwidth = TRUE, width = 0.2) +
  geom_quasirandom(data = pathway_df.anno.filt.bulk %>% 
                     filter(leadingEdge == FALSE & !MAG_ID %in% target_mags), color = "grey", pch = 16, size = 3, groupOnX = FALSE, varwidth = TRUE, width = 0.2) +
  geom_quasirandom(data = pathway_df.anno.filt.bulk %>% 
                     filter(leadingEdge == FALSE & MAG_ID %in% target_mags), aes(color = MAG_ID), pch = 16, size = 3, groupOnX = FALSE, varwidth = TRUE, width = 0.2) +
  geom_quasirandom(data = pathway_df.anno.filt.bulk %>% 
                     filter(leadingEdge == TRUE & MAG_ID %in% target_mags), aes(fill = MAG_ID), pch = 21, size = 5, groupOnX = FALSE, varwidth = TRUE, width = 0.2) +
  scale_fill_manual(breaks = target_mags, 
                    values = cbbPalette_mod) +
  scale_color_manual(breaks = target_mags, 
                     values = cbbPalette_mod) +
  scale_x_continuous(limits = c(0, length(stats))) +
  theme_classic() +
  facet_grid(. ~ contrast)
print(p)
dev.off()

# a variation of GSEA result plotting
ticksSize = 0.2
gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway_df$rank, returnAllExtremes = TRUE)
bottoms <- gseaRes$bottoms
tops <- gseaRes$tops
n <- length(statsAdj)
xs <- as.vector(rbind(pathway_df$rank - 1, pathway_df$rank))
ys <- as.vector(rbind(bottoms, tops))
toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
diff <- (max(tops) - min(bottoms))/8
x = y = NULL
g <- ggplot(toPlot, aes(x = x, y = y)) + 
#  geom_point(color = "darkblue", size = 0.1) + 
  geom_hline(yintercept = max(tops), colour = "red", linetype = "dashed") + 
  geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") + 
  geom_hline(yintercept = 0, colour = "black") + 
#  geom_line(color = "darkblue") + 
  theme_classic() + 
#  geom_segment(data = pathway_df, mapping = aes(x = rank, y = -diff/2, xend = rank, yend = diff/2), size = ticksSize) + 
  geom_point(data = pathway_df %>% filter(MAG_ID %in% c("MAGBg0018", "MAGBg0019")), mapping = aes(x = rank, y = 0, fill = MAG_ID), size = 4, pch = 24) +
  geom_point(data = pathway_df %>% filter(! MAG_ID %in% c("MAGBg0018", "MAGBg0019")), mapping = aes(x = rank, y = 0, fill = MAG_ID), size = 1, pch = 24) + 
  theme(panel.border = element_blank(), panel.grid.minor = element_blank()) + 
  labs(x = "rank", y = "enrichment score")
g

