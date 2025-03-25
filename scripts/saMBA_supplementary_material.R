# data wraggling & ploting libraries
library(tidyverse)
library(see) # colorblind-friendly palettes

# Compare south america vs latin america and caribbean --------------------

# read metadata
metadata <- read_tsv(file = "HMC/data/v1.1.0/sample_metadata.tsv") %>% 
  filter(grepl(x = region, pattern = "Latin America and the Caribbean"))

# read HMC count table
ct <- read_csv(file = "HMC/data/v1.1.0/taxonomic_table_v1.1.0.csv",
               # we exclude the first column (with numbers of rows)
               col_select = !c(1)) %>% 
  # we keep only latin american samples
  filter(sample %in% paste0(metadata$project, "_", metadata$srr))


# phylum relative abundances in long format
long_ct_phylum <- ct %>% 
  pivot_longer(!sample,
               values_to = "abundance",
               names_to = "full_taxonomy") %>% 
  # remove absent taxa from each sample
  filter(abundance > 0) %>% 
  # only keep phylum information
  mutate(phylum = word(full_taxonomy, start = 2, end = 2, sep = "[.]"), 
         .keep = "unused") %>% 
  # remove phyla labeled as 'NA'
  filter(phylum != "NA") %>% 
  # Aggregate abundance at phylum level
  summarise(abundance = sum(abundance), .by = c(sample, phylum))

# ID the top phyla across samples
hmc_top_7_phyla <- long_ct_phylum %>% 
  # calculate abundance of the phyla across samples
  summarise(abundance = sum(abundance, na.rm = TRUE), .by = phylum) %>% 
  slice_max(order_by = abundance, n = 7) %>% 
  pull(phylum)

# calculate relative abundance of the top phyla per sample
long_ct_phylum <- long_ct_phylum %>% 
  mutate(phylum = ifelse(phylum %in% hmc_top_7_phyla, yes = phylum, no = "Other")) %>% 
  mutate(rel_abundance = abundance/sum(abundance) * 100 , .by = sample)

# order samples by the abundance of the most abundant phyla
sample_order <- long_ct_phylum %>% 
  filter(phylum == hmc_top_7_phyla[1]) %>% 
  arrange(rel_abundance) %>% 
  pull(sample) %>% 
  unique()


# non-south-american countries
non_southamerican <- c("MX|GT")

# we force the 2 most abundant phyla to be on the top and bottom position
# to make the comparison across samples easier.
hmc_top_7_phyla_levels <- c(hmc_top_7_phyla[1],
                            hmc_top_7_phyla[3:7],
                            "Other", 
                            hmc_top_7_phyla[2])

long_ct_phylum <- metadata %>% 
  select(project, srr, iso) %>% 
  mutate(sample = paste(project, srr, sep = "_"),
         region = ifelse(grepl(x = iso, pattern = non_southamerican),
                         yes = "Central America", 
                         no = "South America")) %>% 
  inner_join(x = long_ct_phylum,
             y = ., 
             by = "sample") %>% 

  mutate(sample = factor(x = sample, levels = sample_order),
         phylum = factor(x = phylum, levels = hmc_top_7_phyla_levels))


# Plot to compare overall composition of microbiomes from
# central vs south american subjects
plot_south_vs_central_america <- long_ct_phylum %>% 
  ggplot(aes(x = sample, y = rel_abundance, fill = phylum)) + 
  geom_col() + 
  labs(x = NULL, y = NULL) +
  facet_wrap(region~., scales = "free") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_okabeito(order = c(1, 3:8, 2)) + 
  
  labs(y = "Relative abundances (%)",
       fill = "Phyla") +
  
  theme_bw() +
  theme(axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 16, color = "black"),
        
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.position = "bottom",
        
        strip.text = element_text(size = 18, color = "black"#, face = "bold"
                                  ),
        strip.background = element_blank(),
        
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()); plot_south_vs_central_america

ggsave(filename = "samba/outputs/sfig1.jpg",
       plot = plot_south_vs_central_america,
       height = 6, width = 14, units = "in")





# Qualitative comparison of saMBA vs HMC ----------------------------------

# Make a safe copy
hmc_ct <- ct

# remove taxa with 0 abundance from HMC
taxa_to_keep_hmc <- (column_to_rownames(hmc_ct, "sample") |> colSums()) > 0
hmc_ct <- hmc_ct[, c(TRUE, taxa_to_keep_hmc)]

# check if any sample has no counts across remaining taxa
any(column_to_rownames(hmc_ct, "sample") |> rowSums() == 0)


# get taxa and samples present in the HMC and match their format in saMBA
hmc_taxa <- colnames(hmc_ct)[2:ncol(hmc_ct)]
hmc_taxa <- gsub(x = hmc_taxa, pattern = "[.]", replacement = ";")

hmc_samples <- gsub(x = ct$sample, pattern = ".*_", replacement = "")

hmc_ct <- hmc_ct %>% 
  mutate(sample = gsub(x = sample, pattern = ".*_", replacement = "")) %>% 
  rename_with(.fn = ~gsub(x = .x, pattern = "[.]", replacement = ";")) 


# 1. read saMBA metadata
# samba_metadata <- read_tsv(file = "samba/temp/samba_metadata.tsv") 

# 2. read saMBA count table
samba_ct <- read_tsv(file = "samba/inputs/clean_genus_count_table_v0.1.0.tsv") %>% 
  # reformat to match the structure of the HMC count table
  column_to_rownames("full_taxonomy") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("sample") %>% 
  as_tibble()

samba_taxa <- colnames(samba_ct)[2:ncol(samba_ct)]
samba_samples <- samba_ct$sample


# Identify shared taxa and samples to subset the count tables
shared_taxa <- c("sample", base::intersect(hmc_taxa, samba_taxa))
shared_samples <- base::intersect(hmc_samples, samba_samples)



# Subset both count tables using shared taxa and samples
samba_ct <- samba_ct %>% 
  filter(sample %in% shared_samples) %>% 
  select(shared_taxa)

hmc_ct <- hmc_ct %>% 
  filter(sample %in% shared_samples) %>% 
  select(shared_taxa)

# samba in long format with total abundance per phyla
long_samba_ct_phyla <- samba_ct %>% 
  pivot_longer(!sample, 
               names_to = "taxa",
               values_to = "abundance") %>% 
  filter(abundance > 0) %>% 
  mutate(phylum = word(string = taxa, start = 2, end = 2, sep = ";")) %>% 
  filter(phylum != "NA") %>% 
  summarise(abundance = sum(abundance), .by = c(sample, phylum))

# HMC in long format with total abundance per phyla
long_hmc_ct_phyla <- hmc_ct %>% 
  pivot_longer(!sample, 
               names_to = "taxa",
               values_to = "abundance") %>% 
  filter(abundance > 0) %>% 
  mutate(phylum = word(string = taxa, start = 2, end = 2, sep = ";")) %>% 
  filter(phylum != "NA") %>% 
  summarise(abundance = sum(abundance), .by = c(sample, phylum))


long_joined_worflows_ct <- inner_join(x = long_samba_ct_phyla,
                                      y = long_hmc_ct_phyla,
                                      by = c("sample", "phylum"), 
                                      suffix = c("_samba", "_hmc")) %>% 
  pivot_longer(cols = starts_with("abundance"),
               names_to = "workflow",
               values_to = "abundance", 
               names_prefix = "abundance_")


# most abundant phyla in HMC
across_wf_top_7_phyla <- long_joined_worflows_ct %>%  
  #filter(workflow == "hmc") %>% 
  summarise(abundance = sum(abundance, na.rm = TRUE), .by = c(phylum)) %>% 
  mutate(rel_abund = abundance/sum(abundance, na.rm = TRUE)) %>% 
  slice_max(order_by = rel_abund, n = 7) %>% 
  pull(phylum)


# we force the 2 most abundant phyla to be on the top and bottom position
# to make the comparison across samples easier.
across_wf_top_7_phyla_levels <- c(across_wf_top_7_phyla[1],
                                  across_wf_top_7_phyla[3:7],
                                  "Other", 
                                  across_wf_top_7_phyla[2])
  

# order samples based on HMC using the order created in previous section
new_sample_order <- sample_order[gsub(x = sample_order, pattern = ".*_", replacement = "") %in% long_hmc_ct_phyla$sample]
new_sample_order <- gsub(x = new_sample_order, pattern = ".*_", replacement = "")
  
long_joined_worflows_ct <- long_joined_worflows_ct %>% 
  mutate(phylum = ifelse(phylum %in% across_wf_top_7_phyla,
                         yes = phylum,
                         no = "Other"),
         phylum = factor(x = phylum, 
                         levels = across_wf_top_7_phyla_levels),
         
         sample = factor(x = sample, levels = new_sample_order),
         
         workflow = ifelse(workflow == "hmc", yes = "HMC", no = "saMBA")) %>% 
  
  mutate(rel_abundance = abundance/sum(abundance) * 100 , .by = c(sample, workflow))

# plot
plot_compare_workflows_qualitative <- long_joined_worflows_ct %>% 
  ggplot(aes(x = sample, y = rel_abundance, fill = phylum)) + 
  geom_col() + 
  labs(x = NULL, y = NULL) +
  facet_wrap(workflow~., scales = "free", ncol = 1, strip.position = "right") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_okabeito(order = c(1, 3:8, 2)) + 
  
  labs(y = "Relative abundances (%)",
       fill = "Phyla") +
  
  theme_bw() +
  theme(axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.position = "bottom",
        
        strip.background.y = element_blank(),
        strip.text.y.right = element_text(size = 18, color = "black", face = "bold"),
        
        panel.spacing = unit(1.5, "lines"),
        
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) ; plot_compare_workflows_qualitative

ggsave(filename = "samba/outputs/sfig3.jpg",
       plot = plot_compare_workflows_qualitative,
       height = 6, width = 14, units = "in")


# Quantitative comparison of saMBA vs HMC ---------------------------------

# Join ct of both workflows in long format to calculate distances
samba_ct_t <- samba_ct %>% 
  mutate(sample = paste0(sample , "_samba")) %>% 
  column_to_rownames("sample") %>% 
  t()

hmc_ct_t <- hmc_ct  %>% 
  mutate(sample = paste0(sample , "_hmc")) %>% 
  column_to_rownames("sample") %>% 
  t() 

wide_joined_workflow_ct <- merge(x = samba_ct_t, y = hmc_ct_t, by = "row.names") %>% 
  column_to_rownames("Row.names")

unwanted_rows <- c("Eukaryota;NA;NA;NA;NA;NA", "NA;NA;NA;NA;NA;NA", "Archaea;NA;NA;NA;NA;NA")
rows_to_keep <- !(rownames(wide_joined_workflow_ct) %in% unwanted_rows)
wide_joined_workflow_ct <- wide_joined_workflow_ct[rows_to_keep, ]

# calculate distances between samples
bray_workflows <- vegan::vegdist(t(wide_joined_workflow_ct), method = "bray") |> as.matrix()

long_df_bray_dists <- bray_workflows %>% 
  as.data.frame(row.names = rownames(bray_workflows)) %>% 
  rownames_to_column("origin") %>% 
  pivot_longer(!origin,
               names_to = "target",
               values_to = "distance") %>% 
  # keep only 1 triangle of the dist matrix
  filter(origin >= target) %>% 
  mutate(origin_wf = gsub(x = origin, pattern = ".*_", replacement = ""),
         target_wf = gsub(x = target, pattern = ".*_", replacement = ""),
         same_workflow = ifelse(origin_wf == target_wf, yes = TRUE, no = FALSE),
       
         origin_sample = gsub(x = origin, pattern = "_.*", replacement = ""),
         target_sample = gsub(x = target, pattern = "_.*", replacement = ""),
         same_sample = ifelse(origin_sample == target_sample, yes = TRUE, no = FALSE),
       
         cases = case_when(same_sample & same_workflow ~ "Same sample\non\nsame workflow",
                           same_sample & !same_workflow ~ "Same sample\non\ndifferent workflow",
                           
                          !same_sample & same_workflow ~ "Different sample\non\nsame workflow",
                          !same_sample & !same_workflow ~ "Different sample\non\ndifferent workflow"),
       
         x_axis = factor(x = cases,
                         levels = c("Different sample\non\ndifferent workflow",
                                    "Different sample\non\nsame workflow",
                                    
                                    "Same sample\non\ndifferent workflow",
                                    "Same sample\non\nsame workflow")))
  

plot_bray_dists_across_workflows <- long_df_bray_dists %>% 
  ggplot(aes(x = x_axis,
             y = distance)) + 
  geom_point(size = 2, shape = 21, alpha = 1/2, 
             position = position_jitter(seed = 01021997)) +
  geom_boxplot(alpha = .75, outlier.shape = NA) + 
  labs(y = "Bray-Curtis dissimilarity",
       x = NULL) + 
  theme_bw() + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16)) #; plot_bray_dists_across_workflows

ggsave(filename = "samba/outputs/sfig4.jpg",
       plot = plot_bray_dists_across_workflows,
       height = 4, width = 6, units = "in")




# PERMANOVA ---------------------------------------------------------------

metadata <- data.frame(id = colnames(wide_joined_workflow_ct),
                       workflow = gsub(pattern = ".*_", replacement = "", x = colnames(wide_joined_workflow_ct)))

wide_joined_workflow_ct_t <- t(wide_joined_workflow_ct)

permanova_bray <- vegan::adonis2(formula = wide_joined_workflow_ct_t ~ workflow,
                                 data = metadata, 
                                 method = "bray", 
                                 permutations = 1e3)

permanova_bray |> broom::tidy()


# PCA ---------------------------------------------------------------------

bray_prcomp_workflows <- bray_workflows %>% 
  prcomp()

bray_prcomp_workflows$x

pc1 <- round((bray_prcomp_workflows$sdev[1] ^ 2) / sum(bray_prcomp_workflows$sdev ^ 2), 3) * 100
pc2 <- round((bray_prcomp_workflows$sdev[2] ^ 2) / sum(bray_prcomp_workflows$sdev ^ 2), 3) * 100

df_bray_prcomp_workflows <- bray_prcomp_workflows$x[, 1:4] %>% as.data.frame()

df_bray_prcomp_workflows <- df_bray_prcomp_workflows %>% 
  rownames_to_column("id") %>% 
  mutate(workflow = gsub(x = id, pattern = ".*_", replacement = ""))


workflow_colors <- c("samba" = "dodgerblue2", "hmc" = "#e05800")


# shuffle rows to improve visibility in PCA
row_order <- sample(x = 1:nrow(df_bray_prcomp_workflows),
                    size = nrow(df_bray_prcomp_workflows), 
                    replace = FALSE)

df_bray_prcomp_workflows <- df_bray_prcomp_workflows[row_order,]

plot_bray_pca <- df_bray_prcomp_workflows %>% 
  ggplot(aes(x = PC1, y = PC2, fill = workflow)) + 
  stat_ellipse(geom = "polygon", alpha = .35, show.legend = FALSE) +
  geom_point(shape = 21, alpha = .65, size = 3) +
  
  labs(fill = "Workflow",
       x = paste0("PC 1 (", pc1, "%)"),
       y = paste0("PC 2 (", pc2, "%)")) +
  
  scale_fill_manual(values = workflow_colors,
                    labels = c("HMC", "saMBA")) +
  coord_fixed() +
  
  guides(fill = guide_legend(override.aes = list(size = 5))) +
  
  theme_bw() + 
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.position = "top",
        ); plot_bray_pca
  
  
ggsave(filename = "samba/outputs/sfig5.jpg",
       plot = plot_bray_pca,
       height = 8, width = 8, units = "in")
