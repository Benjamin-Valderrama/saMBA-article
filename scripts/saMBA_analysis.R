#' This is the analysis of saMBA that goes to the paper

# data wraggling & ploting libraries
library(tidyverse)
library(scales)

# geography libraries
library(sf)
library(spData)


# read in saMBA files -----------------------------------------------------

track_reads <- read_tsv(file = "samba/inputs/track_reads_v0.1.0.tsv")
filereports <- read_tsv(file = "samba/inputs/insdc_file_reports_v0.1.0.tsv")

pipeline_metadata <- inner_join(x = filereports, 
                                y = track_reads,
                                by = "run_accession")

readxl::read_xlsx(path = "samba/projects_metadata.xlsx", sheet = "full_table") %>% 
  janitor::clean_names() %>% 
  write_tsv(file = "samba/projects_metadata.tsv")

projects_metadata <- read_tsv(file = "samba/projects_metadata.tsv") %>% 
  select(study_accession = prj_number, country, sequencing_technology, region, industrialised, disease_boolean)


# One project has data from multiple nationalities
multicountry_project <- projects_metadata %>% 
  filter(country == "multiple") %>% 
  pull(study_accession)

multicountry_project_df <- read_tsv(file = paste0("manual_curation_of_samples/processed/filereport_", multicountry_project, ".txt")) %>% 
  # Add capitals to the first letter of country
  mutate(country = stringr::str_to_title(country))


multiindustrial_project <- projects_metadata %>%
  filter(industrialised == "multiple") %>%
  pull(study_accession)


multiindustrial_project_df <- read_tsv(file = paste0("samba/inputs/INSDC_filereport_", multiindustrial_project,".txt")) %>% 
  # select relevant columns
  select(study_accession, run_accession, original_library_layout = library_layout, submitted_ftp) %>% 
  # do cleanup to identify samples with "BEL", comming from subjects living in industrialised settings
  mutate(industrialised = word(string = submitted_ftp, start = -1, end = -1, sep = "/"),
         industrialised = ifelse(grepl(x = industrialised, pattern = "^BEL"),
                                 yes = "Industrialised",
                                 no = "Non-Industrialised")) %>% 
  select(!submitted_ftp)


# add projects metadata and samples metadata
metadata <- inner_join(x = pipeline_metadata,
                       y = projects_metadata,
                       by = "study_accession") %>% 
  
  # add specific country to samples from multi-country projects
  mutate(country = ifelse(country == "multiple", NA, country)) %>% 
  full_join(x = .,
            y = multicountry_project_df,
            by = c("study_accession", "run_accession", "original_library_layout")) %>% 
  mutate(country = coalesce(country.x, country.y)) %>% 
  select(!c(country.x, country.y)) %>% 
  
  # add industrilization level to samples from multiindustrialised projects
  mutate(industrialised = ifelse(industrialised == "multiple", NA, industrialised)) %>% 
  full_join(x = .,
            y = multiindustrial_project_df,
            by = c("study_accession", "run_accession", "original_library_layout")) %>% 
  mutate(industrialised = coalesce(industrialised.x, industrialised.y)) %>% 
  select(!c(industrialised.x, industrialised.y)) %>% 

  # fix name of guyana
  mutate(country = ifelse(grepl(x = country, pattern = "Guiana"),
                          yes = "Guyana",
                          no = country))

# export metadata to be used later when comparing saMBA and the HMC
write_tsv(x = metadata, file = "samba/temp/samba_metadata.tsv")
  
# count table (after cleanup)
ct <- read_tsv(file = "samba/inputs/clean_genus_count_table_v0.1.0.tsv")

clean_samples <- colnames(ct)[2:ncol(ct)]

clean_metadata <- metadata %>% 
  filter(run_accession %in% clean_samples)




# Prepare world data for plotting maps ------------------------------------

south_america <- spData::world %>% 
  filter(continent == "South America" | name_long == "France") %>% 
  filter(name_long != "Falkland Islands") %>% 
  select(1:6, 11)

# Remove europe france from the map and only keep french guiana
south_america[south_america$name_long == "France", ][["geom"]][[1]][[3]] <- NULL
south_america[south_america$name_long == "France", ][["geom"]][[1]][[2]] <- NULL



# Panels figure 1 ---------------------------------------------------------

# number of studies per country
plot_studies_per_country <- clean_metadata %>% 
  summarise(n_studies = n_distinct(study_accession), .by = country) %>% 
  full_join(x = .,
            y = south_america,
            by = join_by(country == name_long)) %>% 
  ggplot() +
  geom_sf(aes(fill = n_studies, geometry = geom), 
          color = "black", linewidth = 0.75) +

  scale_fill_steps(low = "#f5d8c0", high = "#D55E00",
                   limits = c(1, 9),
                   na.value = "grey65", 
                   breaks = c(2, 4, 5, 9), 
                   labels = c("1 or 2",
                              "3 or 4",
                              "5 or 6",
                              "7 or more"),
                   guide = "legend") +
  guides(fill = guide_legend(override.aes = list(linewidth = 0.5))) +
  labs(fill = "Studies\nincluded") +
  theme_void() +
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)); plot_studies_per_country


ggsave(file = "samba/outputs/fig1c.jpg", 
       plot = plot_studies_per_country,
       height = 4, width = 6, units = "in")


# number of sequenced samples per country
geo_samples_per_country <- clean_metadata %>% 
  count(country, name = "n_samples") %>% 
  full_join(x = .,
            y = south_america,
            by = join_by(country == name_long)) %>% 
  ggplot() +
  geom_sf(aes(fill = n_samples, geometry = geom), 
          color = "black", linewidth = 0.75) +
  
  scale_fill_gradient(low = "#befeff", high = "#0072B2",
                      na.value = "grey65",
                      limits = c(0, 1e3),  
                      breaks = seq(250, 999, by = 250)
                      ) +
  
  labs(fill = "Samples\nincluded") +
  theme_void() +
  theme(legend.title = element_text(size = 16),
        legend.text = element_text(size = 14)); geo_samples_per_country 


ggsave(file = "samba/outputs/fig1a.jpg", plot = last_plot(),
       height = 4, width = 6, units = "in")




# Piechart: percentage of projects per study size category
fig1b <- clean_metadata %>% 
  # samples per study
  count(study_accession) %>% 
  mutate(category = case_when(n < 50 ~ "Less than 50 samples",
                              (50 <= n) & (n <= 100) ~ "Between 50 and 100 samples",
                              100 < n ~ "More than 100 samples"),
         category = factor(x = category,
                           levels = c("More than 100 samples",
                                      "Between 50 and 100 samples",
                                      "Less than 50 samples"
                           ))) %>% 
  # number of studies of each size category created above
  count(category) %>% 
  
  # Make percent text to plot and Y for positioning of text
  mutate(percent = round(n/sum(n) * 100 , 0),
         percent_text = paste0(percent, "%"),
         y = cumsum(percent) - 0.2*percent) %>% 
  
  ggplot(aes(x = "", y = percent, fill = category)) +
  
  geom_col() +
  coord_polar("y", start = 0) +
  
  geom_text(aes(y = y - 0.45 * percent, label = rev(percent_text)),
            size = 7) +
  
  scale_fill_manual(breaks = c("Less than 50 samples",
                               "Between 50 and 100 samples",
                               "More than 100 samples"),
                    values = c("#c9261b", "grey", "grey65")) +
  
  
  guides(fill = guide_legend(ncol = 1)) +
  
  theme_void() + 
  theme(legend.title = element_blank(),
        legend.margin = margin(t = 20, b = -20),
        legend.text = element_text(size = 14),
        legend.position = "top"); fig1b

ggsave(filename = "samba/outputs/fig1b.jpg", 
       height = 4, width = 6, units = "in")



# Distribution of non-chimeric reads
fig1d <- clean_metadata %>% 
  ggplot(aes(x = nochim)) + 
  geom_histogram(fill = "grey40") +
  geom_vline(xintercept = median(metadata$nochim), 
             linetype = "dashed", 
             linewidth = 1.25,
             color = "#c9261b") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 3),
                labels = trans_format("log10", math_format(10^.x))) +
  
  labs(x = "Non-chimeric reads",
       y = "Number of samples") +
  
  scale_y_continuous(expand = c(0,0)) +
  cowplot::theme_cowplot() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 15))); fig1d

ggsave(filename = "samba/outputs/fig1d.jpg", 
       height = 4, width = 6, units = "in")




# Panels figure 2 ---------------------------------------------------------

# calculate alpha diversity metrics
alpha_divs <- ct %>% 
  column_to_rownames("full_taxonomy") %>% 
  microbiome::alpha(index = c("chao1", "shannon")) %>% 
  rownames_to_column("run_accession")

# chao1
chao1_distribution <- alpha_divs %>% 
  ggplot(aes(x = chao1)) + 
  geom_histogram(fill = "grey40", binwidth = 10) +
  geom_vline(xintercept = median(alpha_divs$chao1), 
             linetype = "dashed", 
             linewidth = 1.25,
             color = "#c9261b") +
  
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  
  labs(x = "Chao1",
       y = "Number of samples") +
  
  cowplot::theme_cowplot() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 15))); chao1_distribution

ggsave(filename = "samba/outputs/fig2b.jpg", 
       plot = chao1_distribution,
       height = 4, width = 6, units = "in")



# shannon
shannon_distribution <- alpha_divs %>% 
  ggplot(aes(x = diversity_shannon)) + 
  geom_histogram(fill = "grey40") +
  geom_vline(xintercept = median(alpha_divs$diversity_shannon), 
             linetype = "dashed", 
             linewidth = 1.25,
             color = "#c9261b") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  
  labs(x = "Shannon diversity",
       y = "Number of samples") +
  
  cowplot::theme_cowplot() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 15))); shannon_distribution

ggsave(filename = "samba/outputs/fig2c.jpg", 
       plot = shannon_distribution,
       height = 4, width = 6, units = "in")



# Venn diagram ------------------------------------------------------------

#' The venn diagram is a panel of figure 2. However, there are many
#' steps in the generation of the plot so I made an independent section

#' The count table shared in the Abdill et al 2025 paper (PMID: 39848248)
#' that can be accessible through https://zenodo.org/records/13733642 ,
#' which is the same source from which the R package they developed 
#' retrieves the count table from, as explained in the associated
#' github repository (https://github.com/blekhmanlab/MicroBioMap)

#' Thus, since I will work with MicroBioMap's 'raw' count table, I will
#' also use raw saMBA's count table.

# Read the saMBA raw count table in
raw_ct_samba <- read_tsv(file = "samba/inputs/genus_count_table_v0.1.0.tsv")

# List of samples in saMBA 
samba_samples_full_id <- raw_ct_samba %>% 
  pivot_longer(!full_taxonomy, 
               names_to = "run_accession",
               values_to = "abundance") %>% 
  inner_join(x = .,
             y = metadata,
             by = "run_accession") %>% 
  mutate(full_id = paste0(study_accession, "_", run_accession)) %>% 
  pull(full_id) %>% 
  unique()


# Read MicroBioMap's raw count table in
microbiomap <- read_csv(file = "HMC/data/v1.1.0/taxonomic_table_v1.1.0.csv",
                        col_select = !c(1)) %>% 
  # Keep only the samples that are present in saMBA
  filter(sample %in% samba_samples_full_id)

# Remove taxa not present in the remaining samples
to_keep <- colSums(microbiomap[, 2:ncol(microbiomap)]) > 0
microbiomap <- microbiomap[,c(TRUE, to_keep)]


# Make a list of taxa present on each resource
microbiomap_taxa <- colnames(microbiomap)[2:ncol(microbiomap)]
samba_taxa <- gsub(x = raw_ct_samba$full_taxonomy, pattern = ";", replacement = ".")


venn_df <- data.frame(workflow = c("samba", 
                                   "shared",
                                   "microbiomap"),
                      taxa_count = c(sum(!samba_taxa %in% microbiomap_taxa),
                                     length(base::intersect(microbiomap_taxa, samba_taxa)),
                                     sum(!microbiomap_taxa %in% samba_taxa)
                                     )
                      ) %>% 
  mutate(percentage = round(taxa_count/sum(taxa_count) * 100, 1),
         percentage = paste0("[", percentage, "%", "]"))

venn_df


# venn diagram from scratch
venn_samba_microbiomap <- ggplot(data = filter(venn_df, workflow != "shared"),
         aes(x = c(-0.25, 0.25), 
         y = 0, 
         fill = workflow,
         alpha = workflow,
         label = workflow)) +
  
  geom_point(shape = 21, size = c(75, 50), stroke = 1, show.legend = FALSE) + 
  
  # stroke arround the circles
  annotate(geom = "point",
           x = c(-0.25, 0.25), 
           y = 0,
           color = "black",
           shape = 21,
           size = c(75, 50),
           stroke = 1.5
  ) +
  
  # titles
  annotate(geom = "text", 
           x = c(-0.7, 0.8), 
           y = c(0.175, 0.15),
           size = c(7, 5),
           fontface = "bold",
           label = c("saMBA", "Previous largest\narchive of South\nAmerican\nmicrobiomes")) +
  
  # counts
  annotate(geom = "text",
           x = c(-0.35, 0.08, 0.4), 
           y = 0.01,
           size = 6,
           fontface = "bold",
           label = venn_df$taxa_count) +
  
  # percentage
  annotate(geom = "text",
           x = c(-0.33, 0.08, 0.4), 
           y = -0.03,
           size = 3.5,
           hjust = 0.5,
           label = venn_df$percentage) +
  
  
  scale_alpha_discrete(breaks = c("samba", "microbiomap"),
                       range = c(0.7, 1)) +
  
  coord_cartesian(xlim = c(-1.2, 1.2),
                  ylim = c(-0.25, 0.25), expand = FALSE) +
  
  scale_fill_manual(values = rev(c("dodgerblue2", 
                                   "#e05800"))) +
  theme_void() ; venn_samba_microbiomap

ggsave(filename = "samba/outputs/fig2a.jpg", plot = last_plot(),
       height = 4, width = 6, units = "in")






# Subsampling by country curves -------------------------------------------

#' Create a function that takes:
#' 1) a dataframe with all samples and taxa found in a country
#' 2) a vector with the number of samples (depth) at which to perform the subsampling for each country
#' 3) the number of iterations at which the subsampling at each depth has to be performed

iterate_subsampling <- function(df, depths, iterations){
  
  # Move sample name (i.e., run accession to rownames)
  df <- column_to_rownames(.data = df, var = "run_accession")
  
  # Get number of samples available per country
  n_samples_of_the_country <- nrow(df)
  
  # Objects to store final results
  mean_number_of_taxa_for_each_depths <- numeric()
  sd_number_of_taxa_for_each_depths <- numeric()
  
  
  for (depth in depths){
    
    # Store number of unique taxa across all iterations at each subsampling depth
    number_of_taxa_for_each_iterations_at_fixed_depth <- numeric()
    
    
    #' *Subsampling* requires that the number of samples of the country is
    #' bigger than the number of samples used in the subset
    if(depth > n_samples_of_the_country){ next } else {
      
      for (iteration in 1:iterations) {
        #print(paste0("iteration = ", iteration))
        set.seed(iteration)
        
        
        # Subsample the df
        indices <- sample(1:n_samples_of_the_country, size = depth)
        subsampled_df <- df[indices, ]
        
        # Remove taxa that is not present in the subsampled df
        taxa_to_keep <- colSums(subsampled_df) > 0
        clean_df <- subsampled_df[, taxa_to_keep]
        
        
        # Calculate number of unique taxa in the subset 
        number_unique_taxa_this_iteration <- clean_df %>% 
          colnames() %>% 
          n_distinct()
        
        # Add the number of 
        number_of_taxa_for_each_iterations_at_fixed_depth <- 
          c(number_of_taxa_for_each_iterations_at_fixed_depth,
            number_unique_taxa_this_iteration)
      }
      
      #' After all iterations of that depth have been calculated, get 
      #' mean and SD of the number of unique taxa at that depth
      
      mean_number_of_taxa_at_this_depth <- mean(number_of_taxa_for_each_iterations_at_fixed_depth)
      sd_number_of_taxa_at_this_depth <- sd(number_of_taxa_for_each_iterations_at_fixed_depth)
      
    }
    
    # Add the mean and sd of number of taxa at the depth just analysed
    mean_number_of_taxa_for_each_depths <- c(mean_number_of_taxa_for_each_depths,
                                             mean_number_of_taxa_at_this_depth)
    
    sd_number_of_taxa_for_each_depths <- c(sd_number_of_taxa_for_each_depths,
                                           sd_number_of_taxa_at_this_depth)
  }
  
  # Return final dataframe
  results <- data.frame(depth = depths[depths <= n_samples_of_the_country],
                        mean_unique_taxa = mean_number_of_taxa_for_each_depths,
                        sd_unique_taxa = sd_number_of_taxa_for_each_depths)
  
  return(results)
  
}


# Double check that no taxonomy has 0 counts
sum(rowSums(column_to_rownames(ct, "full_taxonomy")) == 0)

#' Check number of samples per country to curate a list of sampling
#' depths used in the rarefaction-like approach
clean_metadata %>% 
  count(country, name = "n_samples")

sampling_depths <- c(5, 10, 15, 20, 30, 50, 75, 
                     seq(from = 100, to = 900, by = 50),
                     seq(from = 1e3, to = 3*1e3, by = 250))



# Make a binary count table with country info 
binary_ct_long <- ct %>% 
  pivot_longer(!full_taxonomy,
               names_to = "run_accession",
               values_to = "presence") %>% 
  # make counts to binary (presence/abscense)
  mutate(presence = ifelse(presence > 0, yes = 1, no = 0)) %>% 
  # keep only present taxa
  filter(as.logical(presence)) %>% 
  inner_join(x = .,
             y = clean_metadata,
             by = "run_accession") %>% 
  select(full_taxonomy, run_accession, presence, country, industrialised, disease_boolean)


# Bind the per-country binary count table to the continent-wide count table
binary_ct_long <- rbind(binary_ct_long,
                        mutate(binary_ct_long, country = "South-america")) 



# Iterate the subsampling process at increasing depths to each country
iterations <- 1e3

subsampled_countries <- binary_ct_long %>% 
  # we don't care about industrialization status here
  select(!industrialised) %>% 
  nest(.by = c(country)) %>% 
  mutate(data = map(.x = data,
                    .f = ~pivot_wider(data = .x, 
                                      names_from = "full_taxonomy", 
                                      values_from = "presence", 
                                      values_fill = 0))) %>% 
  mutate(rarefaction = map(.x = data,
                           .f = ~iterate_subsampling(df = .x,
                                                     depths = sampling_depths, 
                                                     iterations = iterations)))


subsampled_df <- subsampled_countries %>% 
  select(!data) %>% 
  unnest(rarefaction)


write_tsv(x = subsampled_df, file = "samba/temp/subsampled_counties.tsv")

colors_subsample <- c(RColorBrewer::brewer.pal(8, "Dark2")[1:7],
                      "#0072B2",
                      "#000000",
                      "#D55E00")


# Force a 0 into the data
subsampled_df_0_depth <- subsampled_df %>% 
  filter(depth == min(depth), .by = country) %>% 
  mutate(across(where(is.numeric), ~ 0))

subsampled_df <- rbind(subsampled_df_0_depth, subsampled_df)


subsampling_plot <- subsampled_df %>% 
  ggplot(aes(x = depth, y = mean_unique_taxa, group = country, color = country)) +
  
  geom_ribbon(data = filter(.data = subsampled_df, country == "South-america"),
              aes(ymin = mean_unique_taxa + sd_unique_taxa,
                  ymax = mean_unique_taxa - sd_unique_taxa),
              fill = "grey75",
              alpha = 0.25) +
  
  geom_point() +
  geom_line() +
  
  labs(x = "Samples",
       y = "Novel taxa", fill = "Country", color = "Country") +
  
  scale_color_manual(values = colors_subsample) +
  
  scale_y_continuous(limits = c(0, 1.4e3),
                     breaks = seq(from = 0, to = 1600, by = 400),
                     expand = c(0,0)) +
  
  scale_x_continuous(limits = c(-50, 3.1*1e3),
                     breaks = seq(from = 0, to = 3.1*1e3, by = 1e3),
                     expand = c(0,0)) +
  
  guides(colour = guide_legend(override.aes = list(size=3.5,
                                                   linewidth = 1),
                               nrow = 5)) +
  
  theme_bw() + 
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 14, color = "black"),
        legend.position = "inside",
        legend.direction = "horizontal",
        legend.text = element_text(size = 14),
        legend.background = element_rect(fill = NA),
        
        legend.position.inside = c(0.6, 0.3)); subsampling_plot

ggsave(filename = "samba/outputs/fig3a.jpg",
       plot = subsampling_plot,
       height = 4, width = 6, units = "in")



# Subsampling by industrilised curves -------------------------------------

# from the object created in previous section...
binary_ct_long_industrialised <- binary_ct_long %>% 
  # remove continent-wide observations
  filter(country != "South-america") %>% 
  # remove samples where we lack information about industrialized context
  filter(industrialised == "Non-Industrialised") %>% 
  #filter(!grepl(x = industrialised, pattern = "Prob_")) %>% 
  #mutate(industrialised = gsub(x = industrialised, pattern = ".*_", replacement = "")) %>% 
  mutate(disease = ifelse(disease_boolean, "disease", "no disease"), .keep = "unused") %>% 
  mutate(industrial_disease = paste0(industrialised, ", ", disease))


# Proportion of studies studying subjects with diseases either in 
# industrialised or non-industrialised subjects
binary_ct_long_industrialised %>% 
  select(!c(full_taxonomy)) %>% 
  unique() %>% 
  count(industrialised, disease) %>% 
  mutate(prop = n / sum(n) * 100, .by = industrialised)
# two thirds of non-industrialised microbiomes are comming from subjects
# with a disease

# Check same as above by bioproject
binary_ct_long_industrialised %>% 
  select(!c(full_taxonomy)) %>% 
  unique() %>% 
  inner_join(x = .,
            y = select(clean_metadata, run_accession, study_accession)) %>% 
  count(study_accession, industrialised, disease)


  
# Iterate the subsampling process at increasing depths to each country
iterations <- 1e3

subsampled_industrial <- binary_ct_long_industrialised %>% 
  # we don't care about the country samples come from here
  select(!c(country, industrialised, disease)) %>% 
  nest(.by = c(industrial_disease)) %>% 
  mutate(data = map(.x = data,
                    .f = ~pivot_wider(data = .x, 
                                      names_from = "full_taxonomy", 
                                      values_from = "presence", 
                                      values_fill = 0))) %>% 
  mutate(rarefaction = map(.x = data,
                           .f = ~iterate_subsampling(df = .x,
                                                     depths = sampling_depths, 
                                                     iterations = iterations)))


subsampled_industrial_df <- subsampled_industrial %>% 
  select(!data) %>% 
  unnest(rarefaction)

write_tsv(x = subsampled_industrial_df, file = "samba/temp/subsampled_industrial.tsv")


# Force a 0 into the data
subsampled_industrial_df_0_depth <- subsampled_industrial_df %>% 
  filter(depth == min(depth), .by = industrial_disease) %>% 
  mutate(across(where(is.numeric), ~ 0))

subsampled_industrial_df <- rbind(subsampled_industrial_df_0_depth, 
                                  subsampled_industrial_df)


industrial_colors <- c("Non-Industrialised, disease" = "grey45",
                       "Non-Industrialised, no disease" = "#009E73"
                       )

subsampling_plot <- subsampled_industrial_df %>% 
  ggplot(aes(x = depth, y = mean_unique_taxa,
             color = industrial_disease,
             )) +

  geom_point(size = 2.5) +
  geom_line(linewidth = 1.25) +
  
  labs(x = "Samples",
       y = "Novel taxa",
       color = NULL
       ) +
  
  scale_color_manual(values = industrial_colors) +
  
  scale_y_continuous(limits = c(0, 1.1e3),
                     breaks = seq(from = 0, to = 1e3, by = 200),
                     expand = c(0,0)) +
  
  scale_x_continuous(limits = c(-50, 1.6*1e3),
                     breaks = seq(from = 0, to = 1.25*1e3, by = 5e2),
                     expand = c(0,0)) +
  
  guides(colour = guide_legend(override.aes = list(size=3.5,
                                                   linewidth = 1),
                               nrow = 2)
         ) +
  
  theme_bw() + 
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 14, color = "black"),
        
        legend.position = "inside",
        legend.direction = "horizontal",
        legend.text = element_text(size = 13),
        legend.background = element_rect(fill = NA),
        
        legend.position.inside = c(0.55, 0.15)); subsampling_plot

ggsave(filename = "samba/outputs/fig3d.jpg",
       plot = subsampling_plot,
       height = 4, width = 4, units = "in")


# Representation index per country ----------------------------------------

#' Local representation index is a metric built using the population
#' at each country in 2022 and the number of gut microbiome samples
#' from each of those countries.

# Population data for south american countries
south_american_countries <- south_america$name_long

south_american_pop <- read_csv(file = "samba/global_pop_data/world_bank_population_estimates.csv", skip = 3) %>% 
  janitor::clean_names() %>%
  # Venezuela appears as 'Venezuela, BR', change it to just 'Venezuela'
  mutate(country_name = case_when(grepl(x = country_name, pattern = "Venezuela") ~ "Venezuela",
                                  T ~ country_name)) %>%
  # Get population for south american countries
  filter(country_name %in% south_american_countries) %>%
  select(name_long = country_name, pop_2022 = x2022) %>%
  # Force france to pop = NA (wich represents french guyana) %>% 
  mutate(pop_2022 = ifelse(name_long == "France", NA, pop_2022)) %>% 
  mutate(total_sa_pop = sum(pop_2022, na.rm = TRUE),
         pop_share = round((pop_2022/total_sa_pop) * 100, 2)) %>%
  arrange(desc(pop_share)) %>% 
  select(!total_sa_pop)



south_american_samples <- clean_metadata %>% 
  # Count samples per country
  count(country, name = "n_samples") %>%  
  # Calculate the share of samples per country
  mutate(total_samples = sum(n_samples, na.rm = TRUE),
         samples_share = round((n_samples/total_samples) * 100, 2)) %>% 
  select(!total_samples)


# Join data of per-country population and sample share
continent_representation_df <- full_join(x = south_american_samples,
            y = south_american_pop,
            by = join_by(country == name_long)) %>% 

  # Calculate representation index as outlined in PMID: 35167588
  mutate(higher_samples_share = ifelse(samples_share > pop_share, TRUE, FALSE)) %>% 
    
  mutate(representation_index = ifelse(higher_samples_share, 
                                       yes = samples_share/pop_share,
                                       no = pop_share/samples_share * -1)) %>% 
  full_join(x = .,
            y = south_america,
            by = join_by(country == name_long))


geo_representation_index <- continent_representation_df %>% 
  ggplot() +
  geom_sf(aes(fill = representation_index, geometry = geom),
          color = "black", linewidth = 0.75) +
  
  scale_fill_gradient2(low = "tomato3", high = "dodgerblue3", #"#0072B2"
                       breaks = c(3.0, -8.6),
                       labels = c("Overrepresented", "Underrepresented")
                       ) +

  labs(fill = "Local\nrepresentation\nindex") +
  theme_void() +
  theme(text = element_text(size = 18)); geo_representation_index


ggsave(file = "samba/outputs/fig3b.jpg", 
       plot = geo_representation_index,
       height = 4, width = 6, units = "in")







# Jaccard distances within countries --------------------------------------

binary_ct <- ct %>% 
  #filter(!grepl(x = full_taxonomy, pattern = "Archaea")) %>% 
  mutate(across(where(is.double), .fn = ~ifelse(.x > 0, yes = 1, no = 0))) %>% 
  column_to_rownames("full_taxonomy") %>% 
  t()


jaccard_distance <- vegan::vegdist(x = binary_ct, 
                                   method = "jaccard", 
                                   diag = TRUE, upper = TRUE) %>% 
  as.matrix()


distance_df <- jaccard_distance %>% 
  reshape2::melt(varnames = c("from", "to")) %>% 
  as_tibble() %>% 
  # add metada to column `from`
  inner_join(x = .,
             y = select(.data = clean_metadata, 
                        run_accession,
                        study_accession_from = study_accession, 
                        country_from = country),
             by = join_by(from == run_accession)) %>% 
  # add metadata to column `to`
  inner_join(x = .,
             y = select(.data = clean_metadata, 
                        run_accession,
                        study_accession_to = study_accession, 
                        country_to = country),
             by = join_by(to == run_accession)) %>% 
  select(from, study_accession_from, country_from, to, study_accession_to, country_to, dist = value)

clean_distance_df <- distance_df %>% 
  filter(country_from == country_to) %>% 
  filter(dist > 0)


country_sorted <- clean_distance_df %>% 
  summarise(median_dist = median(dist), .by = country_from) %>% 
  arrange(desc(median_dist)) %>% 
  pull(country_from)



layer_df <- layer_data(geo_representation_index)
representation_index_colors <- cbind.data.frame(country = continent_representation_df$country,
                                                fill = layer_df$fill) %>% 
  # remove countries filled with grey as they don't have data
  filter(fill != "grey50")
  
colors_jaccard <- representation_index_colors$fill
names(colors_jaccard) <- representation_index_colors$country

plot_jaccard_per_country <- clean_distance_df %>% 
  mutate(country_from = factor(country_from, 
                               levels = country_sorted)) %>% 
  ggplot(aes(x = country_from, y = dist, fill = country_from)) +
  geom_jitter(alpha = 0.25, size = 2, show.legend = FALSE) + 
  geom_boxplot(outliers = FALSE, show.legend = FALSE,
               #alpha = 1
               ) + 
  
  scale_fill_manual(values = colors_jaccard) +
  
  labs(y = "Jaccard distances", x = NULL) +
  
  theme_bw() +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20))#; plot_jaccard_per_country


ggsave(filename = "samba/outputs/fig3c.jpg",
       plot = plot_jaccard_per_country,
       height = 6, width = 14, units = "in")
