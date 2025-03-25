
library(tidyverse)

#' The following projects are projects only included in saMBA


# (1/4) PRJEB36789 --------------------------------------------------------

read_tsv(file = "manual_curation_of_samples/raw_filereport/filereport_PRJEB36789.txt") %>% 
  select(study_accession, run_accession, original_library_layout = library_layout, sample_alias) %>% 
  filter(grepl(x = sample_alias, pattern = "_AR-|_CH-")) %>% 
  mutate(country = ifelse(grepl(x = sample_alias, pattern = "_AR-"), yes = "argentina", no = "chile")) %>% 
  select(!sample_alias) %>% 
  write_tsv(file = "manual_curation_of_samples/processed/filereport_PRJEB36789.txt")



# (2/4) PRJNA300541 -------------------------------------------------------

read_tsv(file = "manual_curation_of_samples/raw_filereport/filereport_PRJNA300541.txt") %>% 
  select(study_accession, run_accession,
         original_library_layout = library_layout, library_strategy, experiment_title) %>%
  filter(library_strategy == "AMPLICON") %>% 
  filter(!grepl(x = experiment_title, pattern = "Salvador|[Aa]nimal") & grepl(x = experiment_title, pattern = "fecal")) %>% 
  mutate(country = "peru") %>% 
  select(!c(library_strategy, experiment_title)) %>% 
  write_tsv(file = "manual_curation_of_samples/processed/filereport_PRJNA300541.txt")



# (3/4) PRJEB55648 --------------------------------------------------------

read_tsv(file = "manual_curation_of_samples/raw_filereport/filereport_PRJEB55648.txt") %>% 
  select(study_accession, run_accession,
         original_library_layout = library_layout, library_strategy, scientific_name, read_count) %>% 
  filter(grepl(x = scientific_name, pattern = "feces") & read_count > 0 & library_strategy == "AMPLICON") %>% 
  mutate(country = "venezuela") %>% 
  select(!c(library_strategy, scientific_name, read_count)) %>% 
  write_tsv(file = "manual_curation_of_samples/processed/filereport_PRJEB55648.txt")



# (4/4) PRJEB31199 --------------------------------------------------------

read_tsv(file = "manual_curation_of_samples/raw_filereport/filereport_PRJEB31199.txt") %>% 
  select(study_accession, run_accession, original_library_layout = library_layout, scientific_name) %>% 
  filter(grepl(x = scientific_name, pattern = "gut")) %>% 
  mutate(country = "bolivia") %>% 
  select(!scientific_name) %>% 
  write_tsv(file = "manual_curation_of_samples/processed/filereport_PRJEB31199.txt")









#' The following projects are projects also included in microbiomap


# (1/5) PRJNA421371 -------------------------------------------------------

read_tsv(file = "manual_curation_of_samples/raw_filereport/filereport_PRJNA421371.txt") %>%
  select(study_accession, run_accession, original_library_layout = library_layout, scientific_name) %>% 
  # pull(scientific_name) %>% unique()
  mutate(country = "brazil") %>% 
  select(!scientific_name) %>% 
  write_tsv(file = "manual_curation_of_samples/processed/filereport_PRJNA421371.txt")



# (2/5) PRJEB11419 --------------------------------------------------------

#' This is the American Gut Project (AGP), which is stored in Qiita.
#' I made an account in the website and download the metadata of all the samples
#' in the AGP.

# read metadata of AGP in qiita
qiita <- read_tsv(file = "manual_curation_of_samples/raw_filereport/qiita_PRJEB11419.txt")

south_american_countries_in_qiita <- c("Argentina","Bolivia","Brazil","Colombia","Paraguay")

# get the samples id from qiita to extract those samples from ENA
southamerican_samples_in_qiita <- qiita %>% 
  mutate(stable_country = case_when(country_residence == "not provided" ~ NA,
                                    country_of_birth == country_residence ~ TRUE,
                                    TRUE ~ FALSE)) %>% 
  filter(stable_country) %>% 
  filter(country_residence %in% south_american_countries_in_qiita) %>% 
  filter(grepl(x = body_site, pattern = "feces|faeces")) %>% 
  select(sample_name, country_residence)


# read the file report of the AGP in ENA
read_tsv(file = "manual_curation_of_samples/raw_filereport/filereport_PRJEB11419.txt") %>% 
  filter(sample_title %in% southamerican_samples_in_qiita$sample_name) %>% 
  select(study_accession, run_accession, original_library_layout = library_layout, sample_title) %>% 
  # add the country to each sample
  inner_join(x = .,
             y = southamerican_samples_in_qiita,
             by = join_by(sample_title == sample_name)) %>% 
  # country and country_residence have the same info
  mutate(country = tolower(country_residence)) %>% 
  select(!c(sample_title, country_residence)) %>% 
  write_tsv(file = "manual_curation_of_samples/processed/filereport_PRJEB11419.txt")
  



# (3/5) PRJNA486009 -------------------------------------------------------

read_tsv(file = "manual_curation_of_samples/raw_filereport/filereport_PRJNA486009.txt") %>% 
  filter(grepl(x = scientific_name, pattern = "human gut metagenome")) %>% 
  filter(grepl(x = library_strategy, pattern = "AMPLICON")) %>% 
  select(study_accession, run_accession, original_library_layout = library_layout, scientific_name) %>% 
  mutate(country = "ecuador") %>% 
  write_tsv(file = "manual_curation_of_samples/processed/filereport_PRJNA486009.txt")



# (4/5) PRJEB3079 ---------------------------------------------------------

read_tsv(file = "manual_curation_of_samples/raw_filereport/filereport_PRJEB3079.txt") %>% 
  filter(grepl(x = library_name, pattern = "Amz")) %>% 
  select(study_accession, run_accession, original_library_layout = library_layout, scientific_name) %>% 
  mutate(country = "venezuela") %>% 
  write_tsv(file = "manual_curation_of_samples/processed/filereport_PRJEB3079.txt")
  



# (5/5) PRJEB5714 ---------------------------------------------------------

qiita <- read_tsv(file = "manual_curation_of_samples/raw_filereport/qiita_PRJEB5714.txt")

south_american_countries_in_qiita <- c("Venezuela")

# get the samples id from qiita to extract those samples from ENA
southamerican_samples_in_qiita <- qiita %>% 
  filter(host_common_name == "human") %>% #colnames() %>% sort()
  filter(grepl(x = country, pattern = south_american_countries_in_qiita)) %>%
  filter(grepl(x = host_body_site, pattern = "feces|faeces")) %>% 
  select(sample_name, country) %>% 
  mutate(country = gsub(x = country, pattern = ".*:", replacement = ""),
         sample_name = gsub(x = sample_name, pattern = ".*[.]", replacement = ""))


# read the file report of the project in ENA
read_tsv(file = "manual_curation_of_samples/raw_filereport/filereport_PRJEB5714.txt") %>% 
  filter(sample_title %in% southamerican_samples_in_qiita$sample_name) %>% 
  select(study_accession, run_accession, original_library_layout = library_layout, sample_title) %>% 
  # add the country to each sample
  inner_join(x = .,
             y = southamerican_samples_in_qiita,
             by = join_by(sample_title == sample_name)) %>% 
  # country and country_residence have the same info
  mutate(country = tolower(country)) %>% 
  select(!c(sample_title, country)) %>% 
  write_tsv(file = "manual_curation_of_samples/processed/filereport_PRJEB5714.txt")





# read metadata of all projects -------------------------------------------

readxl::read_xlsx(path = "manual_curation_of_samples/raw_filereport/samba_projects_metadata.xlsx",
                  sheet = "full_table") %>%
  filter(`can be used straight away`) %>% 
  select(project = `PRJ NUMBER`) %>% 
  pull(project) %>% 
  write_lines(x = .,
              file = "manual_curation_of_samples/processed/00.accession_codes_for_full_projects.txt")


# merge accession code for projects and runs ------------------------------

processed_projects <- list.files(path = "manual_curation_of_samples/processed/", pattern = "filereport*") %>% 
  paste0("manual_curation_of_samples/processed/", .) %>% 
  purrr::map(.x = .,
             .f = read_tsv) %>% 
  purrr::list_rbind()



processed_projects %>% 
  select(study_accession, run_accession) %>% 
  write_tsv(x = ., file = "manual_curation_of_samples/processed/04.sample_accession_codes_from_not_full_projects.txt")
