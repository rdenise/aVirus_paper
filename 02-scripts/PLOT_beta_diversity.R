library(phyloseq)
library(tidyverse)


taxpasta <- read_tsv('taxpasta.all.full.species.virus.tsv')
# taxpasta <- read_tsv('taxpasta.all.full.genus.virus.tsv')

taxpasta <- taxpasta |>
  select(c(
    "taxonomy_id", "name", "rank",
    "lineage", "id_lineage", "rank_lineage",
    'long - none - 1x - kraken2 - Full',
    'medium - heavy - 1x - kraken2 - Full',
    'short - heavy - 1x - kraken2 - Full',
    'long - none - 1x - centrifuge - Full',
    'medium - heavy - 1x - centrifuge - Full',
    'short - heavy - 1x - centrifuge - Full',
    'long - none - 1x - metabuli - Full',
    'medium - heavy - 1x - metabuli - Full',
    'short - heavy - 1x - metabuli - Full',
    'long - none - 1x - centrifuge - Virus',
    'medium - heavy - 1x - centrifuge - Virus',
    'short - heavy - 1x - centrifuge - Virus',
    'long - none - 1x - kraken2 - Virus',
    'medium - heavy - 1x - kraken2 - Virus',
    'short - heavy - 1x - kraken2 - Virus',
    'long - none - 1x - metabuli - Virus',
    'medium - heavy - 1x - metabuli - Virus',
    'short - heavy - 1x - metabuli - Virus',
    'long - none - 1x - truth - truth',
    'medium - heavy - 1x - truth - truth',
    'short - heavy - 1x - truth - truth',
    'long - none - 1x - krakenuniq - Full',
    'medium - heavy - 1x - krakenuniq - Full',
    'short - heavy - 1x - krakenuniq - Full',
    'long - none - 1x - krakenuniq - Virus',
    'medium - heavy - 1x - krakenuniq - Virus',
    'short - heavy - 1x - krakenuniq - Virus'
    ))

# Assume df_long with columns:
# - 'reads': the read counts,
# - 'reads_size': the identifier for the condition (e.g. "long", "medium", "short"),
# - 'software': a field that equals "truth" for the truth rows

df_long <- taxpasta %>%
  pivot_longer(
    cols = -c(taxonomy_id, name, rank, lineage, id_lineage, rank_lineage),
    names_to = "sample",
    values_to = "reads"
  )

# First, select only the columns you need (here we assume all sample columns are in the table)
df_long <- taxpasta %>%
  pivot_longer(
    cols = -c(taxonomy_id, name, rank, lineage, id_lineage, rank_lineage),
    names_to = "sample",
    values_to = "reads"
  )

# Extract parts of the sample name.
# Adjust the extraction as needed. Here we assume the sample name is something like:
# "long - none - 1x - truth - truth" or "long - none - 1x - kraken2 - Full"
df_long <- df_long %>%
  mutate(
    # The first element is the condition (read length)
    reads_size = word(sample, 1, sep = " - "),
    # The tool is in the 4th or 5th position; here we assume that the 'truth' columns can be detected
    software = if_else(str_detect(sample, "truth"), "truth", word(sample, 5, sep = " - "))
  )

# Extract the truth rows and compute the reference truth value per species
truth_ref <- df_long %>%
  filter(software == "truth") %>%
  group_by(name) %>%
  summarize(ref_truth = mean(reads, na.rm = TRUE))

df_norm <- df_long %>%
  # Join the reference truth per species
  left_join(truth_ref, by = "name") %>%
  group_by(name, reads_size) %>%
  # Get the truth count for this species and condition (should be unique per group)
  mutate(
    condition_truth = reads[software == "truth"][1],  # [1] in case there is only one
    norm_factor = ref_truth / condition_truth,
    normalized_reads = reads * norm_factor
  ) %>%
  ungroup()

df_wide_norm <- df_norm %>%
  select(taxonomy_id, name, rank, lineage, id_lineage, rank_lineage, sample, normalized_reads) %>%
  pivot_wider(
    names_from = sample,
    values_from = normalized_reads
  )

df_wide_norm <- df_wide_norm %>%
  mutate(across(everything(), ~ replace(., is.nan(.), 0)))

write.csv(df_wide_norm, file = "normalized_data.csv", row.names = FALSE)

# Generate OTU table
otu_mat <- df_wide_norm |>
  filter(rank == "species") |>
  # filter(rank == "genus") |>
  select(-c(taxonomy_id, rank, lineage, rank_lineage, id_lineage)) |>
  column_to_rownames(var = "name") |>
  as.matrix()



# Generate taxa table
tax_info <- df_wide_norm |>
  filter(rank == "species") |>
  # filter(rank == "genus") |>
  select(name, lineage, rank_lineage)

tax_mat <- map_dfr(seq(1, nrow(tax_info)), function(i) {
    lineage <- str_split(tax_info$lineage[i], ";")[[1]]
    rank <- str_split(tax_info$rank_lineage[i], ";")[[1]]
    tibble(name = tax_info$name[i],
           lineage = lineage,
           rank = rank)
  }) |>
  filter(rank %in% c("superkingdom", "phylum", "class", "order", "family", "genus", "species")) |>
  # filter(rank %in% c("superkingdom", "phylum", "class", "order", "family", "genus")) |>
  pivot_wider(names_from = "rank", values_from = "lineage") |>
  column_to_rownames(var = "name") |>
  as.matrix()

sample_tbl <- tibble(id = colnames(otu_mat),
                     reads_size = str_match(colnames(otu_mat), "([a-z]+) - ([a-z]+) - ([a-z0-9]+) - ([a-z0-9]+) - ([A-ZA-z]+)")[,2],
                     deamination = str_match(colnames(otu_mat), "([a-z]+) - ([a-z]+) - ([a-z0-9]+) - ([a-z0-9]+) - ([A-ZA-z]+)")[,3],
                     software = str_match(colnames(otu_mat), "([a-z]+) - ([a-z]+) - ([a-z0-9]+) - ([a-z0-9]+) - ([A-ZA-z]+)")[,5],
                     dataset_used = str_match(colnames(otu_mat), "([a-z]+) - ([a-z]+) - ([a-z0-9]+) - ([a-z0-9]+) - ([A-ZA-z]+)")[,6]) |>
   column_to_rownames(var = "id")

sample_tbl$software_deamination = paste(sample_tbl$software,sample_tbl$deamination,sep=" - ")
sample_tbl$condition = paste(sample_tbl$reads_size, sample_tbl$software_deamination, sep=' - ')

# Construct the phyloSeq object
phy <- phyloseq(otu_table(otu_mat, taxa_are_rows = T),
                tax_table(tax_mat),
                sample_data(sample_tbl))

# Calculate beta diversity
beta_div <- ordinate(phy, method = "PCoA", distance = "bray")

# Plot PCoA
p = plot_ordination(phy, beta_div, color="software", shape='reads_size')
p + geom_polygon(aes(fill=condition)) + geom_point(size=3)
