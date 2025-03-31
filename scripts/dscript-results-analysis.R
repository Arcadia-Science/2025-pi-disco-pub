library(tidyverse)

#################
# Analyze results from D-script for PPI predictions, filter based on multiple bait:candidate hits in different proteins
#################

# serine protease results
serine_protease_dscript_results <- read_tsv("results/serine_protease_results/2024-03-20-10:48.predictions.tsv", col_names = c("Amblyomma_protein", "Uniprot_target", "confidence_score"))

serine_protease_uniprot_metadata <- read_tsv("metadata/uniprotkb_UP000005640_serine_protease_2024_03_05.tsv")

serine_protease_results_filtered_metadata <- serine_protease_dscript_results %>%
  filter(confidence_score > 0.5) %>%
  mutate(Entry = str_extract(Uniprot_target, "(?<=\\|)[^|]+(?=\\|)")) %>%
  left_join(serine_protease_uniprot_metadata) %>%
  select(Amblyomma_protein, Entry, confidence_score, `Protein names`, Length)

high_confidence_hits <- serine_protease_results_filtered_metadata %>%
  filter(!(Amblyomma_protein %in% c("Amblyomma-americanum_evm.model.contig-30844-1.1", "Amblyomma-americanum_evm.model.contig-78294-1.1"))) %>%
  filter(confidence_score > .8) %>%
  group_by(Entry, `Protein names`) %>%
  summarize(
    count = n(),
    median_confidence_score = median(confidence_score),
    min_confidence_score = min(confidence_score),
    max_confidence_score = max(confidence_score),
    .groups = 'drop'
  ) %>%
  arrange(desc(count))


frequent_high_confidence_hits <- high_confidence_hits %>%
  filter(count >= 5) %>%
  pull(Entry)

serine_protease_high_confident_hits_table <- serine_protease_results_filtered_metadata %>%
  filter(Entry %in% frequent_high_confidence_hits) %>%
  filter(confidence_score > 0.8) %>%
  filter(!(Amblyomma_protein %in% c("Amblyomma-americanum_evm.model.contig-30844-1.1", "Amblyomma-americanum_evm.model.contig-78294-1.1")))

minimum_confidence_hits <- serine_protease_results_filtered_metadata %>%
  filter(!(Amblyomma_protein %in% c("Amblyomma-americanum_evm.model.contig-30844-1.1", "Amblyomma-americanum_evm.model.contig-78294-1.1"))) %>%
  filter(confidence_score > .45) %>%  # lower confidence score from 0.8 to .45 for a comparison
  group_by(Entry, `Protein names`) %>%
  summarize(
    count = n(),
    median_confidence_score = median(confidence_score),
    min_confidence_score = min(confidence_score),
    max_confidence_score = max(confidence_score),
    .groups = 'drop'
  ) %>%
  arrange(desc(count))

write.table(frequent_high_confidence_hits, "results/serine_protease_results/high_confidence_uniprot_hits.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(serine_protease_high_confident_hits_table, "results/serine_protease_results/high_confidence_serine_protease_hits_table.tsv", quote = FALSE, row.names = FALSE, sep = "\t")
