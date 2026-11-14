# 01_prepare_data.R
# Clean the sequences downloaded from NCBI GenBank.
#
# Input : data/raw_download_all_segments.fasta   (as downloaded; contains 20 records)
# Output: data/ha_sequences.fasta        only the HA (segment 4) sequences
#         data/ha_metadata.csv           accession, subtype, location, year, ...
#         data/selected_accessions.csv   accession list for the report
#
# Run from the project root:  source("code/01_prepare_data.R")
#
# Header format of the download:
#   >KC190041 A/Sri Lanka/2001/2009 2009/11/02 4 (HA)
#   accession | strain name | collection date | segment number | (gene)

library(Biostrings)
library(tidyverse)

raw <- readDNAStringSet("data/raw_download_all_segments.fasta")
cat("Records in raw file:", length(raw), "\n")

hdr <- names(raw)
m   <- str_match(hdr, "^(\\S+) A/(.+?)/([^/]+)/(\\d{4}) (\\d{4}/\\d{2}/\\d{2}) (\\d) \\((\\w+)\\)$")

meta <- tibble(
  header          = hdr,
  accession       = m[, 2],
  location        = m[, 3],
  strain_id       = m[, 4],
  year            = as.integer(m[, 5]),
  collection_date = m[, 6],
  segment         = as.integer(m[, 7]),
  gene            = m[, 8],
  length          = width(raw)
)

# Which genes are in the file?  (the HA project needs only "HA")
print(count(meta, gene))
stopifnot(!any(is.na(meta$accession)))     # every header must have been parsed

meta_ha <- meta |>
  filter(gene == "HA") |>
  mutate(
    subtype = "H1N1",                      # subtype of the downloaded data set
    label   = str_replace_all(paste(accession, location, year, sep = "_"),
                              "[^A-Za-z0-9_.]", "_")
  )

ha <- raw[meta$gene == "HA"]
names(ha) <- meta_ha$accession

cat("HA sequences kept:", length(ha), "\n")
cat("Length range:", range(width(ha)), "\n")

writeXStringSet(ha, "data/ha_sequences.fasta")
write_csv(meta_ha |> select(accession, subtype, location, year, collection_date,
                            length, label),
          "data/ha_metadata.csv")
write_csv(meta_ha |> select(accession, subtype, location, year, collection_date, length),
          "data/selected_accessions.csv")
