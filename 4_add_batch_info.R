frangieh_dir <- .get_config_path("LOCAL_FRANGIEH_2021_DATA_DIR")
library(ondisc)
library(readr)
library(dplyr)
library(tibble)

# read in the cell barcode to label mapping
barcode_to_label_path <- sprintf("%s/raw/correspondence/cell_barcode_label_mapping.csv", frangieh_dir)
barcode_to_label_map <- read_csv(barcode_to_label_path)

# get the batches by splitting the cell barcodes
batch_info <- barcode_to_label_map |>
  rename(barcode = `Cell Barcode`, cell_name = `Cell Name`) |>
  rowwise() |>
  mutate(batch = strsplit(barcode, split = "[.]")[[1]][1]) |>
  ungroup() |>
  select(cell_name, batch)

# loop over the three experimental conditions
exper_names <- c("co_culture", "control", "ifn_gamma")
for(exper_name in exper_names) {
  processed_gene_dir <- sprintf(
    "%sprocessed/%s/gene",
    frangieh_dir, exper_name
  )

  # get expression ODM
  odm_fp <- sprintf("%s/gene_expression_matrix.odm", processed_gene_dir)
  metadata_fp <- sprintf("%s/gene_expression_metadata.rds", processed_gene_dir)
  odm <- read_odm(odm_fp, metadata_fp)

  # join the batch information into the cell covariates
  batches <- odm |>
    get_cell_covariates() |>
    rownames_to_column(var = "cell_name") |>
    left_join(batch_info, by = "cell_name") |>
    pull(batch)
  odm_with_batch <- odm |>
    mutate_cell_covariates(batch = batches)

  # save ODM with batch information
  save_odm(odm_with_batch, metadata_fp)
}
