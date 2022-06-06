frangieh_dir <- .get_config_path("LOCAL_FRANGIEH_2021_DATA_DIR")
library(ondisc)
library(readr)
library(dplyr)
library(tibble)

# path to backing odm
backing_odm <- paste0(frangieh_dir, "processed/perturb-cite-seq/gene/gene_expression_matrix.odm")

# read gRNA assignment info
grna_assignments_path <- paste0(frangieh_dir, "raw/single-cell-portal/documentation/all_sgRNA_assignments.txt")
grna_assignments <- read_csv(grna_assignments_path) |>
  column_to_rownames(var = "Cell")

# list of three experiments
experiments <- c("control", "co-culture", "ifn-gamma")

# loop over the three experiments
for(experiment in experiments){
  metadata_fp <- sprintf("%sprocessed/perturb-cite-seq/gene/gene_expression_metadata_%s.rds", frangieh_dir, experiment)
  odm <- read_odm(odm_fp = backing_odm, metadata_fp = metadata_fp)
  perturbations <- grna_assignments[get_cell_barcodes(odm), "sgRNAs"] |>
    lapply(function(sgRNA){
      if(is.na(sgRNA)) character(0)
      else strsplit(sgRNA, split = ",")[[1]]
    })
  odm |> 
    mutate_cell_covariates(perturbation = perturbations) |>
    save_odm(metadata_fp = metadata_fp)
}
