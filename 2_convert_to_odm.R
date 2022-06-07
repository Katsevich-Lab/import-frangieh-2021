###############################################################################
#
# Import Frangieh et al (2020) data.
#
#
# Notes: This script took 70 minutes and 225 GB of RAM to run. 
###############################################################################

### retrieve top-level data directory ###
frangieh_dir <- .get_config_path("LOCAL_FRANGIEH_2021_DATA_DIR")

exper_names <- c("co_culture", "control", "ifn_gamma")
modalities <- c("gene", "grna_assignment", "protein")

processed_dir_names <- matrix(NA, 3, 3, dimnames = list(exper_names, modalities))

for(exper_name in exper_names){
  for(modality in modalities){
    processed_dir_name <- sprintf(
      "%sprocessed/%s/%s",
      frangieh_dir, exper_name, modality
    )
    processed_dir_names[exper_name, modality] <- processed_dir_name
    if (!dir.exists(processed_dir_name)) dir.create(processed_dir_name, recursive = TRUE)
  }
}

# read gene expression metadata
gene_expr_metadata_filename <- sprintf("%sraw/single-cell-portal/metadata/RNA_metadata.csv", frangieh_dir)
gene_expr_metadata <- readr::read_csv(gene_expr_metadata_filename,
                                      col_types = list(
                                        MOI = readr::col_integer(),
                                        UMI_count = readr::col_double(),
                                        .default = readr::col_character()
                                      ),
                                      comment = "TYPE"
) |>
  dplyr::rowwise() |>
  dplyr::mutate(condition = tolower(gsub("γ", "-gamma", condition)),
                condition = gsub("-", "_", condition)) |>
  dplyr::ungroup()

### import protein data ###
cat("Reading protein expression matrix from file...\n")
prot_expr_filename <- sprintf("%sraw/single-cell-portal/other/raw_CITE_expression.csv", frangieh_dir)
prot_expr_data <- readr::read_csv(prot_expr_filename,
                                  # n_max = 3,
                                  col_types = list(
                                    `...1` = readr::col_character(),
                                    .default = readr::col_integer()
                                  )
) |>
  dplyr::rename(Protein = `...1`)

# extract cell barcodes and protein names
cell_barcodes_protein <- colnames(prot_expr_data)[-1]
protein_names <- prot_expr_data |>
  dplyr::select(Protein) |>
  dplyr::rename(protein_name = Protein)

# convert to matrix
cat("Converting protein expression data to matrix...\n")
prot_expr_data <- prot_expr_data |>
  dplyr::select(-Protein) |>
  as.matrix()

# split by experimental condition and save to disk
for(exper_name in exper_names){
  cat(sprintf("Creating ODM for %s protein expression matrix...\n", exper_name))
  processed_protein_dir <- processed_dir_names[exper_name, "protein"]
  odm_fp <- sprintf("%s/protein_expression_matrix.odm", processed_protein_dir)
  metadata_fp <- sprintf("%s/protein_expression_metadata.rds", processed_protein_dir)
  # find cells in this experimental condition
  cells_to_keep <- gene_expr_metadata |> 
    dplyr::filter(condition == exper_name) |>
    dplyr::pull(NAME)
  # create ondisc matrix
  ondisc::create_ondisc_matrix_from_R_matrix(
    r_matrix = prot_expr_data[,cells_to_keep],
    barcodes = cells_to_keep,
    features_df = protein_names,
    odm_fp = odm_fp,
    metadata_fp = metadata_fp
  ) |>
    # add condition information to cell covariates
    ondisc::mutate_cell_covariates(condition = exper_name) |>
    # save to disk
    ondisc::save_odm(metadata_fp = metadata_fp)
}

# remove prot_expr_data from workspace to save memory
rm(prot_expr_data)

### import gRNA data ###

cat("Reading gRNA assignments from file...\n")
gRNA_assignments_filename <- sprintf("%sraw/single-cell-portal/documentation/all_sgRNA_assignments.txt", frangieh_dir)
gRNA_list_filename <- sprintf("%sraw/supp_tables/41588_2021_779_MOESM3_ESM.xlsx", frangieh_dir)

gRNA_assignments <- readr::read_csv(gRNA_assignments_filename)
gRNA_list <- readxl::read_excel(gRNA_list_filename,
                                sheet = 1,  # first sheet corresponds to gRNA list
                                skip = 2,   # two first lines are header
                                n_max = 818 # only the first 818 gRNAs used for perturb-CITE-seq
) 

# extract the experimental design information
experimental_design <- gRNA_list |>
  dplyr::rowwise() |>
  dplyr::mutate(
    # those gRNAs with "SITE" in their names are non-targeting
    target_type = ifelse(grepl("SITE", `Guide Name`),
                         "non-targeting",
                         "gene"
    ),
    target = ifelse(target_type == "non-targeting",
                    "non-targeting",
                    strsplit(`Guide Name`, split = "_")[[1]][1]
    )
  ) |>
  dplyr::ungroup()

# extract the gRNA barcodes
gRNA_barcodes <- gRNA_list |>
  dplyr::rename(
    gRNA_barcode = `sgRNA Sequence`,
    gRNA_name = `Guide Name`
  ) |>
  dplyr::select(gRNA_barcode, gRNA_name)

# split by experimental condition and save to disk
for(exper_name in exper_names){
  cat(sprintf("Creating ODM for %s gRNA assignment matrix...\n", exper_name))
  processed_gRNA_dir <- processed_dir_names[exper_name, "grna_assignment"]
  odm_fp <- sprintf("%s/grna_assignments_ungrouped.odm", processed_gRNA_dir)
  metadata_fp <- sprintf("%s/grna_assignments_ungrouped_metadata.rds", processed_gRNA_dir)
  # find cells in this experimental condition
  cells_to_keep <- gene_expr_metadata |> 
    dplyr::filter(condition == exper_name) |>
    dplyr::pull(NAME) |>
    intersect(cell_barcodes_protein)   # there are a few cells for which gRNA
  # assignment data are available but not
  # protein expression data, so subset cells
  # to only those in protein data
  # get gRNA assignment list
  gRNA_assignment_list <- gRNA_assignments[match(cells_to_keep, gRNA_assignments$Cell), ] |>
    dplyr::pull("sgRNAs") |>
    lapply(function(sgRNA){
      if(is.na(sgRNA)) ""
      else strsplit(sgRNA, split = ",")[[1]]
    })
  # create odm 
  ondisc::convert_assign_list_to_sparse_odm(
    cell_barcodes = cells_to_keep,
    gRNA_ids = gRNA_barcodes$gRNA_name,
    gRNA_assignment_list = gRNA_assignment_list,
    odm_fp = odm_fp,
    metadata_fp = metadata_fp
  ) |>
    # add gRNA metadata to feature covariates
    ondisc::mutate_feature_covariates(
      target = experimental_design$target,
      target_type = experimental_design$target_type
    ) |>
    # save to disk
    ondisc::save_odm(metadata_fp = metadata_fp)
}

### import gene data ###
cat("Reading gene expression matrix from file...\n")
gene_expr_filename <- sprintf("%sraw/single-cell-portal/other/RNA_expression.csv", frangieh_dir)
gene_expr_cluster_filename <- sprintf("%sraw/single-cell-portal/cluster/5fd0e449771a5b0db7207711/RNA_UMAP_cluster.csv", frangieh_dir)

# read gene expression data
gene_expr_data <- readr::read_csv(gene_expr_filename,
  n_max = 3,
  col_types = list(
    GENE = readr::col_character(),
    .default = readr::col_double()
  )
)

# read gene expression clusters
gene_expr_clusters <- readr::read_csv(gene_expr_cluster_filename,
  col_types = list(
    NAME = readr::col_character(),
    X = readr::col_double(),
    Y = readr::col_double()
  ),
  comment = "TYPE"
)

# extract cell barcodes and gene names from gene expression data
cell_barcodes_gene <- colnames(gene_expr_data)[-1]
gene_names <- gene_expr_data |>
  dplyr::select(GENE) |>
  dplyr::rename(gene_name = GENE)

# convert to matrix
cat("Converting gene expression data to matrix...\n")
gene_expr_data_norm <- gene_expr_data |>
  dplyr::select(-GENE) |>
  as.matrix()

# remove gene_expr_data from workspace to save memory
rm(gene_expr_data)

# undo normalization
cat("Undoing normalization...\n")
gene_expr_data_unnorm <- gene_expr_data_norm |>
  exp() |>
  sweep(2, gene_expr_metadata$UMI_count / 1e6, "*") |>
  round()

# remove gene_expr_data_norm from workspace to save memory
rm(gene_expr_data_norm)

# split by experimental condition and save to disk
for(exper_name in exper_names){
  cat(sprintf("Creating ODM for %s gene expression matrix...\n", exper_name))
  processed_gene_dir <- processed_dir_names[exper_name, "gene"]
  odm_fp <- sprintf("%s/gene_expression_matrix.odm", processed_gene_dir)
  metadata_fp <- sprintf("%s/gene_expression_metadata.rds", processed_gene_dir)
  # find cells in this experimental condition
  cells_to_keep <- gene_expr_metadata |> 
    dplyr::filter(condition == exper_name) |>
    dplyr::pull(NAME) |>
    intersect(cell_barcodes_protein)   # there are a few cells for which gene
                                       # expression data are available but not
                                       # protein expression data, so subset cells
                                       # to only those in protein data
  # create ondisc matrix
  ondisc::create_ondisc_matrix_from_R_matrix(
    r_matrix = gene_expr_data_unnorm[,cells_to_keep],
    barcodes = cells_to_keep,
    features_df = gene_names,
    odm_fp = odm_fp,
    metadata_fp = metadata_fp
  ) |>
    # add condition information to cell covariates
    ondisc::mutate_cell_covariates(
      condition = gene_expr_metadata$condition[match(cells_to_keep, gene_expr_metadata$NAME)],
      cluster_x = gene_expr_clusters$X[match(cells_to_keep, gene_expr_clusters$NAME)],
      cluster_y = gene_expr_clusters$Y[match(cells_to_keep, gene_expr_clusters$NAME)]
    ) |>
    # save to disk
    ondisc::save_odm(metadata_fp = metadata_fp) 
}

# remove gene_expr_data_unnorm from workspace to save memory
rm(gene_expr_data_unnorm)

cat("Done.\n")
