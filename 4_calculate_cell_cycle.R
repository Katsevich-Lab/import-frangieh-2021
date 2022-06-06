frangieh_dir <- .get_config_path("LOCAL_FRANGIEH_2021_DATA_DIR")
library(ondisc)
library(Seurat)

##########################################################
#  Load the three odms: control, co culture, and ifn gamma
##########################################################

backing_odm <- paste0(frangieh_dir, "processed/perturb-cite-seq/gene/gene_expression_matrix.odm")
control_metadata <- paste0(frangieh_dir, "processed/perturb-cite-seq/gene/gene_expression_metadata_control.rds")
co_culture_metadata <- paste0(frangieh_dir, "processed/perturb-cite-seq/gene/gene_expression_metadata_co-culture.rds")
ifn_gamma_metadata <- paste0(frangieh_dir, "processed/perturb-cite-seq/gene/gene_expression_metadata_ifn-gamma.rds")
dataset_list <- c(control_metadata, co_culture_metadata, ifn_gamma_metadata)

# load the s genes and g2m genes
avail_genes <- read_odm(backing_odm, control_metadata) |> get_feature_ids()
# get_feature_ids(control_odm)
s_genes <- cc.genes$s.genes[cc.genes$s.genes %in% avail_genes]
g2m_genes <- cc.genes$g2m.genes[cc.genes$g2m.genes %in% avail_genes]

# write function to append cell cycle info to given odm
append_cell_cycle_info <- function(odm, s_genes, g2m_genes) {
  # load data and put into seurat object
  exp <- lowmoi::load_whole_odm(odm = odm)
  s_obj <- CreateSeuratObject(exp)
  rm(exp)
  s_obj <- NormalizeData(s_obj)
  s_obj <- CellCycleScoring(s_obj, s_genes, g2m_genes)
  odm_plus_cell_cycle <- odm |> mutate_cell_covariates(s_score = s_obj$S.Score,
                                                       g2m_score = s_obj$G2M.Score,
                                                       phase = factor(s_obj$Phase))
}

# loop and apply
for (dataset in dataset_list) {
  odm <- read_odm(odm_fp = backing_odm, metadata_fp = dataset)
  odm <- append_cell_cycle_info(odm, s_genes, g2m_genes)
  save_odm(odm, dataset)
}
