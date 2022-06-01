frangieh_dir <- .get_config_path("LOCAL_FRANGIEH_2021_DATA_DIR")
library(ondisc)
library(Seurat)

##########################################################
#  Load the three odms: control, co culture, and ifn gamma
##########################################################
backing_odm <- paste0(frangieh_dir, "processed/perturb-cite-seq/gene/gene_expression_matrix.odm")
control_metadata <- paste0(frangieh_dir, "processed/perturb-cite-seq/gene/gene_expression_metadata_control.rds")
co_culture_metadata <- paste0(frangieh_dir, "processed/perturb-cite-seq/gene/gene_expression_metadata_control.rds")
ifn_gamma_metadata <- paste0(frangieh_dir, "processed/perturb-cite-seq/gene/gene_expression_metadata_ifn-gamma.rds")

# load_odms
control_odm <- read_odm(odm_fp = backing_odm,
                        metadata_fp = control_metadata)
avail_genes <- get_feature_ids(control_odm)
s_genes <- cc.genes$s.genes[cc.genes$s.genes %in% avail_genes]
g2m_genes <- cc.genes$g2m.genes[cc.genes$g2m.genes %in% avail_genes]

# load data and put into seurat object
control_exp <- lowmoi::load_whole_odm(odm = control_odm)
s_obj <- CreateSeuratObject(control_exp)
rm(control_exp)
s_obj <- NormalizeData(s_obj)
s_obj <- CellCycleScoring(s_obj, s_genes, g2m_genes)
control_odm |> mutate_cell_covariates(s_score = s_obj$S.Score,
                                      g2m_score = s_obj$G2M.Score,
                                      phase = factor(s_obj$Phase))

