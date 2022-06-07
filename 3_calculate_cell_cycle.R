frangieh_dir <- .get_config_path("LOCAL_FRANGIEH_2021_DATA_DIR")
library(ondisc)
library(Seurat)

# loop over the three experimental conditions
exper_names <- c("co_culture", "control", "ifn_gamma")
for(exper_name in exper_names){
  processed_gene_dir <- sprintf(
    "%sprocessed/%s/gene",
    frangieh_dir, exper_name
  )
  odm_fp <- sprintf("%s/gene_expression_matrix.odm", processed_gene_dir)
  metadata_fp <- sprintf("%s/gene_expression_metadata.rds", processed_gene_dir)
  odm <- read_odm(odm_fp, metadata_fp)
  
  # load the s genes and g2m genes
  avail_genes <- odm |> get_feature_ids()
  s_genes <- cc.genes$s.genes[cc.genes$s.genes %in% avail_genes]
  g2m_genes <- cc.genes$g2m.genes[cc.genes$g2m.genes %in% avail_genes]
  
  # append cell cycle info to odm
  exp <- lowmoi::load_whole_odm(odm = odm)
  s_obj <- CreateSeuratObject(exp)
  rm(exp)
  s_obj <- NormalizeData(s_obj)
  s_obj <- CellCycleScoring(s_obj, s_genes, g2m_genes)
  odm |> mutate_cell_covariates(
    s_score = s_obj$S.Score,
    g2m_score = s_obj$G2M.Score,
    phase = factor(s_obj$Phase)
  ) |>
    save_odm(metadata_fp)
}