Frangieh (2021) Data Documentation
================
Gene Katsevich;
May 4, 2022

# Overview

This repository contains code to import and process the Frangieh 2021
data. Frangieh et al developed the perturb-CITE-seq protocol, a new
single cell CRISPR screen assay providing both gene and protein
expression readouts. They applied perturb-CITE-seq to study 248 genes
and 20 proteins involved in cancer immunotherapy resistance in a large
screen with over 200,000 patient-derived melanoma cells in total. These
cells represent three experimental conditions (control, treated with
IFN-gamma, and co-cultured with tumor-infiltrating lymphocytes).

The `frangieh-2021` directory structure is as follows:

    ├── processed
    │   └── perturb-cite-seq
    │       ├── gene
    │       │   ├── gene_expression_matrix.odm
    │       │   ├── gene_expression_metadata.rds 
    │       │   ├── gene_expression_metadata_control.rds
    │       │   ├── gene_expression_metadata_ifn-gamma.rds
    │       │   └── gene_expression_metadata_co-culture.rds
    │       ├── gRNA
    │       |   ├── ...
    |       └── protein
    │           ├── ...
    ├── raw
    │   ├── ...

The contents of the `raw` directory are suppressed, as they are
unimportant. The `processed` directory contains the single subdirectory
`perturb-cite-seq`, which in turn contains subdirectories for the three
measured modalities: `gene`, `gRNA`, `protein`. Each of these
directories has five files each: one ODM file and four metadata files.
The ODM file contains the expression information for all cells and
features present in the raw data deposited in the Single Cell Portal.
There is one metadata file corresponding to this full set of cells and
features. The three additional metadata files contain cells from each of
the three conditions.

# Experimental design

248 genes with putative roles in immunotherapy resistance were targeted
by three gRNAs each, giving 744 targeting guide RNAs. There were
additionally 74 negative control gRNAs, for a total of 818 gRNAs. The
following table shows the breakdown of the cells into the three
experimental conditions:

``` r
processed_dir <- sprintf("%s/processed", 
                         .get_config_path("LOCAL_FRANGIEH_2021_DATA_DIR"))
processed_gene_dir <- sprintf("%s/perturb-cite-seq/gene", processed_dir)
gene_odm_fp <- sprintf("%s/gene_expression_matrix.odm", processed_gene_dir)
gene_metadata_fp <- sprintf("%s/gene_expression_metadata.rds", processed_gene_dir)
gene_odm <- ondisc::read_odm(gene_odm_fp, gene_metadata_fp)
gene_odm |> ondisc::get_cell_covariates() |> dplyr::pull(condition) |> table()
```

    ## 
    ## Co-culture    Control       IFNγ 
    ##      73114      57627      87590

Interestingly, this study contains not only negative control gRNAs
(i.e., negative control treatments) but also negative control proteins
(i.e. negative control outcomes). Let’s take a look:

``` r
processed_protein_dir <- sprintf("%s/perturb-cite-seq/protein", processed_dir)
protein_odm_fp <- sprintf("%s/protein_expression_matrix.odm", processed_protein_dir)
protein_metadata_fp <- sprintf("%s/protein_expression_metadata.rds", processed_protein_dir)
protein_odm <- ondisc::read_odm(protein_odm_fp, protein_metadata_fp)
```

Note that there are 24 features. Let’s take a closer look:

``` r
protein_odm |> ondisc::get_feature_ids()
```

    ##  [1] "CD117"       "CD119"       "CD140a"      "CD140b"      "CD172a"     
    ##  [6] "CD184"       "CD202b"      "CD274"       "CD29"        "CD309"      
    ## [11] "CD44"        "CD47"        "CD49f"       "CD58"        "CD59"       
    ## [16] "CD61"        "HLA_A"       "Rat_IgG2a"   "Mouse_IgG1"  "Mouse_IgG2a"
    ## [21] "Mouse_IgG2b" "HLA_E"       "CD9"         "CD279"

There are 20 proteins of interest, along with the following four
controls: `Rat_IgG2a`, `Mouse_IgG1`, `Mouse_IgG2a` `Mouse_IgG2b`. My
understanding is that these are four extra antibodies targeting proteins
that are not actually present in the population of (human) cells under
investigation. So whatever fluctuations of these four protein expression
are detected must be technical rather than biological. Theoretically, we
may pair these protein expressions with any perturbations (not
necessarily just negative control perturbations) to assess calibration.
Note that the expression of each protein was normalized not against
sequencing depth but against the expression of “its corresponding IgG
control.” This implies that each of the 20 proteins of interest has
associated with it one of the four controls, but I’m not sure what this
association is.

# Guide RNA assignments

Unlike the other datasets we are working with, this dataset does not
come with raw gRNA expressions. Instead, it comes with binary gRNA
assignments. Below is a histogram of the number of gRNAs per cell:

``` r
processed_gRNA_dir <- sprintf("%s/perturb-cite-seq/gRNA", processed_dir)
gRNA_odm_fp <- sprintf("%s/gRNA_assignments_ungrouped.odm", processed_gRNA_dir)
gRNA_metadata_fp <- sprintf("%s/gRNA_assignments_ungrouped_metadata.rds", processed_gRNA_dir)
gRNA_odm <- ondisc::read_odm(gRNA_odm_fp, gRNA_metadata_fp)
gRNA_odm |> 
  ondisc::get_cell_covariates() |>
  dplyr::filter(n_nonzero <= 5) |>
  ggplot2::ggplot(ggplot2::aes(x = n_nonzero, y = stat(count) / sum(count))) +
  ggplot2::geom_histogram(binwidth = 1, colour = "black") +
  ggplot2::scale_x_continuous(breaks = 0:5) +
  ggplot2::labs(x = "Number of gRNAs", y = "Frequency") +
  ggplot2::theme_bw() +
  ggplot2::theme(panel.grid.major.x = ggplot2::element_blank(),
                 panel.grid.minor.x = ggplot2::element_blank())
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

Note that only about 60% of cells received exactly one gRNA. About 10%
of cells had no detected gRNAs, and about 30% of cells had more than one
gRNA detected. It is unclear whether the association analysis in the
paper was restricted to cells with exactly one gRNA. We keep all cells
(regardless of gRNA count) in the processed data. Below is the number of
cells in each experimental condition with exactly one gRNA:

``` r
dplyr::left_join(
  gene_odm |>
    ondisc::get_cell_covariates() |>
    tibble::rownames_to_column(var = "cell_barcode") |>
    dplyr::select(cell_barcode, condition),
  gRNA_odm |>
    ondisc::get_cell_covariates() |>
    tibble::rownames_to_column(var = "cell_barcode") |>
    dplyr::select(cell_barcode, n_nonzero),
  by = "cell_barcode"
) |>
  dplyr::filter(n_nonzero == 1) |>
  dplyr::count(condition)
```

    ##    condition     n
    ## 1 Co-culture 46427
    ## 2    Control 30486
    ## 3       IFNγ 50053

# QC’d datasets

The QC’d datasets retain all of the original features, but only those
cells that have a given experimental condition. The gene expression and
gRNA indicator data contain three more cells than the protein expression
data. We remove all cells that are not shared across datasets.

Let’s take the control dataset as an example.

``` r
gRNA_metadata_control_fp <- sprintf("%s/gRNA_assignments_ungrouped_metadata_control.rds", processed_gRNA_dir)
gRNA_control_odm <- ondisc::read_odm(gRNA_odm_fp, gRNA_metadata_control_fp)
gRNA_control_odm
```

    ## A covariate_ondisc_matrix with the following components:
    ##  An ondisc_matrix with 818 features and 57624 cells.
    ##  A cell covariate matrix with columns n_nonzero, n_umis.
    ##  A feature covariate matrix with columns mean_expression, coef_of_variation, n_nonzero, target, target_type.

``` r
gene_metadata_control_fp <- sprintf("%s/gene_expression_metadata_control.rds", processed_gene_dir)
gene_control_odm <- ondisc::read_odm(gene_odm_fp, gene_metadata_control_fp)
gene_control_odm
```

    ## A covariate_ondisc_matrix with the following components:
    ##  An ondisc_matrix with 23712 features and 57624 cells.
    ##  A cell covariate matrix with columns n_nonzero, n_umis, condition, cluster_x, cluster_y, s_score, g2m_score, phase.
    ##  A feature covariate matrix with columns mean_expression, coef_of_variation, n_nonzero.

``` r
protein_metadata_control_fp <- sprintf("%s/protein_expression_metadata_control.rds", processed_protein_dir)
protein_control_odm <- ondisc::read_odm(protein_odm_fp, protein_metadata_control_fp)
protein_control_odm
```

    ## A covariate_ondisc_matrix with the following components:
    ##  An ondisc_matrix with 24 features and 57624 cells.
    ##  A cell covariate matrix with columns n_nonzero, n_umis, condition.
    ##  A feature covariate matrix with columns mean_expression, coef_of_variation, n_nonzero.

The gRNA data have 818 features, which is the total number of gRNAs used
in the experiment.

The gene expression data have 23712 features, which is the total number
of genes measured in the experiment. Frangieh et al did additional
feature QC: “Genes detected in fewer than 200 cells were…removed from
further analysis.” Also, Frangieh et al did not analyze all genes
surviving the above QC. They restricted their attention to the “1,000
highly variable genes were selected using the Scanpy v.1.4.4
implementation of highly variably gene selection…In addition, all
features from each condition’s ten most significant Jackstraw PCA
programs were included.” Luckily, neither the feature QC step nor the
highly variable gene selection has been done on the available data.
We’ll apply our own feature QC.

The protein expression data have 24 features. This includes the 20
proteins of interest as well as four controls, as discussed above.

# Absent information

*Batch.* With this number of cells, it’s almost certain that multiple
sequencing batches were used. This seems to be suggested in the methods
section: “15,000 cells loaded onto each of eight channels per condition
using the 10X Chromium system…” I am guessing that “channel” here means
sequencing batch or lane or something like that. However, batch effects
are not accounted for in the analysis, or discussed at all in the paper.
Unfortunately, batch information is also absent from the data published
in the Single Cell Portal.

# Covariates

The gene data include several key cell-specific covariates. We examine
the “control” ODM as an example:

``` r
gene_control_odm |>
  ondisc::get_cell_covariates() |>
  dplyr::select(n_nonzero, n_umis, s_score, g2m_score, phase) |>
  head()
```

    ##        n_nonzero n_umis     s_score   g2m_score phase
    ## CELL_1      3520  10832 -0.19767236  0.34070730   G2M
    ## CELL_2      3531  10731  0.22472879  0.21481543     S
    ## CELL_3      5541  28821 -0.07671482  0.15859704   G2M
    ## CELL_4      4086  15322  0.13803863 -0.03044829     S
    ## CELL_5      3178  10314 -0.09701530 -0.29331199    G1
    ## CELL_6      3124   8810  0.03581570 -0.35851066     S

Covariates include `n_nonzero`, `n_umis`, `s_score`, `g2m_score`, and
`phase`. Frangieh used the covariates `n_umis` and `phase` in their
MIMOSCA regression. We computed the covariates `s_score`, `g2m_score`,
and `phase` using Seurat’s `CellCycleScoring` function, which is
equivalent (as far as I can tell) to `scanpy`’s cell cycle function,
which Frangieh used to compute cell cycle. The covariates `s_score` and
`g2m_score` are continuous covariates related to cell cycle that are
used to compute `phase`. Satija recommends regressing on `s_score` and
`g2m_score` instead of `phase`.

# Note: data size

The raw gene expression data, a CSV file of normalized gene expressions,
is 25GB. There should be no reason for lab members to download this file
onto their local machine, since it is already processed. Therefore, I
discourage using `hpcc pull FRANGIEH_2021` to get the processed data
onto a local machine (ideally, at some point we should extend the
`hpcc pull` functionality to allow pulling of only the processed data).
If you’d like to get the processed data onto your machine, I recommend
using rsync directly:

    rsync -rltvP $REMOTE_FRANGIEH_2021_DATA_DIR/processed/ $LOCAL_FRANGIEH_2021_DATA_DIR/processed/

Note that even the processed data are about 15GB in size.
