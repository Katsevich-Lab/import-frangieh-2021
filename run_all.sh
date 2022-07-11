# required R packages readr, dplyr, ondisc, readxl, Seurat, lowmoi

#$ -l m_mem_free=250G
#$ -q short.q

hpcc=$(hostname | grep "hpcc" | wc -l)
if [[ hpcc ]]
then
  module load R/R-4.1.2
fi

# source ~/.research_config
# source 1_download_data.sh
Rscript 2_convert_to_odm.R
Rscript 3_calculate_cell_cycle.R
Rscript 4_add_batch_info.R
