# Import library to run Python
library(reticulate)

# Set conda environment
reticulate::use_condaenv("base")

# Read in and run make_master_godata_table.py
reticulate::source_python("~/GitHub/snRNA-seq-pipeline/scripts/04_GO_analysis/make_master_godata_table.py")