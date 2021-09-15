#Run Velocyto on the cluster
##Set up the environment
cd /share/lasallelab/Osman/RettSingleCellRNAseq/01-Cellranger/
module load anaconda3
#conda create -p /share/lasallelab/Osman/RettSingleCellRNAseq/01-Cellranger/velocyto
aklog
conda init bash
source ~/.bashrc
conda activate /share/lasallelab/Osman/RettSingleCellRNAseq/01-Cellranger/velocyto

# install prerequisites
#conda install numpy scipy cython numba matplotlib scikit-learn h5py click
# install velocyto
#pip install velocyto

#Run Velocyto
velocyto run --bcfile /share/lasallelab/Osman/RettSingleCellRNAseq/01-Cellranger/19_0106_C/outs/raw_feature_bc_matrix/barcodes.tsv.gz --mask /share/workshop/adv_scrnaseq/$USER/references/GRCh38_rmsk.gtf --outputfolder /share/lasallelab/Osman/RettSingleCellRNAseq/01-Cellranger/02-Velocyto --samtools-threads 42 --samtools-memory 10000 /share/lasallelab/Osman/RettSingleCellRNAseq/01-Cellranger/19_0106_C/outs/possorted_genome_bam.bam /share/workshop/adv_scrnaseq/$USER/references/refdata-gex-GRCh38-2020-A/genes/genes.gtf > velocyto_sample1.err > velocyto_sample1.out
