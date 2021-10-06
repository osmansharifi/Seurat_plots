#!/bin/bash
#
#SBATCH --time=9-0  # days-hours
#SBATCH --job-name=OScount_all_samples # Job name
#SBATCH --array=1-64 # samples in 'samples', usually 1-4
#SBATCH --nodes=1
#SBATCH --ntasks=32  # Number of cores
#SBATCH --mem=64000  # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --partition=production  # Partition to submit to
#SBATCH --output=countArrayJob_%A_%a.out  # File to which STDOUT will be written
#SBATCH --error=countArrayJob_%A_%a.err  # File to which STDERR will be written
#SBATCH --mail-type=END  # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=osharifi@ucdavis.edu  # Email to which notifications will be sent

start=`date +%s`
hostname
#PATH=/share/biocore/software/bin/:$PATH  # for cellranger binary
module load cellranger/3.1.0
THREADS=${SLURM_NTASKS}
MEM=$(expr ${SLURM_MEM_PER_NODE} / 1000)
echo "Slurm job id: " $SLURM_JOB_ID
echo "Allocated threads: " $THREADS
echo "Allocated memory: " $MEM

project="/share/lasallelab/Osman/RettSingleCellRNAseq"
sample=`sed "${SLURM_ARRAY_TASK_ID}q;d" ${project}/00-RawData/samples`
fastqs="/share/lasallelab/Osman/RettSingleCellRNAseq/00-RawData"
transcriptome="/share/genomes/cellranger_genomes/refdata-cellranger-mm10-3.0.0"

call="cellranger count \
  --id=${sample} \
  --sample=${sample} \
  --transcriptome=${transcriptome} \
  --fastqs=${fastqs}/${sample} \
  --localcores=${THREADS} \
  --localmem=${MEM}"

echo $call
eval $call

end=`date +%s`

runtime=$((end-start))

echo $runtime

