conda activate dataprocessing

cd /your/deduplicated/bam/folder

ls *_1_bismark_bt2_pe.deduplicated.bam > files.txt
files="files.txt"
file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $files)
echo "Start task ${SLURM_ARRAY_TASK_ID}, file name: $file"

raw_name=${file%_1_bismark_bt2_pe.deduplicated.bam}
if [ ! -e extraction ]; do mkdir extraction; fi

mkdir /extraction/$raw_name
bismark_methylation_extractor -p --gzip --bedGraph -o extraction/$raw_name $file
conda deactivate
