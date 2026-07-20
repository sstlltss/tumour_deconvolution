cd /your/fastq/folder

conda activate dataprocessing

ls *_1.fastq.gz > files1.txt
ls *_2.fastq.gz > files2.txt
files1="files1.txt"
files2="files2.txt"
file1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $files1)
file2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $files2)
echo "Start task ${SLURM_ARRAY_TASK_ID}, file name: $file1"

bismark --genome /your/ref_genome/folder -1 $file1 -2 $file2 --rg_tag -o .

conda deactivate
