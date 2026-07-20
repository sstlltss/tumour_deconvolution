cd /your/extraction/folder

rm -f SRR*/*.txt*
mkdir cov
mkdir bedGraph
mv SRR*/*.cov.gz cov
mv SRR*/*.bedGraph.gz bedGraph
rm -r SRR*/

ls cov/*.cov.gz > files.txt
conda activate dataprocessing
python split_chr.py
conda deactivate