
for j in Bacillus Streptococcus
do
mkdir -p ~/HGTnew/multi_comparisons/data/mobilome/$j

for i in `ls ../../data/external/fasta/$j/`
do
echo $j $i
date	

python3 extract_mobilome.py -d ~/HGTnew/multi_comparisons/card-data/nucleotide_fasta_protein_homolog_model.fasta \
			    -q ~/HGTnew/multi_comparisons/data/external/fasta/$j/$i \
			    -o ~/HGTnew/multi_comparisons/data/mobilome/$j/$i

done
done
