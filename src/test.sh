
for j in Bacillus Streptococcus
do
mkdir -p ~/HGTnew/multi_comparisons/data/mobilome/$j

for i in `ls ../../data/external/fasta/$j/ |head -n 100`
do
echo $j $i
date	

python3 extract_mobilome.py -d ~/HGTnew/multi_comparisons/data/external/database/AMRFinderPlus_Card_Megares_Resfinder.fa \
			    -q ~/HGTnew/multi_comparisons/data/external/fasta/$j/$i \
			    -o ~/HGTnew/multi_comparisons/data/mobilome/$j/$i

done
done


for i in `ls ~/HGTnew/multi_comparisons/data/mobilome/Bacillus/ |grep 'seq.fa' `
do
for j in `ls ~/HGTnew/multi_comparisons/data/mobilome/Streptococcus/ |grep 'seq.fa'`
do
python3 align_and_extract_coord_matches.py --sp1 Bacillus --sp2 Streptococcus \
	--st1  ../../data/mobilome/Bacillus/$i \
	--st2  ../../data/mobilome/Streptococcus/$j \
	-o ../../data/processed/test_mobilome/$i\_$j.out
done
done

