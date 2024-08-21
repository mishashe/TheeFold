# -*-coding:utf8 -*

import os
import subprocess
import sys
import numpy as np
import pickle
import random
import re
from collections import defaultdict
from optparse import OptionParser
import pandas as pd
from Bio import SeqIO

################################





parser = OptionParser()
parser.add_option("-c","--clade_dir",dest="clade_dir",action="store",
help="the directory containing all strains of the clade we want to analyse")
#parser.add_option("-q","--query_file",dest="query_file",action="store",
#help="file containing the species full genome")
parser.add_option("-o","--output_file",dest="output_file",action="store",
help="output file name")
parser.add_option("-d", "--db_file",dest="db_file",action="store",
help="file containing the database of genes in the mobilome")
parser.add_option("-n", "--nb_strains",dest="nb_strains",action="store",type = int, default = 20,
help="the number of strains from which we want to extract the mobilome")


(options, args) = parser.parse_args()
patt=re.compile("start=([0-9]+)")

#query_file = options.query_file
clade_dir = options.clade_dir
db_file = options.db_file
output_file=options.output_file
nb_strains = options.nb_strains
pattNewline = re.compile("\\n")

#blastn -subject ../card-data/nucleotide_fasta_protein_homolog_model.fasta -query ../data/e
#xternal/fasta/Bacillus/GCA_001286765.1_Bacillus_JRS7_genomic.fna -outfmt 6 | sort -k9n



def launch_blast(db_file,query_file,clade_dir):

	Matches = {}	
	soft = "/cluster/CBIO/home/fmassip/mambaforge/envs/HGT/bin/blastn"
	query_full = clade_dir + "/"+ query_file
	print(query_full)
	print('blastn','-subject',db_file,'-query',query_full)
	res_blast = subprocess.run(['blastn','-subject',db_file,'-query',query_full,'-outfmt','6'], stdout=subprocess.PIPE).stdout.decode('utf-8')
#	f1.write(res_blast+"\n")
	for line in pattNewline.split(res_blast):
		row = line.split('\t')
		
		if len(row)>2 :
			row[7] = int(row[7])
			row[6] = int(row[6])
			if row[7] > row[6]:
				start = row[6]
				end = row[7]
			else:
				start = row[7]
				end = row[6]
 
			if row[0] in Matches.keys():
				Matches[row[0]].append([start,end])
			else:
				Matches[row[0]]=[[start,end]]
	return(Matches)



def mergeOverlap_and_extract_fasta(Matches,query_file,opFile,clade_dir): 
	os.system("mkdir -p " + output_file + "_tmp")
	for cle in Matches.keys():
		all_pos = np.array(Matches[cle])
		all_pos_sorted = all_pos[all_pos[:, 0].argsort()]

		chrs = np.ones(all_pos_sorted.shape[0],dtype = 'object')
		chrs[:] = cle
		final_array = np.vstack((chrs,all_pos_sorted.transpose().astype('object'))).transpose()
		np.savetxt(output_file+'_tmp/'+query_file,final_array,delimiter = "\t",fmt = '%s')
		

	concatenated_dat = subprocess.run(['mergeBed','-i',output_file+'_tmp/'+query_file,'-d','100'], stdout=subprocess.PIPE).stdout.decode('utf-8')
	toExtract = {}
	concatenated_dat = concatenated_dat.rstrip()

	for line in pattNewline.split(concatenated_dat):
		row = line.split('\t')
		if row[0] in toExtract.keys():
			toExtract[row[0]].append([row[1],row[2]])
		else:
			toExtract[row[0]]=[[row[1],row[2]]]
		
	query_full = clade_dir + "/"+ query_file
	fasta_data = open(query_full,"r")
	records = list(SeqIO.parse(query_full, "fasta"))
	all_contigs = SeqIO.to_dict(SeqIO.parse(query_full, "fasta"))

	for cle in toExtract.keys():
		opFile.write(">"+cle+"\n")
		for index in range(0,len(toExtract[cle])):
			start = int(toExtract[cle][index][0])
			end = int(toExtract[cle][index][1])
			opFile.write(str(all_contigs[cle].seq[start:end]))
			opFile.write("N")

		opFile.write("\n")

#	os.system("rm "+output_file+'_tmp')	


allFiles = os.listdir(clade_dir)

randList = random.sample(allFiles,k = nb_strains)
f2 = open(output_file+'_listFiles.txt',"w")
f2.write("\n".join(randList))
f2.close()


f1= open(output_file+'_seq.fa',"w")

for myFile in randList:
	print(myFile)
	Matches = launch_blast(db_file = db_file,query_file = myFile,clade_dir = clade_dir)
	mergeOverlap_and_extract_fasta(Matches = Matches, query_file = myFile ,opFile = f1,clade_dir = clade_dir)



