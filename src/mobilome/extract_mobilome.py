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

parser.add_option("-p", "--pdb_file",dest="pdb_file",action="store",default = "NO",
help="file containing an other database of genes in the mobilome, but this time with protein AA sequence instead of nucleotides")

parser.add_option("-n", "--nb_strains",dest="nb_strains",action="store",type = int, default = 50,
help="the number of strains from which we want to extract the mobilome")


(options, args) = parser.parse_args()
patt=re.compile("start=([0-9]+)")
pattATCG = re.compile("[^ATCGNatcgn\n]")

#query_file = options.query_file
clade_dir = options.clade_dir
db_file = options.db_file
pdb_file = options.pdb_file
output_file=options.output_file
nb_strains = options.nb_strains
pattNewline = re.compile("\\n")

#blastn -subject ../card-data/nucleotide_fasta_protein_homolog_model.fasta -query ../data/e
#xternal/fasta/Bacillus/GCA_001286765.1_Bacillus_JRS7_genomic.fna -outfmt 6 | sort -k9n



def launch_blast(db_file,pdb_file,subject_file,clade_dir):

	Matches = {}	
	soft = "/cluster/CBIO/home/fmassip/mambaforge/envs/HGT/bin/blastn"
	subject_full = clade_dir + "/"+ subject_file
	print(subject_full)
	print('blastn','-subject',subject_full,'-query',db_file)
	os.system("date")
	res_blast = subprocess.run(['blastn','-subject',subject_full,'-query',db_file,'-outfmt','6'], stdout=subprocess.PIPE).stdout.decode('utf-8')
	os.system("date")
#	f1.write(res_blast+"\n")
#	if db_file:
#		print("yes test")
	print('tblastn','-subject',subject_full,'-query',pdb_file)
	if pdb_file != "NO":
		res_tblast = subprocess.run(['tblastn','-subject',subject_full,'-query',pdb_file,'-outfmt','6'], stdout=subprocess.PIPE).stdout.decode('utf-8')
		os.system("date")
		all_res = res_blast+res_tblast
	else:
		all_res = res_blast
		print("no pdb")

	for line in pattNewline.split(all_res):
		row = line.split('\t')
		
		if len(row)>2 :
			row[8] = int(row[8])
			row[9] = int(row[9])
			if row[8] > row[9]:
				start = max(0,row[9]-2000)
				end = row[8]+2000
			else:
				start = max(0,row[8]-2000)
				end = row[9]+2000
 
			if row[1] in Matches.keys():
				Matches[row[1]].append([start,end])
			else:
				Matches[row[1]]=[[start,end]]


	return(Matches)



def mergeOverlap_and_extract_fasta(Matches,subject_file,opFile,clade_dir): 
	os.system("mkdir -p " + output_file + "_tmp")
	for cle in Matches.keys():
		all_pos = np.array(Matches[cle])
		all_pos_sorted = all_pos[all_pos[:, 0].argsort()]

		chrs = np.ones(all_pos_sorted.shape[0],dtype = 'object')
		chrs[:] = cle
		final_array = np.vstack((chrs,all_pos_sorted.transpose().astype('object'))).transpose()
		np.savetxt(output_file+'_tmp/'+subject_file,final_array,delimiter = "\t",fmt = '%s')
		

	concatenated_dat = subprocess.run(['mergeBed','-i',output_file+'_tmp/'+subject_file,'-d','4000'], stdout=subprocess.PIPE).stdout.decode('utf-8')
	toExtract = {}
	concatenated_dat = concatenated_dat.rstrip()

	for line in pattNewline.split(concatenated_dat):
		row = line.split('\t')
		if row[0] in toExtract.keys():
			toExtract[row[0]].append([row[1],row[2]])
		else:
			toExtract[row[0]]=[[row[1],row[2]]]
		
	subject_full = clade_dir + "/"+ subject_file
	fasta_data = open(subject_full,"r")
	records = list(SeqIO.parse(subject_full, "fasta"))
	all_contigs = SeqIO.to_dict(SeqIO.parse(subject_full, "fasta"))

	for cle in toExtract.keys():
		opFile.write(">"+cle+"\n")
		for index in range(0,len(toExtract[cle])):
			contig_len = len(all_contigs[cle])
			start = int(toExtract[cle][index][0])
			end = min(int(toExtract[cle][index][1]),contig_len-1)
			subseq = str(all_contigs[cle].seq[start:end])
			subseqN = re.sub(pattATCG,"N",subseq) #convert Ambiguous caracters to Ns
			opFile.write(subseqN)
			opFile.write("N")

		opFile.write("\n")

#	os.system("rm "+output_file+'_tmp')	


allFiles = os.listdir(clade_dir)

randList = random.sample(allFiles,k = min(nb_strains,len(allFiles)))
f2 = open(output_file+'_listFiles.txt',"w")
f2.write("\n".join(randList))
f2.close()


f1= open(output_file+'_seq.fa',"w")

if pdb_file != "NO":
	print("GO")
else:
	print("yes I can pbd file")
	print (pdb_file)


for myFile in randList:
	print(myFile)
	print (pdb_file)
	Matches = launch_blast(db_file = db_file,pdb_file = pdb_file ,subject_file = myFile,clade_dir = clade_dir)
	if Matches:
		mergeOverlap_and_extract_fasta(Matches = Matches, subject_file = myFile ,opFile = f1,clade_dir = clade_dir)



