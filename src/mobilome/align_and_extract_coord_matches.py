# -*-coding:utf8 -*

import os
import sys
import pickle
import numpy as np
import random
import re
from collections import defaultdict
from optparse import OptionParser
import subprocess
from datetime import datetime
################################



parser = OptionParser()
parser.add_option("-o","--out", dest="out",action="store",
		  help="output file")
parser.add_option("--st1",dest="strain1",action="store",
		  help="strain1")

parser.add_option("--st2",dest="strain2",action="store",
		  help="strain2")

parser.add_option("--sp1",dest="species1",action="store",
		  help="species1")

parser.add_option("--sp2",dest="species2",action="store",
		  help="species2")
#parser.add_option("-l", "--list",dest="list_file",action="store",
#		  help="list of multifasta file to read")
#
#
(options, args) = parser.parse_args()
out=options.out
strain1 = options.strain1
strain2 = options.strain2
species1 = options.species1
species2 = options.species2
#list_file = options.list_file
#
#listFile=[]
#
#patt=re.compile("plasmid")

patt = re.compile("END")
pattNum = re.compile("[0-9]+$")
pattMiss = re.compile("\^")
align = ""
pattNewline = re.compile("\\n")
pattEq = re.compile("=")
#in_file = "~/HGTnew/dating_long_distance/results/nucmer/EscherichiaColi_KlebsiellaPneumoniae/SAMN02850627_0.fasta-SAMN02141993_0.fasta.delta"

myDir = "/cluster/CBIO/data1/fmassip/HGT/ProjectMisha/HGTnew/multi_comparisons/data/mobilome/"
#species1 = "Acinetobacter"
#species2 = "EscherichiaColi"

def add_MLDs(full_MLD,new_MLD):
	full_MLD[new_MLD[:,0],1] = full_MLD[new_MLD[:,0],1] + new_MLD[:,1]
	return(full_MLD)

def parse_Align(strain1,strain2):	
	
	iden = []
	blocLen = []
	Matches = []
	pos = 0
	
#	strain1_full_path = myDir+species1+"/"+strain1
#	strain2_full_path = myDir+species2+"/"+strain2
	strain1_full_path = strain1+"[multiple]"
	strain2_full_path = strain2+"[multiple]"


	res_lastZ = subprocess.run(['lastz',strain1_full_path,strain2_full_path,"--format=general:name1,start1,length1,name2,start2,length2,cigarx",
	"--traceback=2000.0M"],
	#,"--seed=111010010111"],
	 stdout=subprocess.PIPE).stdout.decode('utf-8')
	
	f2 = open(out+".lastZ","w")
#	f2.write(res_lastZ)
	f2.close()
	f1 = open(out,"w")
	savedSub = ""

	for alignBloc in res_lastZ.split("\n")[1:-1]:

#		print(alignBloc)
		nameRef,current_pos_ref,length1,nameSub,current_pos_sub,length2,cigar = alignBloc.split("\t")
		if savedSub != nameSub:
			f1.write(">"+nameSub+"\n")
			savedSub = nameSub

		#print(alignBloc)
#		print("percId:"+str(percId))
		blocLen.append(str(length1))
#		print(cigar)
		tmp_digit = "0"
		current_pos_ref = int(current_pos_ref)
		current_pos_sub = int(current_pos_sub)
		for  charac in cigar:
			if re.match("\d",charac):
				tmp_digit = tmp_digit+charac
			else:
				if tmp_digit == "0":
					tmp_digit = 1
				else:
					tmp_digit = int(tmp_digit)

				if charac == "=" :
					f1.write(" "+nameRef+"\t"+str(current_pos_ref)+"\t"+str(current_pos_sub)+"\t"+str(tmp_digit)+"\n")
					current_pos_ref = current_pos_ref + tmp_digit
					current_pos_sub = current_pos_sub + tmp_digit
				elif charac == "D" : 
					current_pos_ref = current_pos_ref + tmp_digit				
				elif charac == "I" : 
					current_pos_sub = current_pos_sub + tmp_digit				

				if charac == "X" :
					current_pos_ref = current_pos_ref + tmp_digit
					current_pos_sub = current_pos_sub + tmp_digit
				tmp_digit = "0"
#		print(nameSub+"\t"+str(current_pos_ref)+"\t"+nameRef+"\t"+str(current_pos_sub)+"\t"+str(tmp_digit)+"\tEOF")

	now = datetime.now()
	dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
	print(dt_string+"\n"+strain1+" #### "+strain2)		
	print("Done")


#	return(iden,blocLen,MLD)
#	return(blocLen,MLD)

parse_Align(strain1,strain2)

#identity,blocLen,MLD = compute_MLD(strain1,strain2)
#res = np.transpose(np.vstack((identity,blocLen)))


#np.savetxt(out+".MLD",MLD,delimiter = "\t",fmt = '%d')
#np.savetxt(out+".id",res,delimiter = "\t",fmt = '%s')




