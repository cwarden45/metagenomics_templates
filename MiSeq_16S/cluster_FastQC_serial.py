import sys
import re
import os
import subprocess
from Bio.Seq import Seq

email = "cwarden@coh.org"
threads = 1
parameterFile = "parameters.txt"
finishedSamples = ()

readsFolder = "../Reads"


fastqcFolder = readsFolder + "/QC"
command = "mkdir " + fastqcFolder
os.system(command)
	
fileResults = os.listdir(readsFolder)

submitAll = "run_FastQC_serial.sh"
outHandle = open(submitAll,"w")
text = "#!/bin/bash\n"
text = text + "#$ -M "+email+"\n"
text = text + "#$ -m bea\n"
text = text + "#$ -N fastQC16S\n"
text = text + "#$ -q all.q\n"
text = text + "#$ -pe shared "+str(threads)+"\n"
text = text + "#$ -l vf=4G\n"
text = text + "#$ -j yes\n"
text = text + "#$ -o fastQC16S.log\n"
text = text + "#$ -cwd\n"
text = text + "#$ -V\n"
outHandle.write(text)

for file in fileResults:
	result = re.search("(.*_L\d{3}_R\d_001).fastq.gz$",file)
	
	if result:
		sample = result.group(1)

		fastqcPrefix = re.sub(".fastq.gz","",file)		
		fastQCsubfolder = fastqcFolder + "/" + fastqcPrefix + "_fastqc"
		
		if (sample not in finishedSamples) and not os.path.isdir(fastQCsubfolder):
			print sample

			read1 = readsFolder + "/" + file
			
			text = "/net/isi-dcnl/ifs/user_data/Seq/FastQC/fastqc -o "+fastqcFolder+" " + read1 + "\n"
			outHandle.write(text)
			
			text = "unzip "+fastqcFolder+"/"+sample+"_fastqc.zip -d "+fastqcFolder+"\n"
			outHandle.write(text)

			text = "rm "+fastqcFolder+"/"+sample+"_fastqc.zip\n";
			outHandle.write(text)

			text = "rm "+fastqcFolder+"/"+sample+"_fastqc.html\n";
			outHandle.write(text)