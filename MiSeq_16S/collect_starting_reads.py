import sys
import re
import os
import subprocess

fullReadsFolder = "../Reads"
statFile = "FastQC_total_reads.txt"


statHandle = open(statFile,"w")
text = "SampleID\tTotalReads\n"
statHandle.write(text)
	
fastqcFolder = fullReadsFolder + "/QC"
fileResults = os.listdir(fullReadsFolder)

for file in fileResults:
	result = re.search("(.*)_L\d{3}_R1_001.fastq.gz$",file)
	
	if result:
		sample = result.group(1)
		print sample
		
		
		#get total reads from FastQC
		fastqcPrefix = re.sub(".fastq.gz","",file)
		fastQCtext = fastqcFolder + "/" + fastqcPrefix + "_fastqc/fastqc_data.txt"
		
		inHandle = open(fastQCtext)
		line = inHandle.readline()
		
		lineCount = 0
		
		while line:
			line = re.sub("\n","",line)
			line = re.sub("\r","",line)
			
			lineCount += 1
			
			if lineCount == 7:
				
				totalResult = re.search("Total Sequences\t(\d+)",line)
				if totalResult:
					totalReads = int(totalResult.group(1))
				else:
					print "Problem parsing FastQC file!\n"
					sys.exit()
			
			line = inHandle.readline()
		
		inHandle.close()
		
		text = sample + "\t" + str(totalReads) + "\n"
		statHandle.write(text)