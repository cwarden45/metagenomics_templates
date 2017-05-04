import sys
import re
import os
import subprocess

fullReadsFolder = "../Reads"
cutadaptReadsFolder = "Cutadapt_Filtered_Reads"
statFile = "cutadapt_filter_rate.txt"


statHandle = open(statFile,"w")
text = "SampleID\tuserID\tTotalReads\tCutadapt.Reads\tPercent.Cutadapt.Filtered\n"
statHandle.write(text)
	
fastqcFolder = fullReadsFolder + "/QC"
fileResults = os.listdir(fullReadsFolder)

for file in fileResults:
	result = re.search("(.*_L\d{3}_R1_001).fastq.gz$",file)
	
	if result:
		sample = result.group(1)
		print sample
		
		r2 = re.search("^(\d+)_.*",sample)
		seqID = r2.group(1)
		
		shortID = "S"+seqID
		
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
		
		#get number of filtered forward reads
		cutadaptFQ = cutadaptReadsFolder + "/" + sample + ".fastq"
		
		command = "wc -l " + cutadaptFQ
		wcText = subprocess.check_output(command, shell=True)
		countResult = re.search("^\s+(\d+)",wcText)
		if countResult:
			fqLines = int(countResult.group(1))
		else:
			print "Problem parsing " + wcText
			sys.exit()
		cutadaptReads = fqLines / 4
		percentKept = 100 * float(cutadaptReads)/float(totalReads)
		
		text = sample + "\t" + shortID + "\t" + str(totalReads) + "\t" + str(cutadaptReads) + "\t" + '{0:.2f}'.format(percentKept) + "%\n"
		statHandle.write(text)