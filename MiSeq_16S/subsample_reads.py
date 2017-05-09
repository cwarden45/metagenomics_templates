import os
import re
import sys
#from Bio import SeqIO

fullReadFolder = "../Reads"
subsampleReadFolder = "../Reads/Subsample_Reads"

subsampleCount = 10000

command = "mkdir " + subsampleReadFolder
os.system(command)

fileResults = os.listdir(fullReadFolder)

lineLimit = 4 * subsampleCount

for file in fileResults:
	result = re.search("(.*).fastq$",file)
	
	if result:
		print file
		inFQ = fullReadFolder + "/" + file
		outFQ = subsampleReadFolder + "/" + file
		
		inHandle = open(inFQ)
		line = inHandle.readline()

		outHandle = open(outFQ, "w")

		lineCount = 0
		
		while line:	
			lineCount += 1
			
			if lineCount <= lineLimit:
				text = line
				outHandle.write(text)
			else:	
				break
		

			line = inHandle.readline()