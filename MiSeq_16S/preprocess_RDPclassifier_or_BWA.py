import sys
import re
import os
from Bio import SeqIO

parameterFile = "parameters.txt"
finishedSamples = []

threads = ""
readsFolder = ""

inHandle = open(parameterFile)
lines = inHandle.readlines()
			
for line in lines:
	print line
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	param = lineInfo[0]
	value = lineInfo[1]

	if param == "Reads_Folder":
		readsFolder = value
		
	if param == "Threads":
		threads = value
		
if (readsFolder == "") or (readsFolder == "[required]"):
	print "Need to enter a value for 'Reads_Folder'!"
	sys.exit()
		
if (threads == "") or (threads == "[required]"):
	print "Need to enter a value for 'Threads'!"
	sys.exit()

mergedReadsFolder = readsFolder + "/PEAR_Merged"
if not os.path.exists(mergedReadsFolder):
	command = "mkdir " + mergedReadsFolder
	os.system(command)
	
fileResults = os.listdir(readsFolder)
for file in fileResults:
	result = re.search("(.*)_L\d{3}_R1_001.fastq$",file)
	if result:
		command = "gzip " + file
		os.system(command)
		file = file + ".gz"
	
	resultGZ = re.search("(.*)_L\d{3}_R1_001.fastq.gz$",file)
	
	if resultGZ:
		sample = resultGZ.group(1)
		if sample not in finishedSamples:
			print sample
			
			read1 = readsFolder + "/" + file
			read2 = re.sub("R1_001.fastq","R2_001.fastq",read1)
								
			print "\n\nMerge Reads using PEAR\n\n"
			pearPrefix = mergedReadsFolder + "/" + sample
			command = "pear -f " + read1 + " -r " + read2 + " -o " + pearPrefix + " -j " + threads
			os.system(command)

			print "\n\nCreate read length table\n\n"
			pearAssembly = pearPrefix + ".assembled.fastq"
			readLengthFile = pearPrefix + ".assembled.read_length"
			
			fastq_parser = SeqIO.parse(pearAssembly, "fastq")
			
			outHandle = open(readLengthFile, 'w')
			text = "Merged.Read\tLength\n"
			outHandle.write(text)

			for fastq in fastq_parser:
				readName = fastq.id
				readLength = len(str(fastq.seq))
				
				text = readName + "\t" +str(readLength)+ "\n"
				outHandle.write(text)

