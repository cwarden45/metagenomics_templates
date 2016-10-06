import sys
import re
import os
from Bio import SeqIO

parameterFile = "parameters.txt"
finishedSamples = []

mem = ""
quantFolder = ""
readsFolder = ""
RDPclassifier = ""
classifier = ""
description_file = ""
statsFile = ""

inHandle = open(parameterFile)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	param = lineInfo[0]
	value = lineInfo[1]

	if param == "sample_description_file":
		description_file = value
	
	if param == "total_counts_file":
		statsFile = value
	
	if param == "Classifier":
		classifier = value

	if param == "Reads_Folder":
		readsFolder = value
	
	if param == "Classification_Folder":
		quantFolder = value
		
	if param == "Java_Mem":
		mem = value
		
	if param == "Threads":
		threads = value
		
	if param == "RDPclassifier_Jar":
		RDPclassifier = value

if (description_file== "") or (description_file == "[required]"):
	print "Need to enter a value for 'sample_description_file'!"
	sys.exit()
	
if (statsFile== "") or (statsFile == "[required]"):
	print "Need to enter a value for 'total_counts_file'!"
	sys.exit()

if (classifier== "") or (classifier == "[required]"):
	print "Need to enter a value for 'Classifier'!"
	sys.exit()

if (readsFolder == "") or (readsFolder == "[required]"):
	print "Need to enter a value for 'Reads_Folder'!"
	sys.exit()
	
if (quantFolder == "") or (quantFolder == "[required]"):
	print "Need to enter a value for 'Classification_Folder'!"
	sys.exit()

#define min / max length per target / sample
minHash = {}
maxHash = {}

lineCount = 0
sampleIndex = -1
minIndex = -1
maxIndex = -1

inHandle = open(description_file)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	
	lineCount += 1
	
	if lineCount == 1:
		for i in range(0,len(lineInfo)):
			if lineInfo[i] == "Min.Length":
				minIndex = i
			elif lineInfo[i] == "Max.Length":
				maxIndex = i
			elif lineInfo[i] == "sampleID":
				sampleIndex = i
	else:
		key = lineInfo[sampleIndex]
		min = int(lineInfo[minIndex])
		max = int(lineInfo[maxIndex])
		
		minHash[key] = min
		maxHash[key] = max

if os.path.isfile(statsFile):
	statHandle = open(statsFile, 'a')
else:
	statHandle = open(statsFile, 'w')
	text = "Sample\tTotal.Reads\tLength.Filtered.Reads\tPercent.Length.Filtered\n"
	statHandle.write(text)
		
if classifier == "RDPclassifier":
	if (RDPclassifier== "") or (RDPclassifier == "[required]"):
		print "Need to enter a value for 'RDPclassifier_Jar'!"
		sys.exit()
		
	if (mem == "") or (mem == "[required]"):
		print "Need to enter a value for 'Java_Mem'!"
		sys.exit()

	mergedReadsFolder = readsFolder + "/PEAR_Merged"
	fileResults = os.listdir(readsFolder)
	for file in fileResults:
		resultGZ = re.search("(.*)_L\d{3}_R1_001.fastq.gz$",file)
		
		if resultGZ:
			sample = resultGZ.group(1)
			if sample not in finishedSamples:
				print sample
				
				readCount = 0
				passCount = 0
				
				print "\n\nFASTQ to FASTA (for length-filtered reads)\n\n"
				pearAssembly = mergedReadsFolder + "/" + sample + ".assembled.fastq"
				rdpInput = mergedReadsFolder + "/" + sample + ".assembled.LENGTH_FILTERED.fasta"
				
				minLength = minHash[sample]
				maxLength = maxHash[sample]
				
				fastq_parser = SeqIO.parse(pearAssembly, "fastq")
				
				outHandle = open(rdpInput, 'w')

				for fastq in fastq_parser:
					readName = fastq.id
					readLength = len(str(fastq.seq))
					
					readCount += 1
					
					if (readLength >= minLength) & (readLength <= maxLength):
						passCount += 1
						
						text = ">" + readName + "\n"
						text = text + str(fastq.seq)+ "\n"
						outHandle.write(text)
						
				lengthPercent = 100 * float(passCount) / float(readCount)
				text = sample + "\t" + str(readCount)+ "\t"+ str(passCount) +"\t" + '{0:.2g}'.format(lengthPercent) + "\n"
				statHandle.write(text)
						
				command = "gzip " + pearAssembly
				os.system(command)
				
				print "\n\nApply RDPclassifier\n\n"
				classificationFolder = quantFolder + "/" + sample
				if not os.path.exists(classificationFolder):
					command = "mkdir " + classificationFolder
					os.system(command)

				treeOut = classificationFolder + "/PEAR.RDP.fastahier.txt"
				percentOut = classificationFolder + "/PEAR.RDP.out"
				command = "java -jar -Xmx" + mem + " " + RDPclassifier + " -o " + percentOut + " -h " + treeOut + " " + rdpInput
				os.system(command)
