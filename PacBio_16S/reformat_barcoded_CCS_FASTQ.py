import sys
import re
import os
from Bio import SeqIO

parameterFile = "parameters.txt"
finishedSamples = []

sample_description_file = ""
fastqFolder = ""
minCycles = 3

inHandle = open(parameterFile)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	param = lineInfo[0]
	value = lineInfo[1]

	if param == "Raw_CCS_Folder":
		intermediateFolder = value

	if param == "sample_description_file":
		sample_description_file = value

	if param == "Min_CCS_Cycles":
		minCycles = int(value)

	if param == "Reads_Folder":
		fastqFolder = value
		
inHandle.close()

if (fastqFolder == "") or (fastqFolder == "[required]"):
	print "Need to enter a value for 'Reads_Folder'!"
	sys.exit()

if (minCycles == "") or (minCycles == "[required]"):
	print "Need to enter a value for 'Min_CCS_Cycles'!"
	sys.exit()
	
if (intermediateFolder == "") or (intermediateFolder == "[required]"):
	print "Need to enter a value for 'Raw_CCS_Folder'!"
	sys.exit()
	
if (sample_description_file == "") or (sample_description_file == "[required]"):
	print "Need to enter a value for 'sample_description_file'!"
	sys.exit()
	
sampleIDindex = -1
inFastqIndex = -1

inHandle = open(sample_description_file)
line = inHandle.readline()

lineCount = 0

while line:
	lineCount += 1
	
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	lineInfo = line.split("\t")
	
	if lineCount == 1:
		for i in range(0,len(lineInfo)):
			colID = lineInfo[i]
			
			if colID == "sampleID":
				sampleIDindex = i
			if colID == "baxH5prefix":
				#leave same column name as one sample per run code
				inFastqIndex = i
	elif sampleIDindex != -1 and inFastqIndex != -1:
		sample = lineInfo[sampleIDindex]
		pbbarcodeFastq = lineInfo[inFastqIndex]
		print sample

		fastq = fastqFolder + "/" + sample + ".ccs."+str(minCycles)+"x.fastq"
		command = "cp " + pbbarcodeFastq + " " + fastq
		os.system(command)	

		fasta = fastqFolder + "/" + sample + ".ccs."+str(minCycles)+"x.fasta"
		outHandle = open(fasta,"w")
		fastq_sequences = SeqIO.parse(open(fastq),'fastq')
		for fasta in fastq_sequences:
			readName = fasta.id
			readSequence = str(fasta.seq)
			
			text = ">" + readName + "\n" + readSequence + "\n"
			outHandle.write(text)
	elif inFastqIndex != -1:
		print "Problem with line: " + line
		print "Missing column specifying paths to demultiplexed .fastq files (can be in different folders from different runs)"
		print "Does your sample description file have headers 'baxH5prefix' (uses column label for pipeline without demultiplexing)?"		
	else:
		print "Problem with line: " + line
		print "Does your sample description file have headers 'sampleID'?"
		sys.exit()
	line = inHandle.readline()
