import sys
import re
import os
from Bio import SeqIO

parameterFile = "parameters.txt"
finishedSamples = []

threads = ""
readsFolder = ""
minCycles = 3
sample_description_file = ""

inHandle = open(parameterFile)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	param = lineInfo[0]
	value = lineInfo[1]

	if param == "Reads_Folder":
		readsFolder = value
		
	if param == "Threads":
		threads = value

	if param == "sample_description_file":
		sample_description_file = value

	if param == "Min_CCS_Cycles":
		minCycles = int(value)
		
if (readsFolder == "") or (readsFolder == "[required]"):
	print "Need to enter a value for 'Reads_Folder'!"
	sys.exit()
		
if (threads == "") or (threads == "[required]"):
	print "Need to enter a value for 'Threads'!"
	sys.exit()

if (minCycles == "") or (minCycles == "[required]"):
	print "Need to enter a value for 'Min_CCS_Cycles'!"
	sys.exit()

if (sample_description_file == "") or (sample_description_file == "[required]"):
	print "Need to enter a value for 'sample_description_file'!"
	sys.exit()
	
sampleIDindex = -1

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
	elif sampleIDindex != -1:
		sample = lineInfo[sampleIDindex]
		
		if sample not in finishedSamples:
			print sample
			
			fasta = readsFolder + "/" + sample + ".ccs."+str(minCycles)+"x.fasta"

			print "\n\nRemove Degenerate Nucleotide Reads\n\n"
			command = "mothur \"#screen.seqs(fasta="+fasta+", processors="+threads+", maxambig=0)\""
			os.system(command)
			
			print "\n\nCreate Read Length File\n\n"
			#still do this to get quality-filtered read lengths
			filteredContigFA = re.sub(".fasta$",".good.fasta",fasta)
			readLengthFile = readsFolder + "/" +sample + ".ccs."+str(minCycles)+"x.good.read_length"
			
			fasta_parser = SeqIO.parse(filteredContigFA, "fasta")
			
			outHandle = open(readLengthFile, 'w')
			text = "Merged.Read\tLength\n"
			outHandle.write(text)

			for fasta in fasta_parser:
				readName = fasta.id
				readLength = len(str(fasta.seq))
				
				text = readName + "\t" +str(readLength)+ "\n"
				outHandle.write(text)
	else:
		print "Problem with line: " + line
		print "Does your sample description file have headers 'sampleID'?"
		sys.exit()
	line = inHandle.readline()
				
#cleanup log files
fileResults = os.listdir(".")
for file in fileResults:
	result = re.search(".logfile$",file)
	
	if result:
		print "deleting " + file
		command = "rm " + file
		os.system(command)
