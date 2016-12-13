import sys
import re
import os

parameterFile = "parameters.txt"
finishedSamples = []

intermediateFolder = ""
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
		print sample
		
		unalignedBam = intermediateFolder + "/" + sample + ".subreads.bam"
		ccsBam = intermediateFolder + "/" + sample + ".ccs."+str(minCycles)+"x.bam"
		ccsStats = intermediateFolder + "/" + sample + ".ccs."+str(minCycles)+"x.ccs_report.txt"
		command = "ccs --minPasses=" + str(minCycles) + " --reportFile=" + ccsStats +" " + unalignedBam + " " + ccsBam
		os.system(command)

		fastq = fastqFolder + "/" + sample + ".ccs."+str(minCycles)+"x.fastq"
		command = "/opt/samtools-1.3/samtools bam2fq " + ccsBam + " > " + fastq
		os.system(command)		

	else:
		print "Problem with line: " + line
		print "Does your sample description file have headers 'sampleID'?"
		sys.exit()
	line = inHandle.readline()
