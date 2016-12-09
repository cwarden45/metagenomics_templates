import sys
import re
import os

parameterFile = "parameters.txt"
finishedSamples = []

intermediateFolder = ""
sample_description_file = ""

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

inHandle.close()
		
if (intermediateFolder == "") or (intermediateFolder == "[required]"):
	print "Need to enter a value for 'Raw_CCS_Folder'!"
	sys.exit()
	
if (sample_description_file == "") or (sample_description_file == "[required]"):
	print "Need to enter a value for 'sample_description_file'!"
	sys.exit()
	
sampleIDindex = -1
baxH5index = -1

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
			elif colID == "baxH5prefix":
				baxH5index = i
	elif (sampleIDindex != -1) and (baxH5index != -1):
		sample = lineInfo[sampleIDindex]
		baxH5prefix = lineInfo[baxH5index]
		print sample
		
		baxH5c1 = baxH5prefix + ".1.bax.h5"
		baxH5c2 = baxH5prefix + ".2.bax.h5"
		baxH5c3 = baxH5prefix + ".3.bax.h5"

		bamFilePrefix = intermediateFolder + "/" + sample
		command = "bax2bam -o " + bamFilePrefix + " " + baxH5c1 + " " + baxH5c2 + " " + baxH5c3
		os.system(command)
	else:
		print "Problem with line: " + line
		print "Does your sample description file have headers 'sampleID' and 'baxH5prefix'?"
		sys.exit()
	line = inHandle.readline()
