import sys
import re
import os
from Bio import SeqIO

parameterFile = "parameters.txt"
finishedSamples = []

def rdpClassStats(assignmentFile, outputfile, threshold, taxLevel):
	#class counts
	totalReadCount = 0
	classHash = {}
	classCounts = 0
	
	bacteriaSkippedLines = 0
	
	inHandle = open(assignmentFile)
	line = inHandle.readline()
	
	while line:
		line = line.replace("\n","")
		line = line.replace("\r","")
	
		lineInfo = line.split("\t")
		readName = lineInfo[0]
		totalReadCount +=1
	
		if (lineInfo[5] == "Bacteria"):
			if len(lineInfo) != 23:
				bacteriaSkippedLines += 1
			else:
				nameIndex = -1
				scoreIndex = -1
					
				if taxLevel == "order":
					nameIndex = 14
					scoreIndex = 16
				elif taxLevel == "family":
					nameIndex = 17
					scoreIndex = 19
				elif taxLevel == "genus":
					nameIndex = 20
					scoreIndex = 22
				else:
					print "Need to define indices for " + taxLevel
					sys.exit()
				
				className = lineInfo[nameIndex]
				classScore = float(lineInfo[scoreIndex])
				
				if classScore >= threshold:
					classCounts += 1
					if className in classHash:
						classHash[className] = classHash[className] + 1
					else:
						classHash[className] = 1
		line = inHandle.readline()	
	
	print str(bacteriaSkippedLines) + " read(s) have weird taxonomy annotations..."
	
	percentClass = 100 * float(classCounts) / float(totalReadCount)
	print "Bacterial Classified Reads : " + str(classCounts) + " ( " + str(percentClass) + " )"

	outHandle = open(outputfile, "w")
	text = "Assignment\tRead.Num\tTotal.Percent\tClassified.Percent\n"
	outHandle.write(text)
	
	for assignment in classHash.keys():
		readCount = classHash[assignment]
		totalPercent = 100 * float(readCount)/float(totalReadCount)
		classPercent = 100 * float(readCount)/float(classCounts)
		
		text = assignment + "\t" + str(readCount) + "\t" + '{0:.6g}'.format(totalPercent) + "\t" + '{0:.6g}'.format(classPercent) + "\n"
		outHandle.write(text)
		
	return([str(totalReadCount),str(classCounts), '{0:.3g}'.format(percentClass)])

mem = ""
quantFolder = ""
classifier = ""
description_file = ""
statsFile = ""
combined_counts_file = ""
combined_abundance_file = ""

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
	
	if param == "Classifier":
		classifier = value

	if param == "counts_file":
		combined_counts_file = value

	if param == "abundance_file":
		combined_abundance_file = value
		
	if param == "classified_stats_file":
		statsFile = value
		
	if param == "Java_Mem":
		mem = value
		
	if param == "Threads":
		threads = value
		
	if param == "Classification_Folder":
		quantFolder = value

if (description_file== "") or (description_file == "[required]"):
	print "Need to enter a value for 'sample_description_file'!"
	sys.exit()

if (combined_counts_file== "") or (combined_counts_file == "[required]"):
	print "Need to enter a value for 'counts_file'!"
	sys.exit()
	
if (combined_abundance_file== "") or (combined_abundance_file == "[required]"):
	print "Need to enter a value for 'abundance_file'!"
	sys.exit()
	
if (statsFile== "") or (statsFile == "[required]"):
	print "Need to enter a value for 'classified_stats_file'!"
	sys.exit()

if (classifier== "") or (classifier == "[required]"):
	print "Need to enter a value for 'Classifier'!"
	sys.exit()
	
if (quantFolder == "") or (quantFolder == "[required]"):
	print "Need to enter a value for 'Classification_Folder'!"
	sys.exit()

#define min / max length per target / sample
shortNameHash = {}

lineCount = 0
sampleIndex = -1
nameIndex = -1

inHandle = open(description_file)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	
	lineCount += 1
	
	if lineCount == 1:
		for i in range(0,len(lineInfo)):
			if lineInfo[i] == "userID":
				nameIndex = i
			elif lineInfo[i] == "sampleID":
				sampleIndex = i
	else:
		key = lineInfo[sampleIndex]
		name = lineInfo[nameIndex]
		
		shortNameHash[key] = name

if os.path.isfile(statsFile):
	statHandle = open(statsFile, 'a')
else:
	statHandle = open(statsFile, 'w')
	text = "SampleID\tuserID\tQC.Reads\tClassified.Reads\tPercent.Classified\tSample.Summary.File\n"
	statHandle.write(text)
		
if classifier == "RDPclassifier":
	fileResults = os.listdir(quantFolder)
	
	for folder in fileResults:	
		classificationFolder = quantFolder + "/" + folder
		assignmentFile = classificationFolder  + "/PEAR.RDP.out"
		
		if os.path.exists(assignmentFile):
			sample = shortNameHash[folder]
			print sample
			
			summaryFile = quantFolder + "/" + sample + "_PEAR_RDP_genus_abundance_count.txt"
			
			text = folder + "\t" + sample + "\t"
			classStats = rdpClassStats(assignmentFile, summaryFile, 0.8,"genus")
			text = text + "\t".join(classStats)
			text = text + "\t" + summaryFile + "\n"
			statHandle.write(text)

print "Need to write code for combined counts and abundance files"