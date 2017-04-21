import sys
import re
import os
from Bio import SeqIO

parameterFile = "parameters.txt"
finishedSamples = []

def bwaClassStats(assignmentFile, outputfile):
	totalReadCount = 0
	classifiedCounts = 0
	taxHash = {}

	inHandle = open(assignmentFile)
	line = inHandle.readline()

	lineCount = 0

	while line:
		line = line.replace("\n","")
		line = line.replace("\r","")

		lineCount += 1
		
		if lineCount  > 1:
			totalReadCount +=1
			lineInfo = line.split("\t")
			tax = lineInfo[1] + "\t" + lineInfo[3]
			
			if tax != "NA":
				classifiedCounts += 1
				
				if tax in taxHash:
					taxHash[tax] = taxHash[tax] + 1
				else:
					taxHash[tax] = 1
			
		line = inHandle.readline()

	totalReadCount = lineCount - 1
	outHandle = open(outputfile, "w")
	text = "Accession\tTaxonomy\tRead.Num\tTotal.Percent\tClassified.Percent\n"
	outHandle.write(text)

	for tax in taxHash:
		readCount = taxHash[tax]
		totalPercent = 100 * float(readCount)/float(totalReadCount)
		classPercent = 100 * float(readCount)/float(classifiedCounts)
			
		text = tax + "\t" + str(readCount) + "\t" + '{0:.2g}'.format(totalPercent) + "\t" + '{0:.2g}'.format(classPercent) + "\n"
		outHandle.write(text)
		
	percentClass = 100 * float(classifiedCounts) / float(totalReadCount)
	return([str(totalReadCount),str(classifiedCounts), '{0:.3g}'.format(percentClass)])

def abundanceTable(inputArr, sampleIDs, outputfile, normMethod):
	abDict={}

	for i in xrange(0, len(sampleIDs)):
		sample = sampleIDs[i]
		inputfile = inputArr[i]

		inHandle = open(inputfile)
		line = inHandle.readline()
	
		lineCount = 0
		
		while line:	
			lineCount += 1
			if lineCount > 1:
				line = line.replace("\n","")
				line = line.replace("\r","")
	
				lineInfo = line.split("\t")

				className = lineInfo[0] + "\t" +  lineInfo[1]
			
				abDict[className] = ""
				
			line = inHandle.readline()
	
	for i in xrange(0, len(sampleIDs)):
		sample = sampleIDs[i]
		inputfile = inputArr[i]
		print sample
	
		inHandle = open(inputfile)
		line = inHandle.readline()
	
		classNum = 0
		totalCounts = 0
		
		lineCount = 0
		
		tempDict = {}
		
		while line:	
			lineCount += 1
			if lineCount > 1:
				line = line.replace("\n","")
				line = line.replace("\r","")
	
				lineInfo = line.split("\t")

				className = lineInfo[0] + "\t" +  lineInfo[1]
				abIndex = -1
				if normMethod == "classified_ab":
					abIndex = 4
				elif normMethod == "total_ab":
					abIndex = 3
				elif normMethod == "raw_counts":
					abIndex = 2
				else:
					print "4th param in 'abundanceTable' function must be 'classified_ab','total_ab',or'raw_counts'"
					sys.exit()
				freq = lineInfo[abIndex]
				
				if normMethod == "raw_counts":
					tempDict[className] = freq
				else:
					tempDict[className] = float(freq)
				
			line = inHandle.readline()
	
		for ID in abDict:
			if ID in tempDict:
				percentAb = tempDict[ID]
				if normMethod == "raw_counts":
					abDict[ID] =  abDict[ID] + "\t" + str(percentAb)
				else:
					abDict[ID] =  abDict[ID] + "\t" + '{0:.2g}'.format(percentAb)
			else:
				abDict[ID] = abDict[ID] + "\tNA"
	
	outHandle = open(outputfile, "w")
	text = "Accession\tTaxonomy\t" + "\t".join(sampleIDs) + "\n";
	outHandle.write(text)
	
	for ID in abDict:
		text = ID + abDict[ID] + "\n"
		outHandle.write(text)

	
mem = ""
quantFolder = ""
classifier = ""
description_file = ""
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
longNameHash = {}

lineCount = 0
sampleIndex = -1
nameIndex = -1

processedIDs = []
processedFiles = []

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
		#depending upon how .fastq files are formatted, you can change subsequent line
		key = re.sub("_L001_R1_001","",key)
		name = lineInfo[nameIndex]
		
		shortNameHash[key] = name
		longNameHash[name] = key
		
if classifier == "BWA":
	fileResults = os.listdir(quantFolder)
	
	for folder in fileResults:	
		classificationFolder = quantFolder + "/" + folder
		assignmentFile = classificationFolder  + "/BWA_genus_hits.txt"
		
		if os.path.exists(assignmentFile):
			sample = shortNameHash[folder]
			print sample
			
			summaryFile = classificationFolder + "/" + sample + "_PEAR_BWA_abundance_count_Rep_Ref_Seq.txt"
			bwaClassStats(assignmentFile, summaryFile)
			processedIDs.append(sample)
			processedFiles.append(summaryFile)
elif classifier == "BLAST":
	fileResults = os.listdir(quantFolder)
		
	for folder in fileResults:
		sample = folder
		classificationFolder = quantFolder + "/" + folder
		assignmentFile = classificationFolder  + "/BLAST_genus_hits.txt"
		
		if os.path.exists(assignmentFile):
			sample = shortNameHash[folder]
			print sample
			
			summaryFile = classificationFolder + "/" + sample + "_BLAST_abundance_count_Rep_Ref_Seq.txt"
			bwaClassStats(assignmentFile, summaryFile)
			processedIDs.append(sample)
			processedFiles.append(summaryFile)
else:
	print "Sorry - code only works with custom alignment with 'BWA' or 'BLAST'"
	sys.exit()

combined_counts_file = re.sub(".txt$","_Rep_Ref_Seq.txt",combined_counts_file)
combined_abundance_file = re.sub(".txt$","_Rep_Ref_Seq.txt",combined_abundance_file)
			
abundanceTable(processedFiles, processedIDs, combined_counts_file, "raw_counts")
abundanceTable(processedFiles, processedIDs, combined_abundance_file, "total_ab")
