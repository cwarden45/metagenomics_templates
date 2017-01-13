import sys
import re
import os
from Bio import SeqIO

parameterFile = "parameters.txt"
finishedSamples = []

def mothurClassStats(assignmentFile, outputfile, threshold, taxLevel):
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
		totalReadCount += 1
		assignmentText = lineInfo[1]
		assignmentInfo = assignmentText.split(";")
	
		if (assignmentInfo[0] == "Bacteria(100)"):
			infoIndex = -1
			if taxLevel == "order":
				infoIndex = 3
			elif taxLevel == "family":
				infoIndex = 4
			elif taxLevel == "genus":
				infoIndex = 5
			else:
				print "Need to define indices for " + taxLevel
				sys.exit()
			
			if assignmentInfo[infoIndex] != "unclassified":
				reResult = re.search("(.*)\((\d+)\)",assignmentInfo[infoIndex])
				if reResult:
					className = reResult.group(1)
					emptyResult = re.search("_unclassified$",className)
					
					classScore = int(reResult.group(2))
					rescaleScore = float(classScore)/100

					if (classScore >= threshold) and (not emptyResult):
						classCounts += 1
						if className in classHash:
							classHash[className] = classHash[className] + 1
						else:
							classHash[className] = 1
				else:
					print "problem with " + assignmentInfo[infoIndex]
					sys.exit()
		line = inHandle.readline()	

	outHandle = open(outputfile, "w")
	text = "Assignment\tRead.Num\tTotal.Percent\tClassified.Percent\n"
	outHandle.write(text)
	
	for assignment in classHash.keys():
		readCount = classHash[assignment]
		totalPercent = 100 * float(readCount)/float(totalReadCount)
		classPercent = 100 * float(readCount)/float(classCounts)
		
		text = assignment + "\t" + str(readCount) + "\t" + '{0:.6g}'.format(totalPercent) + "\t" + '{0:.6g}'.format(classPercent) + "\n"
		outHandle.write(text)
	
	percentClass = 100 * float(classCounts) / float(totalReadCount)
	return([str(totalReadCount),str(classCounts), '{0:.3g}'.format(percentClass)])

def rdpClassStats(assignmentFile, outputfile, threshold, taxLevel):
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
			tax = lineInfo[3]
			
			if tax != "NA":
				classifiedCounts += 1
				
				if tax in taxHash:
					taxHash[tax] = taxHash[tax] + 1
				else:
					taxHash[tax] = 1
			
		line = inHandle.readline()

	totalReadCount = lineCount - 1
	outHandle = open(outputfile, "w")
	text = "Assignment\tRead.Num\tTotal.Percent\tClassified.Percent\n"
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

				className = lineInfo[0]
			
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

				className = lineInfo[0]
				abIndex = -1
				if normMethod == "classified_ab":
					abIndex = 3
				elif normMethod == "total_ab":
					abIndex = 2
				elif normMethod == "raw_counts":
					abIndex = 1
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
	text = "Assignment\t" + "\t".join(sampleIDs) + "\n";
	outHandle.write(text)
	
	for ID in abDict:
		text = ID + abDict[ID] + "\n"
		outHandle.write(text)

	
mem = ""
quantFolder = ""
classifier = ""
minCycles = 3
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

	if param == "Min_CCS_Cycles":
		minCycles = int(value)
	
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

if (minCycles == "") or (minCycles == "[required]"):
	print "Need to enter a value for 'Min_CCS_Cycles'!"
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
	
processedIDs = []
processedFiles = []
	
statHandle = open(statsFile, 'w')
text = "Sample\tQC.Reads\tClassified.Reads\tPercent.Classified\tSample.Summary.File\n"
statHandle.write(text)
		
if classifier == "RDPclassifier":
	fileResults = os.listdir(quantFolder)
	
	for folder in fileResults:
		sample = folder
		classificationFolder = quantFolder + "/" + folder
		assignmentFile = classificationFolder  + "/"+sample+".ccs."+str(minCycles)+"x.RDP.out"
		
		if os.path.exists(assignmentFile):
			print sample
			
			summaryFile = quantFolder + "/" + sample + "_RDP_genus_abundance_count.txt"
			processedIDs.append(sample)
			processedFiles.append(summaryFile)
			
			text = folder + "\t" + sample + "\t"
			classStats = rdpClassStats(assignmentFile, summaryFile, 0.8,"genus")
			text = text + "\t".join(classStats)
			text = text + "%\t" + summaryFile + "\n"
			statHandle.write(text)
elif classifier == "BWA":
	fileResults = os.listdir(quantFolder)
	
	for folder in fileResults:
		sample = folder
		classificationFolder = quantFolder + "/" + folder
		assignmentFile = classificationFolder  + "/BWA_genus_hits.txt"
		
		if os.path.exists(assignmentFile):
			print sample
			
			summaryFile = quantFolder + "/" + sample + "_BWA_abundance_count.txt"
			processedIDs.append(sample)
			processedFiles.append(summaryFile)
			
			text = sample + "\t"
			classStats = bwaClassStats(assignmentFile, summaryFile)
			text = text + "\t".join(classStats)
			text = text + "%\t" + summaryFile + "\n"
			statHandle.write(text)
elif classifier == "mothur":
	fileResults = os.listdir(quantFolder)
	
	for folder in fileResults:
		sample = folder
		classificationFolder = quantFolder + "/" + folder
		assignmentFile = classificationFolder  + "/"+sample+".ccs."+str(minCycles)+"x.LENGTH_FILTERED.good.rdp.wang.taxonomy"
		
		if os.path.exists(assignmentFile):
			print sample
			
			summaryFile = quantFolder + "/" + sample + "_mothur_genus_abundance_count.txt"
			processedIDs.append(sample)
			processedFiles.append(summaryFile)
			
			text = sample + "\t"
			classStats = mothurClassStats(assignmentFile, summaryFile, 0.8,"genus")
			text = text + "\t".join(classStats)
			text = text + "%\t" + summaryFile + "\n"
			statHandle.write(text)

print "Creating Count Table"
abundanceTable(processedFiles, processedIDs, combined_counts_file, "raw_counts")
print "Creating Abundance Table"
abundanceTable(processedFiles, processedIDs, combined_abundance_file, "total_ab")
