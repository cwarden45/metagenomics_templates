import sys
import re
import os
from Bio import SeqIO

parameterFile = "parameters.txt"
finishedSamples = []

mem = ""
threads = ""
quantFolder = ""
readsFolder = ""
RDPclassifier = ""
classifier = ""
description_file = ""
statsFile = ""
bwaRef = ""
mothurRef = ""
mothurTax = ""

inHandle = open(parameterFile)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	param = lineInfo[0]
	value = lineInfo[1]

	if param == "Threads":
		threads = value
	
	if param == "mothur_ref":
		mothurRef = value

	if param == "mothur_tax":
		mothurTax = value
		
	if param == "BWA_Ref":
		bwaRef = value
		
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
shortNameHash = {}

lineCount = 0
sampleIndex = -1
minIndex = -1
maxIndex = -1
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
			if lineInfo[i] == "Min.Length":
				minIndex = i
			elif lineInfo[i] == "Max.Length":
				maxIndex = i
			elif lineInfo[i] == "sampleID":
				sampleIndex = i
			elif lineInfo[i] == "userID":
				nameIndex = i
	else:
		key = lineInfo[sampleIndex]
		#might need to uncomment with RDPclassifier code
		#key = re.sub("_L001_R1_001","",key)
		min = int(lineInfo[minIndex])
		max = int(lineInfo[maxIndex])
		name = lineInfo[nameIndex]
		
		shortNameHash[key] = name
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
				
				gunzipFlag = 0
				if not os.path.isfile(pearAssembly):
					if os.path.isfile(pearAssembly + ".gz"):
						command = "gunzip -c " + pearAssembly + ".gz > " + pearAssembly
						os.system(command)
						gunzipFlag = 1
					else:
						print "Can't find " + pearAssembly + " (or .gz version)."
						print "Did you run PEAR via preprocess_RDPclassifier_or_BWA.py?"
						sys.exit()
				
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
				text = sample + "\t" + str(readCount)+ "\t"+ str(passCount) +"\t" + '{0:.2f}'.format(lengthPercent) + "\n"
				statHandle.write(text)
						
				if gunzipFlag == 1:
					command = "rm " + pearAssembly
					os.system(command)
				else:
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
				
elif classifier == "BWA":
	if (bwaRef== "") or (bwaRef == "[required]"):
		print "Need to enter a value for 'BWA_Ref'!"
		sys.exit()
		
	if (threads == "") or (threads == "[required]"):
		print "Need to enter a value for 'Threads'!"
		sys.exit()

	mergedReadsFolder = readsFolder + "/PEAR_Merged"
	fileResults = os.listdir(readsFolder)
	for file in fileResults:
		resultGZ = re.search("(.*)_L\d{3}_R1_001.fastq.gz$",file)
		
		if resultGZ:
			sample = resultGZ.group(1)
			if sample not in finishedSamples:
				print sample

				print "Length-Filter Reads"
				print "############################################################################"
				print "NOTE: This is meant more for validation than original quantification"
				print "\t\tso, length-filter file not provided"
				print "############################################################################"
				pearAssembly = mergedReadsFolder + "/" + sample + ".assembled.fastq"
				bwaRead = mergedReadsFolder + "/" + sample + ".assembled.LENGTH_FILTERED.fastq"
				minLength = minHash[sample]
				maxLength = maxHash[sample]
				
				gunzipFlag = 0
				if not os.path.isfile(pearAssembly):
					if os.path.isfile(pearAssembly + ".gz"):
						command = "gunzip -c " + pearAssembly + ".gz > " + pearAssembly
						os.system(command)
						gunzipFlag = 1
					else:
						print "Can't find " + pearAssembly + " (or .gz version)."
						print "Did you run PEAR via preprocess_RDPclassifier_or_BWA.py?"
						sys.exit()
				
				if not os.path.isfile(bwaRead):
					def length_filter(records, minLen, maxLen):
						for rec in records:
							readLength = len(str(rec.seq))
							
							if (readLength >= minLen) & (readLength <= maxLen):
								yield rec
					
					fastq_parser = SeqIO.parse(pearAssembly, "fastq") 
					SeqIO.write(length_filter(fastq_parser, minLength, maxLength), bwaRead, "fastq")
				
				if gunzipFlag == 1:
					command = "rm " + pearAssembly
					os.system(command)
				else:
					command = "gzip " + pearAssembly
					os.system(command)
					
				print "BWA-Alignment"
				classificationFolder = quantFolder + "/" + sample
				if not os.path.exists(classificationFolder):
					command = "mkdir " + classificationFolder
					os.system(command)

				alnSam = classificationFolder + "/aligned.sam"
				command = "/opt/bwa/bwa mem -t "+ str(threads) + " " + bwaRef+ " " + bwaRead + " > " + alnSam
				os.system(command)

				command = "gzip " + bwaRead
				os.system(command)
				
				alnBam = classificationFolder + "/aligned.bam"
				command = "/opt/samtools-1.3/samtools view -b -F 2048 " + alnSam + " > " + alnBam
				os.system(command)

				command = "rm " + alnSam
				os.system(command)

				sortBam= quantFolder + "/" + sample + ".bam"
				command = "/opt/samtools-1.3/samtools sort " + alnBam + " -o " + sortBam
				os.system(command)

				command = "rm " + alnBam
				os.system(command)

				command = "/opt/samtools-1.3/samtools index " + sortBam
				os.system(command)

				statsFile = classificationFolder + "/alignment_stats.txt"
				command = "/opt/samtools-1.3/samtools flagstat " + sortBam + " > " + statsFile
				os.system(command)

				statsFile = classificationFolder + "/idxstats.txt"
				command = "samtools idxstats " + sortBam + " > " + statsFile
				os.system(command)

				nameSam= classificationFolder + "/name.sort.sam"
				command = "/opt/samtools-1.3/samtools sort -n " + sortBam + " -O sam -o " + nameSam
				os.system(command)
				
				print "Collapse by RDP Genus"
				#annotation hash
				inHandle = open(bwaRef)
				line = inHandle.readline()

				annHash = {}

				while line:
					line = line.replace("\n","")
					line = line.replace("\r","")

					headerResult = re.search("^>",line)
					
					if headerResult:
						lineInfo = line.split("\t")
						id = lineInfo[0]
						id = re.sub(">","",id)
						tax = lineInfo[1]
						tax = re.sub("\"","",tax)
						
						bacResult = re.search("Root;Bacteria",tax)
						
						if bacResult:
							annHash[id] = tax
							
					line = inHandle.readline()

				inHandle.close()

				#classify reads
				genusHits = classificationFolder + "/BWA_genus_hits.txt"
				minMatchLength = 0.8 * minLength
				
				outHandle = open(genusHits, "w")
				text = "Read\tHit\tHit.Length\tTax.String\n"
				outHandle.write(text)

				inHandle = open(nameSam)
				line = inHandle.readline()

				totalCount = 0
				passCount = 0

				while line:
					line = line.replace("\n","")
					line = line.replace("\r","")

					commentResult = re.search("^@",line)
					
					if not commentResult:
						totalCount += 1
						lineInfo = line.split("\t")
						read = lineInfo[0]
						hit = lineInfo[2]
						cigar = lineInfo[5]
						tax = "NA"
									
						matchSum = 0
						for m in re.finditer("(\d+)M",cigar):
							matchSum += int(m.group(1))
						if matchSum > minMatchLength:
							if hit in annHash:
								tax = annHash[hit]
							passCount += 1
						
						text = read + "\t" + hit + "\t" + str(matchSum) + "\t" + tax + "\n"
						outHandle.write(text)
						
					line = inHandle.readline()
				
				command = "rm " + nameSam
				os.system(command)
elif classifier == "mothur":
	if (mothurRef== "") or (mothurRef == "[required]"):
		print "Need to enter a value for 'mothur_ref'!"
		sys.exit()
		
	if (mothurTax == "") or (mothurTax == "[required]"):
		print "Need to enter a value for 'mothur_tax'!"
		sys.exit()

	if (threads == "") or (threads == "[required]"):
		print "Need to enter a value for 'Threads'!"
		sys.exit()
		
	mergedReadsFolder = readsFolder + "/mothur_merged"
	fileResults = os.listdir(readsFolder)
	for file in fileResults:
		resultGZ = re.search("(.*_L\d{3}_R1_001).fastq.gz$",file)
		
		if resultGZ:
			sample = resultGZ.group(1)
			outName = shortNameHash[sample]
			if (sample not in finishedSamples) and (outName not in finishedSamples):
				print sample
				
				readCount = 0
				passCount = 0
				
				print "\n\nLength-Filter mothur-screened FASTA\n\n"
				mothurContigs = mergedReadsFolder + "/" + sample + ".trim.contigs.good.fasta"
				mothurInput = mergedReadsFolder + "/" + outName + ".trim.contigs.good.LENGTH_FILTERED.fasta"
				
				minLength = minHash[sample]
				maxLength = maxHash[sample]
				
				fastq_parser = SeqIO.parse(mothurContigs, "fasta")
				outHandle = open(mothurInput, 'w')

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
				text = sample + "\t" + str(readCount)+ "\t"+ str(passCount) +"\t" + '{0:.2f}'.format(lengthPercent) + "\n"
				statHandle.write(text)
				
				print "\n\nApply mothur classifier\n\n"
				classificationFolder = quantFolder + "/" + outName
				if not os.path.exists(classificationFolder):
					command = "mkdir " + classificationFolder
					os.system(command)

				command = "mothur \"#classify.seqs(fasta="+mothurInput+", reference="+mothurRef+", processors="+threads+", taxonomy="+mothurTax+", cutoff=80)\""
				os.system(command)
				
				taxOut = re.sub(".fasta",".rdp.wang.tax.summary",mothurInput)
				command = "mv " + taxOut + " " + classificationFolder + "/"
				os.system(command)

				taxOut = re.sub(".fasta",".rdp.wang.taxonomy",mothurInput)
				command = "mv " + taxOut + " " + classificationFolder + "/"
				os.system(command)

				taxOut = re.sub(".fasta",".rdp.wang.flip.accnos",mothurInput)
				command = "mv " + taxOut + " " + classificationFolder + "/"
				os.system(command)
				
				#cleanup log files
				fileResults = os.listdir(".")
				for file in fileResults:
					result = re.search(".logfile$",file)
					
					if result:
						print "deleting " + file
						command = "rm " + file
						os.system(command)
else:
	print "classifier must be 'RDPclassifier', 'mothur', or 'BWA'"
	sys.exit()
