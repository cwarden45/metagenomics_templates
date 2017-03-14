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

mergedReadsFolder = readsFolder + "/mothur_merged"
if not os.path.exists(mergedReadsFolder):
	command = "mkdir " + mergedReadsFolder
	os.system(command)
	
fileResults = os.listdir(readsFolder)
for file in fileResults:
	#if starting from .fastq files, early .gz creation doesn't make sense if you're only running mothur, but I've left this since .gz file would already be created from RDP pipeline
	#however, this strategy is slightly better if staring with MiSeq .fastq.gz files
	result = re.search("(.*_L\d{3}_R1_001).fastq$",file)
	if result:
		command = "gzip " + file
		os.system(command)
		file = file + ".gz"
	
	resultGZ = re.search("(.*_L\d{3}_R1_001).fastq.gz$",file)
	
	if resultGZ:
		sample = resultGZ.group(1)
		if sample not in finishedSamples:
			print sample
			
			read1 = re.sub(".gz","",readsFolder + "/" + file)
			read2 = re.sub("R1_001.fastq","R2_001.fastq",read1)

			gunzipFlag = 0
			if not os.path.isfile(read1) or not os.path.isfile(read2):
				gunzipFlag = 1
				
				if  os.path.isfile(read1 + ".gz"):
					command = "gunzip -c " + read1 + ".gz > " + read1
					os.system(command)
				else:
					print "Can't find " + read1 + " (or .gz version)."
					sys.exit()

				if  os.path.isfile(read2 + ".gz"):
					command = "gunzip -c " + read2 + ".gz > " + read2
					os.system(command)
				else:
					print "Can't find " + read2 + " (or .gz version)."
					sys.exit()
					
			print "\n\nCreate Contigs\n\n"
			command = "mothur \"#make.contigs(ffastq="+read1+", rfastq="+read2+", processors="+threads+")\""
			os.system(command)

			if gunzipFlag == 1:
				command = "rm " + read1
				os.system(command)
				command = "rm " + read2
				os.system(command)
			else:
				command = "gzip " + read1
				os.system(command)
				command = "gzip " + read2
				os.system(command)
				
			scrap = re.sub(".fastq$","",read1) + ".scrap.contigs.fasta"
			command = "rm " + scrap
			os.system(command)
			scrap = re.sub(".fastq$","",read1) + ".scrap.contigs.qual"
			command = "rm " + scrap
			os.system(command)
			
			contigFA = re.sub(".fastq$","",read1) + ".trim.contigs.fasta"
			command = "mv " + contigFA + " " + mergedReadsFolder + "/" 
			os.system(command)
			
			contigQual = re.sub(".fastq$","",read1) + ".trim.contigs.qual"
			command = "mv " + contigQual + " " + mergedReadsFolder + "/" 
			os.system(command)
	
			contigReport = re.sub(".fastq$","",read1) + ".contigs.report"
			command = "mv " + contigReport + " " + mergedReadsFolder + "/" 
			os.system(command)
	
			print "\n\nRemove Degenerate Nucleotide Reads\n\n"
			contigFA = mergedReadsFolder + "/" + re.sub(".fastq.gz$","",file)+ ".trim.contigs.fasta"
			command = "mothur \"#screen.seqs(fasta="+contigFA+", processors="+threads+", maxambig=0)\""
			os.system(command)
			
			print "\n\nCreate Read Length File\n\n"
			#still do this to get quality-filtered read lengths
			filteredContigFA = re.sub(".fasta$",".good.fasta",contigFA)
			readLengthFile = mergedReadsFolder + "/" + re.sub(".fastq.gz$","",file) + ".trim.contigs.good.read_length"
			
			fasta_parser = SeqIO.parse(filteredContigFA, "fasta")
			
			outHandle = open(readLengthFile, 'w')
			text = "Merged.Read\tLength\n"
			outHandle.write(text)

			for fasta in fasta_parser:
				readName = fasta.id
				readLength = len(str(fasta.seq))
				
				text = readName + "\t" +str(readLength)+ "\n"
				outHandle.write(text)

#cleanup log files
fileResults = os.listdir(".")
for file in fileResults:
	result = re.search(".logfile$",file)
	
	if result:
		print "deleting " + file
		command = "rm " + file
		os.system(command)
