import os
import sys
import re
from Bio.Seq import Seq

finishedSamples = []

inputFolder = "../Reads"
outputFolder = "../Reads/Cutadapt_Filtered_Reads"
includeAdapter = True

#make sure sample name goes all the way to .fastq.gz (or .fastq, if you modify the code)
#so, you might need a slightly different sample description file with RDPclassifier / BWA / BLAST code
descriptionFile = "sample_description.txt"

#universal adapter sequence - may want to check if there is a more specific sequence (but this seems to work OK)
forwardAdapter = "TACACTCTTTCCCTACACGACGCTCTTCCGATCT"
reverseAdapter = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"

#leave 'includeAdapter = False' to use primer sequence after adapter sequence
#confirm that primers used match those in experiment (for example, forward V19 below had degenerate sites but forward V13 does not)

forwardHash = {}
forwardHash["V19"]="AGRGTTYGATYMTGGCTCAG"
forwardHash["V13"]="AGAGTTTGATCCTGGCTCAG"
forwardHash["V34"]="CCTACGGGNGGCWGCAG"
forwardHash["V45"]="AYTGGGYDTAAAGNG"

reverseHash = {}
reverseHash["V19"]="RGYTACCTTGTTACGACTT"
reverseHash["V13"]="ATTACCGCGGCTGCTGG"
reverseHash["V34"]="GGACTACHVGGGTWTCTAAT"
reverseHash["V45"]="CCGTCAATTYHTTTRAGT"

targetHash = {}

lineCount = 0
sampleIndex = -1
targetIndex = -1

inHandle = open(descriptionFile)
lines = inHandle.readlines()
			
for line in lines:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
	
	lineInfo = line.split("\t")
	
	lineCount += 1
	
	if lineCount == 1:
		for i in range(0,len(lineInfo)):
			if lineInfo[i] == "Target":
				targetIndex = i
			elif lineInfo[i] == "sampleID":
				sampleIndex = i
	else:
		key = lineInfo[sampleIndex]
		target = lineInfo[targetIndex]
		
		targetHash[key] = target
inHandle.close()

fileResults = os.listdir(inputFolder)

command = "mkdir " + outputFolder
os.system(command)

for file in fileResults:
	result = re.search("(.*_L\d{3}_R1_001).fastq.gz$",file)
	
	if result:
		sample = result.group(1)
		if sample not in finishedSamples:
			if sample in targetHash:
				target = targetHash[sample]
				print sample + " : " + target

				read1 = inputFolder + "/" + file
				read2 = re.sub("_R1_001.fastq","_R2_001.fastq",read1)

				trim1 = outputFolder + "/" + re.sub(".gz","",file)
				trim2 = re.sub("_R1_001.fastq","_R2_001.fastq",trim1)
				
				if (target in forwardHash) and (target in reverseHash):
					Fadapter = forwardHash[target]
					Radapter = reverseHash[target]
				else:
					if target not in forwardHash:
						print "Target " + target + " not specified in forward hash"
					if target not in reverseHash:
						print "Target " + target + " not specified in reverse hash"
					sys.exit()
				
				if includeAdapter:
					print "Including Adapter in Cutadapt filter..."
					Fadapter = forwardAdapter + Fadapter
					Radapter = reverseAdapter + Radapter
				
				seqObj = Seq(Radapter)
				revcomR =  seqObj.reverse_complement()
				seqObj = Seq(Fadapter)
				revcomF =  seqObj.reverse_complement()
		
				command = "cutadapt --max-n 0 -a " + str(revcomR) + " -g " + Fadapter + " -A " + str(revcomF) + " -G " + Radapter + " -m 20 -o " + trim1 + " -p " + trim2 + " " + read1 + " " + read2
				os.system(command)
			else:
				print "Target mapping not specified for " + sample
				sys.exit()
			
