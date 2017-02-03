import os
import re
import sys

refFa = "/path/to/RDP/trainset15.rm.partialseq.fa"
refTax = "/path/to/RDP/trainset15_db_taxid.txt"
sintaxFa = "RDP_trainset15.fa"

#taxonomy hash
taxHash = {}

inHandle = open(refTax)
line = inHandle.readline()

while line:
	line = re.sub("\r","",line)
	line = re.sub("\n","",line)
	
	lineInfo = line.split("*")
	
	id = lineInfo[1]
	level = lineInfo[4]
	taxHash[id]=level
	
	line = inHandle.readline()

inHandle.close()

#parse fasta
outHandle = open(sintaxFa,"w")

inHandle = open(refFa)
line = inHandle.readline()

while line:
	line = re.sub("\r","",line)
	line = re.sub("\n","",line)
	
	headerResult = re.search("^>",line)
	
	if headerResult:
		lineInfo = line.split("\t")
		acc = re.search("^(\S*)\|",lineInfo[0])
		acc = acc.group(1)
		text = acc + ";tax="
		
		taxText = lineInfo[1]
		taxInfo = taxText.split(";")
		
		for tax in taxInfo:
			if tax in taxHash:
				level = taxHash[tax]
				
				if level == "domain":
					text = text + "d:" + tax + ","
				elif level == "phylum":
					text = text + "p:" + tax + ","
					#orderTest = re.search("d:",text)
					#if orderTest:
					#	text = text + "p:" + tax + ","
					#else:
					#	print " phylum issue with " + line
					#	break
				elif level == "class":
					text = text + "c:" + tax + ","
					#orderTest = re.search("p:",text)	
					#if orderTest:
					#	text = text + "c:" + tax + ","
					#else:
					#	print " class issue with " + line
					#	break
				elif level == "order":
					text = text + "o:" + tax + ","
					#orderTest = re.search("c:",text)
					#if orderTest:
					#	text = text + "o:" + tax + ","
					#else:
					#	print " order issue with " + line
					#	break
				elif level == "family":
					text = text + "f:" + tax + ","
					#orderTest = re.search("o:",text)
					#if orderTest:
					#	text = text + "f:" + tax + ","
					#else:
					#	print " family issue with " + line
					#	break
				elif level == "genus":
					text = text + "g:" + tax + ","
					#orderTest = re.search("f:",text)	
					#if orderTest:
					#	text = text + "g:" + tax + ","
					#else:
					#	print " genus issue with " + line
					#	break
				elif (level != "rootrank") & (level != "subclass")& (level != "suborder"):
					print "Need to map " + level
					sys.exit()
					
		text = re.sub(",$",";",text)
		text = text + "\n"
		outHandle.write(text)
	else:
		text = line + "\n"
		outHandle.write(text)
	
	line = inHandle.readline()

inHandle.close()