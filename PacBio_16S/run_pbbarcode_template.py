import os

### 		Option #1: Label pre-exising CCS reads (with desired minimum number of cycles)			  ###
### A (at least partially) successful RS_ReadsOfInsert result can be sufficient to create input files ###

#easiest option would be to create on your own with full paths to bc.h5 files (not links, which might be used by SMRT Portal)
#otherwise, code to label on your own is under option #2
barcodeFofn = "barcode.fofn"

outputdir = "pbbarcode"
inputFofn = "ccs.fofn"
command = "pbbarcode emitFastqs --outDir " + outputdir + " " + inputFofn + " " + barcodeFofn
os.system(command)

### Option #2: Can potentially skip any use of SMRT Portal (but you still need access to ConsensusTools)	  ###

#tested with sequences from https://github.com/PacificBiosciences/Bioinformatics-Training/blob/master/barcoding/pacbio_384_barcodes.fasta
#if you have access, you can find barcodes with the standard names in on the SMRT Portal server
barcodeFa = "pacbio_384_barcodes.fasta"

inputFofn = "baxh5.fofn"
threads = 4
outputdir = "bax_pbbarcode_10barcode"
barcodeFofn = outputdir + "/barcode.fofn"
command = "pbbarcode labelZmws --outDir " + outputdir + " --outFofn " + barcodeFofn + " --nProc " + str(threads)+ " " + barcodeFa + " " + inputFofn
os.system(command)

#this step might need to be run separately
inputFofn = "server_baxh5.fofn"
minCycles = 5
command = "sh /opt/smrtanalysis/current/analysis/bin/ConsensusTools.sh CircularConsensus -o " + outputdir + " -n "+str(threads)+" --minFullPasses=" + str(minCycles) + " --fofn=" + inputFofn
os.system(command)

#you can use minimum number of barcodes as a proxy for minimum cycles, if ccs .h5 files with default were already available #
inputFofn = "new_ccs.fofn"
minBarcode = 10
command = "pbbarcode emitFastqs --outDir " + outputdir + " --minNumBarcodes " + str(minBarcode) + " " + inputFofn + " " + barcodeFofn
os.system(command)