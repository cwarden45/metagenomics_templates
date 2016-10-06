import os

rdp = "/mnt/user_data/Seq/cwarden/rdp_classifier_2.11/dist/classifier.jar"

fullFA = "trainset15_092015.fa"
filterFA1 = "trainset15_092015.rm.dupseq.fa"
command = "java -Xmx2g -jar " + rdp + " rm-dupseq -d -l 1200 -i " + fullFA + " -o " + filterFA1
os.system(command)

filterFA2 = "trainset15_092015.rm.partialseq.fa"
command = "java -Xmx8g -jar " + rdp + " rm-partialseq " + filterFA1+ " " + filterFA1 + " " + filterFA2
os.system(command)

command = "/opt/bwa/bwa index " filterFA2
os.system(command)