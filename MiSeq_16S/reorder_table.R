input.file = "abundance.txt"
output.file = "abundance_ordered.txt"

#input.file = "read_counts.txt"
#output.file = "read_counts_ordered.txt"

meta.file = "sample_description.txt"

meta.table = read.table(meta.file, head=T, sep="\t")


input.table = read.table(input.file, head=T, sep="\t")
Assignment = input.table$Assignment
value.mat = input.table[,match(meta.table$userID,names(input.table))]

output.table = data.frame(Assignment, value.mat)
write.table(output.table, output.file, quote=F, sep="\t", row.names=F)