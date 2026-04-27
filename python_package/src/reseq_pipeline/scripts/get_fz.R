#Rscript - group.xls fz.txt

args <- commandArgs(trailingOnly=TRUE)

group_te = read.table(args[1], header=TRUE)
a = t(group_te)
colnames(a) = group_te$species
group_names = a[-1,]
group_names <- gsub("\\+", "_", group_names)
unique_groups <- sort(unique(group_names))
comparisons <- combn(unique_groups, 2, FUN = function(x) paste(x[2], x[1], sep = "_vs_") )
#print(matrix(comparisons,ncol=1))
write.table(matrix(comparisons,ncol=1), file = args[2], row.names = FALSE, quote = FALSE, col.names = FALSE, sep = "\t")
