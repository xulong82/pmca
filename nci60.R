library(xlsx)
library(genefilter)
library(WGCNA)

rm(list = ls())
setwd("~/Dropbox/GitHub/MCA")

nci60 = list(rna = rna, protein = protein)
save(nci60, file = "nci60.rdt") # Saved!

# DATA
rfile <- "NCI60/RNA_Affy_HG_U133_GCRMA.txt"
pfile <- "NCI60/nci60_Protein__Lysate_Array_log2/Protein__Lysate_Array_log2.xls"

rna <- read.table(rfile, stringsAsFactors = F, header = T, sep = "\t")
rna <- rna[! grepl("^X", names(rna))]
rna$LC.NCI_H23 <- NULL
rna <- rna[! rna$Gene.name.c == "", ]
rna <- sapply(rna[10:68], function(x) tapply(x, rna$Gene.name.c, max))

protein <- read.xlsx(pfile, stringsAsFactors = F, sheetIndex = "Results", startRow = 10)
protein <- sapply(protein[10:69], function(x) tapply(x, protein$Gene.name.c, max))
protein <- protein[, colnames(rna)]

summary(rowVars(rna))
rna <- rna[rowVars(rna) > 0.3, ]
