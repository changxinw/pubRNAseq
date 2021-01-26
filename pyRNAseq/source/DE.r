required_Packages = c("tximport", "readr", "DESeq2",
                      "rjson", "dplyr")
if(!all(required_Packages %in% installed.packages())){
	source("https://bioconductor.org/biocLite.R")
	biocLite(setdiff(required_Packages, installed.packages()))
}

require(tximport)
require(readr)
require(DESeq2)
require(rjson)
require(dplyr)

args=commandArgs(T)
design_path = args[1] # Path of design matrix
output = args[2] # Path of output file
ref = args[3] # Directory of reference file

design_mat <- read.table(design_path, sep="\t", header=F)
gsm <- as.vector(design_mat[,1])
condition <- as.vector(design_mat[,2])
# gsm = c('GSM2746716', 'GSM2746717', 'GSM2746718', 'GSM2746719')
# output = '/data5/home/changxin/tmp/'
files <- sapply(gsm, function(i){paste0(output, '/', i, '/quant.sf')})
# condition = c('Ctrl', 'Ctrl', 'Treat', 'Treat')
names <- names(gsm)
sampleTable <- data.frame(sampleName = files, fileName = files, condition = condition)
txi_tx <- tximport(files, type="salmon", ,ignoreTxVersion = TRUE, txOut = TRUE, dropInfReps = TRUE)
dds_tx <- DESeqDataSetFromTximport(txi_tx,colData=sampleTable,design=~condition)
dds_tx <- dds_tx[rowSums(counts(dds_tx)) > 1, ]
dds_tx <- DESeq(dds_tx)
res_tx <- results(dds_tx, alpha = 0.05)
#de_new = na.omit(res_tx)
res_tx$padj[is.na(res_tx$padj)] <- 1
de_new <- res_tx
raw_refseq <- read.table(ref, sep = "\t")
refseq <- raw_refseq[!duplicated(raw_refseq[,1]),]
rownames(refseq) <- refseq[,1]
refseq <- refseq[,-1]

co_index <- intersect(rownames(de_new), rownames(refseq))
result <- cbind(de_new[co_index,], refseq[co_index,"V6"])
colnames(result)[7]  <- "symbol"
result <- result[order(result[,6]),]
result2 <- distinct(as.data.frame(result), symbol, .keep_all = TRUE)
rownames(result2) <- result2$symbol

write.table(result2[,1:6], paste0(output, '/DE_table.txt'), sep='\t', quote=F)
