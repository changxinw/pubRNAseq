library(tximport)
library(readr)
library(DESeq2)
library(rjson)

args=commandArgs(T)
design_path = args[1] # Path of design matrix
output = args[2] # Path of output file

design_mat = read.table(design_path, sep="\t", header=F)
gsm = as.vector(design_mat[,1])
condition = as.vector(design_mat[,2])
# gsm = c('GSM2746716', 'GSM2746717', 'GSM2746718', 'GSM2746719')
# output = '/data5/home/changxin/tmp/'
files = sapply(gsm, function(i){paste0(output, i, '/quant.sf')})
# condition = c('Ctrl', 'Ctrl', 'Treat', 'Treat')
names = names(gsm)
sampleTable <- data.frame(sampleName = files, fileName = files, condition = condition)
txi_tx <- tximport(files, type="salmon", ignoreTxVersion = T, txOut=T)
dds_tx <- DESeqDataSetFromTximport(txi_tx,colData=sampleTable,design=~condition)
dds_tx <- dds_tx[rowSums(counts(dds_tx)) > 1, ]
dds_tx <- DESeq(dds_tx,fitType="local")
res_tx <- results(dds_tx, alpha = 0.01)
write.csv(res_tx, paste0(output, '/DE_table.txt'), sep='\t', row.names=T, col.names=T)