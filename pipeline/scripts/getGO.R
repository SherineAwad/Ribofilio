library(biomaRt)

args <- commandArgs(trailingOnly = TRUE)

genome = useEnsembl(biomart="ensembl")
yeastgenome = useEnsembl(biomart="ensembl",dataset="scerevisiae_gene_ensembl",host = "aug2020.archive.ensembl.org")

for (id in args)
{

output = paste(id,".txt", sep="") 

GO_biomart <- getBM(attributes=c('ensembl_gene_id'), filters = 'go', values =id, mart = yeastgenome)
write.table(GO_biomart, output,row.names=F,quote=F,col.names=F)

}

