library(biomaRt)

genome = useEnsembl(biomart="ensembl")
yeastgenome = useEnsembl(biomart="ensembl",dataset="scerevisiae_gene_ensembl")


GO_0006119 <- getBM(attributes=c('ensembl_gene_id'), filters = 'go', values ="GO:0006119", mart = yeastgenome)
GO_0006406 <- getBM(attributes=c('ensembl_gene_id'), filters = 'go', values ="GO:0006406", mart = yeastgenome)
GO_0006412 <- getBM(attributes=c('ensembl_gene_id'), filters = 'go', values ="GO:0006412", mart = yeastgenome)
GO_0006950 <- getBM(attributes=c('ensembl_gene_id'), filters = 'go', values ="GO:0006950", mart = yeastgenome)
GO_0007049 <- getBM(attributes=c('ensembl_gene_id'), filters = 'go', values ="GO:0007049", mart = yeastgenome)
GO_0009651 <- getBM(attributes=c('ensembl_gene_id'), filters = 'go', values ="GO:0009651", mart = yeastgenome)
GO_0016458 <- getBM(attributes=c('ensembl_gene_id'), filters = 'go', values ="GO:0016458", mart = yeastgenome)
GO_0023052 <- getBM(attributes=c('ensembl_gene_id'), filters = 'go', values ="GO:0023052", mart = yeastgenome)
GO_0031047 <- getBM(attributes=c('ensembl_gene_id'), filters = 'go', values ="GO:0031047", mart = yeastgenome)
GO_0031990 <- getBM(attributes=c('ensembl_gene_id'), filters = 'go', values ="GO:0031990", mart = yeastgenome)
GO_0042254 <- getBM(attributes=c('ensembl_gene_id'), filters = 'go', values ="GO:0042254", mart = yeastgenome)
GO_0051028 <- getBM(attributes=c('ensembl_gene_id'), filters = 'go', values ="GO:0051028", mart = yeastgenome)
GO_0022857 <- getBM(attributes=c('ensembl_gene_id'), filters = 'go', values ="GO:0022857", mart = yeastgenome)
GO_0005840 <- getBM(attributes=c('ensembl_gene_id'), filters = 'go', values ="GO:0005840", mart = yeastgenome)


write.table(GO_0006119, "GO0006119.txt",row.names=F,quote=F,col.names=F)
write.table(GO_0006406, "GO0006406.txt", row.names=F,quote=F,col.names=F) 
write.table(GO_0006412, "GO0006412.txt",row.names=F,quote=F,col.names=F)
write.table(GO_0006950, "GO0006950.txt",row.names=F,quote=F,col.names=F)
write.table(GO_0007049, "GO0007049.txt",row.names=F,quote=F,col.names=F)
write.table(GO_0009651, "GO0009651.txt",row.names=F,quote=F,col.names=F)
write.table(GO_0016458, "GO0016458.txt",row.names=F,quote=F,col.names=F)
write.table(GO_0023052, "GO0023052.txt",row.names=F,quote=F,col.names=F)
write.table(GO_0031047, "GO0031047.txt",row.names=F,quote=F,col.names=F)
write.table(GO_0031990, "GO0031990.txt",row.names=F,quote=F,col.names=F)
write.table(GO_0042254, "GO0042254.txt",row.names=F,quote=F,col.names=F)
write.table(GO_0051028, "GO0051028.txt",row.names=F,quote=F,col.names=F)
write.table(GO_0022857, "GO0022857.txt",row.names=F,quote=F,col.names=F)
write.table(GO_0005840, "GO0005840.txt",row.names=F,quote=F,col.names=F)


