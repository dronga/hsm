beta <- read.table("GSE39672_matrix_processed.txt",header=T,sep="\t")

rownames(beta) <- beta[,1]
samples <- colnames(beta)[seq(2,ncol(beta),by=2)]

beta <- beta[,samples]
beta.all <- beta

EUR <- read.table("EUR.sample.txt",as.is=T)$V1
beta <- beta.all[,EUR]
save(beta,file="EUR.beta.RData")

AFR <- read.table("AFR.sample.txt", as.is=T)$V1
beta <- beta.all[,AFR]
save(beta,file="AFR.beta.RData")


write.table(colnames(beta)[seq(2,ncol(beta),by=2)], file="samples.txt", row.names=F, col.names=F, quote=F)

load("anno.RData")

write.table(anno[,c("IlmnID", "CHR", "MAPINFO")],file="450k.chrpos.txt",row.names=F, col.names=F, quote=F)
