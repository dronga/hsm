p<-read.table("pvalues.120514")
head(p)


p$ad.j<-p.adjust(p[,4])
head(p)
sum(p[,5]<0.05,na.rm=T)

write.table(p[which(p[,5]<0.05),1],file="signif.120514",quote=F,row.names=F,col.names=F)

##find the process number corresponding to cg IDs (to get meth.haps, hap pos, etc
anno.num <- read.table("450k.auto.chrpos.txt")
signif.num <- which(anno.num[,1] %in% ((p[which(p[,5]<0.05),1])))
write.table(paste("^", signif.num, "\\.", sep=""),file="signif.num.120514",quote=F,row.names=F,col.names=F)

source("~/scripts/qq.R")

png("qq.120514.png")
qq.ci(p[,4])
dev.off()

##
signif.num

pos <- read.table("signif.120514.hap.pos")

pos.cg <- cbind(anno.num[pos[,1],1],pos)
rownames(pos.cg) <- pos.cg[,1]
pos.cg$len <- pos.cg[,4]-pos.cg[,3]


head(pos.cg)

load("EUR.beta.RData")

cpg <- read.table("EUR.cpgscores.120514",sep=" ",row.names=1)

head(cpg)

pos.cg <- pos.cg[rownames(cpg),]

require(foreach)

is.numeric(cpg[1,])


##comapre to meQTL
meqtl <- read.csv("moen.meqtl.csv")
sum(rownames(cpg) %in% meqtl[,1])

pdf("plots.120514.pdf")
eff <- foreach(i=1:nrow(cpg),.combine=c) %do%  {
  coef <- coefficients(lm(as.numeric(beta[rownames(cpg)[i],])~as.numeric(cpg[i,])))
  plot(as.numeric(cpg[i,]),as.numeric(beta[rownames(cpg)[i],]),main=rownames(cpg)[i],xlab=paste("intact CpG sites\nstart=",pos.cg[i,4],"hap len=",pos.cg[i,5], "meQTL=",rownames(cpg)[i] %in%  meqtl[,1] ), ylab="Methylation Beta")
  abline(coef)
  coef[2]
}
dev.off()
