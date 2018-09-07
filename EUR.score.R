library(foreach)
load("EUR.beta.RData")

args=(commandArgs(TRUE))

for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}

##cg <- "cg16262201"
##sfile <- "hap.markers/1.3157480.hap.markers"
##rfile <- "temp.1.raw"
##mfile <- "temp.1.methyl.list"

s<-read.table(sfile,as.is=T)
m<-read.table(mfile, as.is=T)
r<-read.table(rfile, header=T)

m <- m[m$V1 %in% s$V1,]
sm <- merge(s,m,by.x=1,by.y=1)

## get the allele associated with CpG creation
sm[,6] <- substr(sm[,3],2,2)
sm[sm[,5]==2,6] <- substr(sm[sm[,5]==2,4],2,2)

snps<-unlist(lapply(strsplit(colnames(r)[-(1:6)],split="_"),function(x){x[1]}))
##if chrpos snps<-substr(unlist(lapply(strsplit(colnames(r)[-(1:6)],split="_"),function(x){x[1]})),2,11)
alleles<-unlist(lapply(strsplit(colnames(r)[-(1:6)],split="_"),function(x){x[2]}))

cpg.score.df <- foreach(i=which( sm[,2] %in% snps), .combine=cbind) %do% {
  m.snp<-sm[i,2]
  r.allele <- alleles[snps==m.snp]
  if(r.allele == sm[i,6] | r.allele == 0) {
    score <- r[,which(snps==m.snp)+6]
  }
  else
    score <- 2 - r[,which(snps==m.snp)+6]
}

colnames(cpg.score.df) <- sm[which( sm[,2] %in% snps),2]
rownames(cpg.score.df)<-r$IID

cpg.score<-apply(cpg.score.df,1,sum)

#beta<-read.table("YRI160.27578probes.beta.GEOformat.txt" , header=T, sep="\t",row.names=1)
