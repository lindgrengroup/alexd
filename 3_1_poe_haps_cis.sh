#!/bin/bash
#$ -cwd -V
#$ -j y
#$ -P lindgren.prja -q short.qa
#$ -t 1-843
#$ -N HM_POE_A
echo "************************************************************"
echo "SGE job ID: "$JOB_ID
echo "SGE task ID: "$SGE_TASK_ID
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo "************************************************************"






module load R

cd /well/lindgren/alexd/miRNA/HAPMAP_ALSPAC_METHOD/SCRIPTS

source 0_config.sh

R --vanilla <<EOF

#j=$SGE_TASK_ID

j=$1

require(data.table)
require(foreach)

load("/well/lindgren/alexd/miRNA/HAPMAP_ALSPAC_METHOD/INPUT/hapmap3_r3_b36_fwd.consensus.qc.poly.CEU.phased.combined.RData")

load("/well/lindgren/alexd/miRNA/expression.RData")
load("/well/lindgren/alexd/miRNA/annotation.RData")

load("${INPUT}.phased.combined.RData")
load("$EXP")
load("$ANNO")

haps <- phased.combined
expression <- mirna

#sample <-  fread("/well/lindgren/alexd/miRNA/HAPMAP_ALSPAC_METHOD/INPUT/hapmap3_r3_b36_fwd.consensus.qc.poly.CEU.fam", data.table=F)

sample <-  fread("${INPUT}.fam", data.table=F)

 bim <- foreach(k=1:22,.combine=rbind) %do%{
   print(k)
   a <- fread(paste("/well/lindgren/alexd/miRNA/HAPMAP_ALSPAC_METHOD/INPUT/hapmap3_r3_b36_fwd.consensus.qc.poly.CEU_R_",k,".bim",sep=""),data.table=F)
   a
 }

dim(bim)

##save(bim,file="/well/lindgren/alexd/miRNA/HAPMAP_ALSPAC_METHOD/INPUT/hapmap3_r3_b36_fwd.consensus.qc.poly.CEU.phased.combined.anno.RData")

##load("/well/lindgren/alexd/miRNA/HAPMAP_ALSPAC_METHOD/INPUT/hapmap3_r3_b36_fwd.consensus.qc.poly.CEU.phased.combined.anno.RData")

#sample <- sample[-1,]

##colnames(haps)[6:ncol(haps)] <- paste(rep(sample[,2],each=2),c("MT","FT"),sep="_")


children <- sample[sample[,4]!=0,2]

children <- unique(substr(colnames(haps),1,7))

#colnames(expression) <-  link[colnames(expression),1]
colnames(haps)
children
sample

chd.exp <- expression[,children]
pat <- haps[,paste(children,"F_A",sep="")]
mat <- haps[,paste(children,"M_A",sep="")]

require(foreach)



dim(chd.exp)
dim(mat)

i=anno[j,"illuminaID"]

summary(lm(unlist(chd.exp[1,]) ~ unlist(mat[1,])))

haps.chr <- as.numeric(bim[,1])
haps.pos <- as.numeric(bim[,4])

chr <- anno[j,"chromo"]
pos <- anno[j,"probe_coords"]

index <- which(haps.chr==chr & haps.pos > (pos - 500000) & (haps.pos < pos + 500000))

pat <- haps[index,paste(children,"F_A",sep="")]
mat <- haps[index,paste(children,"M_A",sep="")]

probes <- rownames(chd.exp)
snps<- bim[index,2]

mat.lm <-   foreach(j=1:nrow(mat),.combine=rbind) %do% {
 if(j%%1000==0) {print(j)}
    out <- try(coefficients(summary(lm(unlist(chd.exp[i,]) ~ unlist(mat[j,]))))[2,])
   if(!inherits(out,"try-error")){
     out
  } else {
    rep(NA,4)
   }
}


save(mat.lm,snps,probes,file=paste("${OUTPUT}_mat_",j,".B.RData",sep=""))

pat.lm <- foreach(j=1:nrow(pat),.combine=rbind) %do% {
 if(j%%100==0) {print(j)}
    out <-try(coefficients(summary(lm(unlist(chd.exp[i,]) ~ unlist(pat[j,]))))[2,])
  if(!inherits(out,"try-error")){
     out
  } else {
    rep(NA,4)
   }

}

save(pat.lm,snps,probes,file=paste("${OUTPUT}_pat_",j,".B.RData",sep=""))


add.lm <-   foreach(j=1:nrow(pat),.combine=rbind) %do% {
 if(j%%100==0) {print(j)}
    out<-try(coefficients(summary(lm(unlist(chd.exp[i,]) ~ (unlist(pat[j,]) + unlist(mat[j,])))))[2,])

  if(!inherits(out,"try-error")){
     out
  } else {
    rep(NA,4)
   }

}

save(add.lm,snps,probes,file=paste("${OUTPUT}_add_",j,".B.RData",sep=""))

EOF
