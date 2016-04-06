load(paste("../INPUT/hapmap3_r3_b36_fwd.consensus.qc.poly.CEU_",22,".phased.RData",sep=""))
cnames <- colnames(phased)

for(i in 1:10) {
  load(paste("../INPUT/hapmap3_r3_b36_fwd.consensus.qc.poly.CEU_",i,".phased.RData",sep=""))
  colnames(phased) <- cnames
  save(phased, file=paste("../INPUT/hapmap3_r3_b36_fwd.consensus.qc.poly.CEU_",i,".phased.RData",sep=""))
}
