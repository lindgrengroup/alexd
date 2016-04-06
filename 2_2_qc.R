load("/well/lindgren/alexd/miRNA/HAPMAP_ALSPAC_METHOD/INPUT/hapmap3_r3_b36_fwd.consensus.qc.poly.CEU.phased.combined.RData")

pass <- apply(phased.combined,1,function(x){!all(is.na(x)) & length(unique(na.omit(x))) > 1})

phased.combined <- phased.combined[pass,]

save(phased.combined,file="/well/lindgren/alexd/miRNA/HAPMAP_ALSPAC_METHOD/INPUT/hapmap3_r3_b36_fwd.consensus.qc.poly.CEU.phased.combined.qc.RData")
