
source 0_config.sh

R --vanilla <<EOF

require(foreach)

phased.combined <- foreach(i=1:22,.combine=rbind) %do% {
  print(i)
  file <- paste("${INPUT}_",i,".phased.RData",sep="")
  load(file)
  phased  
}

save(phased.combined, file="${INPUT}.phased.combined.RData")

EOF