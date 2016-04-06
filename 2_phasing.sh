source 0_config.sh

CHR=$1

R --vanilla << EOF

fam <- read.table("${INPUT}.fam",as.is=T)
rownames(fam) <- fam[,2]

dim(fam[fam[,3]!=0,])

children <- fam[fam[,3]!=0,2]
maternal <- fam[children,4]
paternal <- fam[children,3]

require(data.table)

sample <- read.table("${INPUT}.sample",header=T,as.is=T)
sample <- sample[-1,]

raw <- fread("${INPUT}_R_${CHR}.raw", header=T, data.table=F)
raw <-  raw[,-grep("HET",colnames(raw))]

rownames(raw) <- raw[,"IID"]

#gen <- read.table("/panfs/panasas01/sscm/xt610/TRIOS/IMPUTE/INPUT/TRIO_chr22.gen")
#colnames(gen) <- c("chr","pos","rs",paste(rep(sample[,1],each=2),c("1","2"),sep='_'))

require(foreach)

phased <-  foreach(j=1:length(children),.combine=cbind) %do% {

print(j)

out <- foreach(i=7:ncol(raw),.combine=rbind) %do% {

off <- raw[paste(children[j],"",sep=''),i]
mat <- raw[paste(paternal[j],"",sep=''),i]
pat <- raw[paste(maternal[j],"",sep=''),i]

p<-c(-9,-9,-9,-9)

if(is.na(off) | is.na(mat) | is.na(pat)){
  p<-c(NA,NA,NA,NA)
}
else {

if(off==1){
  if(pat==1 & mat==0) {
    p<-c(1,0,0,0)}
  if(pat==0 & mat==1) {
    p<-c(0,0,1,0)}
  if(pat==1 & mat==1){
    p<-c(0.5,0.5,0.5,0.5)}
  if(pat==0 & mat==0){
    p<-c(NA,NA,NA,NA) }
  
  if(pat==2 & mat==0) {
    p<-c(1,1,0,0) }
  if(pat==0 & mat==2) {
    p<-c(0,0,1,1) }
  if(pat==2 & mat==2) {
    p<-c(NA,NA,NA,NA) }
 
  if(pat==2 & mat==1) {
    p<-c(1,1,0,1) }
  if(pat==1 & mat==2) {
    p<-c(0,1,1,1)}
}

if(off==0){
  if(pat==1 & mat==0) {
    p<-c(0,1,0,0)}
  if(pat==0 & mat==1) {
    p<-c(0,0,0,1)}
  if(pat==1 & mat==1){
    p<-c(0,1,0,1)}
  if(pat==0 & mat==0){
    p<-c(0,0,0,0) }
  
  if(pat+mat>2) {
    p<-c(NA,NA,NA,NA) }
}

if(off==2){
  if(pat==1 & mat==2) {
    p<-c(1,0,1,1)}
  if(pat==2 & mat==1) {
    p<-c(1,1,1,0)}
  if(pat==1 & mat==1){
    p<-c(1,0,1,0)}
  if(pat==2 & mat==2){
    p<-c(2,2,2,2) }
  
  if(pat+mat<2) {
    p<-c(NA,NA,NA,NA) }
}

}

p
}
colnames(out) <- c(paste(rep(children[j],4),c("F_A","F_B","M_A","M_B"),sep=''))
out
}

save(phased,file="${INPUT}_${CHR}.phased.RData")

EOF
