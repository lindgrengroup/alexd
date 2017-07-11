####PLOT ALL BETA
pdf("beta.rep.diff.all.pdf", width=20)

par(mfrow=c(1,6))
for(i in 1:6) {

pre<-top6.beta[i,rep1] 
pre.case<-top6.beta[i,rep1 & phe$grp2=="case"]
pre.ctrl<-top6.beta[i,rep1 & phe$grp2=="ctrl"]

post<-top6.beta[i,rep2]
post.case<-top6.beta[i,rep2 & phe$grp2=="case"]
post.ctrl<-top6.beta[i,rep2 & phe$grp2=="ctrl"]




s<-seq(length(pre))
par(bty="l")
boxplot(pre,post,main=paste(top6[i], "(",top6.names[i],")"),ylab="Beta",names=c("baseline","follow-up"))
stripchart(list(pre,post),vertical=T,pch=16,cex=0.5,add=T)
segments(rep(1,length(pre.case))[s],pre.case[s],rep(2,length(pre.case))[s],post.case[s],col="red",lwd=0.5)
segments(rep(1,length(pre.ctrl))[s],pre.ctrl[s],rep(2,length(pre.ctrl))[s],post.ctrl[s],col="blue",lwd=0.5)

}
dev.off()

####PLOT ALL BETA VERSION 2
pdf("beta.rep.diff.all.2.pdf", width=14)

par(mfrow=c(1,6))
for(i in 1:6) {

pre<-top6.beta[i,rep1] 
pre.case<-top6.beta[i,rep1 & phe$grp2=="case"]
pre.ctrl<-top6.beta[i,rep1 & phe$grp2=="ctrl"]

post<-top6.beta[i,rep2]
post.case<-top6.beta[i,rep2 & phe$grp2=="case"]
post.ctrl<-top6.beta[i,rep2 & phe$grp2=="ctrl"]



 
res.case<-wilcox.test(post.case,pre.case,paired=T,conf.int=T)
res.ctrl<-wilcox.test(post.ctrl,pre.ctrl,paired=T,conf.int=T)
 
stripchart(list(post.case-pre.case,post.ctrl-pre.ctrl),vertical=T,pch=16,method="jitter",main=paste(top6[i], "(",top6.names[i],")"),ylab="Follow-up - Baseline",xlab="Median+/-95%CI", col=c("red","blue"))


points(1,res.case$estimate,col="red",pch=16,cex=2)
arrows(1,res.case$conf.int[1],1,res.case$conf.int[2],col="red",code=3,lwd=3,angle=90)

points(2,res.ctrl$estimate,col="blue",pch=16,cex=2)
arrows(2,res.ctrl$conf.int[1],2,res.ctrl$conf.int[2],col="blue",code=3,lwd=3,angle=90)
abline(h=0,lty=2)#Zero-effectline

}
dev.off()

####PLOT ALL BETA VERSION 3
pdf("beta.rep.diff.all.3.pdf", width=21)

par(mfrow=c(1,6))
for(i in 1:6) {
pre<-top6.beta[i,rep1] 
pre.case<-top6.beta[i,rep1 & phe$grp2=="case"]
pre.ctrl<-top6.beta[i,rep1 & phe$grp2=="ctrl"]

post<-top6.beta[i,rep2]
post.case<-top6.beta[i,rep2 & phe$grp2=="case"]
post.ctrl<-top6.beta[i,rep2 & phe$grp2=="ctrl"]

s<-seq(length(pre))
par(bty="l")
boxplot(pre.case,pre.ctrl,post.case,post.ctrl,xlim = c(0, 8), at=c(1,2,7,8),main=paste(top6[i], "(",top6.names[i],")"),ylab="Beta",names=c("baseline","","follow-up", ""), col=c("palevioletred","lightblue", "palevioletred","lightblue"))
stripchart(list(pre,post),vertical=T,pch=16,cex=0.5,add=T)
segments(rep(1,length(pre.case))[s],pre.case[s],rep(7,length(pre.case))[s],post.case[s],col="red",lwd=0.5)
segments(rep(2,length(pre.ctrl))[s],pre.ctrl[s],rep(8,length(pre.ctrl))[s],post.ctrl[s],col="blue",lwd=0.5)

}
dev.off()

#### PLOT ALL + DIFFERENCES

pdf("beta.rep.diff.all.1.pdf", width=21)
#PLOT
pre<-top6.beta[i,rep1] 
pre.case<-top6.beta[i,rep1 & phe$grp2=="case"]
pre.ctrl<-top6.beta[i,rep1 & phe$grp2=="ctrl"]

post<-top6.beta[i,rep2]
post.case<-top6.beta[i,rep2 & phe$grp2=="case"]
post.ctrl<-top6.beta[i,rep2 & phe$grp2=="ctrl"]



#Settinguptwoscreens
par(mfrow=c(1,2))
 
#FirstGraph
s<-seq(length(pre))
par(bty="l")
boxplot(pre,post,main=paste("Res at",top6[i], "(",top6.names[i],")"),xlab=paste("pre/post p =",signif(p.time[i,3],3),", ","case/control p =",signif(p.time[i,2],3),"\n","interaction p =",signif(p.time[i,1],3)),ylab="Res",names=c("pre","post"),col=c("lightblue","lightgreen"))
stripchart(list(pre,post),vertical=T,pch=16,method="jitter",cex=0.5,add=T)
segments(rep(0.95,length(pre.case))[s],pre.case[s],rep(2,length(pre.case))[s],post.case[s],col="red",lwd=0.5)
segments(rep(0.95,length(pre.ctrl))[s],pre.ctrl[s],rep(2,length(pre.ctrl))[s],post.ctrl[s],col="blue",lwd=0.5)


#Secondgraph
#Confidenceintervals eitherparametric (t.test) or non-parametric (wilcox.text)
#res<-t.test(post,prä,paired=T,conf.int=T)
res.case<-wilcox.test(post.case,pre.case,paired=T,conf.int=T)
res.ctrl<-wilcox.test(post.ctrl,pre.ctrl,paired=T,conf.int=T)
 
stripchart(post-pre,vertical=T,pch=16,method="jitter",main="Difference in Res",ylab="Difference:Post–Pre",xlab="Median+/-95%CI", , col="white")
stripchart(post.case-pre.case,vertical=T,pch=16,method="jitter",main="Difference in Beta",ylab="Difference:Post–Pre",xlab="Median+/-95%CI", , col="red",add=T)
stripchart(post.ctrl-pre.ctrl,vertical=T,pch=16,method="jitter",main="Difference in Beta",ylab="Difference:Post–Pre",xlab="Median+/-95%CI", add=T, col="blue")

points(1,res.case$estimate,col="red",pch=16,cex=2)
arrows(1,res.case$conf.int[1],1,res.case$conf.int[2],col="red",code=3,lwd=3,angle=90)

points(1,res.ctrl$estimate,col="blue",pch=16,cex=2)
arrows(1,res.ctrl$conf.int[1],1,res.ctrl$conf.int[2],col="blue",code=3,lwd=3,angle=90)
abline(h=0,lty=2)#Zero-effectline

dev.off()
