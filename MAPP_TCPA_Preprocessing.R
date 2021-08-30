# Analysis of antagonistic protein pairs in cancer
# antagonistic pairs; low-high and high-low exp in different cancer types
# preprocessing of the TCPA dataset analysis
# by Ertugrul Dalgic, PhD
# 2020-2021

# import TCPA protein exp dataset
# change missing values (NA) to 0
# count 0s before NA assignment to 0 and check percentage of NA
# remove proteins with more than 50% of NA values 
tu = read.delim('TCGA-PANCAN32-L4.csv',sep=',')
nofp1 = c()
for (i in 4:ncol(tu)){
    nofp1 = c(nofp1,length(which(tu[,i]==0)))}
tu[is.na(tu)] = 0
nofp2 = c()
for (i in 4:ncol(tu)){
    nodif = length(which(tu[,i]==0)) - nofp1[i-3]
    nofp2 = c(nofp2,(nodif/nrow(tu)))}
tu = tu[,c(1,2,3,(which(nofp2<=0.5)+3))]

# plot distribution of all values
tuvm = tu[,4:ncol(tu)]
tx = c()
for (i in 1:ncol(tuvm)){tx=c(tx,tuvm[,i])}
svg(file='TCPA_distribution.svg',width=5, height=5)
par(mar=c(4,4,2,2))
hist(tx,breaks = 100,ylim=c(0,350000),main='',xlab='',ylab='',yaxt='n')
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='',ylab='Frequency',
      xlab='Expression Value',font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# plot the number of dif. exp. proteints w.r.t. varying thresholds
mxv = round(max(c(max(tuvm,na.rm=TRUE),abs(min(tuvm,na.rm=TRUE)))))
nofd = c()
for (ithr in (0:(mxv*2))/2){
    nofd = c(nofd,length(tuvm[tuvm>ithr])+length(tuvm[tuvm<(-ithr)]))}
nofd = nofd/(nrow(tuvm)*ncol(tuvm))
svg(file='TCPA_thresholds.svg',width=5,height=4)
par(mar=c(4,4,2,2))
plot(((0:20)/2),nofd[1:21],pch=16,main='',xlab='',ylab='',yaxt='n',xaxt='n')
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='',ylab='Fraction of DE Values',
      xlab='Expression Threshold',font.lab=2,cex.lab=1.2,cex.main=1.3)
lines(((0:20)/2),nofd[1:21],lwd=2)
dev.off()