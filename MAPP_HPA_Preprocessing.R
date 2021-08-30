# Analysis of antagonistic protein pairs in cancer
# antagonistic pairs; low-high and high-low exp in different cancer types
# preprocessing of the HPA dataset analysis
# by Ertugrul Dalgic, PhD
# 2020-2021

# import HPA protein exp dataset for normal and cancer
# import normal-cancer cell and tissue match files
orjtu = read.delim('pathology.tsv')
nrg = read.delim('normal_tissue.tsv')
catisg = read.delim('CancerNormalTissueMatch.txt',colClasses = 'character')
cacesg = read.delim('CancerNormalCellTypeMatch.txt',colClasses = 'character')

# remove NA rows from pathology dataset (tu) and remove unused columns
# some proteins can have very few scoring values which is removed
tug = na.omit(orjtu[,1:7])
sus = apply(tug[,4:7],1,'sum')

# check how size changes w.r.t. sthr
tsizes = c()
gsizes = c()
for (ss in 1:12){
    tugss = tug[which(sus>=ss),]
    tsizes = c(tsizes,nrow(tugss))
    gsizes = c(gsizes,length(unique(tugss$Gene.name)))}
tsizes = tsizes/nrow(tug)
gsizes = gsizes/length(unique(tug$Gene.name))
svg(file='HPA_proteins_expnumbers.svg',width=5, height=4)
par(mar=c(4,4,2,2))
plot(1:12,tsizes,type='l',lwd=2,col='slateblue4',main='',xlab='',ylab='',xaxt='n',yaxt='n')
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='',ylab='Fraction of Remaining Size',
      xlab='Expression Number Threshold',font.lab=2,cex.lab=1.2,cex.main=1.3)
lines(1:12,gsizes,col='coral4',lwd=2)
legend(x = 1,y = 0.6,bty = 'n',lwd = 2,cex=1.2,text.font=2,
       legend=c('All Values','Unique Proteins'),col = c('slateblue4','coral4'))
dev.off()

# remove rows with less than 10 values
sthr = 10
tug = tug[which(sus>=sthr),]

# functions to get appearance counts for antagonist pairs
antcas <- function(tu,nr,wd,dfrathr,catis,caces){
    
    # Cancer Expression parameter
    # define absent or present proteins based on CE thr (dfrathr)
    dfra = ((wd*tu$High+tu$Medium)-(tu$Low+wd*tu$Not.detected))/(wd*tu$High+tu$Medium+tu$Low+wd*tu$Not.detected)
    th = tu[which(dfra >= dfrathr),]
    tl = tu[which(dfra <= -dfrathr),]
    
    # select differentially expressed proteins
    # for each normal tissue cell type select proteins with not detected (absent) or high (present) values
    # for each cancer type select pairs which are absent-present and present-absent in cancer vs normal
    antc = list()
    calst = unique(catis$Cancer)
    for (casample in calst){
        # tumor
        thigh = th[th$Cancer==casample,]
        tno = tl[tl$Cancer==casample,]
        # normal
        nrtiss = catis[catis$Cancer==casample,2:dim(catis)[2]]
        nrcell = caces[caces$Cancer==casample,2:dim(caces)[2]]
        # for combined lists for each cancer
        hnl = c()
        nhl = c()
        for (nrti in nrtiss){
            if (nrti != 'none'){
                nrr = nr[nr$Tissue==nrti,]
                for (ci in nrcell){
                    if (ci != 'none'){
                        ni = nr[nr$Tissue==nrti & nr$Cell.type==ci,]
                        nno = ni[ni$Level=='Not detected',]
                        nhigh = ni[ni$Level=='High',]
                        highno = intersect(thigh$Gene.name,nno$Gene.name)
                        nohigh = intersect(tno$Gene.name,nhigh$Gene.name)
                        hnl = c(hnl,highno)
                        nhl = c(nhl,nohigh) }}}}
        antc[[as.character(casample)]] = c(list(unique(hnl)),list(unique(nhl))) }
    return(antc)}

# find pairs which are ON-OFF in one cancer type and OFF-ON in another cancer type
antcaprcts <- function(antcgiv,catis){
    calst = unique(catis$Cancer)
    aprs = c()
    for (ci in 1:(length(calst)-1)){
        for (cj in (ci+1):length(calst)){
            ca1 = as.character(calst[ci])
            ca2 = as.character(calst[cj])
            cl1 = intersect(antcgiv[[ca1]][[1]],antcgiv[[ca2]][[2]])
            cl2 = intersect(antcgiv[[ca2]][[1]],antcgiv[[ca1]][[2]])
            if (length(cl1)>0 & length(cl2)>0){
                for (cii in cl1){
                    for (cjj in cl2){
                        aprs = rbind(aprs,c(cii,cjj,ca1,ca2))}}}}}
    return(aprs)}    

# analyze how number of pairs change w.r.t. wd (w) and dfrathr (DE thr)
wdrows = c()
for (di in 1:10){
    dii = di/10
    for (wi in 1:5){
       dfra = ((wi*tug$High+tug$Medium)-(tug$Low+wi*tug$Not.detected))/(wi*tug$High+tug$Medium+tug$Low+wi*tug$Not.detected)
       th = tug[which(dfra >= dii),]
       tl = tug[which(dfra <= -dii),]
       thantc = antcas(tug,nrg,wd=wi,dfrathr=dii,catisg,cacesg)
       aap = antcaprcts(thantc,catisg)
       totdr = (nrow(th)+nrow(tl))/nrow(tug)
       wdrows = rbind(wdrows,c(wi,dii,totdr,nrow(aap)))}}
par(mar=c(4,4,2,2))
svg(file='HPAdifexp_CE.svg',width=6, height=5)
boxplot(wdrows[1:5,3],wdrows[6:10,3],wdrows[11:15,3],wdrows[16:20,3],
        wdrows[21:25,3],wdrows[26:30,3],wdrows[31:35,3],wdrows[36:40,3],
        wdrows[41:45,3],wdrows[46:50,3],ylim=c(0.4,1),
        main='',xlab='',ylab='',yaxt='n')
axis(side=1,1:10,seq(0.1,1,0.1),font=2,cex.axis=1.2)
axis(side=2,font=2,cex.axis=1.2)
title(main ='',ylab='Fraction of Remaining Rows',
      xlab='CE Threshold',font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()
svg(file='HPAMAPP_CE.svg',width=6, height=6)
boxplot(wdrows[1:5,4],wdrows[6:10,4],wdrows[11:15,4],wdrows[16:20,4],
        wdrows[21:25,4],wdrows[26:30,4],wdrows[31:35,4],wdrows[36:40,4],
        wdrows[41:45,4],wdrows[46:50,4],
        main='',xlab='',ylab='',yaxt='n')
axis(side=1,1:10,seq(0.1,1,0.1),font=2,cex.axis=1.2)
axis(side=2,font=2,cex.axis=1.2)
title(main ='',ylab='Number of MAPP',
      xlab='CE Threshold',font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()
svg(file='HPAdifexp_w.svg',width=5, height=5)
boxplot(wdrows[seq(1,50,5),3],wdrows[seq(2,50,5),3],wdrows[seq(3,50,5),3],
        wdrows[seq(4,50,5),3],wdrows[seq(5,50,5),3],xaxt='n',ylim=c(0.4,1),
        main='',xlab='',ylab='',yaxt='n')
axis(side=1,1:5,font=2,cex.axis=1.2)
axis(side=2,font=2,cex.axis=1.2)
title(main ='',ylab='Fraction of Remaining Rows',
      xlab='Weight',font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()
svg(file='HPAMAPP_w.svg',width=5, height=6)
boxplot(wdrows[seq(1,50,5),4],wdrows[seq(2,50,5),4],wdrows[seq(3,50,5),4],
        wdrows[seq(4,50,5),4],wdrows[seq(5,50,5),4],xaxt='n',
        main='',xlab='',ylab='',yaxt='n')
axis(side=1,1:5,font=2,cex.axis=1.2)
axis(side=2,font=2,cex.axis=1.2)
title(main ='',ylab='Number of MAPP',
      xlab='Weight',font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()