# Analysis of antagonistic protein pairs in cancer
# based on TCPA
# antagonistic pairs; low-high and high-low exp in different cancer types
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

# define differentially expressed proteins 
# based on the expression value threshold (valthr)
valthr = 2
prots = colnames(tu)[4:ncol(tu)]
sdowns = list()
sups = list()
for (i in 1:nrow(tu)){
    val = as.numeric(tu[i,4:ncol(tu)])
    sups[[i]] = prots[which(val>valthr)]
    sdowns[[i]] = prots[which(val<(-valthr))]}

# save mutually antagonist pairs
antps = c()
for (i in 1:(nrow(tu)-1)){
    for (j in (i+1):nrow(tu)){
        p1 = intersect(sups[[i]],sdowns[[j]])
        p2 = intersect(sups[[j]],sdowns[[i]])
        if (length(p1)>0 & length(p2)>0){
            ct1 = as.character(tu$Cancer_Type[i])
            ct2 = as.character(tu$Cancer_Type[j])
            for (p1i in p1){
                for (p2i in p2){
                    antps = rbind(antps,c(p1i,p2i,ct1,ct2)) }}}}}
colnames(antps)=c('P1','P2','C1','C2')
antps = as.data.frame(antps)

# based on cancer type (avoiding same cancer pairs)
antps2 = c()
for (i in 1:nrow(antps)){
    if (antps[i,]$C1 != antps[i,]$C2){
        antps2 = rbind(antps2,antps[i,])}}
colnames(antps2)=c('P1','P2','C1','C2')
antps2 = as.data.frame(antps2)

# based on unique cancer type (avoiding repetitions and same cancer pairs)
areps = cbind(as.character(antps2$P1),as.character(antps2$P2),
              as.character(antps2$C2),as.character(antps2$C1))
colnames(areps)=colnames(antps2)
antps3 = unique(rbind(antps2,areps))
areps2 = cbind(as.character(antps3$P2),as.character(antps3$P1),
              as.character(antps3$C1),as.character(antps3$C2))
colnames(areps2)=colnames(antps3)
antps3 = unique(rbind(antps3,areps2))

# save(antps,file = 'antps.RData')
# save(antps3,file = 'antps3.RData')
# load(file='antps.RData')
# load(file='antps3.RData')

# function to get appearance counts for antagonist pairs
antcounts <- function(dfant, divideforreps = FALSE){
    proteins = unique(c(as.character(dfant$P1),as.character(dfant$P2)))
    pantcts = c()
    for (ii in 1:(length(proteins)-1)){
        for (jj in (ii+1):length(proteins)){
            pi = proteins[ii]
            pj = proteins[jj]
            paifw = dfant[dfant$P1==pi & dfant$P2==pj,]
            pairv = dfant[dfant$P1==pj & dfant$P2==pi,]
            paict = nrow(paifw)+nrow(pairv)
            if (paict>0){
                pantcts = rbind(pantcts,c(pi,pj,paict)) }}}
    colnames(pantcts)=c('P1','P2','count')
    pantcts = as.data.frame(pantcts)
    if (divideforreps == TRUE){
        div2 = as.numeric(as.character(pantcts$count))/4
        pantcts = cbind(as.character(pantcts$P1),as.character(pantcts$P2),
                        as.character(div2))
        colnames(pantcts)=c('P1','P2','count')
        pantcts = as.data.frame(pantcts)}
    return(pantcts)}

# count mutually antagonist pairs based on samples or cancer types
ctantps = antcounts(antps)
ctantps3 = antcounts(antps3,divideforreps = TRUE)

# randomizations for significance testing
# shuffle each row seperately
# count cases with equal or more scores
rno = 1000
rct1 = c()
rct3 = c()
for (rr in 1:(length(prots)-1)){
    for (qq in (rr+1):length(prots)){
        rct1 = rbind(rct1,c(prots[rr],prots[qq]))
        rct3 = rbind(rct3,c(prots[rr],prots[qq])) }}
colnames(rct1)=c('P1','P2')
colnames(rct3)=c('P1','P2')
rct1 = as.data.frame(rct1)
rct3 = as.data.frame(rct3)

for (ri in 1:rno){
    print(ri)
    rtu = c()
    for (tui in 1:nrow(tu)){
        rtu = rbind(rtu,sample(as.numeric(tu[tui,4:ncol(tu)])))}
    rsdowns = list()
    rsups = list()
    for (i in 1:nrow(rtu)){
        val = rtu[i,]
        rsups[[i]] = prots[which(val>valthr)]
        rsdowns[[i]] = prots[which(val<(-valthr))]}
    rantps = c()
    for (i in 1:(nrow(rtu)-1)){
        for (j in (i+1):nrow(rtu)){
            p1 = intersect(rsups[[i]],rsdowns[[j]])
            p2 = intersect(rsups[[j]],rsdowns[[i]])
            if (length(p1)>0 & length(p2)>0){
                ct1 = as.character(tu$Cancer_Type[i])
                ct2 = as.character(tu$Cancer_Type[j])
                for (p1i in p1){
                    for (p2i in p2){
                        rantps = rbind(rantps,c(p1i,p2i,ct1,ct2)) }}}}}

    colnames(rantps)=c('P1','P2','C1','C2')
    rantps = as.data.frame(rantps)
    # based on samples (cancer type repetition allowed)
    rctantps = antcounts(rantps)
    # based on cancer type (avoiding same cancer pairs)
    rantps2 = c()
    for (i in 1:nrow(rantps)){
        if (rantps[i,]$C1 != rantps[i,]$C2){
            rantps2 = rbind(rantps2,rantps[i,])}}
    colnames(rantps2)=c('P1','P2','C1','C2')
    rantps2 = as.data.frame(rantps2)
    # based on unique cancer type (avoiding repetitions and same cancer pairs)
    rareps = cbind(as.character(rantps2$P1),as.character(rantps2$P2),
                  as.character(rantps2$C2),as.character(rantps2$C1))
    colnames(rareps) = colnames(rantps2)
    rantps3 = unique(rbind(rantps2,rareps))
    rareps2 = cbind(as.character(rantps3$P2),as.character(rantps3$P1),
                   as.character(rantps3$C1),as.character(rantps3$C2))
    colnames(rareps2)=colnames(rantps3)
    rantps3 = unique(rbind(rantps3,rareps2))
    rctantps3 = antcounts(rantps3,divideforreps = TRUE)
    # adding counts
    coadd1 = c()
    coadd3 = c()
    for (rr in 1:(length(prots)-1)){
        for (qq in (rr+1):length(prots)){
            prr = prots[rr]
            pqq = prots[qq]
            if (nrow(rctantps[rctantps$P1==prr & rctantps$P2==pqq,]) > 0){
                countval1 = as.numeric(as.character(rctantps[rctantps$P1==prr & rctantps$P2==pqq,]$count))}
            else if (nrow(rctantps[rctantps$P2==prr & rctantps$P1==pqq,]) > 0){
                countval1 = as.numeric(as.character(rctantps[rctantps$P2==prr & rctantps$P1==pqq,]$count))}
            else {countval1 = 0}
            if (nrow(rctantps3[rctantps3$P1==prr & rctantps3$P2==pqq,]) > 0){
                countval3 = as.numeric(as.character(rctantps3[rctantps3$P1==prr & rctantps3$P2==pqq,]$count))}
            else if (nrow(rctantps3[rctantps3$P2==prr & rctantps3$P1==pqq,]) > 0){
                countval3 = as.numeric(as.character(rctantps3[rctantps3$P2==prr & rctantps3$P1==pqq,]$count))}
            else {countval3 = 0}
            coadd1 = c(coadd1,countval1)
            coadd3 = c(coadd3,countval3)}}
    rct1 = cbind(rct1,coadd1)
    rct3 = cbind(rct3,coadd3)}

# save(rct1,file = 'rct1_1000.RData')
# save(rct3,file = 'rct3_1000.RData')
# load(file='rct1_1000.RData')
# load(file='rct3_1000.RData')

# calculate p values, and FDR corrected p-values
uncorp1 = c()
uncorp3 = c()
indws1 = c()
indws3 = c()
indw = 0
for (rr in 1:(length(prots)-1)){
    for (qq in (rr+1):length(prots)){
        indw = indw + 1
        prr = prots[rr]
        pqq = prots[qq]
        
        if (nrow(ctantps[ctantps$P1==prr & ctantps$P2==pqq,]) > 0) {
            indws1 = c(indws1,indw)
            co1 = as.numeric(as.character(ctantps[ctantps$P1==prr & ctantps$P2==pqq,]$count))}
        else if (nrow(ctantps[ctantps$P2==prr & ctantps$P1==pqq,]) > 0){
            indws1 = c(indws1,indw)
            co1 = as.numeric(as.character(ctantps[ctantps$P2==prr & ctantps$P1==pqq,]$count))}
        else { co1 = 0 }
        
        x1 = rct1[rct1$P1==prr & rct1$P2==pqq,]
        x11 = as.numeric(as.character(x1[3:length(x1)]))
        x111 = length(which(x11 >= co1))/rno
        uncorp1 = c(uncorp1,x111)
        
        if (nrow(ctantps3[ctantps3$P1==prr & ctantps3$P2==pqq,]) > 0) {
            indws3 = c(indws3,indw)
            co3 = as.numeric(as.character(ctantps3[ctantps3$P1==prr & ctantps3$P2==pqq,]$count)) }
        else if (nrow(ctantps3[ctantps3$P2==prr & ctantps3$P1==pqq,]) > 0){
            indws3 = c(indws3,indw)
            co3 = as.numeric(as.character(ctantps3[ctantps3$P2==prr & ctantps3$P1==pqq,]$count))}
        else { co3 = 0 }
        
        x3 = rct3[rct3$P1==prr & rct3$P2==pqq,]
        x33 = as.numeric(as.character(x3[3:length(x3)]))
        x333 = length(which(x33 >= co3))/rno
        uncorp3 = c(uncorp3,x333) }}

uncorp11 = uncorp1[indws1]
uncorp33 = uncorp3[indws3]

# add real values
uncorp11 = uncorp11 + 0.001
uncorp33 = uncorp33 + 0.001

# FDR p-values
corp1 = p.adjust(uncorp11,method = 'fdr')
corp3 = p.adjust(uncorp33,method = 'fdr')

# add p-values to the antagonist protein pairs data frame
ctantps = cbind(ctantps,uncorp11,corp1)
ctantps3 = cbind(ctantps3,uncorp33,corp3)
colnames(ctantps)[4:5] = c('p-value','fdr')
colnames(ctantps3)[4:5] = c('p-value','fdr')

# save(ctantps,file = 'ctantps.RData')
# save(ctantps3,file = 'ctantps3.RData')
# load(file='ctantps.RData')
# load(file='ctantps3.RData')

# select significant pairs based on the p-value threshold
pvalthr = 0.05
sctantps = ctantps[ctantps$fdr<=pvalthr,]
sctantps3 = ctantps3[ctantps3$fdr<=pvalthr,]

# compare ctantps and ctantps3 as well as sctantps and sctantps3
# make a table showing all pairs belonging to any one of these groups 
# so that they can all be viewed
ctantpsmatch = c()
sctantpsmatch = c()
ctantpsonly1 = c()
ctantpsonly3 = c()
sctantpsonly1 = c()
sctantpsonly3 = c()
allctable = c()
for (rr in 1:(length(prots)-1)){
    for (qq in (rr+1):length(prots)){
        ppmatch1 = 0
        ppmatch3 = 0
        sppmatch1 = 0
        sppmatch3 = 0
        if (nrow(ctantps[ctantps$P1==prots[rr] & ctantps$P2==prots[qq],])>0){
            lx = ctantps[ctantps$P1==prots[rr] & ctantps$P2==prots[qq],3:5]
            ctvallist = c(as.numeric(as.character(lx[[1]])),as.numeric(as.character(lx[[2]])),as.numeric(as.character(lx[[3]])))
            ppmatch1 = 1}
        if (nrow(ctantps[ctantps$P2==prots[rr] & ctantps$P1==prots[qq],])>0){
            lx = ctantps[ctantps$P2==prots[rr] & ctantps$P1==prots[qq],3:5]
            ctvallist = c(as.numeric(as.character(lx[[1]])),as.numeric(as.character(lx[[2]])),as.numeric(as.character(lx[[3]])))
            ppmatch1 = 1}
        if (nrow(ctantps3[ctantps3$P1==prots[rr] & ctantps3$P2==prots[qq],])>0){
            lx = ctantps3[ctantps3$P1==prots[rr] & ctantps3$P2==prots[qq],3:5]
            ctvallist3 = c(as.numeric(as.character(lx[[1]])),as.numeric(as.character(lx[[2]])),as.numeric(as.character(lx[[3]])))
            ppmatch3 = 1}
        if (nrow(ctantps3[ctantps3$P2==prots[rr] & ctantps3$P1==prots[qq],])>0){
            lx = ctantps3[ctantps3$P2==prots[rr] & ctantps3$P1==prots[qq],3:5]
            ctvallist3 = c(as.numeric(as.character(lx[[1]])),as.numeric(as.character(lx[[2]])),as.numeric(as.character(lx[[3]])))
            ppmatch3 = 1}
        
        if (nrow(sctantps[sctantps$P1==prots[rr] & sctantps$P2==prots[qq],])>0){
            sppmatch1 = 1}
        if (nrow(sctantps[sctantps$P2==prots[rr] & sctantps$P1==prots[qq],])>0){
            sppmatch1 = 1}
        if (nrow(sctantps3[sctantps3$P1==prots[rr] & sctantps3$P2==prots[qq],])>0){
            sppmatch3 = 1}
        if (nrow(sctantps3[sctantps3$P2==prots[rr] & sctantps3$P1==prots[qq],])>0){
            sppmatch3 = 1}
        
        if (ppmatch1 == 1 & ppmatch3 == 1){
            allctable = rbind(allctable,c(prots[rr],prots[qq],ctvallist,ctvallist3))
            ctantpsmatch = rbind(ctantpsmatch,c(prots[rr],prots[qq]))}
        if (ppmatch1 == 1 & ppmatch3 == 0){
            allctable = rbind(allctable,c(prots[rr],prots[qq],ctvallist,0,1.0,1.0))
            ctantpsonly1 = rbind(ctantpsonly1,c(prots[rr],prots[qq]))}
        if (ppmatch1 == 0 & ppmatch3 == 1){
            allctable = rbind(allctable,c(prots[rr],prots[qq],0,1.0,1.0,ctvallist3))
            ctantpsonly3 = rbind(ctantpsonly3,c(prots[rr],prots[qq]))}
        
        if (sppmatch1 == 1 & sppmatch3 == 1){
            sctantpsmatch = rbind(sctantpsmatch,c(prots[rr],prots[qq]))}
        if (sppmatch1 == 1 & sppmatch3 == 0){
            sctantpsonly1 = rbind(sctantpsonly1,c(prots[rr],prots[qq]))}
        if (sppmatch1 == 0 & sppmatch3 == 1){
            sctantpsonly3 = rbind(sctantpsonly3,c(prots[rr],prots[qq]))}}}
allctable = as.data.frame(allctable)
colnames(allctable) = c('P1','P2','CS count','CS p-value','CS FDR','CT count','CT p-value','CT FDR')
write.table(allctable,file='TCPA_allMAPPtable.txt',
            row.names=FALSE,sep='\t',quote=FALSE)

# compare antagonist pairs of proteins to their coexpression
# compare number of dif.exp., weight ratios for ctantps (CS)
ctprots = unique(c(as.character(ctantps$P1),as.character(ctantps$P2)))
mixvalspre = c()
for (i in 1:length(prots)){
    si = prots[i]
    tusi = tu[,(3+which(prots==si))]
    wchup = which(tusi>valthr)
    wchdown = which(tusi<(-valthr))
    difexpno = length(wchup)+length(wchdown)
    
    if (difexpno != 0){
        psrt = 0.0
        ngrt = 0.0
        for (sj in prots){
            if (si != sj){
                tusj = tu[,(3+which(prots==sj))]
                wchup2 = which(tusj>valthr)
                wchdown2 = which(tusj<(-valthr))
                posint = length(intersect(wchup,wchup2))+length(intersect(wchdown,wchdown2))
                posuni = length(union(wchup,wchup2))+length(union(wchdown,wchdown2))
                posrat = posint/posuni
                negint = length(intersect(wchup,wchdown2))+length(intersect(wchdown,wchup2))
                neguni = length(union(wchup,wchdown2))+length(union(wchdown,wchup2))
                negrat = negint/neguni
                psrt = psrt + posrat
                ngrt = ngrt + negrat}}
        
        if (si %in% ctprots){
            selct1 = ctantps[ctantps$P1==si,]
            selct2 = ctantps[ctantps$P2==si,]
            selct = rbind(selct1,selct2)
            sumsic = sum(as.numeric(selct$count))}
        else {sumsic = 0.0}
        mixvalspre = rbind(mixvalspre,c(difexpno,psrt,ngrt,sumsic))}
        else {mixvalspre = rbind(mixvalspre,c(0.0,0.0,0.0,0.0))}}

# compare antagonist pairs of proteins to their coexpression
# compare number of dif.exp., weight ratios ctantps3 (CT)
# presence of dif exp in at least one sample of a cancer type is required
# ctprots3 = unique(c(as.character(ctantps3$P1),as.character(ctantps3$P2)))
# mixvalspre = c()
# for (i in 1:length(prots)){
#     si = prots[i]
#     tusi = tu[,(3+which(prots==si))]
#     wchup = unique(tu$Cancer_Type[which(tusi>valthr)])
#     wchdown = unique(tu$Cancer_Type[which(tusi<(-valthr))])
#     difexpno = length(wchup)+length(wchdown)
#     
#     if (difexpno != 0){
#         psrt = 0.0
#         ngrt = 0.0
#         for (sj in prots){
#             if (si != sj){
#                 tusj = tu[,(3+which(prots==sj))]
#                 wchup2 = unique(tu$Cancer_Type[which(tusj>valthr)])
#                 wchdown2 = unique(tu$Cancer_Type[which(tusj<(-valthr))])
#                 posint = length(intersect(wchup,wchup2))+length(intersect(wchdown,wchdown2))
#                 posuni = length(union(wchup,wchup2))+length(union(wchdown,wchdown2))
#                 posrat = posint/posuni
#                 negint = length(intersect(wchup,wchdown2))+length(intersect(wchdown,wchup2))
#                 neguni = length(union(wchup,wchdown2))+length(union(wchdown,wchup2))
#                 negrat = negint/neguni
#                 psrt = psrt + posrat
#                 ngrt = ngrt + negrat}}
#         
#         if (si %in% ctprots3){
#             selct1 = ctantps3[ctantps3$P1==si,]
#             selct2 = ctantps3[ctantps3$P2==si,]
#             selct = rbind(selct1,selct2)
#             sumsic = sum(as.numeric(selct$count))}
#         else {sumsic = 0.0}
#         mixvalspre = rbind(mixvalspre,c(difexpno,psrt,ngrt,sumsic))}
#     else {mixvalspre = rbind(mixvalspre,c(0.0,0.0,0.0,0.0))}}

#ranking weights for either ctantps (CS) or ctantps3 (CT) based values
mixvals = c()
for (mi in 1:length(prots)){
    newrts = c(prots[mi],as.character(mixvalspre[mi,1]/sum(mixvalspre[,1])),
               as.character(mixvalspre[mi,2]/sum(mixvalspre[,2])),
               as.character(mixvalspre[mi,3]/sum(mixvalspre[,3])),
               as.character(mixvalspre[mi,4]/sum(mixvalspre[,4])))
    mixvals = rbind(mixvals,newrts)}
mixvals = cbind(mixvals[,1],mixvalspre,mixvals[,2:5])
mixvals = as.data.frame(mixvals)
colnames(mixvals) = c('protein','dif_exp_no','pos_cor','neg_cor','ant_pairs_no','dif_exp_rat','pos_cor_rat','neg_cor_rat','ant_pairs_rat')
# compare differences 
dif_poscor = as.numeric(as.character(mixvals$pos_cor_rat))-as.numeric(as.character(mixvals$dif_exp_rat))
dif_negcor = as.numeric(as.character(mixvals$neg_cor_rat))-as.numeric(as.character(mixvals$dif_exp_rat))
dif_antpairs = as.numeric(as.character(mixvals$ant_pairs_rat))-as.numeric(as.character(mixvals$dif_exp_rat))
alldif = cbind(dif_poscor,dif_negcor,dif_antpairs)
mixvals = cbind(mixvals,alldif)
write.table(mixvals,file='TCPA_mixvals_CS.txt',
            row.names=FALSE,sep='\t',quote=FALSE)
# write.table(mixvals,file='TCPA_mixvals_CT.txt',
#             row.names=FALSE,sep='\t',quote=FALSE)

# plot differential ratio values
svg(file='TCPA_weightedratio_mixvals_CS.svg',width=4.5, height=4)
alldif = alldif[order(alldif[,3]),]
# svg(file='TCPA_weightedratio_mixvals_CT.svg',width=4.5, height=4.6)
limy = max(c(abs(min(alldif)),max(alldif)))
plot(alldif[,3],type='l',ylim=c((-limy-0.002),(limy+0.02)),lwd = 1.5,
     main='',xlab='',ylab='',xaxt='n',yaxt='n')
lines(alldif[,1],col='darkorchid1',lwd = 1.5)
lines(alldif[,2],col='darkorange1',lwd = 1.5)
legend(x = 0.1,y = 0.09,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('MAPP','PCPP','NCPP'),col = c('black','darkorchid1','darkorange1'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='CS Based Analysis of TCPA',ylab ="Differential Ratio",
      xlab="Protein Index",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# import PPI and TF interactions datasets; human interactions are selected
orjppints = read.delim('BIOGRID-MV-Physical-4.3.194.tab3.txt',colClasses = 'character')
hupints = orjppints[orjppints$Organism.Name.Interactor.A == 'Homo sapiens' & 
                        orjppints$Organism.Name.Interactor.B == 'Homo sapiens',]
hupints = hupints[,c(8,9)]
# two TF interaction datasets are combined
orjtfints1 = read.delim('HTRIdb_data.txt',colClasses = 'character')
orjtfints2 = read.delim('trrust_rawdata.human.tsv',header=FALSE,colClasses = 'character')
colnames(orjtfints2)=c('TF','GENE','TYPE','PUBMED_ID')
tfints1 = orjtfints1[,c(3,5)]
tfints2 = orjtfints2[,c(1,2)]
colnames(tfints1) = c('TF','GENE')
tfints = rbind(tfints1,tfints2)
hutfints = unique(tfints)
namesmat = read.delim('TCPA_NCBI_Symbols.txt',colClasses = 'character')

# CS based pairs with NCBI Symbols
prpairsw = c()
for (pi in 1:nrow(sctantps)){
    p1 = as.character(sctantps[pi,1])
    p2 = as.character(sctantps[pi,2])
    np11 = namesmat[namesmat$TCPA_Symbol==p1,]$NCBI_Symbol1
    np12 = namesmat[namesmat$TCPA_Symbol==p1,]$NCBI_Symbol2
    np13 = namesmat[namesmat$TCPA_Symbol==p1,]$NCBI_Symbol3
    np21 = namesmat[namesmat$TCPA_Symbol==p2,]$NCBI_Symbol1
    np22 = namesmat[namesmat$TCPA_Symbol==p2,]$NCBI_Symbol2
    np23 = namesmat[namesmat$TCPA_Symbol==p2,]$NCBI_Symbol3
    np1s = c(np11,np12,np13)
    np2s = c(np21,np22,np23)
    for (ni1 in np1s){
        for (ni2 in np2s){
            if (ni1 != '-' & ni2 != '-'){
                if (ni1 != ni2){
                    prpairsw = rbind(prpairsw,c(ni1,ni2))}}}}}
prpairsw = rbind(prpairsw,cbind(prpairsw[,2],prpairsw[,1]))
colnames(prpairsw) = colnames(sctantps)[1:2]
prpairsw = unique(as.data.frame(prpairsw))
prpairsw = prpairsw[1:(nrow(prpairsw)/2),]

# CT based pairs with NCBI Symbols
prpairsw3 = c()
for (pi in 1:nrow(sctantps3)){
    p1 = as.character(sctantps3[pi,1])
    p2 = as.character(sctantps3[pi,2])
    np11 = namesmat[namesmat$TCPA_Symbol==p1,]$NCBI_Symbol1
    np12 = namesmat[namesmat$TCPA_Symbol==p1,]$NCBI_Symbol2
    np13 = namesmat[namesmat$TCPA_Symbol==p1,]$NCBI_Symbol3
    np21 = namesmat[namesmat$TCPA_Symbol==p2,]$NCBI_Symbol1
    np22 = namesmat[namesmat$TCPA_Symbol==p2,]$NCBI_Symbol2
    np23 = namesmat[namesmat$TCPA_Symbol==p2,]$NCBI_Symbol3
    np1s = c(np11,np12,np13)
    np2s = c(np21,np22,np23)
    for (ni1 in np1s){
        for (ni2 in np2s){
            if (ni1 != '-' & ni2 != '-'){
                if (ni1 != ni2){
                    prpairsw3 = rbind(prpairsw3,c(ni1,ni2))}}}}}
prpairsw3 = rbind(prpairsw3,cbind(prpairsw3[,2],prpairsw3[,1]))
colnames(prpairsw3) = colnames(sctantps3)[1:2]
prpairsw3 = unique(as.data.frame(prpairsw3))
prpairsw3 = prpairsw3[1:(nrow(prpairsw3)/2),]

# check the protein-protein interaction neighborhood of antagonist pairs
# check for direct interactions (first degree neighbors)
neprpairsw = c()
for (pi in 1:nrow(prpairsw)){
    pr1 = as.character(prpairsw[pi,1])
    pr2 = as.character(prpairsw[pi,2])
    if (pr1 %in% hupints[hupints$Official.Symbol.Interactor.A==pr2,]$Official.Symbol.Interactor.B){
        neprpairsw = rbind(neprpairsw,prpairsw[pi,])}}
neprpairsw = unique(neprpairsw)

neprpairsw3 = c()
for (pi in 1:nrow(prpairsw3)){
    pr1 = as.character(prpairsw3[pi,1])
    pr2 = as.character(prpairsw3[pi,2])
    if (pr1 %in% hupints[hupints$Official.Symbol.Interactor.A==pr2,]$Official.Symbol.Interactor.B){
        neprpairsw3 = rbind(neprpairsw3,prpairsw3[pi,])}}
neprpairsw3 = unique(neprpairsw3)

# check for indirect interactions (second degree neighbors)
ne2prpairsw = c()
for (pi in 1:nrow(prpairsw)){
    pr1 = as.character(prpairsw[pi,1])
    pr2 = as.character(prpairsw[pi,2])
    p2ne = hupints[hupints$Official.Symbol.Interactor.A==pr2,]$Official.Symbol.Interactor.B
    for (pp in p2ne){
        pr3 = as.character(pp)
        if (pr1 %in% hupints[hupints$Official.Symbol.Interactor.A==pr3,]$Official.Symbol.Interactor.B){
            ne2prpairsw = rbind(ne2prpairsw,prpairsw[pi,])}}}
ne2prpairsw = unique(ne2prpairsw)

ne2prpairsw3 = c()
for (pi in 1:nrow(prpairsw3)){
    pr1 = as.character(prpairsw3[pi,1])
    pr2 = as.character(prpairsw3[pi,2])
    p2ne = hupints[hupints$Official.Symbol.Interactor.A==pr2,]$Official.Symbol.Interactor.B
    for (pp in p2ne){
        pr3 = as.character(pp)
        if (pr1 %in% hupints[hupints$Official.Symbol.Interactor.A==pr3,]$Official.Symbol.Interactor.B){
            ne2prpairsw3 = rbind(ne2prpairsw3,prpairsw3[pi,])}}}
ne2prpairsw3 = unique(ne2prpairsw3)

# check the neighborhood of antagonist pairs in the transcriptional network
tfprpairsw = c()
for (pi in 1:nrow(prpairsw)){
    pr1 = as.character(prpairsw[pi,1])
    pr2 = as.character(prpairsw[pi,2])
    if (pr1 %in% hutfints[hutfints$TF==pr2,]$GENE){
        tfprpairsw = rbind(tfprpairsw,c(pr2,pr1))}
    if (pr2 %in% hutfints[hutfints$TF==pr1,]$GENE){
        tfprpairsw = rbind(tfprpairsw,c(pr1,pr2))}}
tfprpairsw = unique(as.data.frame(tfprpairsw))
colnames(tfprpairsw)=c('P1','P2')

tfprpairsw3 = c()
for (pi in 1:nrow(prpairsw3)){
    pr1 = as.character(prpairsw3[pi,1])
    pr2 = as.character(prpairsw3[pi,2])
    if (pr1 %in% hutfints[hutfints$TF==pr2,]$GENE){
        tfprpairsw3 = rbind(tfprpairsw3,c(pr2,pr1))}
    if (pr2 %in% hutfints[hutfints$TF==pr1,]$GENE){
        tfprpairsw3 = rbind(tfprpairsw3,c(pr1,pr2))}}
tfprpairsw3 = unique(as.data.frame(tfprpairsw3))
colnames(tfprpairsw3)=c('P1','P2')

# combine all pairs in NCBI symbols and make a table
intvals = c()
prr = unique(c(as.character(prpairsw$P1),as.character(prpairsw$P2),
               as.character(prpairsw3$P1),as.character(prpairsw3$P2)))
for (pii in 1:(length(prr)-1)){
    for (pjj in (pii+1):length(prr)){
        ntgs = c()
        tgs1 = c()
        pi = prr[pii]
        pj = prr[pjj]
        if (pi %in% prpairsw[prpairsw$P1==pj,]$P2 | 
            pj %in% prpairsw[prpairsw$P1==pi,]$P2){
            tgs1 = c(tgs1,1)}
        else {tgs1 = c(tgs1,0)}
        if (pi %in% prpairsw3[prpairsw3$P1==pj,]$P2 | 
            pj %in% prpairsw3[prpairsw3$P1==pi,]$P2){
            tgs1 = c(tgs1,1)}
        else {tgs1 = c(tgs1,0)}
        # ignore pairs that are not present in either of CS or CT pairs
        # continue if the pair is present in at least one of CS or CT pairs
        if (sum(tgs1)>0){
            tgs2 = c()
            if (pi %in% neprpairsw[neprpairsw$P1==pj,]$P2 | 
                pj %in% neprpairsw[neprpairsw$P1==pi,]$P2){
                tgs2 = c(tgs2,1)}
            else {tgs2 = c(tgs2,0)}
            if (pi %in% neprpairsw3[neprpairsw3$P1==pj,]$P2 | 
                pj %in% neprpairsw3[neprpairsw3$P1==pi,]$P2){
                tgs2 = c(tgs2,1)}
            else {tgs2 = c(tgs2,0)}
            tgs3 = c()
            if (pi %in% ne2prpairsw[ne2prpairsw$P1==pj,]$P2 | 
                pj %in% ne2prpairsw[ne2prpairsw$P1==pi,]$P2){
                tgs3 = c(tgs3,1)}
            else {tgs3 = c(tgs3,0)}
            if (pi %in% ne2prpairsw3[ne2prpairsw3$P1==pj,]$P2 | 
                pj %in% ne2prpairsw3[ne2prpairsw3$P1==pi,]$P2){
                tgs3 = c(tgs3,1)}
            else {tgs3 = c(tgs3,0)}
            tgs4 = c()
            if (pi %in% tfprpairsw[tfprpairsw$P1==pj,]$P2 | 
                pj %in% tfprpairsw[tfprpairsw$P1==pi,]$P2){
                tgs4 = c(tgs4,1)}
            else {tgs4 = c(tgs4,0)}
            if (pi %in% tfprpairsw3[tfprpairsw3$P1==pj,]$P2 | 
                pj %in% tfprpairsw3[tfprpairsw3$P1==pi,]$P2){
                tgs4 = c(tgs4,1)}
            else {tgs4 = c(tgs4,0)}
            # print tags based on presence in pairs, ne, ne2, tfints
            # presence in CS, CT pairs
            if (sum(tgs1)>0){
                if (tgs1[1] == 1){
                    if (tgs1[2] == 1){
                        ntgs = c(ntgs,'Both')}
                    else {ntgs = c(ntgs,'CS')}}
                else {ntgs = c(ntgs,'CT')}}
            else {ntgs = c(ntgs,'None')}
            # presence in ne (direct PPI)
            if (sum(tgs2)>0){
                if (tgs2[1] == 1){
                    if (tgs2[2] == 1){
                        ntgs = c(ntgs,'Both')}
                    else {ntgs = c(ntgs,'CS')}}
                else {ntgs = c(ntgs,'CT')}}
            else {ntgs = c(ntgs,'None')}
            # presence in ne2 (indirect PPI)
            if (sum(tgs3)>0){
                if (tgs3[1] == 1){
                    if (tgs3[2] == 1){
                        ntgs = c(ntgs,'Both')}
                    else {ntgs = c(ntgs,'CS')}}
                else {ntgs = c(ntgs,'CT')}}
            else {ntgs = c(ntgs,'None')}
            # presence in tfints
            if (sum(tgs4)>0){
                if (tgs4[1] == 1){
                    if (tgs4[2] == 1){
                        ntgs = c(ntgs,'Both')}
                    else {ntgs = c(ntgs,'CS')}}
                else {ntgs = c(ntgs,'CT')}}
            else {ntgs = c(ntgs,'None')}
            intvals = rbind(intvals,c(pi,pj,ntgs))}}}
intvals = as.data.frame(intvals)
colnames(intvals) = c('P1','P2','MAPP_type','Direct_PPI_type','Indirect_PPI_type','TFint_type')
write.table(intvals,file='TCPA_intvals.txt',
            row.names=FALSE,sep='\t',quote=FALSE)            
        
# plot degree distributions
# for combined network
pis = c(as.character(intvals$P1),as.character(intvals$P2))
pit = table(pis)
pid = as.data.frame(pit)
svg(file='TCPA_combinedCSCT_MAPPnet_degdist.svg',width=4, height=4)
hist(pid$Freq,breaks=20,ylim=c(0,20),main='',xlab='',ylab='',xaxt='n',yaxt='n')
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='Combined TCPA MAPP Network',ylab ="Frequency",xlab="Degree",
      font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# for CS network
selintvals = intvals[intvals$MAPP_type=='CS' | intvals$MAPP_type=='Both',]
pis = c(as.character(selintvals$P1),as.character(selintvals$P2))
pit = table(pis)
pid = as.data.frame(pit)
svg(file='TCPA_CS_MAPPnet_degdist.svg',width=4, height=4.3)
hist(pid$Freq,breaks=20,ylim=c(0,25),main='',xlab='',ylab='',xaxt='n',yaxt='n')
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='TCPA CS Based MAPP Network',ylab ="Frequency",xlab="Degree",
      font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# for CT network
selintvals = intvals[intvals$MAPP_type=='CT' | intvals$MAPP_type=='Both',]
pis = c(as.character(selintvals$P1),as.character(selintvals$P2))
pit = table(pis)
pid = as.data.frame(pit)
svg(file='TCPA_CT_MAPPnet_degdist.svg',width=4, height=4.3)
hist(pid$Freq,breaks=15,ylim=c(0,25),xlim=c(0,20),main='',xlab='',ylab='',xaxt='n',yaxt='n')
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='TCPA CT Based MAPP Network',ylab ="Frequency",xlab="Degree",
      font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()
        
# enlist the cancer types of the significant MAPP
# get more frequently appearing ones
scct = c()
for (ri in 1:nrow(sctantps3)){
  rip1 = as.character(sctantps3[ri,1])
  rip2 = as.character(sctantps3[ri,2])
  a1 = antps3[antps3$P1==rip1 & antps3$P2==rip2, ]
  a2 = antps3[antps3$P1==rip2 & antps3$P2==rip1, ]
  aa = rbind(a1,a2)
  scct = rbind(scct,aa[,3:4])}

# get unique pairs
scct2 = unique(scct)
unscct = scct2[1,]
for (ri in 1:nrow(scct2)){
  ric1 = as.character(scct2[ri,1])
  ric2 = as.character(scct2[ri,2])
  a1 = unscct[unscct$C1==ric1 & unscct$C2==ric2, ]
  a2 = unscct[unscct$C1==ric2 & unscct$C2==ric1, ]
  aa = rbind(a1,a2)
  if (nrow(aa) == 0){
    unscct = rbind(unscct,c(ric1,ric2))}}
write.table(unscct,file='TCPA_uniqueCTpairs_sctMAPP.txt',
            row.names=FALSE,sep='\t',quote=FALSE)

# get appearance ratio and compare to dataset ratio
tucl = tu$Cancer_Type
tuct = table(tucl)/length(tucl)
scs = unique(c(as.character(scct$C1),as.character(scct$C2)))
ctcount = c()
for (ci in scs){
  c1c = scct[scct$C1==ci,]
  c2c = scct[scct$C2==ci,]
  rc = (nrow(c1c)+nrow(c2c))/nrow(scct)
  dc = tuct[ci]
  ctcount = rbind(ctcount,c(ci,rc,dc))}
ctcount = as.data.frame(ctcount)
colnames(ctcount)=c('CT','MAPP_ratio','Dataset_ratio')
ctcount = ctcount[order(ctcount$MAPP_ratio),]
write.table(ctcount,file='TCPA_CT_freqs.txt',
            row.names=FALSE,sep='\t',quote=FALSE)

# plot histogram cancer type frequencies
allrts = as.numeric(as.character(ctcount$ratio))
svg(file='TCPA_CT_freqs.svg',width=4.5, height=4)
hist(allrts,breaks=15,ylim=c(0,5),xlim=c(0,0.26),
     main='',xlab='',ylab='',xaxt='n',yaxt='n')
text(0.21, 1.5, labels = 'LGG', cex = 1.2, font=2)
text(0.16, 2.0, labels = 'SARC', cex = 1.2, font=2)
text(0.14, 1.4, labels = 'CESC', cex = 1.2, font=2)
text(0.125, 2.5, labels = 'SKCM', cex = 1.2, font=2)
text(0.125, 3.0, labels = 'BLCA', cex = 1.2, font=2)
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='Cancer Frequencies of TCPA',ylab ="Frequency",xlab="Fraction of Count",
      font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()