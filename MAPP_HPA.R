# Analysis of antagonistic protein pairs in cancer
# based on HPA
# antagonistic pairs; low-high and high-low exp in different cancer types
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

# get real data pairs based on the decided values of w and DE thr (dfraset)
wdset = 2
dfraset = 0.5
thantc = antcas(tug,nrg,wd=wdset,dfrathr=dfraset,catisg,cacesg)
aap = antcaprcts(thantc,catisg)
aap = as.data.frame(aap)
colnames(aap)=c('P1','P2','C1','C2')

# randomize each gene (shuffle each row seperately)
rno = 1000
rpall = c()
susg = apply(tug[,4:7],1,'sum')
for (ri in 1:rno){
    print(ri)
    rtug = c()

    for (si in 1:length(susg)){
        vv = 0
        xv = susg[si]
        tv = c()
        for (jj in 1:3){
            xv = xv - vv
            vv = sample(0:xv,1)
            tv = c(tv,vv)}
        tv = c(tv,susg[si]-sum(tv))
        tv = sample(tv)
        rtug = rbind(rtug,tv)}
    
    colnames(rtug) = colnames(tug)[4:7]
    rtug = as.data.frame(rtug)
    rtug = cbind(tug[,1:3],rtug)
    rthantc = antcas(rtug,nrg,wd=wdset,dfrathr=dfraset,catisg,cacesg)
    raap = antcaprcts(rthantc,catisg)
    rpall = rbind(rpall,cbind(raap[,1],raap[,2],raap[,3],raap[,4],ri)) }

rpall = as.data.frame(rpall)
colnames(rpall)=c('P1','P2','C1','C2','iter')

# save(rpall1,file='rpall_1000.RData')
# load('rpall_1000.RData')

# make a table with counts and p-values and export
appv = c()
tprts = unique(c(as.character(aap$P1),as.character(aap$P2)))
for (rr in 1:(length(tprts)-1)){
    for (qq in (rr+1):length(tprts)){
        up1 = as.character(tprts[rr])
        up2 = as.character(tprts[qq])
        cnt1 = nrow(aap[aap$P1==up1 & aap$P2==up2,])
        cnt2 = nrow(aap[aap$P2==up1 & aap$P1==up2,])
        cnt = cnt1+cnt2
        if (cnt > 0){
            pvv = 1 #adding the real result
            xx1 = rpall[rpall$P1==up1 & rpall$P2==up2,]
            xx2 = rpall[rpall$P2==up1 & rpall$P1==up2,]
            xx = rbind(xx1,xx2)
            for (iti in unique(xx$iter)){
                rcnt = nrow(xx[xx$iter == iti,])
                if (rcnt >= cnt){
                    pvv = pvv+1 }}
            if (pvv > rno){pvv = rno}
            pvv = pvv/rno
            appv = rbind(appv,c(up1,up2,cnt,pvv))}}}
        
appv = as.data.frame(appv)
colnames(appv)=c('P1','P2','count','p-value')
write.table(appv,file='HPA_allMAPPtable.txt',
            row.names=FALSE,sep='\t',quote=FALSE)

# get only the significant MAPP based on FDR p-val. thr.
pvalthr = 0.05
apvs = as.numeric(as.character(appv$`p-value`))
fdrapvs = p.adjust(apvs,method = 'fdr')
appv = cbind(appv,fdrapvs)
sigappv = appv[which(fdrapvs <= pvalthr),]

# compare antagonist pairs of proteins to their coexpression
# compare number of dif.exp., weight ratios
ctprots = unique(c(as.character(aap$P1),as.character(aap$P2)))
catypes = unique(as.character(tug$Cancer))
alldifprots = c()
for (ci in catypes){
    alldifprots = c(alldifprots,thantc[[ci]][[1]])
    alldifprots = c(alldifprots,thantc[[ci]][[2]])
}
prots = unique(alldifprots)
freqdifprots = table(alldifprots)
# make a list of proteins matching cancers
thantp = list()
for (pi in prots){
    pih = c()
    pil = c()
    for (ci in catypes){
        if (pi %in% thantc[[ci]][[1]]){
            pih = c(pih,ci)}
        if (pi %in% thantc[[ci]][[2]]){
            pil = c(pil,ci)}}
    thantp[[pi]] = c(list(unique(pih)),list(unique(pil)))}

mixvalspre = c()
for (mi in 1:length(prots)){
    print(mi)
    pi = prots[mi]
    difexpno = freqdifprots[[pi]]
    psrt = 0
    ngrt = 0
    for (pj in prots){
        if (pi != pj){
            posint = length(intersect(thantp[[pi]][[1]],thantp[[pj]][[1]]))+
                length(intersect(thantp[[pi]][[2]],thantp[[pj]][[2]]))
            posuni = length(union(thantp[[pi]][[1]],thantp[[pj]][[1]]))+
                length(union(thantp[[pi]][[2]],thantp[[pj]][[2]]))
            posrat = posint/posuni
            negint = length(intersect(thantp[[pi]][[1]],thantp[[pj]][[2]]))+
                length(intersect(thantp[[pi]][[2]],thantp[[pj]][[1]]))
            neguni = length(union(thantp[[pi]][[1]],thantp[[pj]][[2]]))+
                length(union(thantp[[pi]][[2]],thantp[[pj]][[1]]))
            negrat = negint/neguni
            psrt = psrt + posrat
            ngrt = ngrt + negrat}}
        
    if (pi %in% ctprots){
        selct1 = aap[aap$P1==pi,]
        selct2 = aap[aap$P2==pi,]
        selct = rbind(selct1,selct2)
        sumsic = nrow(selct)}
    else {sumsic = 0.0}
    mixvalspre = rbind(mixvalspre,c(difexpno,psrt,ngrt,sumsic))}

#ranking weights
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
write.table(mixvals,file='HPA_mixvals_CT.txt',
            row.names=FALSE,sep='\t',quote=FALSE)

# plot differential ratio values
svg(file='HPA_weightedratio_mixvals_CT.svg',width=4.5, height=4)
alldif = alldif[order(alldif[,3]),]
limy = max(c(abs(min(alldif)),max(alldif)))
plot(alldif[,3],type='l',ylim=c((-limy-0.005),(limy+0.06)),lwd = 1.5,
     main='',xlab='',ylab='',xaxt='n',yaxt='n')
lines(alldif[,1],col='darkorchid1',lwd = 1.5)
lines(alldif[,2],col='darkorange1',lwd = 1.5)
legend(x = 0.1,y = 0.11,bty = 'n',lwd = 1.5,text.font=2,cex=1.2,
       legend=c('MAPP','PCPP','NCPP'),col = c('black','darkorchid1','darkorange1'))
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='CT Based Analysis of HPA',ylab ="Differential Ratio",
      xlab="Protein Index",font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# import PPI and TF interactions datasets; human interactions are selected
orjppints = read.delim('BIOGRID-MV-Physical-4.3.194.tab3.txt',colClasses = 'character')
orjtfints1 = read.delim('HTRIdb_data.txt',colClasses = 'character')
orjtfints2 = read.delim('trrust_rawdata.human.tsv',header=FALSE,colClasses = 'character')
colnames(orjtfints2)=c('TF','GENE','TYPE','PUBMED_ID')
hupints = orjppints[orjppints$Organism.Name.Interactor.A == 'Homo sapiens' & 
                     orjppints$Organism.Name.Interactor.B == 'Homo sapiens',]
hupints = hupints[,c(8,9)]

# check the neighborhood of antagonist pairs in human protein interaction network
# direct interactions (first degree)
prpairsw = sigappv
neprpairsw = c()
for (pi in 1:nrow(prpairsw)){
    pr1 = as.character(prpairsw[pi,1])
    pr2 = as.character(prpairsw[pi,2])
    if (pr1 %in% hupints[hupints$Official.Symbol.Interactor.A==pr2,]$Official.Symbol.Interactor.B){
        neprpairsw = rbind(neprpairsw,prpairsw[pi,])}}
neprpairsw = unique(neprpairsw)

# indirect interactions (second degree)
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

# check the neighborhood of antagonist pairs in the transcriptional network
# two datasets are combined
tfints1 = orjtfints1[,c(3,5)]
tfints2 = orjtfints2[,c(1,2)]
colnames(tfints1) = c('TF','GENE')
tfints = rbind(tfints1,tfints2)
hutfints = unique(tfints)

tfprpairsw = c()
for (pi in 1:nrow(prpairsw)){
    pr1 = as.character(prpairsw[pi,1])
    pr2 = as.character(prpairsw[pi,2])
    if (pr1 %in% hutfints[hutfints$TF==pr2,]$GENE){
        tfprpairsw = rbind(tfprpairsw,c(pr2,pr1))}
    if (pr2 %in% hutfints[hutfints$TF==pr1,]$GENE){
        tfprpairsw = rbind(tfprpairsw,c(pr1,pr2))}}
tfprpairsw = unique(tfprpairsw)

# combine all pairs in NCBI symbols and make a table
intvals = c()
prr = unique(c(as.character(prpairsw$P1),as.character(prpairsw$P2)))
for (pii in 1:(length(prr)-1)){
    for (pjj in (pii+1):length(prr)){
        ntgs = c()
        pi = prr[pii]
        pj = prr[pjj]
        if (pi %in% prpairsw[prpairsw$P1==pj,]$P2 | 
            pj %in% prpairsw[prpairsw$P1==pi,]$P2){
            if (pi %in% neprpairsw[neprpairsw$P1==pj,]$P2 | 
                pj %in% neprpairsw[neprpairsw$P1==pi,]$P2){
                ntgs = c(ntgs,'Direct')}
            else if (pi %in% ne2prpairsw[ne2prpairsw$P1==pj,]$P2 | 
                pj %in% ne2prpairsw[ne2prpairsw$P1==pi,]$P2){
                ntgs = c(ntgs,'Indirect')}
            else {ntgs = c(ntgs,'None')}
            if (pi %in% tfprpairsw[tfprpairsw$P1==pj,]$P2 | 
                pj %in% tfprpairsw[tfprpairsw$P1==pi,]$P2){
                ntgs = c(ntgs,'Direct')}
            else {ntgs = c(ntgs,'None')}
            intvals = rbind(intvals,c(pi,pj,ntgs))}}}
intvals = as.data.frame(intvals)
colnames(intvals) = c('P1','P2','PPI_type','TFint_type')
write.table(intvals,file='HPA_intvals.txt',
            row.names=FALSE,sep='\t',quote=FALSE)            

# plot degree distributions
pis = c(as.character(intvals$P1),as.character(intvals$P2))
pit = table(pis)
pid = as.data.frame(pit)
svg(file='HPA_MAPPnet_degdist.svg',width=4, height=4)
hist(pid$Freq,breaks=15,ylim=c(0,40),xlim=c(0,15),
     main='',xlab='',ylab='',xaxt='n',yaxt='n')
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='HPA MAPP network',ylab ="Frequency",xlab="Degree",
      font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()

# enlist the cancer types of the significant MAPP
# get more frequently appearing ones
# significant pairs with cancer type pairs
aapsig = c()
for (si in 1:nrow(aap)){
  if (nrow(sigappv[sigappv$P1==as.character(aap[si,1]) & sigappv$P2==as.character(aap[si,2]),])>0 | 
      nrow(sigappv[sigappv$P2==as.character(aap[si,1]) & sigappv$P1==as.character(aap[si,2]),])>0){
    aapsig = rbind(aapsig,aap[si,])}}
scct = aapsig[,3:4]

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
write.table(unscct,file='HPA_uniqueCTpairs_sctMAPP.txt',
            row.names=FALSE,sep='\t',quote=FALSE)

# get appearance ratio and compare to the frequency in the dataset
tuc = tug$Cancer
tuct = table(tuc)/length(tuc)
scs = unique(c(as.character(scct$C1),as.character(scct$C2)))
ctcount = c()
for (ci in scs){
  c1c = scct[scct$C1==ci,]
  c2c = scct[scct$C2==ci,]
  rc = (nrow(c1c)+nrow(c2c))/nrow(scct)
  dc = tuct[ci]
  ctcount = rbind(ctcount,c(ci,rc,dc))}
ctcount = as.data.frame(ctcount)
colnames(ctcount) = c('CT','MAPP_ratio','Dataset_ratio')
ctcount = ctcount[order(ctcount$MAPP_ratio),]
write.table(ctcount,file='HPA_CT_freqs.txt',
            row.names=FALSE,sep='\t',quote=FALSE)

# plot histogram of cancer type frequencies
allrts = as.numeric(as.character(ctcount$ratio))
svg(file='HPA_CT_freqs.svg',width=4, height=3)
hist(allrts,breaks=10,ylim=c(0,6),xlim=c(0,0.65),
     main='',xlab='',ylab='',xaxt='n',yaxt='n')
text(0.58, 1.7, labels = 'BC', cex = 1.2, font=2)
text(0.29, 1.7, labels = 'TC', cex = 1.2, font=2)
text(0.23, 3.7, labels = 'M,G,CC', cex = 1.2, font=2)
text(0.165, 1.7, labels = 'EC', cex = 1.2, font=2)
axis(side=2,font=2,cex.axis=1.2)
axis(side=1,font=2,cex.axis=1.2)
title(main ='Cancer Frequencies of HPA',ylab ="Frequency",xlab="Fraction of Count",
      font.lab=2,cex.lab=1.2,cex.main=1.3)
dev.off()
