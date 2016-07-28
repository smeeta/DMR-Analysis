### R code from vignette source 'DMRcate.Rnw'
### Tim Peters
###################################################
### code chunk number 1: bioconductor (eval = FALSE)
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite("DMRcate")


###################################################
### code chunk number 2: libr
###################################################
library(DMRcate)


###################################################
### code chunk number 3: loaddata
###################################################
data(dmrcatedata)
myMs <- logit2(myBetas)


###################################################
### code chunk number 4: filter
###################################################
nrow(illuminaSNPs)
nrow(myMs)
myMs.noSNPs <- rmSNPandCH(myMs, dist=2, mafcut=0.05)
nrow(myMs.noSNPs)


###################################################
### code chunk number 5: annotate
###################################################
patient <- factor(sub("-.*", "", colnames(myMs)))
type <- factor(sub(".*-", "", colnames(myMs)))
design <- model.matrix(~patient + type) 
myannotation <- cpg.annotate("array", myMs.noSNPs, analysis.type="differential",
                             design=design, coef=39)


###################################################
### code chunk number 6: dmrcate
###################################################
dmrcoutput <- dmrcate(myannotation, lambda=1000, C=2)


###################################################
### code chunk number 7: ranges
###################################################
results.ranges <- extractRanges(dmrcoutput, genome = "hg19")
results.ranges


###################################################
### code chunk number 8: plotting
###################################################
groups <- c(Tumour="magenta", Normal="forestgreen")
cols <- groups[as.character(type)]
samps <- c(1:6, 38+(1:6))
DMR.plot(ranges=results.ranges, dmr=1, CpGs=myBetas, phen.col=cols, genome="hg19", 
         samps=samps)


###################################################
### code chunk number 9: wgbssuitecpgs
###################################################
CpGs


###################################################
### code chunk number 10: prepareDSS
###################################################
meth <- as.data.frame(CpGs)[,c(1:2, grep(".C$", colnames(as.data.frame(CpGs))))]
coverage <- as.data.frame(CpGs)[,c(1:2, grep(".cov$", colnames(as.data.frame(CpGs))))]

treat1 <- data.frame(chr=coverage$seqnames, pos=coverage$start, 
                     N=coverage$Treatment1.cov, X=meth$Treatment1.C)

treat2 <- data.frame(chr=coverage$seqnames, pos=coverage$start, 
                     N=coverage$Treatment2.cov, X=meth$Treatment2.C)

treat3 <- data.frame(chr=coverage$seqnames, pos=coverage$start, 
                     N=coverage$Treatment3.cov, X=meth$Treatment3.C)

ctrl1 <- data.frame(chr=coverage$seqnames, pos=coverage$start, 
                    N=coverage$Control1.cov, X=meth$Control1.C)

ctrl2 <- data.frame(chr=coverage$seqnames, pos=coverage$start, 
                    N=coverage$Control2.cov, X=meth$Control2.C)

ctrl3 <- data.frame(chr=coverage$seqnames, pos=coverage$start, 
                    N=coverage$Control3.cov, X=meth$Control3.C)

samples <- list(treat1, treat2, treat3, ctrl1, ctrl2, ctrl3)
sampnames <- sub("\\..*", "", colnames(meth))[-c(1:2)]

obj_bsseq <- makeBSseqData(samples, sampnames)
DSSres <- DMLtest(obj_bsseq, group1=sampnames[1:3], group2=sampnames[4:6], smoothing=FALSE) 



###################################################
### code chunk number 11: wgbsDMRcate
###################################################
wgbsannot <- cpg.annotate("sequencing", DSSres)
wgbs.DMRs <- dmrcate(wgbsannot, lambda = 1000, C = 50, pcutoff = 0.05, mc.cores = 1)
wgbs.ranges <- extractRanges(wgbs.DMRs, genome = "hg19")
groups <- c(Treatment="darkorange", Control="blue")
cols <- groups[sub("[0-9]", "", sampnames)]
DMR.plot(ranges=wgbs.ranges, dmr=1, CpGs=CpGs, phen.col=cols, genome="hg19")


###################################################
### code chunk number 12: sessionInfo
###################################################
sessionInfo()

##################################################
