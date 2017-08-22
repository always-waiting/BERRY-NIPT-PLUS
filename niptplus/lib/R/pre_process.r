#!/usr/bin/Rscript
#===============================================================================
#  Copyright (c)   berrygenomics 2015
#
#  DESCRIPTION: 包括几个处理过程：1，合并bin  2，去除异常GC的bin 3，去除低复杂度的bin
#               4，去除低覆盖度的bin  5，如果需要，则进行GC校正
#
#      OPTIONS: ---
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: sunhy
#         MAIL: sunhuaiyu@berrygnomics.com
#      VERSION: 2.0
#      CREATED: 
#===============================================================================
library(getopt)
Spec <- matrix(c(
    'rcfile',    'r', 1, 'character', 'the reads count file path,            must',
    'gcfile',    'g', 1, 'character', 'the GC count file path,               must',
    'outdir',    'o', 1, 'character', 'the output dir path,                  must',
    'complexity','c', 2, 'character', 'the bins low-complexity matrix path,  option', 
    'Merge',     'M', 2, 'integer',   'merge some bins to one,               default 5',
    'gccorrect', 'G', 0, 'integer',   'whether do GC correction,             default no',
    'help',      'h', 0, 'logical',   'This help'
), byrow=TRUE, ncol=5)
opt <- getopt(Spec)
# if help was asked for print a friendly message
# and exit with a non-zero error code
if ( !is.null(opt$help) | is.null(opt$rcfile) | is.null(opt$gcfile) | is.null(opt$outdir)) {
    cat(getopt(Spec, usage=TRUE));
    q(status=1);
}
# set options
if (is.null(opt$Merge)) { opt$Merge <- 5 }
if (is.null(opt$gccorrect)) {opt$gccorrect <- FALSE } else { opt$gccorrect <- TRUE }
BIN_LENGTH <- 20000   # the bins of the RD matrix are 20K


if ( !file.exists(opt$outdir) ) dir.create( path=opt$outdir )
## file name : /file/path/14C01700_L6_I333.R1.clean.fastq.gz.20K.txt
prefix <- strsplit(basename(opt$rcfile), '[.]', perl=T)
outprefix = paste(opt$outdir, '/', prefix[[1]][1], sep='')

logfile = paste(outprefix, '.filter.log', sep='')
if (file.exists(logfile)) { file.remove(logfile) }


############## processing start
## read
RC <- as.matrix(read.table(opt$rcfile, header=F, stringsAsFactors=F, row.names=1, sep='\t'))
GC <- as.matrix(read.table(opt$gcfile, header=F, stringsAsFactors=F, row.names=1, sep='\t'))


## only analysis autosome
RC <- RC[1:22, ]
GC <- GC[1:22, ]

## constructing a new matrix according to the merge number, the final left bins will be persist
library(zoo)
RC <- t(apply(RC, 1, FUN=function(x){ rollapply(x, width=opt$Merge, by=opt$Merge, FUN=sum, partial=T, align='left')} ))
GC <- t(apply(GC, 1, FUN=function(x){ rollapply(x, width=opt$Merge, by=opt$Merge, FUN=sum, partial=T, align='left')} ))


## transform GCs to GC rate (percentage), every reads length are 36.
GC <- 100 * GC / ( RC * 36 )
GC[is.nan(GC)] <- 0
GC <- round(GC, 2)

cat('The matrix dim', dim(RC), end='\n', file=logfile)
cat('usefull bins:', sum(RC != 0), end='\n', file=logfile, append=T)

## filter some extreme GC value and extreme RC, remove bins GC percent > 80% or < 20%
cat('remove extreme GC bins:', sum(GC != 0 & (GC<20 | GC>80)), end='\n', file=logfile, append=T)

GC[GC < 20 | GC > 80] <- NA
RC[is.na(GC)] <- NA




## filter low-complexity bins, remove bins which have more 80% repeat regions(considered as low-complexity DNA bin)

if (!is.null(opt$complexity)) {
    Complex_data <- as.matrix(read.table(opt$complexity, header=F, row.names=1, stringsAsFactors=F))
    Complex_data <- Complex_data[1:22, ]
    Complex_data <- t(apply(Complex_data, 1, FUN=function(x){ 
        rollapply(x, width=opt$Merge, by=opt$Merge, FUN=sum, partial=T, align='left')
        } ))
    
    if (nrow(RC) != nrow(Complex_data) | ncol(RC) != ncol(Complex_data)) {
        cat('Warning: the complexity file has not same dim', end='\n', file=stderr())
        quit(status=1)
    }
    Complex_data <- Complex_data/(BIN_LENGTH * opt$Merge)
    OutLier <- Complex_data > 0.8   # 80%
    cat('remove low-complexity bins:', sum(OutLier), end='\n', file=logfile, append=T)
	cat('low-complexity bins RD:', RC[OutLier], end='\n', file=logfile, append=T)
	cat('low-complexity bins RD mean:', mean(RC[OutLier]), end='\n', file=logfile, append=T)
    RC[OutLier] <- NA
    GC[OutLier] <- NA
}


## plot RCs density, and remove low-coverage bins, remove x < 1/3 of peak value 
density_plot <- function(data, rate, main='RC density distribution') {
	Den <- density(data, na.rm=T)
	plot(Den, type='o', xlab='after RC', ylab='density', col='#4682B4', axes=T, xpd=F, pch=1, font.lab=1,2,
		 font.main=2, font.axis=1, lwd=2, cex.lab=1.5, cex.main=1.5, cex.axis=1, cex.sub=1, cex=0.5,
		 main=main)
	den_y_peak <- max(Den$y)
	den_x_peak <- Den$x[Den$y == den_y_peak]
	segments(den_x_peak, den_y_peak, den_x_peak, 0, col='black', lwd=1, lty=2)
	rm_cutoff <- den_x_peak * rate
	segments(rm_cutoff, den_y_peak, rm_cutoff, 0, col='red', lwd=2)
	legend('topright', legend=c(sprintf('Mean: %.4f', mean(data, na.rm=T)),
		 sprintf('Max density: %.4f', den_x_peak), sprintf('%.3f max density: %.4f', rate, rm_cutoff)),
		 cex=0.8,inset = 0.01)
	box()
	rm_cutoff
}


## GC cerrection: get loess curve
loessFit <- function(outprefix, rc, gc){
    rc <- as.vector(rc)
    rc <- rc[!is.na(rc)]
    gc <- as.vector(gc)
    gc <- gc[!is.na(gc)]

    png(file=paste(outprefix, '.GCcorrelation_before.png', sep=''), height=600, width=900)
    par(mar=c(5.1,4.5,4.1,2.1), mgp=c(3,1,0), ps=15)
    plot(gc, rc, main='Correlation between GC% and Read count', xlab='CG%', ylab='Read count', col='#A0522D',
        cex.lab=1.5, cex.main=1.5, cex.axis=1, cex.sub=1, cex=0.5, font.lab=1.2, font.main=2, font.axis=1, lwd=1,
        ylim=c(0, max(rc, na.rm=T)), xlim=c(10, 90), pch=21
    )

#	## fast loess
#	RC.bins.median <- tapply(rc, gc, function(x) median(x, na.rm=T))
#	GC.Levels <- as.numeric(names(RC.bins.median))
#	Loess <- loess(RC.bins.median ~ GC.Levels)
#   loess.fittedRC  <- predict( Loess, GC.Levels)
#	lines(GC.Levels, loess.fittedRC, col=3, lwd=2)

    ## normal loess
    Loess <- loess(rc ~ gc)
    loess.fittedRC  <- predict( Loess, gc[order(gc)])
    lines(gc[order(gc)], loess.fittedRC, col=3, lwd=2)

    # Lowess<-lowess(DF,f=0.05,delta=0.001)
    # lines(Lowess,col=4,lwd=2)

    M <- median(rc)
    segments(min(gc), M, max(gc), M, col=1, lwd=2, lty=2)
    legend('topright', c('loess'), col=c(3,4), lty = 1, pch = 20, inset = .02)
    dev.off()
    cat('Done loessFit\n',file=stderr())
    Loess
}


## rmove low coverage bins

pdf(file=paste(outprefix, 'RC_density.pdf', sep='_'), width=12, height=9)
par(mfrow=c(3,1), lend=1)
cutOff <- density_plot(RC, 1/3)
for (i in c(1:22)) {
	density_plot(RC[i, ], 1/3, paste('Rc density by chr', i, sep=''))
}
dev.off()
cat('Done plot \n', file=stderr())
cat('remove low-coverage bins:', sum( !is.na(RC) & RC < cutOff ), end='\n', file=logfile, append=T)
RC[!is.na(RC) & RC<cutOff] <- NA
GC[is.na(RC)] <- NA

cat('finally, persist bins:', sum(!is.na(RC)), end='\n', file=logfile, append=T) 

if (opt$gccorrect) {
    MedianRC <- median(RC, na.rm=T)
    Loess <- loessFit(outprefix, RC, GC)
    loessFited <- predict( Loess, as.vector(GC))

    # RCgcCorrection <- ( MedianRC / loessFited ) * RC
    RCgcCorrection <- RC + ( MedianRC - loessFited)

    RCgcCorrection[!is.na(RCgcCorrection) & RCgcCorrection <= 0] <- NA

    ## after GC correction, RC distribution
    png(file=paste(outprefix, '.GCcorrelation_after.png', sep=''), height=600, width=900)
    plot(GC, as.vector(RCgcCorrection), main='after Correlation between GC% and Read count',
        xlab='CG%', ylab='Read count', col='#A0522D', cex.lab=1.5, cex.main=1.5, cex.axis=1, cex.sub=1, cex=0.5,
        font.lab=1.2, font.main=2, font.axis=1, lwd=1, ylim=c(0, max(RC, na.rm=T)), xlim=c(10, 90), pch=21
    )
    dev.off()

    cat('Done GC Correction\n', file=stderr())
    
    ## out after GC correction, RC matrix 
    RCgcCorrection.tatal.reads <- sum(RCgcCorrection, na.rm=T)
    Matrixoutfile <- paste(outprefix, '.gc.filter.', opt$Merge * 20, 'K.txt', sep='')
    if (file.exists(Matrixoutfile)) { unlink(Matrixoutfile) }
    cat(paste('#uniqMapped', round(RCgcCorrection.tatal.reads, 2), sep='\t'), end='\n', file=Matrixoutfile)
    rownames(RCgcCorrection) <- paste('chr', c(1:22), sep='')
    write.table(RCgcCorrection, file=Matrixoutfile, col.names=F, row.names=T, append=T, quote=F, sep='\t')

    ## every chr's dis
    pdf(file=paste(outprefix, '.chromosome_dis.pdf', sep=''), height=9, width=12)
    par(mfrow=c(3, 1), lend=1)
    for (i in c(1:22)) {
       Chrs <- RCgcCorrection[i, ]
       last <- max(which(!is.na(Chrs)))
       x <- seq(1, last)
       plot(x * opt$Merge * BIN_LENGTH / 1000000, Chrs[x], type='h', xlab='pos(M)', ylab='reads count', col='gray', main=rownames(RCgcCorrection)[i])
    }
    dev.off()
}else{
    ## out after filtered RC matrix 
    RCgcCorrection.tatal.reads <- sum(RC, na.rm=T)
    Matrixoutfile <- paste(outprefix, '.filter.', opt$Merge * 20, 'K.txt', sep='')
    if (file.exists(Matrixoutfile)) { unlink(Matrixoutfile) }
    cat(paste('#uniqMapped', round(RCgcCorrection.tatal.reads, 2), sep='\t'), end='\n', file=Matrixoutfile)
    rownames(RC) <- paste('chr', c(1:22), sep='')
    write.table(RC, file=Matrixoutfile, col.names=F, row.names=T, append=T, quote=F, sep='\t')

    ## every chr's dis
    pdf(file=paste(outprefix, '.chromosome_dis.pdf', sep=''), height=9, width=12)
    par(mfrow=c(3,1), lend=1)
    for (i in c(1:22)) {
        Chrs <- RC[i, ]
        last <- max(which(!is.na(Chrs)))
        x <- seq(1, last)
        plot(x * opt$Merge * BIN_LENGTH / 1000000, Chrs[x], type='h', xlab='pos(M)', ylab='reads count', col='gray', main=rownames(RC)[i])
    }
    dev.off()
}




