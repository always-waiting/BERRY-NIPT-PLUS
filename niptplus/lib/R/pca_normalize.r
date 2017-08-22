#!/usr/bin/Rscript
#===============================================================================
#  Copyright (c)   berrygenomics 2015
#
#  DESCRIPTION: PCA
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

spec <- matrix(c(
    'help',         'h', 0, 'logical',   'this help',
    'input',        'i', 1, 'character', 'the input file path, read depth by chrs                                   must',
    'outdir',       'o', 1, 'character', 'output dir                                                                must',
    'Chr',          'c', 1, 'integer',   'the chr which the data come from, integer                                 must',
    'redisseq',     'R', 2, 'character', 'the redis seqs matrix file                                                option',
    'rmPCs',        'n', 2, 'integer',   'the remove PCs number                                                     option',
    'rm_mean_pve',  'r', 2, 'numeric',   'the remove PCs PVE mean, this is useless when set rmPCs, default 0.95     option',
    'Merge',        'M', 2, 'integer',   'merge some bins to one, default 5                                         option',
    'nocentered',   't', 0, 'logical',   'whether to centered, logical, default mean centering                      option'
), byrow=T, ncol=5 )
## usage
opt <- getopt(spec)
options(scipen=99)
if ( !is.null(opt$help) | is.null(opt$input)  | is.null(opt$outdir) | is.null(opt$Chr) ) {
    cat(getopt(spec, usage=T))
    q(status=1)
}

centered <- TRUE
if ( is.null(opt$Merge) ) { opt$Merge <- 5 }
if ( is.null(opt$rm_mean_pve) ) { opt$rm_mean_pve <- 0.95 }
if ( !is.null(opt$nocentered) ) { centered <- FALSE }
BIN_LENGTH <- 20000
bins_size <- opt$Merge * BIN_LENGTH
prefix <- gsub('.matrix', '', basename(opt$input))
prefix <- paste(opt$outdir, '/', prefix, sep='')
###### processing start
library(reshape2)
library(ggplot2)

pca.matrix <- as.matrix(read.table(opt$input, header=F, row.names=1, stringsAsFactors=F, sep='\t'))
Ncol <- seq(1, ncol(pca.matrix))
Nrow <- nrow(pca.matrix)

colnames(pca.matrix) <- paste(opt$Chr, ':', (Ncol - 1) * bins_size, '-', Ncol * bins_size - 1, sep='')

x <- apply(pca.matrix, 2, FUN=function(x){ sum(is.na(x)) })

pca.use.matrix <- pca.matrix[ , x == 0]        ## must no NA

out_filter_na <- paste(prefix, '.filter_na.matrix', sep='')
write.table(pca.use.matrix, file=out_filter_na, row.names=T, col.names=NA, quote=F, sep='\t')

## use redis seq number normalize
if ( !is.null(opt$redisseq) ) {
    library(zoo)
    Redisseq <- t(as.matrix(read.table(opt$redisseq, header=F, stringsAsFactors=F)))
    Redisseq <- Redisseq[opt$Chr, ]
    Redisseq <- rollapply(Redisseq, width=opt$Merge, by=opt$Merge, FUN=sum, partial=T, align='left')
    Redisseq <- Redisseq[x==0]
    pca.use.matrix <- t(apply(pca.use.matrix, 1, FUN=function(x) { x / Redisseq }))
    out_scale_matrix <- paste(prefix, '.redis_scale.txt', sep='')
    write.table(pca.use.matrix, file=out_scale_matrix, row.names=T, col.names=NA, quote=F, sep='\t')
}
## mean centered or median centered ?
#	Mean <- apply(pca.use.matrix, 2, median)
Mean <- apply(pca.use.matrix, 2, mean)
Sd   <- apply(pca.use.matrix, 2, sd)
mean_sd_matrix <- t(data.frame(Mean, Sd))
colnames(mean_sd_matrix) <- colnames(pca.use.matrix)
out_mean_sd <- paste(prefix, '.mean_sd.txt', sep='')
write.table(mean_sd_matrix, file=out_mean_sd, row.names=T, col.names=NA, quote=F, sep='\t')
# output foldchange file
pca.use.matrix.fd <- apply(pca.use.matrix, 2, FUN=function(x) { x / mean(x) })
out_fd_matrix <- paste(prefix, '.fd.matrix', sep='')
write.table(pca.use.matrix.fd, file=out_fd_matrix, row.names=T,col.names=NA, quote=F, sep='\t')
fd_constancy <- apply(pca.use.matrix.fd, 1, FUN=function(x) { c(mean(x), sd(x)) })
colnames(fd_constancy) <- rownames(pca.use.matrix)
rownames(fd_constancy) <- c('Mean', 'Sd')
out_fd_mean_sd <- paste(prefix, '.fd.mean_sd.matrix', sep='')
write.table(fd_constancy, file=out_fd_mean_sd, row.names=T, col.names=NA, quote=F, sep='\t')

#	pca.use.matrix <- apply(pca.use.matrix, 2, FUN=function(x){ x - median(x) })
pca.use.matrix <- apply(pca.use.matrix, 2, FUN=function(x){ x - mean(x) })
out_centered_matrix <- paste(prefix, '.centered.matrix', sep='')
write.table(pca.use.matrix, file=out_centered_matrix, row.names=T, col.names=NA, quote=F, sep='\t')


## singular value decomposition , plot singular value distribution, output PCs 
svd_data <- svd(pca.use.matrix)
D <- diag(svd_data$d)
V <- t(svd_data$v)
rownames(V) <- paste('PC', seq(1, nrow(V)), sep='')
colnames(V) <- colnames(pca.use.matrix)
U <- svd_data$u

out_sd_dis_pdf <- paste(prefix, '.SD_DIS.png', sep='')
png(file=out_sd_dis_pdf, height=800, width=600)
par(mfrow=c(2, 1), lend=1)
svd_data_sd_x <- seq(1, length(svd_data$d))
plot(svd_data_sd_x, svd_data$d, type='h', col='blue', xlab='singular values', ylab='value', main='singular value dis')
plot(svd_data_sd_x[1:100], svd_data$d[1:100], type='h', col='blue', xlab='singular values', ylab='value', main='100 singular value dis', xlim=c(0,100))
dev.off()

singular_value <- svd_data$d
contrib <- singular_value ^ 2
contrib_average <- sum(contrib) / Nrow
contrib <- contrib / contrib_average 
singular_value_matrix <- data.frame(ID=seq(1, length(singular_value)), singluarvalue=singular_value, contrib=contrib) 
out_sd_dis_txt <- paste(prefix, '.SD_DIS.txt', sep='')
write.table(singular_value_matrix, file=out_sd_dis_txt, row.names=F, col.names=T, quote=F, sep='\t')

out_PCs_txt <- paste(prefix, '.PCs.txt', sep='')
write.table(V, file=out_PCs_txt, row.names=T, col.names=NA, quote=F, sep='\t')


## get PVE mean
if (is.null(opt$rmPCs)) {
    Remove.pc.number <- which(contrib >= contrib_average * opt$rm_mean_pve)
    if ( length(Remove.pc.number) == 0 ) {
        Remove.pc.number <- c(0)
    }
    out_sd_remove_number <- paste(prefix, '.SD_REMOVE_NUMBER.txt', sep='')
    cat(paste(Remove.pc.number, sep='\t'), end='\n', file=stderr())
    cat(max(Remove.pc.number), end='\n', file=out_sd_remove_number)
}else{
    Remove.pc.number <- seq(1, opt$rmPCs)
    out_sd_remove_number <- paste(prefix, '.SD_REMOVE_NUMBER.txt', sep='')
    cat(paste(Remove.pc.number, sep='\t'), end='\n', file=stderr())
    cat(max(Remove.pc.number), end='\n', file=out_sd_remove_number)
}

## set singular value to 0
if (Remove.pc.number[1] != 0) {
    for ( i in Remove.pc.number) {
        D[i, i] <- 0
    }
}
## get normalize result
pca.normal.matrix <- svd_data$u %*% D %*% t(svd_data$v)

rownames(pca.normal.matrix) <- rownames(pca.use.matrix)
colnames(pca.normal.matrix) <- colnames(pca.use.matrix)

out_normal_matrix <- paste(prefix, '.normalized.matrix', sep='')
write.table(pca.normal.matrix, file=out_normal_matrix, row.names=T, col.names=NA, quote=F, sep='\t')

## the distribution of mean&sd from the final matrix, and cal z-score for HMM call cnv 

pca.normal.matrix.MeanSd <- apply(pca.normal.matrix, 1, FUN=function(x){c( mean(x), sd(x))})
out_normal_matrix_mean_sd_dis <- paste(prefix, '.normal.matrix.BySampleDis.pdf', sep='')
pdf(file=out_normal_matrix_mean_sd_dis)
cur_mean <- pca.normal.matrix.MeanSd[1,]
hist(cur_mean, breaks=50, freq=F, main='mean distribution by sample', xlab='mean', ylab='density')
lines(density(cur_mean), col='blue', lwd=3)

cur_sd <- pca.normal.matrix.MeanSd[2,]
hist(cur_sd, breaks=50, freq=F, main='sd disribution by sample', xlab='sd', ylab='density')
lines(density(cur_sd), col='blue', lwd=3)
dev.off()


pca.normal.matrix.zscore <- t(apply(pca.normal.matrix, 1, FUN=function(x){ (x - mean(x)) / sd(x) }))
out_normal_matrix_zscore <- paste(prefix, '.normalized.zscore.matrix', sep='')
write.table(pca.normal.matrix.zscore, file=out_normal_matrix_zscore, row.names=T, col.names=NA, quote=F, sep='\t')

## plotgraph, the colnames will be change

x.coordinate <- which(x==0) * opt$size

colnames(pca.use.matrix) <- x.coordinate

Temp1 <- melt(pca.use.matrix)
colnames(Temp1) <- c('sample', 'pos', 'value')
P1 <- qplot(data=Temp1, x=pos , y=value, geom='line', colour=sample, ylim=c(min(Temp1$value), max(Temp1$value)), main=paste('chr', opt$Chr, '  before pca', sep=''), xlab='pos', ylab='value')+theme(legend.position = 'none')


colnames(pca.normal.matrix) <- x.coordinate
Temp2 <- melt(pca.normal.matrix)
colnames(Temp2) <- c('sample', 'pos', 'value')
P2 <- qplot(data=Temp2, x=pos , y=value, geom='line', colour=sample, ylim=c(min(Temp1$value), max(Temp1$value)), main=paste('chr', opt$Chr, '  after pca', sep=''), xlab='pos', ylab='value')+theme(legend.position = 'none')


 
colnames(pca.normal.matrix.zscore) <- x.coordinate
Temp3 <- melt(pca.normal.matrix.zscore)
colnames(Temp3) <- c('sample', 'pos', 'value')
P3 <- qplot(data=Temp3, x=pos , y=value, geom='line', colour=sample, ylim=c(-5,5), main=paste('chr', opt$Chr, '  after pca', sep=''), xlab='pos', ylab='zscore')+theme(legend.position = 'none')

library(gridExtra)
png(file=paste(prefix, '.png', sep=''), width=800, height=1200)

grid.arrange(P1, P2, P3, nrow=3)

dev.off()
