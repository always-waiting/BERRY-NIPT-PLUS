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
Spec <- matrix(c(
    'input',       'i', 1, 'character',   'input file path                                                           must',
    'Chr',         'c', 1, 'integer',     'which Chr                                                                 must',
    'outdir',      'o', 1, 'character',   'output dir                                                                must',
    'prefile',     'p', 1, 'character',   'the prepare PCA files prefix                                              must',
    'rmPCs',       'R', 2, 'integer',     'need remove PCs number                                                    option',
    'rm_mean_pve', 'r', 2, 'numeric',     'the remove PCs PVE mean, this is useless when set rmPCs, default 0.7      option',
    'Merge',       'M', 2, 'integer',     'merge some bins to one, default 5                                         option',
    'redisseq',    'd', 2, 'character',   'the redis seqs matrix file, if set, will do redis seq correction          option'
), byrow=T, ncol=5)

opt <- getopt(Spec)

if (is.null(opt$input) | is.null(opt$prefile) | is.null(opt$Chr) ) {
    cat(getopt(Spec, usage=T))
    q(status=1)
}

options(scipen=200)

if (is.null(opt$Merge) ) { opt$Merge <- 5 }
if (is.null(opt$rm_mean_pve) ) { opt$rm_mean_pve <- 0.7 }
BIN_LENGTH <- 20000
bins_size <- opt$Merge * BIN_LENGTH
prefix <- gsub('.matrix', '', basename(opt$input))
prefix <- paste(opt$outdir, '/', prefix, sep='')



PCs_file <- paste(opt$prefile, '.PCs.txt', sep='')
singular_value_file <- paste(opt$prefile, '.SD_DIS.txt', sep='')
mean_sd_file <- paste(opt$prefile, '.mean_sd.txt', sep='')

input <- as.matrix(read.table(opt$input, header=F, row.names=1, stringsAsFactors=F, check.name=F, sep='\t'))
pre_PCs <- as.matrix(read.table(PCs_file, header=T, row.names=1, stringsAsFactors=F, sep='\t', check.name=F))
mean_sd <- as.matrix(read.table(mean_sd_file, header=T, row.names=1, stringsAsFactors=F, sep='\t', check.name=F))
singular_value <- as.matrix(read.table(singular_value_file, stringsAsFactors=F))


Ncol <- seq(1, ncol(input))

colnames(input) <- paste(opt$Chr, ':', (Ncol-1) * bins_size, '-', Ncol * bins_size - 1, sep="")
# some sample maybe some bins has NA, unexcept
x <- apply(input, 2, FUN=function(x){ sum(is.na(x)) == 0 })
input <- subset(input, select=x)


use.index <- which( colnames(input) %in% colnames(pre_PCs))

# x <- apply(input, 2, FUN=function(x){ sum(is.na(x)) })

# pca.use.matrix <- input[ , use.index]
pca.use.matrix <- subset(input, select=use.index)

use.index <- which (colnames(pre_PCs) %in% colnames(pca.use.matrix))
pre_PCs <- subset(pre_PCs, select=use.index)
mean_sd <- subset(mean_sd, select=use.index)

if ( ncol(pca.use.matrix) != ncol(pre_PCs) ) {
	cat('FATAL ERROR: pcas error ncols', file=stderr(), end='\n')
	q(status=1)
}
if ( ncol(pca.use.matrix) != ncol(mean_sd) ) {
	cat('FATAL ERROR: rdfile error ncols', file=stderr(), end='\n')
	q(status=1)
}

write.table(pca.use.matrix, file=paste(prefix, '.filter_na.matrix', sep=''), row.names=T, col.names=NA, quote=F, sep='\t')

## use redis seq number normalize
if ( !is.null(opt$redisseq) ) {
    library(zoo)
    Redisseq <- t(as.matrix(read.table(opt$redisseq, header=F, stringsAsFactors=F)))
    Redisseq <- Redisseq[opt$Chr, ]
    Redisseq <- rollapply(Redisseq, width=opt$Merge, by=opt$Merge, FUN=sum, partial=T, align='left')
    Redisseq <- Redisseq[use.index]
    pca.use.matrix <- t(apply(pca.use.matrix, 1, FUN=function(x) { x / Redisseq }))
    out_scale_matrix <- paste( prefix, '.redis_scale.txt', sep='')
    write.table(pca.use.matrix, file=out_scale_matrix, row.names=T, col.names=NA, quote=F, sep='\t')
}

## mean centered and output foldchange matrix

pca.use.matrix.fd <- t(apply(pca.use.matrix, 1, FUN=function(x) {x/mean_sd[1,]}))
write.table(pca.use.matrix.fd, file=paste(prefix, '.fd.matrix', sep=''), row.names=T, col.names=NA, quote=F, sep='\t')
fd_constancy <- apply(pca.use.matrix.fd, 1, FUN=function(x) { c(mean(x), sd(x)) })
colnames(fd_constancy) <- rownames(pca.use.matrix)
rownames(fd_constancy) <- c('Mean', 'Sd')
out_fd_mean_sd <- paste(prefix, '.fd.mean_sd.matrix', sep='')
write.table(fd_constancy, file=out_fd_mean_sd, row.names=T, col.names=NA, quote=F, sep='\t')

pca.use.matrix.center <- t(apply(pca.use.matrix, 1, FUN=function(x) { x - mean_sd[1, ] })) 

write.table(pca.use.matrix.center, file=paste(prefix, '.centered.matrix', sep=''), row.names=T, col.names=NA, quote=F, sep='\t')


## remove PCs
if (is.null(opt$rmPCs)) {
    Remove.pc.number <- which(singular_value[, 2] ^ 2 >= (sum(singular_value[, 2] ^ 2) / nrow(pre_PCs)) * opt$rm_mean_pve)
    if ( length(Remove.pc.number) == 0 ) {
        Remove.pc.number <- c()
    }
    out_sd_remove_number <- paste( prefix, '.SD_REMOVE_NUMBER.txt', sep='')
    cat(paste(Remove.pc.number, sep='\t'), end='\n', file=stderr())
    cat(max(Remove.pc.number), end='\n', file=out_sd_remove_number)
}else{
    Remove.pc.number <- seq(1, opt$rmPCs)
    out_sd_remove_number <- paste(prefix, '.SD_REMOVE_NUMBER.txt', sep='')
    cat(paste(Remove.pc.number, sep='\t'), end='\n', file=stderr())
    cat(max(Remove.pc.number), end='\n', file=out_sd_remove_number)
}

#### use a formula to remove PCs

for (i in seq(1, nrow(pca.use.matrix.center)) ) {
	for (j in Remove.pc.number ) {
		principal <- sum(pca.use.matrix.center[i,]*pre_PCs[j])
		# cat(paste(principal,i,j,sep='\t'),end='\n', file=stderr())
		pca.use.matrix.center[i, ] <- pca.use.matrix.center[i, ] - pre_PCs[j, ] * sum( pca.use.matrix.center[i, ] * pre_PCs[j,])
	}
}

write.table(pca.use.matrix.center, file=paste(prefix, '.normalized.matrix', sep=''), row.names=T, col.names=NA, quote=F, sep='\t')
 
pca.use.matrix.center.zscore <- t(apply(pca.use.matrix.center, 1, FUN=function(x) { (x - mean(x))/sd(x) }))

write.table(pca.use.matrix.center.zscore, file=paste(prefix, '.normalized.zscore.matrix', sep=''), row.names=T, col.names=NA, quote=F, sep='\t')
#### plot every sample's zscore distribution

coor_x <- colnames(pca.use.matrix.center.zscore)
coor_x <- do.call(rbind, strsplit(coor_x, "[:-]", perl=T))
coor_x <- as.numeric(coor_x[, 2])
pdf(file=paste(prefix, '.zscore.dis.pdf', sep=''), width=12, height=9)
par(mfrow=c(3,1))
for (i in seq(1, nrow(pca.use.matrix.center.zscore))) {
	plot(coor_x/1000000, pca.use.matrix.center.zscore[i,],xlab='Pos(M)',ylab='zscore',main=rownames(pca.use.matrix.center)[i],pch=20,col='blue')
}
dev.off()

