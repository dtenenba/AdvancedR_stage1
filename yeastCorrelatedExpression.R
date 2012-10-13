correlationFinder <- function() {
  dataFile="sub_combined_complete_dataset_526G_198E.txt"
  cor.threshold <- 0.85

  tbl <- read.table(dataFile, sep='\t', header=T, quote='',
    comment.char='', fill=T, stringsAsFactors=FALSE)
  rownames(tbl) <- tbl$X
  exclude.these.columns <-  which(sapply(1:ncol(tbl),
    function(col) class(tbl [,col])) != 'numeric')
  if(length(exclude.these.columns) > 0)
    tbl <- tbl [, -exclude.these.columns]
  mtx.cor <- cor(t(as.matrix(tbl)), use='pairwise.complete.obs')
  mtx.cor <- upper.tri(mtx.cor) * mtx.cor
  max = nrow(mtx.cor)

  correlated.genes <- list()
  ret <- list()
  for(r in 1:max) {
    zz = as.integer(which(mtx.cor [r,] > cor.threshold))
    if(length(zz) > 0) {
      gene.a = rownames(mtx.cor) [r]
      genes.b = rownames(mtx.cor) [zz]
      correlated.genes[[gene.a]] <- genes.b
      ret[[ rownames(mtx.cor)[r] ]] <-
        rownames(mtx.cor)[zz]
      } # if length
     } # for r
  ret
}
