#' Manhattan plot
#' 
#' This function takes a dataframe of GWAS results as input, which
#' requires columns named: "snp", "chr", "bp", and "p", and returns 
#' a Manhattan plot. This function assumes snp position is in 
#' base pairs "bp". Options are: genome-wide significance threshold, 
#' chromosome colors, and snps to highlight.
#' 
#' @param x Dataframe of GWAS results. Requires columns named: "snp", "chr", "bp", "p".
#' @param chrom.colors Chromosome colors to alternate
#' @param threshold Genome-wide significance threshold
#' @param highlight Logical. Highlight character vector of snps in x
#' @param p p-value 
#' @return A Manhattan plot of GWAS results
#' @export

manhattan_plot <- function(x, chrom.colors, chrLabels, genomewideline=5, highlight=NULL, is.logp=TRUE, ...) {
  
  CHR = BP = P = index = NULL
  
  if (!is.null(x[[snp]])) {
    d <- data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]], pos=NA, index=NA, SNP=x[[snp]], stringsAsFactors=FALSE)
  } else {
    d <- data.frame(CHR=x[[chr]], BP=x[[bp]], P=x[[p]], pos=NA, index=NA)
  }
  
  d <- d[order(d$CHR, d$BP), ]
  
  if (is.logp) {
    d$logp <- d$P
  } else {
    d$logp <- -log10(d$P)
  }
  
  d$index <- rep.int(seq_along(unique(d$CHR)), times=tapply(d$SNP, d$CHR, length))
  lastbase <- 0
  ticks <- NULL
  
  for (i in unique(d$index)) {
    if (i == 1) {
      d[d$index == i, ]$pos <- d[d$index == i, ]$BP
    } else {
      lastbase <- lastbase + max(d[d$index == (i - 1), "BP"])
      d[d$index == i, "BP"] <- d[d$index == i, "BP"] - min(d[d$index == i, "BP"]) + 1
      d[d$index == i, "pos"] <- d[d$index == i, "BP"] + lastbase
    }
  }
  
  ticks <- tapply(d$pos, d$index, quantile, probs=0.5)
  labs <- unique(d$CHR)
  
  def_args <- list(
    xaxt="n", 
    bty="n", 
    xaxs="i", 
    axes=FALSE, 
    cex.lab=1,
    yaxs="i", 
    las=1, 
    pch=20, 
    xlim=c(floor(max(d$pos) * -0.03), ceiling(max(d$pos) * 1.03)),
    ylim=c(0, ceiling(max(d$logp))*1.1), 
    xlab="", 
    ylab=""
  )
  
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
  
  if (genomewideline) {
    segments(x0=0, x1=max(d$pos), y0=genomewideline, y1=genomewideline, col="cyan3", lty=1, lwd=1, lend=2)
  } else {
    segments(x0=0, x1=max(d$pos), y0=5, y1=5, col="cyan4", lty=1, lwd=1, lend=2)
  }
  
  chrom.colors <- rep_len(chrom.colors, max(d$index))
  icol <- 1
  for (i in unique(d$index)) {
    points(d[d$index == i, "pos"], d[d$index == i, "logp"], col=chrom.colors[icol],  cex=1, pch=20, ...)
    icol <- icol + 1
  }
  
  chrList <- unique(d$CHR)
  meanrunList <- list()
  d$runid <- c(1:nrow(d))
  
  for (i in 1:length(chrList)) {
    sdat <- d[d$CHR %in% chrList[i], ]
    meanrunList[[i]] <- mean(c(max(sdat$pos), min(sdat$pos)))
  }
  
  axis(1, tck=-0.06, cex.axis=0.85, lty=1, lend=1, col=NA, col.ticks=NA, lwd.ticks=0, tck=-0.02, at=unlist(meanrunList), 
       labels=chrLabels, font=1, mgp=c(0.5, 0.5, 0), ...)
  axis(2, line=-0.85, tck=-0.02, cex.axis=0.85, lend=1, lwd.ticks=1, mgp=c(0.25, 0.5, 0), ...)
  title(ylab=expression('-log'[10]*' ('*italic('p')*')'), mgp=c(1, 0.5, 0), cex.lab=1)
  title(xlab="Chromosome", mgp=c(1.25, 0.5, 0), cex.lab=1)
  
  if (!is.null(highlight)) {
    d.highlight <- d[which(d$SNP %in% highlight), ]
    with(d.highlight, points(pos, logp, col="firebrick2", cex=1, pch=20, ...))
  }
  
  par(xpd=FALSE)
  
}
