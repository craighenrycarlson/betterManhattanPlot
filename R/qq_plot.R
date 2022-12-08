#' QQ-plot
#' 
#' This function takes a vector of pvalues from GWAS results as input,
#' and returns a QQ-plot 
#' 
#' @param pvector Vector of p-values
#' @param traitName Name of trait to plot. 
#' @return A QQ plot of GWAS results
#' @export

qq_plot <- function (pvector, traitName="", ...){
  
  pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & is.finite(pvector) & pvector < 1 & pvector > 0]
  o <- -log10(sort(pvector, decreasing=FALSE))
  e <- -log10(ppoints(length(pvector)))
  
  def_args <- list(
    pch=20, 
    xlim=c(0, max(e)*1.1), 
    bty='n', 
    ylim=c(0, ceiling(max(o))*1.1), 
    xlab="", 
    ylab="",
    axes=FALSE, 
    frame=FALSE, 
    cex.lab=1, 
    cex=0.85
  )
  
  dotargs <- list(...)
  tryCatch(do.call("plot", c(list(x=e, y=o), bty='n', main=traitName, cex.main=0.85, def_args[!names(def_args) %in% names(dotargs)], dotargs)), warn=stop)
  
  axis(1, tck=-0.02, cex.axis=0.85, lend=2, mgp=c(0.5, 0.5, 0), line=0)
  axis(2, tck=-0.02, line=0, cex.axis=0.85, lend=2, mgp=c(0.5, 0.5, 0), las=1)
  segments(x0=0, y0=0, x1=max(e), y1=max(e), col="firebrick2", lwd=1, lend=2)
  
  title(ylab="Observed", xlab="", sub="", mgp=c(2.5, 0.5, 0), cex.lab=0.85)
  title(ylab=expression('-log'[10]*' ('*italic(p)*')'), xlab="", sub="", mgp=c(1.5, 0.5, 0), cex.lab=0.85)
  title(xlab="Expected", ylab="", sub="", mgp=c(1.75, 0.5, 0), cex.lab=0.85)
  title(xlab=expression('-log'[10]*' ('*italic(p)*')'), ylab="", sub="", mgp=c(2.75, 0.5, 0), cex.lab=0.85)
  
}