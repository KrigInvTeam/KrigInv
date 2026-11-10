vorob_threshold <- function(pn){
  if(is.null(pn)) {return(0.5)}
  
  M <- length(pn)
  integ.pn <- mean(pn)
  pn.sort <- sort(pn)
  index <- (1-integ.pn)*M
  # index < 1 when mean(pn) > 1 - 1/M: the threshold is then the smallest pn
  index.floor <- max(1, floor(index))
  index.cap <- min(M,index.floor+1)
  wbig <- max(0, index - index.floor)
  wsmall <- 1- wbig
  alpha <- pn.sort[index.floor]*wsmall + pn.sort[index.cap]*wbig
  return(alpha)
}
