
print_uncertainty_1d <- function(model,T,type="pn",
      lower=0,upper=1,resolution=500,new.points=0,
			xscale=c(0,1),show.points=TRUE,
			cex.points=1,cex.axis=1,pch.points.init=17,pch.points.end=17,
			col.points.init="black",col.points.end="red",
			xaxislab=NULL,yaxislab=NULL,xaxispoint=NULL,yaxispoint=NULL,
      vorobmean=FALSE,krigmeanplot=FALSE,Tplot=FALSE,consQuantile=NULL,...){

	n <- model@n
	initSize <- n -new.points
	obs.X <- as.numeric(model@X)
	alpha <- NULL

	s <- matrix(seq(from=lower,to=upper,length=resolution),ncol=1)

	pred <- predict_nobias_km(object=model,newdata=s,type="UK")
	pn <- excursion_probability(mn = pred$mean,sn = pred$sd,T = T)

	pred_obs <- predict_nobias_km(object=model,newdata=model@X,type="UK")
	pn_obs <- excursion_probability(mn = pred_obs$mean,sn = pred_obs$sd,T = T)

	if(type=="pn"){
		myvect <- pn
		obs.Y <- pn_obs
	}else if(type=="sur"){
		myvect <- pn * (1-pn)
		obs.Y <- rep(0,times=n)
	}else if(type=="timse"){
	  if(length(T)==1){
	    Wn <- 1/(sqrt(2*pi)*pred$sd)*exp(-1/2*((pred$mean-T)/pred$sd)^2)
	    myvect <- (pred$sd)^2 * Wn
	  }else{
	    Wn0 <- 1/(sqrt(2*pi)*pred$sd)
	    Wn <- 0
	    for(i in 1:length(T)){
	      Ti <- T[i]
	      Wn <- Wn + Wn0 * exp(-1/2*((pred$mean-Ti)/pred$sd)^2)
	    }
	    myvect <- (pred$sd)^2 * Wn
	  }
		obs.Y <- rep(0,times=n)
	}else if(type=="imse"){
		myvect <- (pred$sd)^2
		obs.Y <- rep(0,times=n)
	}else if(type=="vorob"){
    alpha <- vorob_threshold(pn=pn)
    pn_bigger_than_alpha <- (pn>alpha)+0
    pn_lower_than_alpha <- 1-pn_bigger_than_alpha
    myvect <- pn*pn_lower_than_alpha + (1-pn)*pn_bigger_than_alpha
		obs.Y <- rep(0,times=n)
	}else{
		print("unknown 'type' value, we set type = pn")
		myvect <- pn
		obs.Y <- pn_obs
	}
	res <- mean(myvect)


	scale.x <- xscale[1] + seq(from=lower,to=upper,length=resolution) * (xscale[2]-xscale[1])
	obs.X <- xscale[1] + obs.X * (xscale[2]-xscale[1])
	axes <- is.null(xaxislab)

	plot(x = scale.x, y = myvect, type = "l",axes=axes,...)

	if(!is.null(xaxislab)) {
		axis(side=1,at=xaxispoint,labels=xaxislab,cex.axis=cex.axis)
		axis(side=2,at=yaxispoint,labels=yaxislab,cex.axis=cex.axis)
	}

	if(show.points){
		points(x=obs.X[1:initSize],y=obs.Y[1:initSize], col=col.points.init, pch=pch.points.init,cex=cex.points)

		if (new.points!=0){
			indices <- c((initSize+1):n)
			points(x=obs.X[indices],y=obs.Y[indices], col=col.points.end, pch=pch.points.end,cex=cex.points)
		}
	}

	if(vorobmean){
		if(is.null(alpha)) alpha <- vorob_threshold(pn=pn)
		if(!is.numeric(alpha)) alpha <- 0.5
		lines(scale.x,rep(alpha,times=resolution),col="blue",lty=2,lwd=3)
	}

	if(!is.null(consQuantile)){
	  if(is.numeric(consQuantile)){
	    lines(scale.x,rep(consQuantile,times=resolution),col="green",lty=2,lwd=3)
	  }else{
	    if(!is.null(consQuantile$consLevel) && length(T) > 1){
	      print("Conservative estimates are not available with multiple thresholds.")
	    }else if(!is.null(consQuantile$consLevel)){
	      pred <- predict_nobias_km(object=model,newdata=s,type="UK",cov.compute = TRUE)
	      pred$cov <- pred$cov +1e-7*diag(nrow = nrow(pred$cov),ncol = ncol(pred$cov))
	      # conservativeEstimate expects pn = P(Y < T) for excursions below T;
	      # its level lvs refers to that probability, so on the pn = P(Y > T)
	      # scale used here the level is 1 - lvs.
	      ceType <- if (is.null(consQuantile$typeEx)) ">" else consQuantile$typeEx
	      current.CE <- conservativeEstimate(alpha = consQuantile$consLevel, pred=pred,
	                                         design=s, threshold=T, pn = if (ceType == ">") pn else 1 - pn,
	                                         type = ceType, verb = 0,
	                                         lightReturn = TRUE, algo = "GANMC")
	      ceLevel <- if (ceType == ">") current.CE$lvs else 1 - current.CE$lvs
	      lines(scale.x,rep(ceLevel,times=resolution),col="green",lty=2,lwd=3)
	      res <- c(res, current.CE$lvs)
	      names(res) <- c("integrated pn","CEquantile")
	    }
	  }

	}

	if(krigmeanplot){
	  par(new = TRUE)
	  plot(x=scale.x, y = pred$mean, type = "l", axes = FALSE, xlab = "", ylab = "",col="blue",lwd=2)
	  axis(side=4, at = pretty(range(pred$mean)))
	  if(Tplot){
	    abline(h = T,col="red",lwd=2)
	  }
	}

	return(res)

}
