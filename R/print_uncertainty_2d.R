print_uncertainty_2d <-
function(model,T,type="pn",lower=c(0,0),upper=c(1,1),
			resolution=200,new.points=0,
			xscale=c(0,1),yscale=c(0,1),show.points=TRUE,
			cex.contourlab=1,cex.points=1,cex.axis=1,pch.points.init=17,pch.points.end=17,
			col.points.init="black",col.points.end="red",nlevels=10,levels=NULL,xaxislab=NULL,
			yaxislab=NULL,xaxispoint=NULL,yaxispoint=NULL,krigmeanplot=FALSE,vorobmean=FALSE,consQuantile=NULL,...){

  cQcalcul<-FALSE
  if(!is.numeric(consQuantile) && !is.null(consQuantile$consLevel)){
    resolution=min(resolution,100)
    cQcalcul<-TRUE
  }

  alpha <- NULL
	s1 <- seq(from=lower[1],to=upper[1],length=resolution)
	s2 <- seq(from=lower[2],to=upper[2],length=resolution)

	img_design <- expand.grid(s1,s2)

  pred <- predict_nobias_km(object=model,newdata=img_design,type="UK",cov.compute = cQcalcul)

	pn <- excursion_probability(mn = pred$mean,sn = pred$sd,T = T)

	if(type=="pn"){
		#we print pn(x)
		myvect <- pn
	}else if(type=="sur"){
		myvect <- pn * (1-pn)
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
	}
	else if(type=="imse"){
		myvect <- (pred$sd)^2
	}else if(type=="vorob"){
	  alpha <- vorob_threshold(pn=pn)
	  pn_bigger_than_alpha <- (pn>alpha)+0
	  pn_lower_than_alpha <- 1-pn_bigger_than_alpha
	  myvect <- pn*pn_lower_than_alpha + (1-pn)*pn_bigger_than_alpha
	}else{
		print("unknown value for the type argument, we take type=pn")
		myvect <- pn
	}
	res <- mean(myvect)

	mymatrix <- matrix(myvect, nrow=resolution,ncol=resolution)
	krigmean <- matrix(pred$mean, nrow=resolution,ncol=resolution)
  if(krigmeanplot) contourmatrix <- krigmean
  if(!krigmeanplot) contourmatrix <- mymatrix

	scale.x <- xscale[1] + s1 * (xscale[2]-xscale[1])
	scale.y <- yscale[1] + s2 * (yscale[2]-yscale[1])

	axes <- is.null(xaxislab)
	image(x=scale.x,y=scale.y,z=mymatrix,col=grey.colors(10),cex.axis=cex.axis,axes=axes,...)

	if(!is.null(xaxislab)){
		axis(side=1,at=xaxispoint,labels=xaxislab,cex.axis=cex.axis)
		axis(side=2,at=yaxispoint,labels=yaxislab,cex.axis=cex.axis)
	}
	if(!is.null(levels)){
		contour(x=scale.x,y=scale.y,z=contourmatrix,add=TRUE,labcex=cex.contourlab,levels=levels)
	} else {
		contour(x=scale.x,y=scale.y,z=contourmatrix,add=TRUE,labcex=cex.contourlab,nlevels=nlevels)
	}
  if(vorobmean){
    if(is.null(alpha)) alpha <- vorob_threshold(pn=pn)
    pnmatrix <- matrix(pn, nrow=resolution,ncol=resolution)
    contour(x=scale.x,y=scale.y,z=pnmatrix,add=TRUE,levels=c(alpha),lwd=5,col="blue",labcex=0.01)
  }

	if(!is.null(consQuantile)){
	  pnmatrix <- matrix(pn, nrow=resolution,ncol=resolution)
	  if(is.numeric(consQuantile)){
	    contour(x=scale.x,y=scale.y,z=pnmatrix,add=TRUE,levels=c(consQuantile),lwd=5,col="green",labcex=0.01)
	  }else{
	    if(!is.null(consQuantile$consLevel) && length(T) > 1){
	      print("Conservative estimates are not available with multiple thresholds.")
	    }else if(!is.null(consQuantile$consLevel)){
	      pred$cov <- pred$cov +1e-7*diag(nrow = nrow(pred$cov),ncol = ncol(pred$cov))
	      # conservativeEstimate expects pn = P(Y < T) for excursions below T;
	      # its level lvs refers to that probability, so on the pn = P(Y > T)
	      # scale used here the level is 1 - lvs.
	      ceType <- if (is.null(consQuantile$typeEx)) ">" else consQuantile$typeEx
	      current.CE <- conservativeEstimate(alpha = consQuantile$consLevel, pred=pred,
	                                         design=img_design, threshold=T, pn = if (ceType == ">") pn else 1 - pn,
	                                         type = ceType, verb = 0,
	                                         lightReturn = TRUE, algo = "GANMC")
	      ceLevel <- if (ceType == ">") current.CE$lvs else 1 - current.CE$lvs
	      contour(x=scale.x,y=scale.y,z=pnmatrix,add=TRUE,levels=c(ceLevel),lwd=5,col="green",labcex=0.01)
	      res <- c(res, current.CE$lvs)
	      names(res) <- c("integrated pn","CEquantile")
	    }
	  }
	}

	if(show.points){
		initialpoints <- model@X[1:(model@n-new.points),]
		initialpoints[,1]<- xscale[1] + initialpoints[,1] * (xscale[2]-xscale[1])
		initialpoints[,2]<- yscale[1] + initialpoints[,2] * (yscale[2]-yscale[1])

		points(initialpoints, col=col.points.init, pch=pch.points.init, lwd=4,cex=cex.points)
		if (new.points!=0){

			finalpoints <- model@X[(model@n-new.points+1):model@n,]
			if (new.points==1) {
				finalpoints[1]<- xscale[1] + finalpoints[1] * (xscale[2]-xscale[1])
				finalpoints[2]<- yscale[1] + finalpoints[2] * (yscale[2]-yscale[1])
			}else{
				finalpoints[,1]<- xscale[1] + finalpoints[,1] * (xscale[2]-xscale[1])
				finalpoints[,2]<- yscale[1] + finalpoints[,2] * (yscale[2]-yscale[1])
			}

			if (new.points==1) {points(t(finalpoints), col=col.points.end, pch=pch.points.end, lwd=4,cex=cex.points)
			} else points(finalpoints, col=col.points.end, pch=pch.points.end, lwd=4,cex=cex.points)
		}
	}

	return(res)
}
