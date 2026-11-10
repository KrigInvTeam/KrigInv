tmse_optim <- function(x, model, T, method.param=NULL){

	y <- t(x)
	if(ncol(y) == model@d) z <- y
	if(ncol(y) != model@d) z <- x
	
	krig <- predict_nobias_km(object=model,newdata=as.data.frame(z),type="UK",se.compute=TRUE)
	mk    <- krig$mean
	sk    <- krig$sd
	
	if(is.null(method.param)) method.param <- 0
	epsilon <- method.param
	
	if(length(T)==1){
	  W <- 1/sqrt(2*pi*(sk^2+epsilon^2))*exp(-0.5*((mk-T)/sqrt(sk^2+epsilon^2))^2)
	}else{
	  W0 <- 1/sqrt(2*pi*(sk^2+epsilon^2))
	  W <- 0
	  for(i in 1:length(T)){
	    Ti <- T[i]
      W <- W + W0 *  exp(-0.5*((mk-Ti)/sqrt(sk^2+epsilon^2))^2)
	  }
	}
	
	tmse <- W*sk^2
	tmse[is.nan(tmse)] <- 0 # sk == 0 gives Inf * 0; the limit of W*sk^2 as sk -> 0 is 0
	return(tmse)
}
