EGI <- function(T, model, method=NULL, method.param=NULL,
            fun, iter, lower, upper, new.noise.var=0,
            optimcontrol=NULL, kmcontrol=NULL,integcontrol=NULL,...) {

  if (is.null(method)) method <- "ranjan"
  if (method == "timse" || method == "imse" || method == "sur" || method == "jn" || method == "vorob" || method == "vorobCons" || method == "vorobVol") return(EGIparallel(T=T,model=model,method = method,method.param = method.param,fun = fun,iter = iter,
                                                                lower=lower,upper=upper,batchsize = 1,new.noise.var = new.noise.var,
                                                                optimcontrol = optimcontrol,integcontrol = integcontrol,kmcontrol = kmcontrol,...))

	n <- nrow(model@X); d <- model@d
	if (is.null(kmcontrol$penalty)) kmcontrol$penalty <- model@penalty
	if (length(model@penalty) == 0) kmcontrol$penalty <- NULL
	if (is.null(kmcontrol$optim.method)) kmcontrol$optim.method <- model@optim.method
	if (is.null(kmcontrol$parinit)) kmcontrol$parinit <- model@parinit
	if (is.null(kmcontrol$control)) kmcontrol$control <- model@control
	if (is.null(kmcontrol$CovReEstimate)) kmcontrol$CovReEstimate <- model@param.estim

	if(length(T)>1 & method!= "tmse"){
	  print("Multiple thresholds are supported only for the following criteria: tmse,timse,sur,jn")
	  return(0)
	}

	for (i in 1:iter) {
		oEGI <- max_infill_criterion(lower=lower, upper=upper,T=T, method=method,
								method.param=method.param,model=model,optimcontrol=optimcontrol)

	  print("New point"); print(oEGI$par)
	  X.new <- oEGI$par; y.new <- fun(oEGI$par,...)
	  X.new <- as.numeric(X.new); X.new <- matrix(X.new,nrow=1,ncol=d)

	  #model <- update_km(model=model,NewX=X.new,NewY=y.new,CovReEstimate=kmcontrol$CovReEstimate,new.noise.var=new.noise.var,kmcontrol=kmcontrol)
	  model <-  update(object = model,newX=X.new,newy=y.new,cov.reestim  =kmcontrol$CovReEstimate,newnoise.var= new.noise.var,kmcontrol=kmcontrol)
	}

	return(list(
				par=model@X[(n+1):(n+iter),, drop=FALSE],
				value=model@y[(n+1):(n+iter),, drop=FALSE],
				npoints=1,
				nsteps=iter,
				lastmodel=model,
				lastvalue=oEGI$value,
				allvalues=oEGI$allvalues,
				variance.volume=oEGI$variance.volume
				))
}
