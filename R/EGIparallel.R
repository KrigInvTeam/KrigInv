EGIparallel <- function(T, model, method=NULL, method.param=NULL,
                        fun, iter, batchsize = 1, lower, upper, new.noise.var=0,
                        optimcontrol=NULL, kmcontrol=NULL,integcontrol=NULL,...) {

  if (is.null(method)) method <- "sur"
  if (method == "tmse" || method == "ranjan" || method == "bichon" || method == "tsee"){
    if(batchsize > 1) print("For this criterion, batchsize needs to be set to 1. Switching to bachsize = 1...")
    if(!is.numeric(method.param)) method.param=NULL #print("For this criterion, method.param is a scalar. Switching to default")
    return(EGI(T = T,model = model,method=method,method.param = method.param,fun = fun,iter = iter,lower = lower,upper = upper,new.noise.var = new.noise.var,optimcontrol = optimcontrol,kmcontrol = kmcontrol,integcontrol = integcontrol,...))
  }

  n <- nrow(model@X); d <- model@d
  if (is.null(kmcontrol$penalty)) kmcontrol$penalty <- model@penalty
  if (length(model@penalty) == 0) kmcontrol$penalty <- NULL
  if (is.null(kmcontrol$optim.method)) kmcontrol$optim.method <- model@optim.method
  if (is.null(kmcontrol$parinit)) kmcontrol$parinit <- model@parinit
  if (is.null(kmcontrol$control)) kmcontrol$control <- model@control
  if (is.null(kmcontrol$CovReEstimate)) kmcontrol$CovReEstimate <- model@param.estim

  if (is.null(optimcontrol$optim.option)) optimcontrol$optim.option <- 2

  # methods vorobCons and vorobVol require a call to the conservativeEstimate function of the anMC package, 
  # with a prediction on a fine design.
  # Here we initialize the parameters for such call
  if(method=="vorobCons" || method=="vorobVol"){
    if(!is.list(method.param)) method.param <- list(penalization=1, typeEx=">")
    if(is.null(method.param$consLevel)) method.param$consLevel = 0.95
    if(is.null(method.param$n_discrete_design)) method.param$n_discrete_design=500*model@d
    if(is.null(method.param$design)) method.param$design=as.matrix (sobol (
                n = method.param$n_discrete_design,dim = model@d))

    colnames(method.param$design) <- colnames(model@X)
    current.pred <- current.CE <- list()
    allCE_lvs <- rep(NA,iter+1)
  }

  for (i in 1:iter) {
    if (method == "sur" || method == "jn"){

      if(length(T) > 1) T <- sort(T)
      integration.param <- integration_design(integcontrol,d,lower,upper,model,T)
      real.volume.variance <- (method == "jn")

      oEGI <- max_sur_parallel(lower=lower, upper=upper,optimcontrol=optimcontrol,T=T,model=model,
                               integration.param=integration.param,batchsize=batchsize,
                               new.noise.var=new.noise.var,real.volume.variance=real.volume.variance)
    }
    else if(method == "timse" || method == "imse"){

      if(length(T) > 1) T <- sort(T)
      integration.param <- integration_design(integcontrol,d,lower,upper,model,T)
      if(!is.numeric(method.param)) method.param=NULL

      oEGI <- max_timse_parallel(lower=lower, upper=upper,optimcontrol=optimcontrol,T=T,model=model,
                                 integration.param=integration.param,batchsize=batchsize,
                                 new.noise.var=new.noise.var,epsilon=method.param,imse=(method == "imse"))
    }
    else if(method == "vorob"){

      if(length(T) > 1){
        print("vorob criterion not available with multiple thresholds.")
        return(0)
      }
      if(!is.list(method.param)) method.param <- list(penalization=1, typeEx=">")

      integration.param <- integration_design(integcontrol,d,lower,upper,model,T)

      oEGI <- max_vorob_parallel(lower=lower, upper=upper,optimcontrol=optimcontrol,T=T,model=model,
                                 integration.param=integration.param,batchsize=batchsize,
                                 new.noise.var=new.noise.var,penalisation = method.param$penalization,
                                 typeEx = method.param$typeEx)
    }
    else if(method=="vorobVol" || method=="vorobCons"){
      if(length(T) > 1){
        print("Conservative Vorobev criteria not available with multiple thresholds.")
        return(0)
      }
      if(!is.list(method.param)) method.param <- list(penalization=1, typeEx=">")
      integration.param <- integration_design(integcontrol,d,lower,upper,model,T)

      if(is.null(method.param$pred) || i>1){
        current.pred <- predict.km(object = model,
                                  newdata = method.param$design,
                                  type = "UK",cov.compute = TRUE)
        current.pred$cov <- current.pred$cov +1e-7*diag(nrow = nrow(current.pred$cov),ncol = ncol(current.pred$cov))
      }else{
        current.pred <- method.param$pred
      }

      if(is.null(method.param$consVorbLevel)){
        current.CE <- conservativeEstimate(alpha = method.param$consLevel, pred=current.pred,
                                          design=method.param$design, threshold=T, pn = NULL,
                                          type = method.param$typeEx, verb = 1,
                                          lightReturn = TRUE, algo = "GANMC")
        integration.param$alpha <- current.CE$lvs
      }else{
        integration.param$alpha <- method.param$consVorbLevel
      }
      allCE_lvs[i] <- integration.param$alpha

      if(method=="vorobCons"){
        oEGI <- max_vorob_parallel(lower=lower, upper=upper,optimcontrol=optimcontrol,T=T,model=model,
                                   integration.param=integration.param,batchsize=batchsize,
                                   new.noise.var=new.noise.var,penalisation = method.param$penalization,
                                   typeEx = method.param$typeEx)
      }
      else if(method=="vorobVol"){

        oEGI <- max_futureVol_parallel(lower=lower, upper=upper,optimcontrol=optimcontrol, batchsize=batchsize,
                               integration.param=integration.param, T=T, model=model,
                               new.noise.var =new.noise.var, typeEx = method.param$typeEx)
      }
    }

    print("New points"); print(oEGI$par)
    X.new <- oEGI$par; y.new <- rep(0,times=nrow(X.new))
    for (i in 1:nrow(X.new)) y.new[i] <- fun((oEGI$par)[i,],...)

    #model <- update_km(model=model,NewX=X.new,NewY=y.new,CovReEstimate=kmcontrol$CovReEstimate,new.noise.var=rep(new.noise.var,times=batchsize),kmcontrol=kmcontrol)
    model <- update(object = model,newX=X.new,newy=y.new,cov.reestim  =kmcontrol$CovReEstimate,newnoise.var= rep(new.noise.var,times=batchsize),kmcontrol=kmcontrol)
  }
  
  
  if(method=="vorobCons" || method=="vorobVol"){
    current.pred <- predict.km(object = model,
                            newdata = method.param$design,
                            type = "UK",cov.compute = TRUE)
    current.pred$cov <- current.pred$cov +1e-7*diag(nrow = nrow(current.pred$cov),ncol = ncol(current.pred$cov))
    current.CE <- conservativeEstimate(alpha = method.param$consLevel, pred=current.pred,
                                   design=method.param$design, threshold=T, pn = NULL,
                                   type = method.param$typeEx, verb = 1,
                                    lightReturn = TRUE, algo = "GANMC")
    allCE_lvs[iter+1] <- current.CE$lvs
    return(list(
      par=model@X[(n+1):(n+iter*batchsize),, drop=FALSE],
      value=model@y[(n+1):(n+iter*batchsize),, drop=FALSE],
      npoints=batchsize,
      nsteps=iter*batchsize,
      lastmodel=model,
      lastvalue=oEGI$value,
      allvalues=oEGI$allvalues,
      current.CE=current.CE,
      allCE_lvs = allCE_lvs
    ))
  }
  
  return(list(
    par=model@X[(n+1):(n+iter*batchsize),, drop=FALSE],
    value=model@y[(n+1):(n+iter*batchsize),, drop=FALSE],
    npoints=batchsize,
    nsteps=iter*batchsize,
    lastmodel=model,
    lastvalue=oEGI$value,
    allvalues=oEGI$allvalues
  ))
}
