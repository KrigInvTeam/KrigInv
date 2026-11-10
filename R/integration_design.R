
integration_design <- function(integcontrol=NULL,d=NULL,lower,upper,model=NULL,T=NULL,min.prob=0.001){
	
	#Generic function to build integration points for some criterion
	#available important sampling schemes: for SUR (and Jn) AND tIMSE
	result<-NULL
	
	if(is.null(d)) d <- length(lower)
	if (length(lower) != length(upper) ){
		print("Error in integration_design: 'lower' and 'upper' must have the same length")
		return(NULL)
	}
	
	#Trivial case 1
	if(is.null(integcontrol)){
		#nothing has been specified, thus we use default values
		n.int.points<-d*100
		integration.points <- t(lower+t(sobol(n=n.int.points,dim=d))*(upper-lower))
    if(d==1) integration.points <- matrix(integration.points,ncol=1)
    if(!is.null(model)) colnames(integration.points)<- colnames(model@X)
		
		result$integration.points <- integration.points
		result$integration.weights <- NULL
		return(result)
	}
	
	#Trivial case 2
	if(!is.null(integcontrol$integration.points)){
		#integration points pre-specified
		if(!is.null(model)) colnames(integcontrol$integration.points) <- colnames(model@X)
		result$integration.points <- integcontrol$integration.points
		result$integration.weights <- integcontrol$integration.weights
		return(result)
	}
	
	#non trivial cases:
	if(is.null(integcontrol$n.points)) integcontrol$n.points <- d*100
	if(is.null(integcontrol$distrib)) integcontrol$distrib <- "sobol"
	
	if(integcontrol$distrib=="sobol"){
		integration.points <- t(lower+t(sobol(n=integcontrol$n.points,dim=d))*(upper-lower))
    if(d==1) integration.points <- matrix(integration.points,ncol=1)
		if(!is.null(model)) colnames(integration.points)<- colnames(model@X)
		result$integration.points <- integration.points
		result$integration.weights<-NULL
		return(result)
	}
	
	if(integcontrol$distrib=="MC"){
		integration.points <- t(lower+t(matrix(runif(d*integcontrol$n.points),ncol=d))*(upper-lower))
    if(d==1) integration.points <- matrix(integration.points,ncol=1)
		if(!is.null(model)) colnames(integration.points)<- colnames(model@X)
		result$integration.points <- integration.points
		result$integration.weights<-NULL
		return(result)
	}
	
	if(integcontrol$distrib=="sur"){
		if(is.null(integcontrol$n.candidates)) integcontrol$n.candidates <- integcontrol$n.points*10
		if(is.null(integcontrol$init.distrib)) integcontrol$init.distrib <- "sobol"
		
		#generation of the initial candidates points:
		if(integcontrol$init.distrib=="sobol") initial.integration.points <- t(lower+t(sobol(n=integcontrol$n.candidates,dim=d))*(upper-lower))
		if(integcontrol$init.distrib=="MC") initial.integration.points <- t(lower+t(matrix(runif(d*integcontrol$n.candidates),ncol=d))*(upper-lower))
		if(integcontrol$init.distrib=="spec") initial.integration.points <- integcontrol$init.distrib.spec
		
		if(d==1) initial.integration.points <- matrix(initial.integration.points,ncol=1)
		
		#prediction on these initial candidate points
		if(is.null(model)){
			print("Error in integration_design: for 'sur', 'jn', 'imse' or 'timse' importance sampling distribution you must set the argument 'model'")
			return(NULL)
		}
		predictions <- predict_nobias_km(object=model,newdata=initial.integration.points,type="UK")
		
		pn <- excursion_probability(mn = predictions$mean,sn = predictions$sd,T = T)
		Tau.n <- pn*(1-pn)
		
		Tau.n.sum <- sum(Tau.n)
    if(Tau.n.sum==0) Tau.n.sum <- 1
		prob.n <- pmax(Tau.n/Tau.n.sum,min.prob/integcontrol$n.candidates)
		prob.n <- prob.n/sum(prob.n)
    weight.n <- 1/(prob.n*integcontrol$n.candidates*integcontrol$n.points)
    
    prob.n.copy <-c(0,prob.n)
    prob.n.cum <-cumsum(prob.n.copy)
    
		my.indices <- findInterval(runif(integcontrol$n.points),prob.n.cum,all.inside=TRUE)
		integration.points <- initial.integration.points[my.indices,]
		integration.weights <- weight.n[my.indices]
		
		if(d==1) integration.points <- matrix(integration.points,ncol=1)
		if(integcontrol$n.points==1) integration.points <- matrix(integration.points,ncol=d)
    
		if(!is.null(model)) colnames(integration.points)<- colnames(model@X)
		result$integration.points <- integration.points
		result$integration.weights <- integration.weights
		return(result)
	}
	if(integcontrol$distrib=="vorob"){
	  if(is.null(integcontrol$n.candidates)) integcontrol$n.candidates <- integcontrol$n.points*10
	  if(is.null(integcontrol$init.distrib)) integcontrol$init.distrib <- "sobol"
	  
	  #generation of the initial candidates points:
	  if(integcontrol$init.distrib=="sobol") initial.integration.points <- t(lower+t(sobol(n=integcontrol$n.candidates,dim=d))*(upper-lower))
	  if(integcontrol$init.distrib=="MC") initial.integration.points <- t(lower+t(matrix(runif(d*integcontrol$n.candidates),ncol=d))*(upper-lower))
	  if(integcontrol$init.distrib=="spec") initial.integration.points <- integcontrol$init.distrib.spec
	  
	  if(d==1) initial.integration.points<-matrix(initial.integration.points,ncol=1)
    
	  #prediction on these initial candidate points
	  if(is.null(model)){
	    print("Error in integration_design: for 'sur', 'jn', 'imse' or 'timse' importance sampling distribution you must set the argument 'model'")
	    return(NULL)
	  }
	  predictions <- predict_nobias_km(object=model,newdata=initial.integration.points,type="UK")
	  
	  pn <- excursion_probability(mn = predictions$mean,sn = predictions$sd,T = T)
	  alpha <- vorob_threshold(pn)

	  pn_bigger_than_alpha <- (pn>alpha)+0
	  pn_lower_than_alpha <- 1-pn_bigger_than_alpha
    
	  Tau.n <- pn*pn_lower_than_alpha+(1-pn)*pn_bigger_than_alpha
	  
	  Tau.n.sum <- sum(Tau.n)
	  if(Tau.n.sum==0) Tau.n.sum <- 1
	  prob.n <- pmax(Tau.n/Tau.n.sum,min.prob/integcontrol$n.candidates)
	  prob.n <- prob.n/sum(prob.n)
	  weight.n <- 1/(prob.n*integcontrol$n.candidates*integcontrol$n.points)
	  
	  prob.n.copy <-c(0,prob.n)
	  prob.n.cum <-cumsum(prob.n.copy)
	  
	  my.indices <- findInterval(runif(integcontrol$n.points),prob.n.cum,all.inside=TRUE)
	  integration.points <- initial.integration.points[my.indices,]
	  integration.weights <- weight.n[my.indices]
	  
	  if(d==1) integration.points <- matrix(integration.points,ncol=1)
	  if(integcontrol$n.points==1) integration.points <- matrix(integration.points,ncol=d)
	  
	  if(!is.null(model)) colnames(integration.points)<- colnames(model@X)
	  result$integration.points <- integration.points
	  result$integration.weights <- integration.weights
    result$alpha <- alpha
	  return(result)
	}
	if(integcontrol$distrib=="jn"){
	
		if(is.null(integcontrol$n.candidates)) integcontrol$n.candidates <- integcontrol$n.points*10
		if(is.null(integcontrol$init.distrib)) integcontrol$init.distrib <- "sobol"
		
		#generation of the initial candidates points:
		if(integcontrol$init.distrib=="sobol") initial.integration.points <- t(lower+t(sobol(n=integcontrol$n.candidates,dim=d))*(upper-lower))
		if(integcontrol$init.distrib=="MC") initial.integration.points <- t(lower+t(matrix(runif(d*integcontrol$n.candidates),ncol=d))*(upper-lower))
		if(integcontrol$init.distrib=="spec") initial.integration.points <- integcontrol$init.distrib.spec
		
		if(d==1) initial.integration.points<-matrix(initial.integration.points,ncol=1)
		
		if(is.null(model)){
			print("Error in integration_design: for 'sur', 'jn', 'imse' or 'timse' importance sampling distribution you must set the argument 'model'")
			return(NULL)
		}		
		predictions <- predict_nobias_km(object=model,newdata = initial.integration.points,type="UK",cov.compute=FALSE)
		
    M <- integcontrol$n.candidates
    pn <- excursion_probability(mn = predictions$mean,sn = predictions$sd,T = T)
		pn.sum <- sum(pn)
		prob.n <- pmax(pn/pn.sum,min.prob/M)
		prob.n <- prob.n/sum(prob.n)
    weight.n <- 1/(prob.n*M*integcontrol$n.points)
    		     
    prob.n.copy <-c(0,prob.n)
    prob.n.cum <-cumsum(prob.n.copy)
    my.indices.1 <- findInterval(runif(integcontrol$n.points),prob.n.cum,all.inside=TRUE)
    my.indices.2 <- findInterval(runif(integcontrol$n.points),prob.n.cum,all.inside=TRUE)
    
		integration.points.1 <- initial.integration.points[my.indices.1,]
		integration.points.2 <- initial.integration.points[my.indices.2,]
		
		integration.points <- rbind(integration.points.1,integration.points.2)
		integration.weights <- weight.n[my.indices.1]*weight.n[my.indices.2]*integcontrol$n.points
				
		if(d==1) integration.points <- matrix(integration.points,ncol=1)
		if(!is.null(model)) colnames(integration.points)<- colnames(model@X)
		result$integration.points <- integration.points
		result$integration.weights <- integration.weights
		
		return(result)
	}
	if(integcontrol$distrib=="timse" || integcontrol$distrib=="imse"){
		if(is.null(integcontrol$n.candidates)) integcontrol$n.candidates <- integcontrol$n.points*10
		if(is.null(integcontrol$init.distrib)) integcontrol$init.distrib <- "sobol"
		
		if(integcontrol$init.distrib=="sobol") initial.integration.points <- t(lower+t(sobol(n=integcontrol$n.candidates,dim=d))*(upper-lower))
		if(integcontrol$init.distrib=="MC") initial.integration.points <- t(lower+t(matrix(runif(d*integcontrol$n.candidates),ncol=d))*(upper-lower))
		if(integcontrol$init.distrib=="spec") initial.integration.points <- integcontrol$init.distrib.spec
		if(d==1) initial.integration.points <- matrix(initial.integration.points,ncol=1)
		
		#prediction on these initial candidate points
		if(is.null(model)){
			print("Error in integration_design: for 'sur', 'jn', 'imse' or 'timse' importance sampling distribution you must set the argument 'model'")
			return(NULL)
		}
		
		predictions <- predict_nobias_km(object=model,newdata=initial.integration.points,type="UK")
		
		mk <- predictions$mean
		sk <- predictions$sd
		
		if(length(T)==1){
		  weight <- 1/sqrt(2*pi*sk^2) * exp(-0.5*((mk-T)/sk)^2)
		}else{
		  weight0 <- 1/sqrt(2*pi*sk^2)
		  weight <- 0
		  for(i in 1:length(T)){
		    Ti <- T[i]
		    weight <- weight + weight0 * exp(-0.5*((mk-Ti)/sk)^2)
		  }
		}
		
		weight[is.nan(weight)] <- 0
		
		if(integcontrol$distrib=="timse") timse <- weight * sk*sk
		if(integcontrol$distrib=="imse") timse <- sk*sk
		
		timse.sum <- sum(timse)
    if(timse.sum==0) timse.sum <- 1
		prob.n <- pmax(timse/timse.sum,min.prob/integcontrol$n.candidates)
		prob.n <- prob.n/sum(prob.n)
    weight.n <- 1/(prob.n*integcontrol$n.candidates*integcontrol$n.points)
		
		prob.n.copy <-c(0,prob.n)
		prob.n.cum <- cumsum(prob.n.copy)
    
		my.indices <- findInterval(runif(integcontrol$n.points),prob.n.cum,all.inside=TRUE)
		integration.points <- initial.integration.points[my.indices,]
		integration.weights <- weight.n[my.indices]
		
		if(d==1) integration.points <- matrix(integration.points,ncol=1)
		if(integcontrol$n.points==1) integration.points <- matrix(integration.points,ncol=d)
    
		if(!is.null(model)) colnames(integration.points)<- colnames(model@X)
		result$integration.points <- integration.points
		result$integration.weights <- integration.weights
		return(result)
	}
	
}

