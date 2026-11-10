\name{max_sur_parallel}
\Rdversion{1.4.1}
\alias{max_sur_parallel}

\title{ Minimizer of the parallel \code{"sur"} or \code{"jn"} criterion}

\description{
Minimization, based on the package rgenoud (or on exhaustive search on a discrete set), of the \code{"sur"} or \code{"jn"} criterion for a batch of candidate sampling points.  
}

\usage{
max_sur_parallel(lower, upper, optimcontrol = NULL, 
batchsize, integration.param, T,
model, new.noise.var = 0,real.volume.variance=FALSE)
}

\arguments{
  \item{lower}{
Vector containing the lower bounds of the design space.
}
  \item{upper}{
Vector containing the upper bounds of the design space.
}
  \item{optimcontrol}{
Optional list of control parameters for the optimization of the sampling criterion. The field \code{method} defines which optimization method is used: it can be either \code{"genoud"} (default) for an optimisation using the genoud algorithm, or \code{"discrete"} for an optimisation over a specified discrete set. 
If the field \code{method} is set to \code{"genoud"}, one can set some parameters of this algorithm: 
\code{pop.size}  (default : 50*d),  \code{max.generations} (10*d), \code{wait.generations} (2),  \code{BFGSburnin} (2) and the mutations \code{P1}, \code{P2}, up to \code{P9} (see \code{\link[rgenoud]{genoud}}). Numbers into brackets are the default values.
If the field \code{method} is set to \code{"discrete"}, one can set the field \code{optim.points}: p * d matrix corresponding to the p points where the criterion will be evaluated. If nothing is specified, 100*d points are chosen randomly. 
Finally, one can control the field \code{optim.option} in order to decide how to optimize the sampling criterion. 
If \code{optim.option} is set to 2 (default), batchsize sequential optimizations in dimension d are performed to find the optimum.
If \code{optim.option} is set to 1, only one optimization in dimension batchsize*d is performed. This option is only available with \code{"genoud"}. This option might provide more global and accurate solutions, but is a lot more expensive.
}
  \item{batchsize}{
Number of points to sample simultaneously. The sampling criterion will return batchsize points at a time for sampling.
}
  \item{integration.param}{
Optional list of control parameter for the computation of integrals, containing the fields \code{integration.points}: a p*d matrix corresponding to p integrations points and \code{integration.weights}: a vector of size p corresponding to the weights of these integration points. 
If nothing is specified, default values are used (see: function \code{\link{integration_design}} for more details). 
}
  \item{T}{
Target value (scalar).
}
  \item{model}{
A Kriging model of \code{\link[DiceKriging]{km}} class.
}
  \item{new.noise.var}{
Optional scalar value of the noise variance of the new observations.
}
\item{real.volume.variance}{Optional argument to use the \code{"jn"} criterion instead of \code{"sur"}. Default: FALSE}
}

\value{A list with components:
\item{par}{the best set of points found.}
\item{value}{the value of the sur criterion at par.}
\item{allvalues}{If an optimization on a discrete set of points is chosen, the value of the criterion at all these points.}
}

\references{
Chevalier C., Bect J., Ginsbourger D., Vazquez E., Picheny V., Richet Y. (2014), \emph{Fast parallel kriging-based stepwise uncertainty reduction with application to the identification of an excursion set}, Technometrics, vol. 56(4), pp 455-465. \doi{10.1080/00401706.2013.860918}

Chevalier C., Ginsbourger D. (2014), \emph{Corrected Kriging update formulae for batch-sequential data assimilation}, in Pardo-Iguzquiza, E., et al. (Eds.) Mathematics of Planet Earth, pp 119-122. \doi{10.1007/978-3-642-32408-6_29}
}

\author{
Clement Chevalier (University of Neuchatel, Switzerland)
}

\seealso{\code{\link{EGIparallel}},\code{\link{sur_optim_parallel}},\code{\link{jn_optim_parallel}}}

\examples{
#max_sur_parallel


set.seed(9)
N <- 20 #number of observations
T <- c(40,80) #thresholds
testfun <- branin
lower <- c(0,0)
upper <- c(1,1)

#a 20 points initial design
design <- data.frame( matrix(runif(2*N),ncol=2) )
response <- testfun(design)

#km object with matern3_2 covariance
#params estimated by ML from the observations
model <- km(formula=~., design = design, 
	response = response,covtype="matern3_2")

optimcontrol <- list(method="genoud",pop.size=50,optim.option=1)
integcontrol <- list(distrib="sur",n.points=50,init.distrib="MC")
integration.param <- integration_design(integcontrol=integcontrol,d=2,
                                            lower=lower,upper=upper,model=model,
                                            T=T)
batchsize <- 5 #number of new points

\donttest{
obj <- max_sur_parallel(lower=lower,upper=upper,optimcontrol=optimcontrol,
                        batchsize=batchsize,T=T,model=model,
                        integration.param=integration.param)
                        #one (hard) optim in dimension 5*2 !

obj$par;obj$value #optimum in 5 new points
new.model <- update(object=model,newX=obj$par,newy=apply(obj$par,1,testfun),
                       cov.reestim=TRUE)


par(mfrow=c(1,2))
print_uncertainty(model=model,T=T,type="pn",lower=lower,upper=upper,
cex.points=2.5,main="probability of excursion")

print_uncertainty(model=new.model,T=T,type="pn",lower=lower,upper=upper,
new.points=batchsize,col.points.end="red",cex.points=2.5,
main="updated probability of excursion")
}
}
