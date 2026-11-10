\name{max_timse_parallel}
\Rdversion{1.4.1}
\alias{max_timse_parallel}

\title{ Minimizer of the parallel timse criterion}

\description{
Minimization, based on the package rgenoud (or on exhaustive search on a discrete set), of the timse criterion for a batch of candidate sampling points.  
}

\usage{
max_timse_parallel(lower, upper, optimcontrol = NULL, 
batchsize, integration.param, T,
model, new.noise.var = 0,
epsilon=0,imse=FALSE)
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
Array containing one or several thresholds.
}
  \item{model}{
A Kriging model of \code{\link[DiceKriging]{km}} class.
}
  \item{new.noise.var}{
Optional scalar value of the noise variance of the new observations.
}
  \item{epsilon}{
Optional tolerance value (a real positive number). Default value is 0.
}
  \item{imse}{
Optional boolean to decide if the \code{"imse"} criterion should be used instead of \code{"timse"}. default: \code{FALSE}.
}
}

\value{A list with components:
\item{par}{the best set of parameters found.}
\item{value}{the value of the sur criterion at par.}
\item{allvalues}{If an optimization on a discrete set of points is chosen, the value of the criterion at all these points.}
}

\references{
Picheny V., Ginsbourger D., Roustant O., Haftka R.T., (2010) \emph{Adaptive designs of experiments for accurate approximation of a target region}, J. Mech. Des. vol. 132(7). \doi{10.1115/1.4001873}

Picheny V. (2009) \emph{Improving accuracy and compensating for uncertainty in surrogate modeling}, Ph.D. thesis, University of Florida and Ecole Nationale Superieure des Mines de Saint-Etienne

Chevalier C., Bect J., Ginsbourger D., Vazquez E., Picheny V., Richet Y. (2014), \emph{Fast parallel kriging-based stepwise uncertainty reduction with application to the identification of an excursion set}, Technometrics, vol. 56(4), pp 455-465. \doi{10.1080/00401706.2013.860918}
}

\author{
Victor Picheny (INRA, Toulouse, France)

Clement Chevalier (University of Neuchatel, Switzerland)
}

\seealso{\code{\link{EGIparallel}},\code{\link{max_sur_parallel}}}

\examples{
#max_timse_parallel


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

optimcontrol <- list(method="genoud",pop.size=200,optim.option=2)
integcontrol <- list(distrib="timse",n.points=400,init.distrib="MC")
integration.param <- integration_design(integcontrol=integcontrol,d=2,
                                            lower=lower,upper=upper,model=model,
                                            T=T)
batchsize <- 5 #number of new points

\donttest{
obj <- max_timse_parallel(lower=lower,upper=upper,optimcontrol=optimcontrol,
                        batchsize=batchsize,T=T,model=model,
                        integration.param=integration.param,epsilon=0,imse=FALSE)
                        #5 optims in dimension 2 !

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
