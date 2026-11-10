\name{max_vorob_parallel}
\Rdversion{1.4.1}
\alias{max_vorob_parallel}

\title{ Minimizer of the parallel vorob criterion}

\description{
Minimization, based on the package rgenoud (or on exhaustive search on a discrete set), of the Vorob'ev criterion for a batch of candidate sampling points.
}

\usage{
max_vorob_parallel(lower, upper, optimcontrol = NULL,
batchsize, integration.param, T,
model, new.noise.var = 0,
penalisation = NULL, typeEx = ">")
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
Target value (scalar). The criterion CANNOT be used with multiple thresholds.
}
  \item{model}{
A Kriging model of \code{\link[DiceKriging]{km}} class.
}
  \item{new.noise.var}{
Optional scalar value of the noise variance of the new observations.
}
	\item{penalisation}{Optional penalization constant for type I errors. If equal to zero,  computes the Type II criterion.}
    \item{typeEx}{A character (">" or "<") identifying the type of excursion}
}

\value{A list with components:
\item{par}{the best set of parameters found.}
\item{value}{the value of the Vorob'ev criterion at par.}
\item{allvalues}{If an optimization on a discrete set of points is chosen, the value of the criterion at all these points.}
}

\references{
Chevalier C., Ginsbourger D., Bect J., Molchanov I. (2013) \emph{Estimating and quantifying uncertainties on level sets using the Vorob'ev expectation and deviation with gaussian process models} mODa 10, Advances in Model-Oriented Design and Analysis, Contributions to Statistics, pp 35-43. \doi{10.1007/978-3-319-00218-7_5}

Chevalier C. (2013) \emph{Fast uncertainty reduction strategies relying on Gaussian process models} Ph.D Thesis, University of Bern

Azzimonti, D., Ginsbourger, D., Chevalier, C., Bect, J., and Richet, Y. (2021). \emph{Adaptive design of experiments for conservative estimation of excursion sets}. Technometrics, 63(1):13-26. \doi{10.1080/00401706.2019.1693427}
}

\author{
Clement Chevalier (University of Neuchatel, Switzerland)

Dario Azzimonti (IDSIA, Switzerland)
}

\seealso{\code{\link{EGIparallel}},\code{\link{max_sur_parallel}}}

\examples{
#max_vorob_parallel

set.seed(9)
N <- 20 #number of observations
T <- 80 #threshold
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
obj <- max_vorob_parallel(lower=lower,upper=upper,optimcontrol=optimcontrol,
                        batchsize=batchsize,T=T,model=model,
                        integration.param=integration.param)
                        #5 optims in dimension 2 !

obj$par;obj$value #optimum in 5 new points
new.model <- update(object=model,newX=obj$par,newy=apply(obj$par,1,testfun),
                       cov.reestim=TRUE)

par(mfrow=c(1,2))
print_uncertainty(model=model,T=T,type="pn",lower=lower,upper=upper,vorobmean=TRUE,
cex.points=2.5,main="probability of excursion")

print_uncertainty(model=new.model,T=T,type="pn",lower=lower,upper=upper,vorobmean=TRUE,
new.points=batchsize,col.points.end="red",cex.points=2.5,
main="updated probability of excursion")

}
}
