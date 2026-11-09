\name{max_infill_criterion}
\alias{max_infill_criterion}
\title{Optimizer for the infill criteria}
\Rdversion{1.4.1}

\description{
Optimization, of the chosen infill criterion (maximization or minimization, depending on the case)
}
\usage{
max_infill_criterion(lower, upper, optimcontrol = NULL, 
method, T, model, method.param = NULL)
}

\arguments{
  \item{lower}{
Vector containing the lower bounds of the design space.
}
  \item{upper}{
Vector containing the upper bounds of the design space.
}
  \item{optimcontrol}{
Optional list of control parameters for the optimization of the sampling criterion. The field \code{method} defines which optimization method is used. It can be either \code{"genoud"} (default) for an optimization using the genoud algorithm, or \code{"discrete"} for an optimization over a specified discrete set. 
If the field \code{method} is set to \code{"genoud"}, one can set some parameters of this algorithm: 
\code{pop.size}  (default : 50*d),  \code{max.generations} (10*d), \code{wait.generations} (2),  \code{BFGSburnin} (2) and the mutations \code{P1}, \code{P2}, up to \code{P9} (see \code{\link[rgenoud]{genoud}}). Numbers into brackets are the default values.
If the field \code{method} is set to \code{"discrete"}, one can set the field \code{optim.points}: p * d matrix corresponding to the p points where the criterion will be evaluated. If nothing is specified, 100*d points are chosen randomly. 
}
  \item{method}{
Criterion used for choosing observations: \code{"ranjan"} (default) , \code{"bichon"}, \code{"tsee"}, or \code{"tmse"}.
}
  \item{T}{
Array containing one or several thresholds. The \code{"tmse"} criterion can be used with multiple thresholds. The \code{"ranjan"}, \code{"bichon"}, \code{"tsee"} criteria can be used with only one threshold.
}
  \item{model}{
A Kriging model of \code{\link[DiceKriging]{km}} class.
}
  \item{method.param}{
Optional tolerance value (scalar). Default value is 1 for \code{"ranjan"} and \code{"bichon"}, and 0 for \code{"tmse"}.
}
}

\value{A list with components:
\item{par}{The best set of parameters found.}
\item{value}{The value of the chosen criterion at par.}
\item{allvalues}{If an optimization on a discrete set of points is chosen, the value of the criterion at all these points.}
}

\references{
Bect J., Ginsbourger D., Li L., Picheny V., Vazquez E. (2012), \emph{Sequential design of computer experiments for the estimation of a probability of failure}, Statistics and Computing vol. 22(3), pp 773-793. \doi{10.1007/s11222-011-9241-4}

Picheny V., Ginsbourger D., Roustant O., Haftka R.T., (2010) \emph{Adaptive designs of experiments for accurate approximation of a target region}, J. Mech. Des. vol. 132(7). \doi{10.1115/1.4001873}

Bichon B.J., Eldred M.S., Swiler L.P., Mahadevan S., McFarland J.M. (2008) \emph{Efficient global reliability analysis for nonlinear implicit performance functions}, AIAA Journal 46(10), pp 2459-2468. \doi{10.2514/1.34321}

Ranjan P., Bingham D., Michailidis G. (2008) \emph{Sequential experiment design for contour estimation from complex computer codes} Technometrics 50(4), pp 527-541. \doi{10.1198/004017008000000541}
}
\author{
Victor Picheny (INRA, Toulouse, France)

David Ginsbourger (IDIAP Martigny and University of Bern, Switzerland)

Clement Chevalier (University of Neuchatel, Switzerland)
}

\seealso{\code{\link{EGI}},\code{\link{ranjan_optim}},\code{\link{tmse_optim}},\code{\link{bichon_optim}},\code{\link{tsee_optim}}}

\examples{
#max_infill_criterion

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

optimcontrol <- list(method="genoud",pop.size=50)

\donttest{
obj <- max_infill_criterion(lower=lower,upper=upper,optimcontrol=optimcontrol,
                             method="bichon",T=T,model=model)

obj$par;obj$value
new.model <- update(object=model,newX=obj$par,newy=testfun(obj$par),cov.reestim=TRUE)

par(mfrow=c(1,2))
print_uncertainty(model=model,T=T,type="pn",lower=lower,upper=upper,
cex.points=2.5,main="probability of excursion")

print_uncertainty(model=new.model,T=T,type="pn",lower=lower,upper=upper,
new.points=1,col.points.end="red",cex.points=2.5,main="updated probability of excursion")
}
}
