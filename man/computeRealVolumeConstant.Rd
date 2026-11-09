
\name{computeRealVolumeConstant}
\Rdversion{1.4.1}
\alias{computeRealVolumeConstant}

\title{A constant used to calculate the expected excursion set's volume variance}

\description{
This function computes a constant used to calculate exactly the value of the \code{"jn"} criterion. 
Computing this constant does NOT change the optimum of the \code{"jn"} criterion. 
Therefore, its calculation is indicative only and is only necessary to know exactly (in expectation) the excursion set's volume variance. 
}

\usage{
computeRealVolumeConstant(model,integration.points,
integration.weights=NULL,T)
}

\arguments{

  \item{model}{A Kriging model of \code{\link[DiceKriging]{km}} class.}
	\item{integration.points}{p*d matrix of points for numerical integration in the X space.}
	\item{integration.weights}{(Optional) Vector of size p corresponding to the weights of these integration points. If not provided, all weights are set to 1.}
 	\item{T}{Target threshold.}
 	
}

\details{
Note that, even if the \code{"jn"} criterion can be used with more than one threshold, the computation of this constant is implemented only 
when the number of threshold is equal to 1.
}

\value{a scalar
}

\references{
Chevalier C., Bect J., Ginsbourger D., Vazquez E., Picheny V., Richet Y. (2014), \emph{Fast parallel kriging-based stepwise uncertainty reduction with application to the identification of an excursion set}, Technometrics, vol. 56(4), pp 455-465. \doi{10.1080/00401706.2013.860918}
}

\author{
Clement Chevalier (University of Neuchatel, Switzerland)
}

\seealso{\code{\link{precomputeUpdateData}}, \code{\link{predict_nobias_km}}}

\examples{
#computeRealVolumeConstant

set.seed(9)
N <- 20 #number of observations
testfun <- branin
T <- 80

#a 20 points initial design
design <- data.frame( matrix(runif(2*N),ncol=2) )
response <- testfun(design)

#km object with matern3_2 covariance
#params estimated by ML from the observations
model <- km(formula=~., design = design, 
            response = response,covtype="matern3_2")

integcontrol <- list(n.points=500,distrib="jn",init.distrib="MC")
obj <- integration_design(integcontrol=integcontrol,
	lower=c(0,0),upper=c(1,1),model=model,T=T)
	
integration.points <- obj$integration.points
integration.weights <- obj$integration.weights

\donttest{
computeRealVolumeConstant(model=model,
	integration.points=integration.points,
        integration.weights=integration.weights,T=T)
}
}
