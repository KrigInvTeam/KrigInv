\name{vorob_threshold}
\Rdversion{1.4.1}
\alias{vorob_threshold}

\title{ Calculation of the Vorob'ev threshold }
\description{
Evaluation of the Vorob'ev threshold given an excursion probability vector. This threshold is such that the volume of the set (x : pn(x) > threshold) is equal to the integral of pn. 
}
\usage{
vorob_threshold(pn)
}

\arguments{
  	\item{pn}{Input vector of arbitrary size containing the excursion probabilities pn(x).}
}
\details{
In this function, all the points x are supposed to be equaly weighted.
}

\value{a scalar: the Vorob'ev thresold}
\references{
Chevalier C., Ginsbourger D., Bect J., Molchanov I. (2013) \emph{Estimating and quantifying uncertainties on level sets using the Vorob'ev expectation and deviation with gaussian process models} mODa 10, Advances in Model-Oriented Design and Analysis, Contributions to Statistics, pp 35-43. \doi{10.1007/978-3-319-00218-7_5}

Chevalier C. (2013) \emph{Fast uncertainty reduction strategies relying on Gaussian process models} Ph.D Thesis, University of Bern
}

\author{
Clement Chevalier (University of Neuchatel, Switzerland)
}

\seealso{\code{\link{max_vorob_parallel}}, \code{\link{vorob_optim_parallel}} }

\examples{
#vorob_threshold

set.seed(9)
N <- 20 #number of observations
T <- 80 #threshold
testfun <- branin

#a 20 points initial design
design <- data.frame( matrix(runif(2*N),ncol=2) )
response <- testfun(design)

#km object with matern3_2 covariance
#params estimated by ML from the observations
model <- km(formula=~., design = design, 
	response = response,covtype="matern3_2")

\donttest{
###we need to compute some additional arguments:
#integration points, and current kriging means and variances at these points
integcontrol <- list(n.points=50,distrib="sobol")
obj <- integration_design(integcontrol=integcontrol,
lower=c(0,0),upper=c(1,1),model=model,T=T)

integration.points <- obj$integration.points

pred <- predict_nobias_km(object=model,newdata=integration.points,
type="UK",se.compute=TRUE)
pn <- pnorm((pred$mean-T)/pred$sd)

vorob_threshold(pn)
}
}
