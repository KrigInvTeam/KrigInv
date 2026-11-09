\name{computeQuickKrigcov}
\Rdversion{1.4.1}
\alias{computeQuickKrigcov}

\title{Quick computation of kriging covariances}

\description{
Computes kriging covariances between some new points and many integration points, using precomputed data.
}

\usage{
computeQuickKrigcov(model,integration.points,X.new,
precalc.data, F.newdata , c.newdata)
}

\arguments{
 
  \item{model}{A Kriging model of \code{\link[DiceKriging]{km}} class.}
	\item{integration.points}{p*d matrix of fixed integration points in the X space.}
 	\item{X.new}{q*d matrix of new points. The calculated covariances are the covariances between these new point and the integration points.}
 	\item{precalc.data}{List containing precalculated data. This list is generated using the function \code{\link{precomputeUpdateData}}.}
 	\item{F.newdata}{The value of the kriging trend basis function at point X.new.}
 	\item{c.newdata}{The (unconditional) covariance between X.new and the design points.}

}

\details{
This function requires to use another function in order to generate the proper arguments. 
The argument \code{precalc.data} can be generated using \code{\link{precomputeUpdateData}}.
The arguments \code{F.newdata} and \code{c.newdata} can be obtained using \code{\link{predict_nobias_km}}.
}

\value{Matrix of size p*q containing kriging covariances
}

\references{
Chevalier C., Bect J., Ginsbourger D., Vazquez E., Picheny V., Richet Y. (2014), \emph{Fast parallel kriging-based stepwise uncertainty reduction with application to the identification of an excursion set}, Technometrics, vol. 56(4), pp 455-465. \doi{10.1080/00401706.2013.860918}

Chevalier C., Ginsbourger D. (2014), \emph{Corrected Kriging update formulae for batch-sequential data assimilation}, in Pardo-Iguzquiza, E., et al. (Eds.) Mathematics of Planet Earth, pp 119-122. \doi{10.1007/978-3-642-32408-6_29}
}

\author{
Clement Chevalier (University of Neuchatel, Switzerland)
}

\seealso{\code{\link{precomputeUpdateData}}, \code{\link{predict_nobias_km}}}


\examples{
#computeQuickKrigcov

set.seed(9)
N <- 20 #number of observations
testfun <- branin

#a 20 points initial design
design <- data.frame( matrix(runif(2*N),ncol=2) )
response <- testfun(design)

#km object with matern3_2 covariance
#params estimated by ML from the observations
model <- km(formula=~., design = design, 
            response = response,covtype="matern3_2")

#the integration.points are the points where we want to 
#compute predictions/covariances if a point new.x is added 
#to the DOE
x.grid <- seq(0,1,length=20)
integration.points <- expand.grid(x.grid,x.grid)
integration.points <- as.matrix(integration.points)

#precalculation
precalc.data <- precomputeUpdateData(model=model,
                     integration.points=integration.points)

#now we can compute quickly kriging covariances 
#between these data and any other points.
#example if 5 new points are added:
X.new <- matrix(runif(10),ncol=2)
pred <- predict_nobias_km(object=model,
                          newdata=X.new,type="UK",se.compute=TRUE)

kn <- computeQuickKrigcov(model=model,
                    integration.points=integration.points,X.new=X.new,
                    precalc.data=precalc.data,
                    F.newdata=pred$F.newdata,
                    c.newdata=pred$c)
}
