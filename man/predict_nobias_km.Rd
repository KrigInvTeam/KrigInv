\name{predict_nobias_km}
\Rdversion{1.4.1}
\alias{predict_nobias_km}

\title{ Kriging predictions}
\description{
This function is similar to the predict.km function from the DiceKriging package. 
The only change is the additionnal F.newdata output.
}
\usage{
predict_nobias_km(object, newdata, type = "UK", 
se.compute = TRUE, cov.compute = FALSE, low.memory=FALSE,...)
}

\arguments{
  \item{object}{
A Kriging model of \code{\link[DiceKriging]{km}} class.
}
  \item{newdata}{
Vector, matrix or data frame containing the points where to perform predictions.
}
  \item{type}{
Character string corresponding to the kriging family, to be chosen between simple kriging ("SK"), or universal kriging ("UK").
}
  \item{se.compute}{
Optional boolean. If \code{FALSE}, only the kriging mean is computed. If \code{TRUE}, the kriging standard deviation and confidence intervals are computed too.
}
  \item{cov.compute}{
Optional boolean. If \code{TRUE} the conditional covariance matrix is computed.
}
  \item{low.memory}{
Optional boolean. If set to \code{TRUE} the function will only return kriging means and standard deviations.
}
  \item{\dots}{
No other arguments.
}
}

\value{
  \item{mean }{ kriging mean (including the trend) computed at \code{newdata}. }
  \item{sd }{ kriging standard deviation computed at \code{newdata}. 
  Not computed if \code{se.compute=FALSE}. }
  \item{cov }{ kriging conditional covariance matrix. Not computed if \code{cov.compute=FALSE} (default).}
  \item{lower95, }{ }
  \item{upper95 }{ bounds of the 95 \% confidence interval computed at \code{newdata} (to be interpreted with special care when parameters are estimated, see description above). Not computed if \code{se.compute=FALSE}.}
  \item{c }{ an auxiliary matrix, containing all the covariances between newdata and the initial design points. }
  \item{Tinv.c }{ an auxiliary vector, equal to \code{T^(-1)*c}.}
  \item{F.newdata}{value of the trend function at \code{newdata}.}
}

\section{Warning}{
Beware that the only consistency check between \code{newdata} and the experimental design is to test whether they have same number of columns. In that case, the columns of \code{newdata} are interpreted in the same order as the initial design.
}

\references{

N.A.C. Cressie (1993), \emph{Statistics for spatial data}, Wiley series in probability and mathematical statistics.

A.G. Journel and C.J. Huijbregts (1978), \emph{Mining Geostatistics}, Academic Press, London.

D.G. Krige (1951), A statistical approach to some basic mine valuation problems on the witwatersrand, \emph{J. of the Chem., Metal. and Mining Soc. of South Africa}, \bold{52} no. 6, 119-139.

J.D. Martin and T.W. Simpson (2005), Use of kriging models to approximate deterministic computer models, \emph{AIAA Journal}, \bold{43} no. 4, 853-863.

G. Matheron (1963), Principles of geostatistics, \emph{Economic Geology}, \bold{58},
1246-1266.

G. Matheron (1969), Le krigeage universel, \emph{Les Cahiers du Centre de Morphologie Mathematique de Fontainebleau}, \bold{1}.

J.-S. Park and J. Baek (2001), Efficient computation of maximum likelihood estimators in a spatial linear model with power exponential covariogram, \emph{Computer Geosciences}, \bold{27} no. 1, 1-7.

C.E. Rasmussen and C.K.I. Williams (2006), \emph{Gaussian Processes for Machine Learning}, the MIT Press, \url{https://gaussianprocess.org/gpml/}

J. Sacks, W.J. Welch, T.J. Mitchell, and H.P. Wynn (1989), Design and analysis of computer experiments, \emph{Statistical Science}, \bold{4}, 409-435.

}

\author{
O. Roustant (Ecole des Mines de St-Etienne, France)

David Ginsbourger (IDIAP Martigny and University of Bern, Switzerland)
}

\seealso{  \code{\link[DiceKriging]{predict.km}},  \code{\link[DiceKriging]{km}} }

\keyword{models}
\keyword{methods}

\examples{
#predict_nobias_km
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

n.grid <- 100
x.grid <- y.grid <- seq(0,1,length=n.grid)

newdata <- expand.grid(x.grid,y.grid)
pred <- predict_nobias_km(object=model,newdata=newdata,type="UK",se.compute=TRUE)

z.grid1 <- matrix(pred$mean, n.grid, n.grid)
z.grid2 <- matrix(pred$sd, n.grid, n.grid)

par(mfrow=c(1,2))

#plots: contour of the kriging mean and stdev
image(x=x.grid,y=y.grid,z=z.grid1,col=grey.colors(10))
contour(x=x.grid,y=y.grid,z=z.grid1,15,add=TRUE)
points(design, col="black", pch=17, lwd=4,cex=2)
title("Kriging mean")

image(x=x.grid,y=y.grid,z=z.grid2,col=grey.colors(10))
contour(x=x.grid,y=y.grid,z=z.grid2,15,add=TRUE)
points(design, col="black", pch=17, lwd=4,cex=2)
title("Kriging standard deviation")
}
