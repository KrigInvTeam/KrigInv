\name{print_uncertainty_1d}
\Rdversion{1.4.1}
\alias{print_uncertainty_1d}

\title{
Prints a measure of uncertainty for 1d function
}
\description{
This function draws the value of a given measure of uncertainty over the whole input domain (1D).
The function can be used to print relevant outputs after having used the function \code{\link{EGI}} or \code{\link{EGIparallel}}.
}
\usage{
print_uncertainty_1d(model, T, type = "pn",
lower = 0, upper = 1, resolution = 500, new.points = 0,
xscale = c(0, 1), show.points = TRUE, cex.points = 1,
cex.axis = 1, pch.points.init = 17, pch.points.end = 17,
col.points.init = "black", col.points.end = "red", xaxislab = NULL,
yaxislab = NULL, xaxispoint = NULL, yaxispoint = NULL,
vorobmean=FALSE,krigmeanplot=FALSE,Tplot=FALSE,consQuantile=NULL,...)
}

\arguments{
  \item{model}{
Kriging model of \code{\link[DiceKriging]{km}} class.
}
  \item{T}{
Array containing one or several thresholds.
}
  \item{type}{
Type of uncertainty that the user wants to print.
	Possible values are \code{"pn"} (probability of excursion), or
	\code{"sur"}, \code{"imse"}, \code{"timse"}, \code{"vorob"} if we print a measure of uncertainty corresponding to one criterion.
}
  \item{lower}{
Lower bound for the input domain.
}
  \item{upper}{
Upper bound for the input domain.
}
  \item{resolution}{
Number of points to discretize the interval (lower,upper).
}
  \item{new.points}{
Number of new observations.
These observations are the last new.points observations and can be printed in another color and the initial observations (see argument: \code{col.points.end}).
}
  \item{xscale}{
If one wants to rescale the input domain on another interval it is possible to set this vector of size 2. The new interval will be translated by \code{xscale[1]} and expanded by a factor \code{xscale[2] - xscale[1]}.
}
  \item{show.points}{
Boolean: should we show the observations on the graph ?
}
  \item{cex.points}{
Multiplicative factor for the size of the points.
}
  \item{cex.axis}{
Multiplicative factor for the size of the axis graduations.
}
  \item{pch.points.init}{
Symbol for the \code{n-new.points} first observations.
}
  \item{pch.points.end}{
Symbol for the \code{new.points} last observations.
}
  \item{col.points.init}{
Color for the \code{n-new.points} first observations.
}
  \item{col.points.end}{
Color for the \code{new.points} last observations.
}
  \item{xaxislab}{
Optional new labels that will replace the normal levels on x axis.
}
  \item{yaxislab}{
Optional new labels that will replace the normal levels on y axis.
}
  \item{xaxispoint}{
Position of these new labels on x axis.
}
  \item{yaxispoint}{
Position of these new labels on y axis.
}
  \item{vorobmean}{
Optional boolean. When it is set to \code{TRUE} the Vorob'ev expectation is plotted. It corresponds to the averaged excursion set, using the definition of Vorob'ev. Here, the estimated set is the set above the Vorob'ev threshold (plotted in blue).
}
\item{krigmeanplot}{
When set to \code{TRUE} a kriging mean is plotted on a second y axis.
}
\item{Tplot}{
When set to \code{TRUE}, and if \code{krigmeanplot} is also set to \code{TRUE}, draws horizontal lines corresponding to the different thresholds \code{T}.
}
\item{consQuantile}{
Optional value for plotting conservative quantiles. In order to plot \itemize{
\item Conservative estimates: \code{consQuantile} is a list containing at least \code{consLevel} (scalar), with the option \code{typeEx} (character, default = ">").
\item Generic Vorob'ev quantiles: \code{consQuantile} is a scalar corresponding to the Vorob'ev quantile level.
}
}
\item{...}{
Additional arguments to the \code{plot} function.
}
}

\value{
The integrated uncertainty. If the conservative estimate is computed, it also returns the conservative quantile level.
}

\references{
Bect J., Ginsbourger D., Li L., Picheny V., Vazquez E. (2012), \emph{Sequential design of computer experiments for the estimation of a probability of failure}, Statistics and Computing vol. 22(3), pp 773-793. \doi{10.1007/s11222-011-9241-4}
}

\author{
Clement Chevalier (University of Neuchatel, Switzerland)

Dario Azzimonti (IDSIA, Switzerland)
}

\seealso{\code{\link{print_uncertainty_2d}},\code{\link{print_uncertainty_nd}}}

\examples{
#print_uncertainty_1d

set.seed(9)
N <- 9 #number of observations
T <- c(-0.2,0.2) #thresholds
testfun <- sin
lower <- c(0)
upper <- c(6)

#a 20 points initial design
design <- data.frame( lower+(upper-lower)*matrix(runif(N),ncol=1) )
response <- testfun(design)

#km object with matern3_2 covariance
#params estimated by ML from the observations
model <- km(formula=~., design = design,
	response = response,covtype="matern3_2")

print_uncertainty_1d(model=model,T=T,lower=lower,upper=upper,
  main="probability of excursion",xlab="x",ylab="pn",
  cex.points=1.5,col.points.init="red",
  krigmeanplot=TRUE,Tplot=TRUE)

\donttest{
uq1d <- print_uncertainty_1d(model=model,T=T,lower=lower,upper=upper,
  main="probability of excursion",xlab="x",ylab="pn",
  cex.points=1.5,col.points.init="red",
  krigmeanplot=TRUE,Tplot=TRUE,consQuantile =list(consLevel=0.95))

  print_uncertainty_1d(model=model,T=T,lower=lower,upper=upper,
  main="probability of excursion",xlab="x",ylab="pn",
  cex.points=1.5,col.points.init="red",
  krigmeanplot=TRUE,Tplot=TRUE,consQuantile =uq1d[2])
}

}
