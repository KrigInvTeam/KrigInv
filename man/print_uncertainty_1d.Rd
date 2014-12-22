\name{print_uncertainty_1d}
\Rdversion{1.3}
\alias{print_uncertainty_1d}

\title{
Prints a measure of uncertainty for 1d function.
}
\description{
This function draws the value of a given measure of uncertainty over the whole input domain (1D). 
Possible measures are \code{"pn"} (beeing the probability of excursion) and measures specific to a sampling criterion: \code{"sur"}, \code{"timse"} and \code{"imse"}.
This function can be used to print relevant outputs after having used the function \code{\link{EGI}}.
}
\usage{
print_uncertainty_1d(model, T, type = "pn", 
lower = 0, upper = 1, resolution = 500, 
new.points = 0, xlab = "", ylab = "", main = "", 
xscale = c(0, 1), show.points = TRUE, cex.main = 1, cex.lab = 1, 
cex.points = 1, cex.axis = 1, pch.points.init = 17, pch.points.end = 17, 
col.points.init = "black", col.points.end = "red", xaxislab = NULL, 
yaxislab = NULL, xaxispoint = NULL, yaxispoint = NULL, 
xdecal = 3, ydecal = 3,DiceViewplot=FALSE,vorobmean=FALSE)
}

\arguments{
  \item{model}{
Kriging model of \code{\link[DiceKriging]{km}} class.
}
  \item{T}{
Target value (scalar).
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
  \item{xlab}{
Label for the x axis.
}
  \item{ylab}{
Label for the y axis.
}
  \item{main}{
Title of the graph.
}
  \item{xscale}{
If one wants to rescale the input domain on another interval it is possible to set this vector of size 2. The new interval will be translated by \code{xscale[1]} and expanded by a factor \code{xscale[2] - xscale[1]}.
}
  \item{show.points}{
Boolean: should we show the observations on the graph ?
}
  \item{cex.main}{
Multiplicative factor for the size of the title.
}
  \item{cex.lab}{
Multiplicative factor for the size of titles of the axis.
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
  \item{xdecal}{
Optional position shifting of the titles of the x axis.
}
  \item{ydecal}{
Optional position shifting of the titles of the y axis.
}
  \item{DiceViewplot}{
Optional boolean. When it is set to \code{TRUE} (default) a second plot is added, generated with the DiceView package. This plot shows the kriging mean and confidence intervals on the whole input domain.
}
\item{vorobmean}{
Optional boolean. When it is set to \code{TRUE} the Vorob'ev expectation is plotted. It corresponds to the averaged excursion set, using the definition of Vorob'ev. Here, the estimated set is the set above the Vorob'ev threshold (plotted in blue).
}
}

\value{
the integrated uncertainty
}

\references{
Bect J., Ginsbourger D., Li L., Picheny V., Vazquez E. (2010), \emph{Sequential design of computer experiments for the estimation of a probability of failure}, Statistics and Computing, pp.1-21, 2011, \url{http://arxiv.org/abs/1009.5177}
}

\author{ 
Clement Chevalier (IMSV, Switzerland, and IRSN, France)
}

\seealso{\code{\link{EGI}}}

\examples{
#print_uncertainty_1d

set.seed(8)
N <- 9 #number of observations
T <- 1 #threshold
testfun <- fundet
lower <- c(0)
upper <- c(1)

#a 9 points initial design
design <- data.frame( matrix(runif(N),ncol=1) )
response <- testfun(design)

#km object with matern3_2 covariance
#params estimated by ML from the observations
model <- km(formula=~., design = design, 
	response = response,covtype="matern3_2")

print_uncertainty_1d(DiceViewplot=FALSE,model=model,T=T,
main="probability of excursion",cex.points=1.5,col.points.init="red",vorobmean=TRUE)


}