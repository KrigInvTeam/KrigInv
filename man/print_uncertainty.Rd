\name{print_uncertainty}
\Rdversion{1.4.1}
\alias{print_uncertainty}

\title{
Prints a measure of uncertainty for a function of any dimension
}
\description{
This function prints the value of a given measure of uncertainty. 
The function can be used to print relevant outputs after having used the function \code{\link{EGI}} or \code{\link{EGIparallel}}.
}
\usage{
print_uncertainty(model, T, type = "pn", ...)
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
  \item{\dots}{
Other arguments of the functions \code{print_uncertainty_1d}, \code{2d} or \code{nd}.
}
}

\value{
the integrated uncertainty
}

\references{
Bect J., Ginsbourger D., Li L., Picheny V., Vazquez E. (2012), \emph{Sequential design of computer experiments for the estimation of a probability of failure}, Statistics and Computing vol. 22(3), pp 773-793. \doi{10.1007/s11222-011-9241-4}
}

\author{ 
Clement Chevalier (University of Neuchatel, Switzerland)
}

\seealso{\code{\link{print_uncertainty_1d}},\code{\link{print_uncertainty_2d}},\code{\link{print_uncertainty_nd}}}

\examples{
#print_uncertainty

set.seed(9)
N <- 20 #number of observations
T <- c(80,100) #threshold
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

#you could do many plots, but only one is run here
print_uncertainty(model=model,T=T,main="probability of excursion",type="pn")
#print_uncertainty(model=model,T=T,main="Vorob'ev uncertainty",type="vorob")
#print_uncertainty(model=model,T=T,main="imse uncertainty",type="imse")
#print_uncertainty(model=model,T=T,main="timse uncertainty",type="timse")
#print_uncertainty(model=model,T=T,main="sur uncertainty",type="sur")
#print_uncertainty(model=model,T=T,main="probability of excursion",type="pn",
#vorobmean=TRUE)
}
