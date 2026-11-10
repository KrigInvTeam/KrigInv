\name{print_uncertainty_nd}
\Rdversion{1.4.1}
\alias{print_uncertainty_nd}

\title{
Print a measure of uncertainty for functions with dimension d strictly larger than 2
}
\description{
This function draws projections on various plans of a given measure of uncertainty. 
The function can be used to print relevant outputs after having used the function \code{\link{EGI}} or \code{\link{EGIparallel}}.
}

\usage{
print_uncertainty_nd(model,T,type="pn",lower=NULL,upper=NULL,
  		resolution=20, nintegpoints=400,
  		cex.lab=1,cex.contourlab=1,cex.axis=1,
			nlevels=10,levels=NULL,
			xdecal=3,ydecal=3, option="mean", pairs=NULL,...)
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
	\code{"sur"}, \code{"imse"}, \code{"timse"} if we print a measure of uncertainty corresponding to one criterion. 
	The \code{"vorob"} option is not available due to the difficulty of computing a Vorob'ev threshold.
}
  \item{lower}{
Vector containing the lower bounds of the input domain. If nothing is set we use a vector of 0.
}
  \item{upper}{
Vector containing the upper bounds of the input domain. If nothing is set we use a vector of 1.
}
  \item{resolution}{
Number of points to discretize a plan included in the domain. For the moment, we cannot use values higher than 40 do to 
computation time, except when the argument \code{"pairs"} is not set to its default value.
}
  \item{nintegpoints}{
to do
}
  \item{cex.lab}{
Multiplicative factor for the size of titles of the axis.
}
  \item{cex.contourlab}{
Multiplicative factor for the size of labels of the contour plot.
}
  \item{cex.axis}{
Multiplicative factor for the size of the axis graduations.
}
  \item{nlevels}{
Integer corresponding to the number of levels of the contour plot.
}
  \item{levels}{
Array: one can directly set the levels of the contour plot.
}
  \item{xdecal}{
Optional position shifting of the titles of the x axis.
}
  \item{ydecal}{
Optional position shifting of the titles of the y axis.
}
  \item{option}{
Optional argument (a string). The 3 possible values are \code{"mean"} (default), \code{"max"} and \code{"min"}.
}
  \item{pairs}{
Optional argument. When set to \code{NULL} (default) the function performs the projections on plans spanned by each pair (i,j) of dimension. Otherwise, the argument is an array of size 2 corresponding to the dimensions spanning the (only) plan on which the projection is performed.
}
\item{...}{
Additional arguments to the \code{image} function.
}
}

\value{
The integrated uncertainty
}

\references{
Bect J., Ginsbourger D., Li L., Picheny V., Vazquez E. (2012), \emph{Sequential design of computer experiments for the estimation of a probability of failure}, Statistics and Computing vol. 22(3), pp 773-793. \doi{10.1007/s11222-011-9241-4}
}

\author{ 
Clement Chevalier (University of Neuchatel, Switzerland)
}

\seealso{\code{\link{print_uncertainty_1d}},\code{\link{print_uncertainty_2d}}}

\examples{
#print_uncertainty_nd

set.seed(9)
N <- 30 #number of observations
T <- -1 #threshold
testfun <- hartman3
#The hartman3 function is defined over the domain [0,1]^3. 

lower <- rep(0,times=3)
upper <- rep(1,times=3)

#a 30 points initial design
design <- data.frame( matrix(runif(3*N),ncol=3) )
response <- apply(design,1,testfun)

#km object with matern3_2 covariance
#params estimated by ML from the observations
model <- km(formula=~., design = design, 
	response = response,covtype="matern3_2")

\donttest{
print_uncertainty_nd(model=model,T=T,main="average probability of excursion",type="pn",
                    option="mean")

print_uncertainty_nd(model=model,T=T,main="maximum probability of excursion",type="pn",
                     option="max")

}
}
