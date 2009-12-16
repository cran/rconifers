
\name{print.sample.data}
\alias{print.sample.data}
\title{Print a CONIFERS sample.data object}

\description{
   These functions print a few summary values from a sample.data object.
}

\usage{
   \method{print}{sample.data}( x,
      digits = max( 3, getOption("digits") - 1 ), \dots )
}

\arguments{
   \item{x}{an object of class \code{sample.data}.}
   \item{digits}{number of digits to print.}
   \item{\dots}{other arguments.}
}



\references{

Ritchie, M.W. 2008. User's Guide and Help System for CONIFERS: A Simulator for Young Conifer Plantations Version
4.10. See \url{http://www.fs.fed.us/psw/programs/ecology_of_western_forests/projects/conifers/}

}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com},\cr
	     Martin W. Ritchie \email{mritchie@fs.fed.us} }

\seealso{    \code{\link{calc.max.sdi}},
	     \code{\link{sample.data}},
 	     \code{\link{summary.sample.data}}
}

\examples{
library( rconifers )

## set the variant to the SWO variant
# and set the species map
set.species.map( set.variant(0) )

# load and display CONIFERS SWO example plots
data( plots.swo )
print( plots.swo )

# load and display CONIFERS SWO example plants
data( plants.swo )
print( plants.swo )

# create the sample.data list object
sample.swo.3 <- list( plots=plots.swo, plants=plants.swo, age=3, x0=0.0 )
class(sample.swo.3)  <- "sample.data"
print( sample.swo.3)

}

\keyword{models}
