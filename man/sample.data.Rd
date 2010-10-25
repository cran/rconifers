
%%	$Id $	

\name{sample.data}
\alias{sample.data}

\title{CONIFERS forest growth model sample data}

\description{ A list object of type \code{sample.data} stores all of
the basic information about a \code{\link{data.frame}} object
representing a sample of plants.}

\usage{
sample.swo.3 <- list( plots=plots.swo, plants=plants.swo, age=3, x0=0.0 )
class(sample.swo.3)  <- "sample.data"
}

\details{
	To create the basic data type used in rconifers, you create a
	list object with the following elements (order is not
	important):
\describe{
\item{plots.swo}{is a \link{data.frame} with the the same elements as \code{\link{plots.swo}}.}
\item{plants.swo}{is a \link{data.frame} with the the same elements as \code{\link{plants.swo}}.}
\item{age}{is an integer value that represents the age of the plants,
in years.}
\item{x0}{is the $x_{0}$ coefficient for the Hann and Wang (1990) mortality model.}
}
}





\references{

Ritchie, M.W. 2008. User's Guide and Help System for CONIFERS: A Simulator for Young Conifer Plantations Version
4.10. See \url{http://www.fs.fed.us/psw/programs/ecology_of_western_forests/projects/conifers/}

}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com},\cr
	     Martin W. Ritchie \email{mritchie@fs.fed.us} }

\seealso{    \code{\link{plants.smc}},
	     \code{\link{plots.smc}}
	     \code{\link{plants.swo}},
	     \code{\link{plots.swo}},
	     \code{\link{sample.data}},
	     \code{\link{set.species.map}},		  
	     \code{\link{set.variant}},		  
 	     \code{\link{smc}}, 
 	     \code{\link{summary.sample.data}},
	     \code{\link{swo}}, 
	     \code{\link{thin}}
}

\examples{
library( rconifers )
## Example for SWO variant
## set the variant to the SWO variant and set species map
##set.species.map( set.variant(0) )
set.variant(0)

## grow the data that was originally swo in the smc variant
# load and display CONIFERS example plots
data( plots.swo )
print( plots.swo )

# load and display CONIFERS example plants
data( plants.swo )
print( plants.swo )

# create the sample.data list object
sample.swo.3 <- list( plots=plots.swo, plants=plants.swo, age=3, x0=0.0 )
class(sample.swo.3)  <- "sample.data"


}

\keyword{models}


