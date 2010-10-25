
\keyword{methods}
\name{set.species.map}
\alias{set.species.map}

\title{Set and update the species mapping table in the CONIFERS simulator}

\description{ This function assigns the mapping of user defined
  species codes to the functional species within the CONIFERS growth
  model.}

\usage{set.species.map( x, verbose=FALSE )}

		   
\arguments{
  \item{x}{ is a data.frame object with the same columns as \code{\link{smc}} and \code{\link{swo}}}
  \item{verbose}{ is switch for controlling the verbosity of the operation.}
}


%at least the zero-offset species
%    index in the species data frame, the functional species code (See
%    conifers manual, a user defined species code as a vector of
%    characters, a vector that contains the endemic mortality
%    coefficients, a vector that contains the maximum stand density
%    index value, a vector that contains the browse damage
%    coefficients, a vector that contains the mechanical damage
%    coefficients, a vector that contains genetic worth coefficients
%    for height growth, and a vector that contains genetic worth
%    coefficients for diameter growth.}

\details{The set.species.map is how the user controls the behavior of
the functional species that is used to project individual plants
forward. The idx *must* be a zero-offset index vector which is
assigned the element \code{SPECIES_RECORD[idx]} in the model (see \code{\link{smc}}). This
function is used to control the mapping between the user species code
and the functional species codes in the growth model.
}

\value{
	no value is returned.
}



\references{

Ritchie, M.W. 2008. User's Guide and Help System for CONIFERS: A Simulator for Young Conifer Plantations Version
4.10. See \url{http://www.fs.fed.us/psw/programs/ecology_of_western_forests/projects/conifers/}

}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com},\cr
	     Martin W. Ritchie \email{mritchie@fs.fed.us} }

\seealso{    \code{\link{calc.max.sdi}},
	     \code{\link{impute}},
	     \code{\link{plants.smc}},
	     \code{\link{plots.smc}}
	     \code{\link{plants.swo}},
	     \code{\link{plots.swo}},
	     \code{\link{rconifers}},
	     \code{\link{sample.data}},	     		  
	     \code{\link{set.variant}},
 	     \code{\link{smc}}, 
 	     \code{\link{summary.sample.data}},
	     \code{\link{swo}}
}

\examples{
library( rconifers )

# Example 1
## set the variant to the SWO variant and set species map
set.variant(0)
data( swo )
set.species.map( swo )

# Example 2
## set the variant to the SMC variant
set.variant(1)

## and manually set the genetic worth values in the functional species coefficients
data( smc )
smc[smc$code == 'DF',]$genetic.worth.d <- 20
smc[smc$code == 'DF',]$genetic.worth.h <- 20
set.species.map( smc, verbose=TRUE )


}

\keyword{models}


