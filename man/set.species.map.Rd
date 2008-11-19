
\keyword{methods}
\name{set.species.map}
\alias{set.species.map}

\title{Set and update the species mapping table in the CONIFERS simulator}

\description{ This function assigns the mapping of user defined
  species codes to the functional species within the CONIFERS growth
  model.}

\usage{
set.species.map( sp.map=list(idx,fsp,code,em,msdi,b,m,gwh,gwd) )
}
		   
\arguments{
  \item{sp.map}{a list that contains the following elements:
  \describe{
  	    \item{idx}{a vector of index values (zero-offset)}
  	    \item{fsp}{a vector of functional species codes}
  	    \item{code}{a vector of strings representing the species
  code used during field data collection}
       	    \item{em}{endemic (background) mortality}
	    \item{msdi}{a vector containing the maximum stand density
  index values}
  	\item{b}{a vector that contains the browse damage}
  	\item{m}{a vector that contains mechanical damage}
  	\item{gwh}{genetic worth height multiplier}
  	\item{gwd}{genetic worth diameter multiplier}
  }
}


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
assigned the element \code{SPECIES_RECORD[idx]} in the model. This
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
	     \code{\link{plants}},	
	     \code{\link{plots}} 
	     \code{\link{project}},
	     \code{\link{rand.seed}},
	     \code{\link{rconifers}},
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

## set the variant to the SWO variant
set.variant( 0 )

# load the Southwest-Oregon species coefficients into R as a data.frame object
data( swo )

# set the species map 
sp.map <- list(idx=swo$idx,
                 fsp=swo$fsp,
                 code=as.character(swo$code),
                 em=swo$endemic.mort,
		 msdi=swo$max.sdi,
                 b=swo$browse.damage,
                 m=swo$mechanical.damage,
                 gwh=swo$genetic.worth.h,
                 gwd=swo$genetic.worth.d)
set.species.map( sp.map )

## grow the data that was originally swo in the smc variant
# load and display CONIFERS example plots
data( plots )
print( plots )

# load and display CONIFERS example plants
data( plants )
print( plants )

# randomly remove 10 crown.width observations to test
# the impute function
blanks <- sample( 1:nrow( plants ), 10, replace=FALSE )
plants[blanks,]$crown.width <- NA

# create the sample.data list object
sample.3 <- list( plots=plots, plants=plants, age=3, x0=0.0 )
class(sample.3)  <- "sample.data"

# fill in missing values
sample.3.imp <- impute( sample.3 )

# print the maximum stand density index for the current settings
print( calc.max.sdi( sample.3.imp ) )

# print a summary of the sample
print( sample.3.imp )

# now, project the sample forward for 20 years
# with all of the options turned off
sample.23 <- project( sample.3.imp, 
	  20, 
	  control=list(rand.err=0,rand.seed=0,endemic.mort=0,sdi.mort=0))

## print the projected summaries
print( sample.23 )

## plot the diagnostics from the fit a linear dbh-tht model 
## before thinning
opar <- par( mfcol=c(2,2 ) )
plot( lm( sample.23$plants$tht ~ sample.23$plants$dbh ) )
par( opar )

## thin the stand to capture mortality, redistribute growth, 
## and possibly generate revenue
## Proportional thin for selected tree species, does not remove shrubs
sample.23.t1 <- thin( sample.23, 
	     control=list(type=1, target=50.0, target.sp="DF" ) )
print( sample.23.t1 )

## Proportional thin across all tree species
sample.23.t2 <- thin( sample.23, 
	     control=list(type=2, target=50.0 ) )
print( sample.23.t2 )

## Thin from below, by dbh, all species
sample.23.t3 <- thin( sample.23, 
	     control=list(type=3, target=50.0 ) )
print( sample.23.t3 )

## Thin from below, by dbh for species "PM"
sample.23.t4 <- thin( sample.23, 
	     control=list(type=4, target=50.0, target.sp="PM" ) )
print( sample.23.t4 )

## plot the diagnostics from the fit a linear dbh-tht model
## after proportional thinning
opar <- par( mfcol=c(2,2 ) )
plot( lm( sample.23.t2$plants$tht ~ sample.23.t2$plants$dbh ) )
par( opar )

## print the differences, by species
print( sp.sums( sample.23.t4 ) - sp.sums( sample.23 ) )

}

\keyword{models}


