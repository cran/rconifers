
\name{sp.sums}
\alias{sp.sums}
\title{A summary of a CONIFERS sample.data object by species}

\description{ This function returns a \code{\link{data.frame}} object
   that contain common summaries by species for
   \code{\link{sample.data}} object.  
}

\usage{
   sp.sums( x )
}

\arguments{
   \item{x}{an object of class \code{\link{sample.data}}.}
}

\format{

The \code{\link{data.frame}} object returned from sp.sums contains the
following columns for each species:

\describe{
\item{expf}{The total plants per acre.}
\item{tht}{The mean total height.}
\item{ba}{The total basal area for all plants.}
\item{qmd}{The quadratic mean diameter for all trees. See Details.}
}
}

\details{ This function returns a \code{\link{data.frame}} object that
	contains species level summary information. It is intended
	demonstration only and users are encouraged to examine the
	source code and modify. All results include all stems and not
	only those above breast height.
}

\references{

Ritchie, M.W. 2008. User's Guide and Help System for CONIFERS: A Simulator for Young Conifer Plantations Version
4.10. See \url{http://www.fs.fed.us/psw/programs/ecology_of_western_forests/projects/conifers/}

}

\author{Jeff D. Hamann \email{jeff.hamann@forestinformatics.com},\cr
	     Martin W. Ritchie \email{mritchie@fs.fed.us} }

\seealso{    \code{\link{calc.max.sdi}},
	     \code{\link{sample.data}},
	     \code{\link{set.species.map}},		  
	     \code{\link{set.variant}},		  
 	     \code{\link{smc}}, 
 	     \code{\link{summary.sample.data}},
	     \code{\link{swo}}
}

\examples{
library( rconifers )

#Example 1 Using SMC
## set the variant to the SMC variant
set.species.map(set.variant( 1 ))

## grow the data that was originally swo in the smc variant
# load and display CONIFERS example plots
data( plots.smc )
print( plots.smc )

# load and display CONIFERS example plants
data( plants.smc )
print( plants.smc )

# create the sample.data list object
sample.smc.3 <- list( plots=plots.smc, plants=plants.smc, age=3, x0=0.0 )
class(sample.smc.3)  <- "sample.data"

print( sp.sums( sample.smc.3 ) )
smry<-sp.sums( sample.smc.3 )
dfqmd<-smry["DF","qmd"]

#Example 2 Using SWO
## set the variant to the SWO variant
set.species.map(set.variant( 0 ) )

## grow the swo data
# load and display CONIFERS example plots
data( plots.swo )
print( plots.swo )

# load and display CONIFERS example plants
data( plants.swo )
print( plants.swo )

# create the sample.data list object
sample.swo.3 <- list( plots=plots.swo, plants=plants.swo, age=3, x0=0.0 )
class(sample.swo.3)  <- "sample.data"

# print the maximum stand density index for the current settings
print( calc.max.sdi( sample.swo.3 ) )

# print a summary of the sample
print( sample.swo.3 )

# now, project the sample forward for 20 years
# with all of the options turned off
sample.swo.23 <- project( sample.swo.3,
	  20, 
	  control=list(rand.err=0,rand.seed=0,endemic.mort=0,sdi.mort=0))

## print the projected summaries
print( sample.swo.23 )

## thin the stand
## Proportional thin for selected tree species, does not remove shrubs
sample.swo.23.t1 <- thin( sample.swo.23,
	     control=list(type=1, target=50.0, target.sp="DF" ) )
print( sample.swo.23.t1 )

## Proportional thin across all tree species
sample.swo.23.t2 <- thin( sample.swo.23,
	     control=list(type=2, target=50.0 ) )
print( sample.swo.23.t2 )

## Thin from below, by dbh, all species
sample.swo.23.t3 <- thin( sample.swo.23,
	     control=list(type=3, target=50.0 ) )
print( sample.swo.23.t3 )

## Thin from below, by dbh for species "PM"
sample.swo.23.t4 <- thin( sample.swo.23,
	     control=list(type=4, target=50.0, target.sp="PM" ) )
print( sample.swo.23.t4 )

## print the differences, by species
print( sp.sums( sample.swo.23.t4 ) - sp.sums( sample.swo.23 ) )

## generate a plot of the species summaries
opar <- par(mfrow=c(2,2))

##hist( sample.swo.23$plants$dbh, xlim=c(0,50), ylim=c(0,50), main="Diameter Distribution" )
##hist( sample.swo.23$plants$tht, xlim=c(0,100), ylim=c(0,110), main="Height Distribution" )

hist( sample.swo.23$plants$dbh, main="Diameter Distribution" )
hist( sample.swo.23$plants$tht, main="Height Distribution" )

ba <- sp.sums(sample.swo.23)$ba
names(ba) <- rownames( sp.sums(sample.swo.23) )
pie( ba[ba>0], main="Basal Area" )
plot( sample.swo.23$plants$tht ~ sample.swo.23$plants$dbh, ylim=c(0,220), xlim=c(0,40), main="Height vs. Basal Diameter" )

par(opar)


}

\keyword{models}



