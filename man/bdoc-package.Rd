\name{bdoc-package}
\alias{bdoc-package}
\alias{bdoc-packages}
\docType{package}
\title{
Bayesian Discrete Ordered Classification of DNA Barcodes
}
\description{
This package contains the "bdoc" function that will classify DNA barcodes in a test data set
to a species in the reference data set of DNA barcodes.  This function will produce an assignment probability
together with plots of the posterior probabilities of belonging to any of the species in the reference data set.  These plots can be used to determine if a test barcode comes 
from a species not contained in the reference data set. 
}
\details{
\tabular{ll}{
Package: \tab bdoc\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2009-08-31\cr
License: \tab GPL (version 2 or later)\cr
LazyLoad: \tab yes\cr
}
The object "traindata" should be of type data.frame and contain the species-level identification in the second column.  This column should be named "species" in order for the function to construct the correct conditional probabilities.  The object "testdata" should be of type data.frame with only the barcodes of the DNA sequences to be classified.  All the rest of the options have default values that are strongly recommended.  Plots of the posterior probabilities for each of the barcodes in the test data set are constructed and saved with format plot.file to the current R directory.  See example below. 
}
\author{

Michael Anderson and Suzanne Dubnicka

Maintainer: Michael Anderson <mpa@ksu.edu>

}
\references{
Hebert, P., A. Cywinska, S. Ball, and J. deWaard (2003). Biological identifications through DNA barcodes. Proc. R. Soc. Lond. (B) 270, 313-322.
}

\keyword{ classif }
\seealso{
\code{\link{data.frame}}
}
\examples{
data(battraindata1)
data(battestdata1)

traindata<-battraindata1
#battraindata1 contains the genus (column 1) and species (column 2) 
#barcode information for 758 bats representing 96 unique species.
#The length of each barcode is 659 nucleotides long.
 
testdata<-battestdata1
#battetdata1 contains the genus (column 1) and species (column 2) 
#barcode information for 82 bats that were held out of battraindata1.
#The length of each barcode is 659 nucleotides long and to classify,
#the first two columns need to be removed as these will usually not 
#be known.

\donttest{
result<-bdoc(traindata,testdata[,-c(1:2)])  #after this executes, plots of type
                                            #plot.file names "seq1", "seq2", 
                                            #and so on can be found in the 
                                            #folder identified by getwd().

result$priors

result$species.class           #gives the matrix of species assignments
                               #and probabilities.

result$posteriors[[1]]$post }  #gives the matrix of posterior probabilities
                               #at each position for barcode 1.  Change 
                               #posteriors[[1]] to posteriors[[2]] for the 
                               #posteriors for barcode 2, etc.



}
