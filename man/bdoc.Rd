\name{bdoc}
\alias{bdoc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Bayesian Discrete Ordered Classification of DNA Barcodes }
\description{
  This package contains the "bdoc" function that will classify DNA barcodes in a test data set
to a species in the reference data set of DNA barcodes.  This function will produce an assignment probability
together with plots of the posterior probabilities of belonging to any of the species in the reference data set.  These plots can be used to determine if a test barcode comes 
from a species not contained in the reference data set. 
}
\usage{
bdoc(traindata, testdata, delta = 9.7e-08, epsilon = 0.2, priors = "equal", stoppingrule = TRUE, impute = 1, plot.file = "pdf")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{traindata}{ Contains a training dataset of type data.frame with the first column reserved for an ID (possibly genus), the second column reserved for the actual species name, and the remainder of the columns containing the nucleotide sequence of the DNA barcode.}
  \item{testdata}{ Contains a test data set of type data.frame of the DNA barcodes to be classified.  Note, column 1 should contain the first nucleotide, position 2 the second, and so on.}
  \item{delta}{Scalar value between 0 and 0.1 used to adjust the conditional probabilities.}
  \item{epsilon}{ Scalar value between 0 and 1 used to adjust the posterior probability calculations.}
  \item{priors}{The prior probabilities to be used.  This can be a vector of probabilities for each species in the reference data set (should sum to 1) or any of the following options: "equal" - to use prior probabilities all eqaul to 1/s if s is the number of species in the reference data set; "data" - to use prior probabilities equal to the prevalence of each species in the reference data set; "dir" - to use unequal, arbitrary probabilities generated from a Dirichlet(1,1,...,1) distribution.}
  \item{stoppingrule}{Logical.  By defalut stoppingrule=TRUE, which will terminate the sequential calculation when the posterior probability for a species in the reference data set equals 1.  If set to FALSE, the calculation continues until the end of the barcode is reached.}
  \item{impute}{Imputation method.  If impute=1, the the proportional allocation method will be used.  If impute=2 the majority rule impuation will be performed.}
  \item{plot.file}{ Type of posterior probability plot to be saved to the current directory.  By default, plot.file="pdf", which will save a PDF of the posterior proability plot(s).  Other options include: "jpg" - save a JPEG file of the plot(s); "png" - save a PNG file of the plot(s); "wmf" - save a WMF file of the plot(s); "ps" - save a PS file of the plot(s).}
}
\details{
  The object "traindata" should be of type data.frame and contain the species-level identification in the second column.  This column should be named "species" in order for the function to construct the correct conditional probabilities.  The object "testdata" should be of type data.frame with only the barcodes of the DNA sequences to be classified.  All the rest of the options have default values that are strongly recommended.  Plots of the posterior probabilities for each of the barcodes in the test data set are constructed and saved with format plot.file to the current R directory.  See example below. 
}
\value{
  \item{k }{The total number of barcodes in the test data set.}
  \item{totaltime }{The total time used for all of the barcodes in the test data set.}
  \item{imp }{The time used to impute the missing values.}
  \item{like }{The time used construct and adjust the conditional probabilities.}
  \item{class }{The time used to compute the posteriors and make the species level assignment.}
  \item{delta }{The value used to adjust the conditional probabilities.}
  \item{species.class }{A matrix of the species assignment as well as the probability of assignment for each barcode in the test data set.}
  \item{priors }{A vector of the initial prior probabilities.}
  \item{posteriors }{A list containing: 1.  the species-level assignment for each barcode in the test data set; 2. the matrix of posterior probabilities at each position for each barcode in the test data set.  See the example below.}
  \item{Posterior Probability Plots }{Posterior probability plots are constructed and saved in the format of plot.file to the current R directory named "seq1", "seq2", and so on.}
}
\references{ Hebert, P., A. Cywinska, S. Ball, and J. deWaard (2003). Biological identifications
through DNA barcodes. Proc. R. Soc. Lond. (B) 270, 313-322. }
\author{ Michael Anderson and Suzanne Dubnicka }

\seealso{  \code{\link{data.frame}} }
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

result$posteriors[[1]]$post }   #gives the matrix of posterior probabilities
                               #at each position for barcode 1.  Change 
                               #posteriors[[1]] to posteriors[[2]] for the 
                               #posteriors for barcode 2, etc.




}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ classif }
% __ONLY ONE__ keyword per line
