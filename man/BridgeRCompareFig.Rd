\name{BridgeRCompareFig}
\alias{BridgeRCompareFig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw scattered plot to compare RNA half-life between two conditions.
}
\description{
Draw scattered plot to compare RNA half-life between two conditions.
}
\usage{
BridgeRCompareFig(InputFile="BridgeR_6_HalfLife_Pvalue_estimation.txt",
                              InforColumn = 4,
                              group, 
                              hour, 
                              ComparisonFile,
                              ModelMode="R2_selection")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{InputFile}{
File path/name
}
  \item{InforColumn}{
Integer
}
  \item{group}{
Vector(string)
}
  \item{hour}{
Vector(number)
}
  \item{ComparisonFile}{
Vector(string)
}
  \item{ModelMode}{
Mode
}
}
\details{
Draw scattered plot to compare RNA half-life between two conditions.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
fig file
}
\references{
https://github.com/Naoto-Imamachi/BRIC-seq_data_analysis/tree/master/BridgeR
}
\author{
Naoto Imamachi
}
\note{
2015-11-25
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (InputFile="BridgeR_6_HalfLife_Pvalue_estimation.txt",
                              InforColumn = 4,
                              group, 
                              hour, 
                              ComparisonFile,
                              ModelMode="R2_selection")
{
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line
