\name{BridgeRCustom}
\alias{BridgeRCustom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Custom BRIC-seq data analysis pipeline
}
\description{
RNA half-life for each gene is calculated using your custom normalization factor such as the RPKM of Luciferase.
}
\usage{
BridgeRCustom(YourNormFactor, 
              SelectNormFactor=F,
              InputFiles, 
              InforColumn=4,
              group, 
              hour, 
              RPKMcutoff=0.1,
              CutoffDataPointNumber = 4,
              CutoffDataPoint1 = c(1,2),
              CutoffDataPoint2 = c(8,12),
              ThresholdHalfLife = c(8,12),
              CutoffRelExp=0.001,
              ModelMode="R2_selection")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{YourNormFactor}{
File pass: tab-delimited text about normalization factor for each sample
}
  \item{SelectNormFactor}{
Bool(True or False)
}
  \item{InputFiles}{
Vector
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
  \item{RPKMcutoff}{
Float
}
  \item{CutoffDataPointNumber}{
Integer
}
  \item{CutoffDataPoint1}{
Vector(Integer or float)
}
  \item{CutoffDataPoint2}{
Vector(Integer or float)
}
  \item{ThresholdHalfLife}{
Vector(Integer or float)
}
  \item{CutoffRelExp}{
Float
}
  \item{ModelMode}{
Bool(True or False)
}
}
\details{
RNA half-life for each gene is calculated Using your custom normalization factor such as the RPKM of Luciferase.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
text files and fig files
}
\references{
https://github.com/Naoto-Imamachi/BRIC-seq_data_analysis/tree/master/BridgeR
}
\author{
Naoto Imamachi
}
\note{
2015-11-05
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
#mynf <- "Luciferase2_Rel_RPKM.txt"
#files <- c("Control_genes_RefSeq_result_mRNA.fpkm_table",
#           "Knockdown1_genes_RefSeq_result_mRNA.fpkm_table",
#           "Knockdown2_genes_RefSeq_result_mRNA.fpkm_table")
#group <- c("Control","knockdown1","knockdown2")
#hour <- c(0,1,2,4,8,12)

#BridgeRCustom(YourNormFactor=mynf, InputFiles=files, group=group, hour=hour)

## The function is currently defined as
function (YourNormFactor, 
                          InputFiles, 
                          SelectNormFactor=F,
                          InforColumn=4, 
                          group, 
                          hour, 
                          RPKMcutoff=0.1,
                          CutoffDataPointNumber = 4,
                          CutoffDataPoint1 = c(1,2),
                          CutoffDataPoint2 = c(8,12),
                          ThresholdHalfLife = c(8,12),
                          CutoffRelExp=0.001,
                          ModelMode="R2_selection") 
{
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line
