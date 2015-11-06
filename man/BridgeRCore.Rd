\name{BridgeRCore}
\alias{BridgeRCore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Default BRIC-seq data analysis pipeline
}
\description{
RNA half-life for each gene is calculated using 99-97.5 percentile normalization.
}
\usage{
BridgeRCore(InputFiles, InforColumn = 4, group, hour, RPKMcutoff = 0.1, RelRPKMFig = F, SelectNormFactor = T, CutoffRelExp = 0.1, CutoffDataPoint = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{InputFiles}{
Vector
}
  \item{InforColumn}{
Integer
}
  \item{group}{
Vector
}
  \item{hour}{
Vector(number)
}
  \item{RPKMcutoff}{
Float
}
  \item{RelRPKMFig}{
Bool(True or False)
}
  \item{SelectNormFactor}{
Bool(True or False)
}
  \item{CutoffRelExp}{
Float
}
  \item{CutoffDataPoint}{
Integer
}
}
\details{
RNA half-life for each gene is calculated Using 99-97.5 percentile normalization.
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
#files <- c("Control_genes_RefSeq_result_mRNA.fpkm_table",
#           "Knockdown1_genes_RefSeq_result_mRNA.fpkm_table",
#           "Knockdown2_genes_RefSeq_result_mRNA.fpkm_table")
#group <- c("Control","knockdown1","knockdown2")
#hour <- c(0,1,2,4,8,12)

#BridgeRCore(InputFiles=files, group=group, hour=hour)

## The function is currently defined as
function (InputFiles, InforColumn = 4, group, hour, RPKMcutoff = 0.1, 
    RelRPKMFig = F, SelectNormFactor = T, CutoffRelExp = 0.1, 
    CutoffDataPoint = 3) 
{
    library(data.table)
    library(ggplot2)
    BridgeRDataSetFromCuffnorm(CuffnormFiles = InputFiles, group = group, 
        hour = hour, cutoff = RPKMcutoff, InforColumn = InforColumn, 
        OutputFile = "BridgeR_1_Relative_expression_dataset.txt")
    BridgeRDatasetChecker(InputFile = "BridgeR_1_Relative_expression_dataset.txt", 
        group = group, hour = hour, InforColumn = InforColumn, 
        OutputFile = "BridgeR_2_Relative_RPKM_distribution")
    BridgeRNormalizationFactors(InputFile = "BridgeR_1_Relative_expression_dataset.txt", 
        group = group, hour = hour, InforColumn = InforColumn, 
        YMin = -2, YMax = 2, MakeFig = RelRPKMFig, figname = "BridgeR_3_Normalizaion_factor", 
        nfname = "BridgeR_3_Normalizaion_factor")
    BridgeRNormalization(filename = "BridgeR_1_Relative_expression_dataset.txt", 
        group = group, hour = hour, InforColumn = InforColumn, 
        SelectNormFactor = SelectNormFactor, NormFactor = "BridgeR_3_Normalizaion_factor", 
        OutputFile = "BridgeR_4_Normalized_expression_dataset.txt")
    BridgeRDatasetChecker(InputFile = "BridgeR_4_Normalized_expression_dataset.txt", 
        group = group, hour = hour, InforColumn = InforColumn, 
        OutputFile = "BridgeR_4B_Normalized_RPKM_distribution")
    BridgeRHalfLifeCalculation(filename = "BridgeR_4_Normalized_expression_dataset.txt", 
        group = group, hour = hour, InforColumn = InforColumn, 
        CutoffRelExp = CutoffRelExp, CutoffDataPoint = CutoffDataPoint, 
        OutputFile = "BridgeR_5_HalfLife_calculation.txt")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line