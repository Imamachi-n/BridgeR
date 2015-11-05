\name{BridgeRNormalizationFactorsHK}
\alias{BridgeRNormalizationFactorsHK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate normalization factor using house-keeping genes
}
\description{
Estimate normalization factor using house-keeping genes
}
\usage{
BridgeRNormalizationFactorsHK(InputFile, group, hour, InforColumn = 4, InforHKGenes = 2, HKGenes = c("GAPDH", "PGK1", "PPIA", "ENO1", "ATP5B", "ALDOA"), OutputFile = "BridgeR_3B_Normalizaion_factor_from_House_keeping_genes.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{InputFile}{
File path/name
}
  \item{group}{
Vector(string)
}
  \item{hour}{
Vector(number)
}
  \item{InforColumn}{
Integer
}
  \item{InforHKGenes}{
Integer
}
  \item{HKGenes}{
Vector(String)
}
  \item{OutputFile}{
File path/name
}
}
\details{
Estimate normalization factor using house-keeping genes
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
text file
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
#inputfile <- "BridgeR_1_Relative_expression_dataset.txt"
#outputfile <- "BridgeR_3_Normalizaion_factor_HouseKeepingGenes.txt"
#group <- c("Control","knockdown1","knockdown2")
#hour <- c(0,1,2,4,8,12)

#BridgeRNormalizationFactorsHK(InputFile=inputfile, group=group, hour=hour, OutputFile=outputfile)

## The function is currently defined as
function (InputFile, group, hour, InforColumn = 4, InforHKGenes = 2, 
    HKGenes = c("GAPDH", "PGK1", "PPIA", "ENO1", "ATP5B", "ALDOA"), 
    OutputFile = "BridgeR_3B_Normalizaion_factor_from_House_keeping_genes.txt") 
{
    group_number <- length(group)
    time_points <- length(hour)
    input_file <- fread(InputFile, header = T)
    hour_label <- NULL
    for (x in hour) {
        label <- x
        if (x < 10) {
            label <- paste("0", x, sep = "")
        }
        hour_label <- append(hour_label, paste("T", label, "_", 
            a, sep = ""))
    }
    cat("Sample", hour_label, sep = "\t", file = OutputFile)
    cat("\n", file = OutputFile, append = T)
    for (a in 1:group_number) {
        HKgenes_infor <- input_file[[InforHKGenes]]
        HKgenes_infor_index <- NULL
        for (x in 1:length(HKGenes)) {
            HKgenes_infor_index <- append(HKgenes_infor_index, 
                which(HKgenes_infor == HKGenes[x]))
        }
        infor_st <- 1 + (a - 1) * (time_points + InforColumn)
        infor_ed <- (InforColumn) * a + (a - 1) * time_points
        exp_st <- infor_ed + 1
        exp_ed <- infor_ed + time_points
        HKGenes_raw_data <- NULL
        for (x in 1:length(HKgenes_infor_index)) {
            if (x == 1) {
                HKGenes_raw_data <- input_file[HKgenes_infor_index[x], 
                  exp_st:exp_ed, with = F]
            }
            else {
                HKGenes_raw_data <- rbind(HKGenes_raw_data, input_file[HKgenes_infor_index[x], 
                  exp_st:exp_ed, with = F])
            }
        }
        nf <- NULL
        for (x in 1:length(HKGenes_raw_data)) {
            if (x == 1) {
                nf <- as.vector(as.matrix(HKGenes_raw_data[x, 
                  ]))
            }
            else {
                nf <- HKGenes_raw_data[x, ] * nf
            }
        }
        nf <- nf^(1/length(HKGenes))
        cat(group[a], nf, sep = "\t", file = OutputFile, append = T)
        cat("\n", file = OutputFile, append = T)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line
