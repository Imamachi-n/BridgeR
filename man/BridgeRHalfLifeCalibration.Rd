\name{BridgeRHalfLifeCalibration}
\alias{BridgeRHalfLifeCalibration}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calibrate RNA half-life.
}
\description{
Calibrate RNA half-life.
}
\usage{
BridgeRHalfLifeCalibration(InputFile = "BridgeR_5_HalfLife_calculation.txt", group, hour, ComparisonFile, InforColumn = 4, OutputFile = "BridgeR_6_HalfLife_Calibration")
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
  \item{ComparisonFile}{
Vector(string)
}
  \item{InforColumn}{
Integer
}
  \item{OutputFile}{
File path/name
}
}
\details{
Calibrate RNA half-life.
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
2015-11-06
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

#BridgeRHalfLifeComparison(filename=inputfile, group=group, hour=hour, ComparisonFile=compfile, OutputFig=outputfig)

## The function is currently defined as
function(filename = "BridgeR_5_HalfLife_calculation.txt", group,
    hour, ComparisonFile, InforColumn = 4, OutputFile = "BridgeR_7_HalfLife_Calibration")
{
    #ComparisonFile: The length of vector must be 2 => c("Control","Knockdown")
    ###Prepare_file_infor###
    time_points <- length(hour)
    group_number <- length(group)
    comp_file_number <- NULL
    for(a in 1:length(ComparisonFile)){
        comp_file_number <- append(comp_file_number, which(group == ComparisonFile[a]))
    }
    
    input_file <- fread(filename, header=T)
    #figfile <- paste(OutputFig,"_",group[comp_file_number[1]],"_vs_",group[comp_file_number[2]],".png", sep="")
    #png(filename=figfile,width = 1200, height = 1200)
    
    ###Plot_Half-life_comparison###
    gene_infor <- input_file[,1:InforColumn,with=F]
    half_life_column_1 <- comp_file_number[1]*(time_points + InforColumn + 8) #number
    half_1 <- as.numeric(input_file[[half_life_column_1]])
    half_life_column_2 <- comp_file_number[2]*(time_points + InforColumn + 8) #number
    half_2 <- as.numeric(input_file[[half_life_column_2]])
    half_data <- data.table(half1=half_1,half2=half_2)
    
    test <- lm(half2 ~ half1 + 0, data=half_data)
    coef <- as.numeric(test$coefficients)
    adjusted_half_2 <- half_2/coef
    
    adjusted_half_data <- data.table(half1=half_1,half2=adjusted_half_2)
    adjusted_half_data <- cbind(gene_infor, adjusted_half_data)
    
    result_file <- paste(OutputFile,"_",group[comp_file_number[1]],"_vs_",group[comp_file_number[2]],".txt",sep="")
    write.table(adjusted_half_data, file=result_file,sep="\t",row.names=FALSE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line
