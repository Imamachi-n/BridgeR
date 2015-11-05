\name{BridgeRHalfLifeDifferenceHist}
\alias{BridgeRHalfLifeDifferenceHist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw histgram to compare RNA half-life difference between two conditions.
}
\description{
Draw histgram to compare RNA half-life difference between two conditions.
}
\usage{
BridgeRHalfLifeDifferenceHist(filename = "BridgeR_4_half-life_calculation.txt", group, hour, ComparisonFile, InforColumn = 4, BinwidthFig = 0.01, OutputFig = "BridgeR_5_Half-life_difference_Histgram")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
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
  \item{BinwidthFig}{
Number(Integer or Float)
}
  \item{OutputFig}{
File path/name
}
}
\details{
Draw histgram to compare RNA half-life difference between two conditions.
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
#inputfile <- "BridgeR_5_HalfLife_calculation.txt"
#outputfig <- "BridgeR_5_HalfLife_Difference_LineGraph"
#group <- c("Control","knockdown1","knockdown2")
#hour <- c(0,1,2,4,8,12)
#compfile <- c("Control","Knockdown1")

#BridgeRHalfLifeDifferenceHist(filename=inputfile, group=group, hour=hour, ComparisonFile=compfile, OutputFig=outputfig)

## The function is currently defined as
function (filename = "BridgeR_4_half-life_calculation.txt", group, 
    hour, ComparisonFile, InforColumn = 4, BinwidthFig = 0.01, 
    OutputFig = "BridgeR_5_Half-life_difference_Histgram") 
{
    time_points <- length(hour)
    group_number <- length(group)
    comp_file_number <- NULL
    for (a in 1:length(ComparisonFile)) {
        comp_file_number <- append(comp_file_number, which(group == 
            ComparisonFile[a]))
    }
    input_file <- fread(filename, header = T)
    figfile <- paste(OutputFig, "_", group[comp_file_number[1]], 
        "_vs_", group[comp_file_number[2]], ".png", sep = "")
    png(filename = figfile, width = 1200, height = 1200)
    half_life_column_1 <- comp_file_number[1] * (time_points + 
        InforColumn + 8)
    half_1 <- input_file[[half_life_column_1]]
    half_life_column_2 <- comp_file_number[2] * (time_points + 
        InforColumn + 8)
    half_2 <- input_file[[half_life_column_2]]
    div_half <- log2(half_2/half_1)
    print(summary(half_1))
    print(summary(half_2))
    plot_data <- data.frame(div_half)
    p.scatter <- ggplot()
    p.scatter <- p.scatter + layer(data = plot_data, mapping = aes(x = div_half), 
        geom = "freqpoly", binwidth = BinwidthFig, size = 1.2)
    p.scatter <- p.scatter + xlim(min(plot_data$div_half), max(plot_data$div_half))
    p.scatter <- p.scatter + ggtitle("Half-life difference")
    name_xlab <- paste("Relative half-life(", group[comp_file_number[1]], 
        "_vs_", group[comp_file_number[2]], ")", sep = "")
    p.scatter <- p.scatter + xlab(name_xlab)
    p.scatter <- p.scatter + ylab("Transcripts #")
    plot(p.scatter)
    dev.off()
    plot.new()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
