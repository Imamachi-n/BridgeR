\name{BridgeRHalfLifeDistribution}
\alias{BridgeRHalfLifeDistribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw boxplot to compare RNA half-life distribution between two conditions.
}
\description{
Draw boxplot to compare RNA half-life distribution between two conditions.
}
\usage{
BridgeRHalfLifeDistribution(filename = "BridgeR_4_half-life_calculation.txt", group, hour, ComparisonFile, InforColumn = 4, OutputFig = "BridgeR_5_Half-life_distribution")
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
  \item{OutputFig}{
File path/name
}
}
\details{
Draw boxplot to compare RNA half-life distribution between two conditions.
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
#outputfig <- "BridgeR_5_HalfLife_Distribution_LineGraph"
#group <- c("Control","knockdown1","knockdown2")
#hour <- c(0,1,2,4,8,12)
#compfile <- c("Control","Knockdown1")

#BridgeRHalfLifeDistribution(filename=inputfile, group=group, hour=hour, ComparisonFile=compfile, OutputFig=outputfig)

## The function is currently defined as
function (filename = "BridgeR_4_half-life_calculation.txt", group, 
    hour, ComparisonFile, InforColumn = 4, OutputFig = "BridgeR_5_Half-life_distribution") 
{
    time_points <- length(hour)
    group_number <- length(group)
    sample_number <- length(ComparisonFile)
    comp_file_number <- NULL
    figfile <- OutputFig
    for (a in 1:sample_number) {
        comp_file_number <- append(comp_file_number, which(group == 
            ComparisonFile[a]))
        figfile <- paste(figfile, "_", ComparisonFile[a], sep = "")
    }
    figfile <- paste(figfile, ".png", sep = "")
    input_file <- fread(filename, header = T)
    png(filename = figfile, width = 1200, height = 1200)
    half_life_fig <- NULL
    for (x in 1:sample_number) {
        data_file <- NULL
        half_life_column <- comp_file_number[x] * (time_points + 
            InforColumn + 8)
        half_life_data <- input_file[[half_life_column]]
        half_life_data <- half_life_data[!is.na(half_life_data)]
        if (x == 1) {
            Sample <- as.factor(rep(ComparisonFile[x], length(half_life_data)))
            half_life_fig <- data.frame(half_life_data, Sample)
        }
        else {
            Sample <- as.factor(rep(ComparisonFile[x], length(half_life_data)))
            half_life_fig <- rbind(half_life_fig, data.frame(half_life_data, 
                Sample))
        }
    }
    p.scatter <- ggplot()
    p.scatter <- p.scatter + layer(data = half_life_fig, mapping = aes(x = half_life_data, 
        colour = Sample), geom = "freqpoly", binwidth = 0.1, 
        size = 1.2, alpha = 0.5)
    p.scatter <- p.scatter + xlim(0, 25)
    p.scatter <- p.scatter + ggtitle("Half-life distribution")
    p.scatter <- p.scatter + xlab("half-life")
    p.scatter <- p.scatter + ylab("Transcripts #")
    plot(p.scatter)
    dev.off()
    plot.new()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line
