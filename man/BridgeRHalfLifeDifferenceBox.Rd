\name{BridgeRHalfLifeDifferenceBox}
\alias{BridgeRHalfLifeDifferenceBox}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw boxplot to compare RNA half-life distribution between two conditions.
}
\description{
Draw boxplot to compare RNA half-life distribution between two conditions.
}
\usage{
BridgeRHalfLifeDifferenceBox(filename = "BridgeR_4_half-life_calculation.txt", group, hour, ComparisonFile, InforColumn = 4, OutputFig = "BridgeR_5_Half-life_difference_Boxplot")
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
#outputfig <- "BridgeR_5_HalfLife_Comparison_BoxPlot"
#group <- c("Control","knockdown1","knockdown2")
#hour <- c(0,1,2,4,8,12)
#compfile <- c("Control","Knockdown1")

#BridgeRHalfLifeDifferenceBox(filename=inputfile, group=group, hour=hour, ComparisonFile=compfile, OutputFig=outputfig)

## The function is currently defined as
function (filename = "BridgeR_4_half-life_calculation.txt", group, 
    hour, ComparisonFile, InforColumn = 4, OutputFig = "BridgeR_5_Half-life_difference_Boxplot") 
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
    fig_width <- 200 * length(ComparisonFile)
    png(filename = figfile, width = fig_width, height = 1200)
    half_life_column_1 <- comp_file_number[1] * (time_points + 
        InforColumn + 8)
    half_1 <- input_file[[half_life_column_1]]
    half_1 <- half_1[!is.na(half_1)]
    half_1_number <- length(half_1)
    half_1_label <- rep(ComparisonFile[1], half_1_number)
    plot_data_1 <- data.frame(half_1, half_1_label)
    colnames(plot_data_1) <- c("half_data", "label")
    half_life_column_2 <- comp_file_number[2] * (time_points + 
        InforColumn + 8)
    half_2 <- input_file[[half_life_column_2]]
    half_2 <- half_2[!is.na(half_2)]
    half_2_number <- length(half_2)
    half_2_label <- rep(ComparisonFile[2], half_2_number)
    plot_data_2 <- data.frame(half_2, half_2_label)
    colnames(plot_data_2) <- c("half_data", "label")
    print(summary(half_1))
    print(summary(half_2))
    plot_data <- rbind(plot_data_1, plot_data_2)
    p.boxplot <- ggplot()
    p.boxplot <- p.boxplot + layer(data = plot_data, mapping = aes(x = label, 
        y = half_data), geom = "boxplot")
    p.boxplot <- p.boxplot + ylim(0, 24)
    p.boxplot <- p.boxplot + ggtitle("Half-life difference")
    p.boxplot <- p.boxplot + ylab("Half-life")
    plot(p.boxplot)
    dev.off()
    plot.new()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line
