\name{BridgeRHalfLifeComparison}
\alias{BridgeRHalfLifeComparison}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw scattered plot to compare RNA half-life between two conditions.
}
\description{
Draw scattered plot to compare RNA half-life between two conditions.
}
\usage{
BridgeRHalfLifeComparison(filename = "BridgeR_5_half-life_calculation.txt", group, hour, ComparisonFile, InforColumn = 4, LogScale = F, Calibration=F, OutputFig = "BridgeR_6_Half-life_comparison")
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
  \item{LogScale}{
Bool(True or False)
}
  \item{Calibration}{
Bool(True or False)
}
  \item{OutputFig}{
File path/name
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
#outputfig <- "BridgeR_5_HalfLife_Comparison_ScatteredPlot"
#group <- c("Control","knockdown1","knockdown2")
#hour <- c(0,1,2,4,8,12)
#compfile <- c("Control","Knockdown1")

#BridgeRHalfLifeComparison(filename=inputfile, group=group, hour=hour, ComparisonFile=compfile, OutputFig=outputfig)

## The function is currently defined as
function (filename = "BridgeR_5_half-life_calculation.txt", group, 
    hour, ComparisonFile, InforColumn = 4, LogScale = F, OutputFig = "BridgeR_6_Half-life_comparison") 
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
    gene_number <- length(half_1)
    half_1_fig <- NULL
    half_2_fig <- NULL
    factor_fig <- NULL
    for (x in 1:gene_number) {
        if (is.na(half_1[x]) || is.na(half_2[x])) {
            next
        }
        half_1_fig <- append(half_1_fig, half_1[x])
        half_2_fig <- append(half_2_fig, half_2[x])
        div <- half_1[x]/half_2[x]
        if (div <= 0.5) {
            factor_fig <- append(factor_fig, 2)
        }
        else if (div >= 2) {
            factor_fig <- append(factor_fig, 1)
        }
        else {
            factor_fig <- append(factor_fig, 0)
        }
    }
    if (LogScale) {
        half_1_fig <- log2(half_1_fig)
        half_2_fig <- log2(half_2_fig)
    }
    up_genes <- length(which(factor_fig == 1))
    down_genes <- length(which(factor_fig == 2))
    plot_data <- data.frame(half_1_fig, half_2_fig, factor_fig)
    print_out <- paste("Plotted: ", length(plot_data[, 1]), " genes", 
        sep = "")
    print_out2 <- paste("At least 2-fold upregulated: ", down_genes, 
        " genes", sep = "")
    print_out3 <- paste("At least 2-fold downregulated: ", up_genes, 
        " genes", sep = "")
    print(print_out)
    print(print_out2)
    print(print_out3)
    p.scatter <- ggplot()
    p.scatter <- p.scatter + layer(data = plot_data, mapping = aes(x = half_1_fig, 
        y = half_2_fig, colour = factor(factor_fig)), geom = "point", 
        size = 2.5, alpha = 0.3)
    p.scatter <- p.scatter + layer(data = plot_data, mapping = aes(x = half_1_fig, 
        y = half_2_fig), geom = "smooth", geom_params = list(color = "blue", 
        size = 1.2), stat = "smooth", stat_params = list(method = "lm", 
        se = F))
    p.scatter <- p.scatter + xlim(0, max(plot_data$half_1_fig)) + 
        ylim(0, max(plot_data$half_2_fig))
    p.scatter <- p.scatter + ggtitle("Half-life comparison")
    name_xlab <- paste(group[comp_file_number[1]], " (Time)", 
        sep = "")
    name_ylab <- paste(group[comp_file_number[2]], " (Time)", 
        sep = "")
    p.scatter <- p.scatter + xlab(name_xlab)
    p.scatter <- p.scatter + ylab(name_ylab)
    p.scatter <- p.scatter + theme(legend.position = "none")
    p.scatter <- p.scatter + scale_colour_manual(values = c("black", 
        "blue", "red"))
    plot(p.scatter)
    dev.off()
    plot.new()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line
