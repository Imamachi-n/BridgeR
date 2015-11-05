\name{test_q}
\alias{test_q}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
estimate percentile data
}
\description{
estimate percentile data
}
\usage{
test_q(x, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
exp_data
}
  \item{y}{
hour_label
}
}
\details{
estimate percentile data
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
data_table
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
#q_data <- test_q(log10(exp_data[[x]]),hour_label[x])

## The function is currently defined as
function (x, y) 
{
    q_99 <- as.vector(quantile(x, prob = 0.99, na.rm = T))
    q_95 <- as.vector(quantile(x, prob = 0.95, na.rm = T))
    q_90 <- as.vector(quantile(x, prob = 0.9, na.rm = T))
    q_80 <- as.vector(quantile(x, prob = 0.8, na.rm = T))
    q_70 <- as.vector(quantile(x, prob = 0.7, na.rm = T))
    q_60 <- as.vector(quantile(x, prob = 0.6, na.rm = T))
    q_50 <- as.vector(quantile(x, prob = 0.5, na.rm = T))
    q_40 <- as.vector(quantile(x, prob = 0.4, na.rm = T))
    q_30 <- as.vector(quantile(x, prob = 0.3, na.rm = T))
    q_20 <- as.vector(quantile(x, prob = 0.2, na.rm = T))
    q_10 <- as.vector(quantile(x, prob = 0.1, na.rm = T))
    q_5 <- as.vector(quantile(x, prob = 0.05, na.rm = T))
    q_1 <- as.vector(quantile(x, prob = 0.01, na.rm = T))
    vec <- c(q_99, q_95, q_90, q_80, q_70, q_60, q_50, q_40, 
        q_30, q_20, q_10, q_5, q_1)
    factor_label <- c("99\%", "95\%", "90\%", "80\%", "70\%", "60\%", 
        "50\%", "40\%", "30\%", "20\%", "10\%", "05\%", "01\%")
    label <- rep(y, 13)
    q_table <- data.frame(name = label, q = vec, factor = factor_label)
    return(q_table)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line
