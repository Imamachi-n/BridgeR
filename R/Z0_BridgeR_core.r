#Title: BridgeR_Core
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-11-05

BridgeRCore <- function(InputFiles, InforColumn=4, group, hour, RPKMcutoff=0.1, RelRPKMFig=F, SelectNormFactor=T, CutoffRelExp=0.1, CutoffDataPoint=3){
    library(data.table)
    library(ggplot2)
    BridgeRDataSetFromCuffnorm(CuffnormFiles=InputFiles, group=group, hour=hour, cutoff=RPKMcutoff, InforColumn=InforColumn, OutputFile="BridgeR_1_Relative_expression_dataset.txt")
    BridgeRDatasetChecker(InputFile="BridgeR_1_Relative_expression_dataset.txt", group=group, hour=hour, InforColumn=InforColumn, OutputFile="BridgeR_2_Relative_RPKM_distribution")
    BridgeRNormalizationFactors(InputFile="BridgeR_1_Relative_expression_dataset.txt",group=group, hour=hour, InforColumn=InforColumn, YMin=-2, YMax=2, MakeFig=RelRPKMFig, figname="BridgeR_3_Normalizaion_factor", nfname="BridgeR_3_Normalizaion_factor")
    BridgeRNormalization(filename="BridgeR_1_Relative_expression_dataset.txt", group=group, hour=hour, InforColumn=InforColumn, SelectNormFactor=SelectNormFactor, NormFactor="BridgeR_3_Normalizaion_factor", OutputFile="BridgeR_4_Normalized_expression_dataset.txt")
    BridgeRDatasetChecker(InputFile="BridgeR_4_Normalized_expression_dataset.txt", group=group, hour=hour, InforColumn=InforColumn, OutputFile="BridgeR_4B_Normalized_RPKM_distribution")
    BridgeRHalfLifeCalculation(filename="BridgeR_4_Normalized_expression_dataset.txt", group=group, hour=hour, InforColumn=InforColumn, CutoffRelExp=CutoffRelExp, CutoffDataPoint=CutoffDataPoint, OutputFile="BridgeR_5_HalfLife_calculation.txt")
}

BridgeRHKGenes <- function(InputFiles, InforColumn=4, InforHKGenes=2, HKGenes=c("GAPDH","PGK1","PPIA","ENO1","ATP5B","ALDOA"), group, hour, RPKMcutoff=0.1, SelectNormFactor=T, CutoffRelExp=0.1, CutoffDataPoint=3){
    library(data.table)
    library(ggplot2)
    BridgeRDataSetFromCuffnorm(CuffnormFiles=InputFiles, group=group, hour=hour, cutoff=RPKMcutoff, InforColumn=InforColumn, OutputFile="BridgeR_1_Relative_expression_dataset.txt")
    BridgeRDatasetChecker(InputFile="BridgeR_1_Relative_expression_dataset.txt", group=group, hour=hour, InforColumn=InforColumn, OutputFile="BridgeR_2_Relative_RPKM_distribution")
    BridgeRNormalizationFactorsHK(InputFile="BridgeR_1_Relative_expression_dataset.txt", group=group, hour=hour, InforColumn=InforColumn, InforHKGenes=InforHKGenes, HKGenes=HKGenes, OutputFile="BridgeR_3_Normalizaion_factor_HouseKeepingGenes.txt")
    BridgeRNormalizationForLuc2(filename="BridgeR_1_Relative_expression_dataset.txt", group=group, hour=hour, InforColumn=InforColumn, NormFactor="BridgeR_3_Normalizaion_factor_HouseKeepingGenes.txt", OutputFile="BridgeR_4_Normalized_expression_dataset.txt")
    BridgeRDatasetChecker(InputFile="BridgeR_4_Normalized_expression_dataset.txt", group=group, hour=hour, InforColumn=InforColumn, OutputFile="BridgeR_4B_Normalized_RPKM_distribution")
    BridgeRHalfLifeCalculation(filename="BridgeR_4_Normalized_expression_dataset.txt", group=group, hour=hour, InforColumn=InforColumn, CutoffRelExp=CutoffRelExp, CutoffDataPoint=CutoffDataPoint, OutputFile="BridgeR_5_HalfLife_calculation.txt")
}

BridgeRCustom <- function(YourNormFactor, InputFiles, InforColumn=4, group, hour, RPKMcutoff=0.1, SelectNormFactor=T, CutoffRelExp=0.1, CutoffDataPoint=3){
    library(data.table)
    library(ggplot2)
    BridgeRDataSetFromCuffnorm(CuffnormFiles=InputFiles, group=group, hour=hour, cutoff=RPKMcutoff, InforColumn=InforColumn, OutputFile="BridgeR_1_Relative_expression_dataset.txt")
    BridgeRDatasetChecker(InputFile="BridgeR_1_Relative_expression_dataset.txt", group=group, hour=hour, InforColumn=InforColumn, OutputFile="BridgeR_2_Relative_RPKM_distribution")
    BridgeRNormalizationForLuc2(filename="BridgeR_1_Relative_expression_dataset.txt", group=group, hour=hour, InforColumn=InforColumn, NormFactor=YourNormFactor, OutputFile="BridgeR_4_Normalized_expression_dataset.txt")
    BridgeRDatasetChecker(InputFile="BridgeR_4_Normalized_expression_dataset.txt", group=group, hour=hour, InforColumn=InforColumn, OutputFile="BridgeR_4B_Normalized_RPKM_distribution")
    BridgeRHalfLifeCalculation(filename="BridgeR_4_Normalized_expression_dataset.txt", group=group, hour=hour, InforColumn=InforColumn, CutoffRelExp=CutoffRelExp, CutoffDataPoint=CutoffDataPoint, OutputFile="BridgeR_5_HalfLife_calculation.txt")
}

BridgeRCompare <- function(InputFile="BridgeR_5_HalfLife_calculation.txt", InforColumn=4, group, hour, ComparisonFile){
    library(data.table)
    library(ggplot2)
    BridgeRHalfLifeComparison(filename=InputFile, InforColumn=4, group=group, hour=hour, ComparisonFile=ComparisonFile, LogScale=F, OutputFig="BridgeR_5_HalfLife_Comparison_ScatteredPlot")
    BridgeRHalfLifeDistribution(filename=InputFile, InforColumn=4, group=group, hour=hour, ComparisonFile=ComparisonFile, OutputFig="BridgeR_5_HalfLife_Distribution_LineGraph")
    BridgeRHalfLifeDifferenceHist(filename=InputFile, InforColumn=4, group=group, hour=hour, ComparisonFile=ComparisonFile, BinwidthFig=0.01, OutputFig="BridgeR_5_HalfLife_Difference_LineGraph")
    BridgeRHalfLifeDifferenceBox(filename=InputFile, InforColumn=4, group=group, hour=hour, ComparisonFile=ComparisonFile, OutputFig= "BridgeR_5_HalfLife_Comparison_BoxPlot")
}

BridgeRCalibration <- function(InputFile="BridgeR_5_HalfLife_calculation.txt", InforColumn=4, group, hour, ComparisonFile){
    library(data.table)
    library(ggplot2)
    BridgeRHalfLifeCalibration(InputFile=InputFile, group=group, hour=hour, ComparisonFile=ComparisonFile, InforColumn=4, OutputFile="BridgeR_6_HalfLife_Calibration")
    BridgeRHalfLifeComparison(filename=InputFile, InforColumn=4, group=group, hour=hour, ComparisonFile=ComparisonFile, LogScale=F, Calibration=T, OutputFig="BridgeR_6_Adjusted_HalfLife_Comparison_ScatteredPlot")
    BridgeRHalfLifeDifferenceHist(filename=InputFile, InforColumn=4, group=group, hour=hour, ComparisonFile=ComparisonFile, BinwidthFig=0.01, Calibration=T, OutputFig="BridgeR_6_Adjusted_HalfLife_Difference_LineGraph")
    BridgeRHalfLifeDifferenceBox(filename=InputFile, InforColumn=4, group=group, hour=hour, ComparisonFile=ComparisonFile, Calibration=T, OutputFig= "BridgeR_6_Adjusted_HalfLife_Comparison_BoxPlot")
}
