#Title: BridgeR_Core
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-11-05

BridgeRCore <- function(InputFiles,
                        InforColumn=4,
                        group,
                        hour,
                        RPKMcutoff=0.1,
                        RelRPKMFig=F,
                        SelectNormFactor=T,
                        CutoffDataPointNumber = 4,
                        CutoffDataPoint1 = c(1,2),
                        CutoffDataPoint2 = c(8,12),
                        ThresholdHalfLife = c(8,12),
                        CutoffRelExp=0.001,
                        ModelMode="R2_selection"
                        ){
    BridgeRDataSetFromCuffnorm(CuffnormFiles=InputFiles, 
                               group=group, 
                               hour=hour, 
                               cutoff=RPKMcutoff, 
                               InforColumn=InforColumn, 
                               OutputFile="BridgeR_1_Relative_expression_dataset.txt")
    BridgeRDatasetChecker(InputFile="BridgeR_1_Relative_expression_dataset.txt", 
                          group=group, 
                          hour=hour, 
                          InforColumn=InforColumn, 
                          OutputFile="BridgeR_2_Relative_RPKM_distribution")
    BridgeRNormalizationFactors(InputFile="BridgeR_1_Relative_expression_dataset.txt",
                                group=group, 
                                hour=hour, 
                                InforColumn=InforColumn, 
                                YMin=-2, 
                                YMax=2, 
                                MakeFig=RelRPKMFig, 
                                figname="BridgeR_3_Normalizaion_factor", 
                                nfname="BridgeR_3_Normalizaion_factor")
    BridgeRNormalization(filename="BridgeR_1_Relative_expression_dataset.txt", 
                         group=group, 
                         hour=hour, 
                         InforColumn=InforColumn, 
                         SelectNormFactor=SelectNormFactor, 
                         NormFactor="BridgeR_3_Normalizaion_factor", 
                         OutputFile="BridgeR_4_Normalized_expression_dataset.txt")
    BridgeRDatasetChecker(InputFile="BridgeR_4_Normalized_expression_dataset.txt", 
                          group=group, 
                          hour=hour, 
                          InforColumn=InforColumn, 
                          OutputFile="BridgeR_4B_Normalized_RPKM_distribution")
    BridgeRHalfLifeCalc(InputFiles = "BridgeR_4_Normalized_expression_dataset.txt",
                        InforColumn = InforColumn,
                        group = group, 
                        hour = hour,
                        CutoffDataPointNumber = CutoffDataPointNumber,
                        CutoffDataPoint1 = CutoffDataPoint1,
                        CutoffDataPoint2 = CutoffDataPoint2,
                        ThresholdHalfLife = ThresholdHalfLife,
                        CutoffRelExp = CutoffRelExp,
                        ModelMode = ModelMode)
}

BridgeRHKGenes <- function(InputFiles,
                           InforColumn=4, 
                           InforHKGenes=2, 
                           HKGenes=c("GAPDH","PGK1","PPIA","ENO1","ATP5B","ALDOA"), 
                           group, 
                           hour, 
                           RPKMcutoff=0.1,
                           CutoffDataPointNumber = 4,
                           CutoffDataPoint1 = c(1,2),
                           CutoffDataPoint2 = c(8,12),
                           ThresholdHalfLife = c(8,12),
                           CutoffRelExp=0.001,
                           ModelMode="R2_selection"){
    BridgeRDataSetFromCuffnorm(CuffnormFiles=InputFiles, 
                               group=group, 
                               hour=hour, 
                               cutoff=RPKMcutoff, 
                               InforColumn=InforColumn, 
                               OutputFile="BridgeR_1_Relative_expression_dataset.txt")
    BridgeRDatasetChecker(InputFile="BridgeR_1_Relative_expression_dataset.txt", 
                          group=group, 
                          hour=hour, 
                          InforColumn=InforColumn, 
                          OutputFile="BridgeR_2_Relative_RPKM_distribution")
    BridgeRNormalizationFactorsHK(InputFile="BridgeR_1_Relative_expression_dataset.txt", 
                                  group=group, 
                                  hour=hour, 
                                  InforColumn=InforColumn, 
                                  InforHKGenes=InforHKGenes, 
                                  HKGenes=HKGenes, 
                                  OutputFile="BridgeR_3_Normalizaion_factor_HouseKeepingGenes.txt")
    BridgeRNormalizationForLuc2(filename="BridgeR_1_Relative_expression_dataset.txt", 
                                group=group, 
                                hour=hour, 
                                InforColumn=InforColumn, 
                                NormFactor="BridgeR_3_Normalizaion_factor_HouseKeepingGenes.txt", OutputFile="BridgeR_4_Normalized_expression_dataset.txt")
    BridgeRDatasetChecker(InputFile="BridgeR_4_Normalized_expression_dataset.txt", 
                          group=group, 
                          hour=hour, 
                          InforColumn=InforColumn, 
                          OutputFile="BridgeR_4B_Normalized_RPKM_distribution")
    BridgeRHalfLifeCalc(InputFiles = "BridgeR_4_Normalized_expression_dataset.txt",
                        InforColumn = InforColumn,
                        group = group, 
                        hour = hour,
                        CutoffDataPointNumber = CutoffDataPointNumber,
                        CutoffDataPoint1 = CutoffDataPoint1,
                        CutoffDataPoint2 = CutoffDataPoint2,
                        ThresholdHalfLife = ThresholdHalfLife,
                        CutoffRelExp = CutoffRelExp,
                        ModelMode = ModelMode)
}

BridgeRCustom <- function(YourNormFactor, 
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
                          ModelMode="R2_selection"){
    BridgeRDataSetFromCuffnorm(CuffnormFiles=InputFiles, 
                               group=group, 
                               hour=hour, 
                               cutoff=RPKMcutoff, 
                               InforColumn=InforColumn, 
                               OutputFile="BridgeR_1_Relative_expression_dataset.txt")
    BridgeRDatasetChecker(InputFile="BridgeR_1_Relative_expression_dataset.txt", 
                          group=group, 
                          hour=hour, 
                          InforColumn=InforColumn, 
                          OutputFile="BridgeR_2_Relative_RPKM_distribution")
    BridgeRNormalizationForLuc2(filename="BridgeR_1_Relative_expression_dataset.txt", 
                                group=group, 
                                hour=hour, 
                                InforColumn=InforColumn, 
                                NormFactor=YourNormFactor, 
                                OutputFile="BridgeR_4_Normalized_expression_dataset.txt")
    BridgeRDatasetChecker(InputFile="BridgeR_4_Normalized_expression_dataset.txt", 
                          group=group, 
                          hour=hour, 
                          InforColumn=InforColumn, 
                          OutputFile="BridgeR_4B_Normalized_RPKM_distribution")
    BridgeRHalfLifeCalc(InputFiles = "BridgeR_4_Normalized_expression_dataset.txt",
                        InforColumn = InforColumn,
                        group = group, 
                        hour = hour,
                        CutoffDataPointNumber = CutoffDataPointNumber,
                        CutoffDataPoint1 = CutoffDataPoint1,
                        CutoffDataPoint2 = CutoffDataPoint2,
                        ThresholdHalfLife = ThresholdHalfLife,
                        CutoffRelExp = CutoffRelExp,
                        ModelMode = "R2_selection")
}

BridgeRHalfLifeCalc <- function(InputFiles = "BridgeR_4_Normalized_expression_dataset.txt", 
                                InforColumn = 4, 
                                group, 
                                hour,
                                CutoffDataPointNumber = 4,
                                CutoffDataPoint1 = c(1,2),
                                CutoffDataPoint2 = c(8,12),
                                ThresholdHalfLife = c(8,12),
                                CutoffRelExp = 0.001,
                                ModelMode = "R2_selection"){
    if(ModelMode == "Raw_model"){
        BridgeRHalfLifeCalcRaw(filename=InputFiles, 
                                   group=group, 
                                   hour=hour, 
                                   InforColumn=InforColumn, 
                                   CutoffRelExp=CutoffRelExp, 
                                   CutoffDataPoint=CutoffDataPointNumber, 
                                   OutputFile="BridgeR_5_HalfLife_calculation.txt")
    }else if(ModelMode == "Three_model"){
        BridgeRHalfLifeCalcModel3(filename = InputFiles, 
                                  group = group, 
                                  hour = hour, 
                                  InforColumn = InforColumn, 
                                  CutoffRelExp = CutoffRelExp, 
                                  CutoffDataPoint = CutoffDataPointNumber, 
                                  OutputFile = "BridgeR_5B_HalfLife_calculation_3model.txt")
    }else if(ModelMode == "R2_selection"){
        BridgeRHalfLifeCalcR2Select(InputFile = InputFiles,
                                    group = group, 
                                    hour = hour, 
                                    InforColumn = InforColumn, 
                                    CutoffDataPointNumber = CutoffDataPointNumber,
                                    CutoffDataPoint1 = CutoffDataPoint1,
                                    CutoffDataPoint2 = CutoffDataPoint2,
                                    ThresholdHalfLife = ThresholdHalfLife,
                                    OutputFile = "BridgeR_5C_HalfLife_calculation_R2_selection.txt")
    }else{
        print("ERROR: Defined wrong ModelMode...")
        print("Choose the following ModelMode: Raw_model, Three_model, R2_selection.")
    }
}

BridgeRCompare <- function(InputFile="BridgeR_5C_HalfLife_calculation_R2_selection.txt",
                           InforColumn = 4,
                           group, 
                           hour, 
                           ComparisonFile,
                           CutoffDataPointNumber = 4,
                           ModelMode="R2_selection",
                           Calibration=F){
    outputfile_name <- NULL
    if(ModelMode == "Three_model"){
        if(Calibration == T){
            outputfile_name <- "BridgeR_6_Calibrated_HalfLife_calculation_3model.txt"
            BridgeRHalfLifeCalibration(InputFile = "BridgeR_5B_HalfLife_calculation_3model.txt",
                                       group = group, 
                                       hour = hour, 
                                       ComparisonFile = ComparisonFile,
                                       InforColumn = InforColumn, 
                                       OutputFile = "BridgeR_6_Calibrated_HalfLife_calculation_3model.txt")
            
        }else if(Calibration == F){
            outputfile_name <- InputFile
        }
    }else if(ModelMode == "Raw_model" || ModelMode == "R2_selection"){
        if(Calibration == T){
            outputfile_name <- "BridgeR_6B_Calibrated_HalfLife_Pvalue_estimation.txt"
            BridgeRPvalueEvaluation(InputFile=InputFile,
                                    group,
                                    hour,
                                    ComparisonFile,
                                    InforColumn = InforColumn,
                                    CutoffDataPointNumber = CutoffDataPointNumber,
                                    OutputFile=outputfile_name,
                                    Calibration=Calibration)
        }else if(Calibration == F){
            outputfile_name <- "BridgeR_6_HalfLife_Pvalue_estimation.txt"
            BridgeRPvalueEvaluation(InputFile=InputFile,
                                    group,
                                    hour,
                                    ComparisonFile,
                                    InforColumn = InforColumn,
                                    CutoffDataPointNumber = CutoffDataPointNumber,
                                    OutputFile=outputfile_name,
                                    Calibration=Calibration)
        }
    }
    BridgeRCompareFig(InputFile=outputfile_name,
                      InforColumn = InforColumn,
                      group = group, 
                      hour = hour, 
                      ComparisonFile = ComparisonFile,
                      ModelMode="R2_selection")
}

BridgeRCompareFig <- function(InputFile="BridgeR_6_HalfLife_Pvalue_estimation.txt",
                              InforColumn = 4,
                              group, 
                              hour, 
                              ComparisonFile,
                              ModelMode="R2_selection"){
        BridgeRHalfLifeComparison(filename=InputFile, 
                                  InforColumn=4, 
                                  group=group, 
                                  hour=hour, 
                                  ComparisonFile=ComparisonFile, 
                                  LogScale=F, 
                                  OutputFig="BridgeR_7_HalfLife_Comparison_ScatteredPlot",
                                  ModelMode = ModelMode)
        BridgeRHalfLifeDistribution(filename=InputFile, 
                                    InforColumn=4, 
                                    group=group, 
                                    hour=hour, 
                                    ComparisonFile=ComparisonFile, 
                                    OutputFig="BridgeR_7_HalfLife_Distribution_LineGraph",
                                    ModelMode = ModelMode)
        BridgeRHalfLifeDifferenceHist(filename=InputFile,
                                      InforColumn=4, 
                                      group=group, 
                                      hour=hour, 
                                      ComparisonFile=ComparisonFile, 
                                      BinwidthFig=0.01, 
                                      OutputFig="BridgeR_7_HalfLife_Difference_LineGraph",
                                      ModelMode = ModelMode)
        BridgeRHalfLifeDifferenceBox(filename=InputFile, 
                                     InforColumn=4, 
                                     group=group, 
                                     hour=hour, 
                                     ComparisonFile=ComparisonFile, 
                                     OutputFig= "BridgeR_7_HalfLife_Comparison_BoxPlot",
                                     ModelMode = ModelMode)
}
