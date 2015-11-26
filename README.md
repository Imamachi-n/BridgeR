# BridgeR
An R library for comprehensive BRIC-seq data analysis tool.  
BridgeR library was made under R version 3.2.2.
***
##Dependencies (Required R libraries)
data.table  
ggplot2  
shiny  
BSDA  

##Install
###Linux
```
sudo R CMD INSTALL BridgeR
```
###Windows
(1) Update the path variables on your Windows computer(64-bit)
```
C:\Program Files\R\R-3.2.2\bin
```
(2) Install BridgeR library on Command Prompt
```
R CMD INSTALL BridgeR
```
##Example
(1) Calculate RNA half-life from your BRIC-seq dataset and compare RNA half-life bwtween two conditions.
```
library(BridgeR)

files <- c("Control.fpkm_table",
           "Knockdown.fpkm_table")
hour <- c(0,1,2,4,8,12)
group <- c("Control","Knockdown")

BridgeRCore(InputFiles=files,
            InforColumn=4,
            group=group,
            hour=hour,
            SelectNormFactor=T,
            CutoffDataPointNumber = 4,
            CutoffDataPoint1 = c(1,2),
            CutoffDataPoint2 = c(8,12),
            ThresholdHalfLife = c(8,12),
            ModelMode="R2_selection")

BridgeRCompare(InputFile="BridgeR_5_HalfLife_calculation.txt",
               InforColumn=4,
               group=group,
               hour=hour,
               ComparisonFile=group,
               CutoffDataPointNumber = 4,
               ModelMode="R2_selection",
               Calibration=F)
```
(2) Draw RNA decay curve predicted from your BRIC-seq dataset.
```
library(BridgeR)

RPKM_ctrl <- "Control.fpkm_table"
RPKM_kd <- "Knockdown.fpkm_table"
Normalized_data <- "BridgeR_4_Normalized_expression_dataset.txt"
group <- c("Control","Knockdown")
hour <- c(0,1,2,4,8,12)

BridgeReport6P(filename1=RPKM_ctrl,
               filename2=RPKM_kd,
               filename3=Normalized_data,
               group=group,
               hour=hour, 
               ComparisonFile=group,
               SearchRow = "symbol",
               InforColumn = 4,
               Color = c("black","red"),
               CutoffDataPoint1 = c(1,2),
               CutoffDataPoint2 = c(8,12))
```
