# BridgeR
An R library for comprehensive BRIC-seq data analysis tool.  
BridgeR library was made under R version 3.2.2.
***
##Dependencies(Required R libraries)
data.table  
ggplot2  

##Install
###Linux
```
sudo R CMD INSTALL BridgeR
```
###Windows
(1)Update the path variables on your Windows computer(64-bit)
```
C:\Program Files\R\R-3.2.2\bin
```
(2)Install BridgeR library on Command Prompt
```
R CMD INSTALL BridgeR
```
##Example
```
library(BridgeR)

files <- c("Control.fpkm_table",
           "Knockdown.fpkm_table")
hour <- c(0,1,2,4,8,12)
group <- c("Control","Knockdown")

BridgeRCore(InputFiles=files, group=group, hour=hour, RelRPKMFig=T)
BridgeRCompare(InputFile="BridgeR_5_HalfLife_calculation.txt", group=group, hour=hour, ComparisonFile=group)

```