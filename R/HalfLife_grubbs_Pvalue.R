#Title: Halflife_SD_grubbs_test
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-12-01

BridgeRHalfSDGrubbsTest <- function(HalfLifeFile = "HalfLife_RPKM_mean_SD.txt",
                                    HalfLifeCompFile = "BridgeR_6_HalfLife_Pvalue_estimation_PUM2_study.txt",
                                    hour = c(0,1,2,4,8,12),
                                    InforColumn = 4,
                                    CompFile = 2,
                                    CompName = "siPUM2_PUM2_study",
                                    NormNames = c("siStealth_PUM1_study","siCTRL_PUM2_study","siCTRL_PUM1_2_study"),
                                    OutputFile = "HalfLife_Grubbs_test.txt"){
    ###Prepare_file_infor###
    TimeDataSize <- length(hour)
    
    Half_index <- NULL
    for(SampleIndex in CompFile){
        infor_st <- 1 + (SampleIndex - 1)*(TimeDataSize + InforColumn + 13) #13 => Infor: Model, R2, half_life, etc...
        infor_ed <- (InforColumn)*SampleIndex + (SampleIndex - 1)*(TimeDataSize + 13) #13 => Infor: Model, R2, half_life, etc...
        Half_index <- infor_ed + TimeDataSize + 8
    }
    
    half_st <- 1 + InforColumn
    half_ed <- InforColumn + length(NormNames)
    input_file_CTRL_half <- fread(HalfLifeFile, header=T)[,half_st:half_ed,with=F]
    input_file_KD_half <- fread(HalfLifeCompFile, header=T)[,Half_index,with=F]
    
    input_file_all_half <- cbind(input_file_CTRL_half, input_file_KD_half)
    
    header1 <- paste("HalfLife_",NormNames,sep="")
    header2 <- paste("HalfLife_",CompName,sep="")
    header3 <- c("P-value(Grubbs test)","Status")
    
    cat(header1, file=OutputFile, sep="\t")
    cat("\t", file=OutputFile, append=T)
    cat(header2, file=OutputFile, sep="\t", append=T)
    cat("\t", file=OutputFile, sep="\t", append=T)
    cat(header3, file=OutputFile, sep="\t", append=T)
    cat("\n", file=OutputFile, sep="\t", append=T)
    
    for(x in 1:length(input_file_all_half[[1]])){
        half_list <- as.numeric(as.vector(as.matrix(input_file_all_half[x,])))
        half_list_test <- half_list[!is.na(half_list)]
        if(length(half_list_test) == 0 || is.na(half_list[4])){
            cat(half_list,sep="\t",file=OutputFile,append=T)
            cat("\t",file=OutputFile,append=T)
            
            cat("NA",file=OutputFile,append=T)
            cat("\t",file=OutputFile,append=T)
            cat("Notest",file=OutputFile,append=T)
            cat("\n",file=OutputFile,append=T)
            next
        }else{
            grubbs_result <- grubbs.test(half_list)
            grubbs_alternative <- grubbs_result$alternative
            grubbs_alternative <- gsub("highest value ","",grubbs_alternative)
            grubbs_alternative <- gsub("lowest value ","",grubbs_alternative)
            grubbs_alternative <- as.numeric(gsub(" is an outlier","",grubbs_alternative))
            
            if(grubbs_alternative == half_list[4]){
                grubbs_pvalue <- grubbs_result$p.value
                
                cat(half_list,sep="\t",file=OutputFile,append=T)
                cat("\t",file=OutputFile,append=T)
                
                cat(grubbs_pvalue,file=OutputFile,append=T)
                cat("\t",file=OutputFile,append=T)
                cat("Grubbs",file=OutputFile,append=T)
                cat("\n",file=OutputFile,append=T)
            }else{
                cat(half_list,sep="\t",file=OutputFile,append=T)
                cat("\t",file=OutputFile,append=T)
                
                cat("NA",file=OutputFile,append=T)
                cat("\t",file=OutputFile,append=T)
                cat("Notest",file=OutputFile,append=T)
                cat("\n",file=OutputFile,append=T)
            }
        }
    }
}

###TEST###
setwd("C:/Users/Naoto/Documents/github/test_R_libraries/BridgeR_test/2015-11-30_PUM1_study_CTRL_comparison/test")
BridgeRHalfSDGrubbsTest()
