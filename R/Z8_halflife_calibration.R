#Title: Half-life calibration
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-11-06

BridgeRHalfLifeCalibration <- function(InputFile = "BridgeR_5B_HalfLife_calculation_3model.txt", 
                                       group, 
                                       hour, 
                                       ComparisonFile, 
                                       InforColumn = 4, 
                                       OutputFile = "BridgeR_6_Calibrated_HalfLife_calculation_3model.txt"){
    #ComparisonFile: The length of vector must be 2 => c("Control","Knockdown")
    ###Prepare_file_infor###
    time_points <- length(hour)
    group_number <- length(group)
    comp_file_number <- NULL
    for(a in 1:length(ComparisonFile)){
        comp_file_number <- append(comp_file_number, which(group == ComparisonFile[a]))
    }
    
    input_file <- fread(InputFile, header=T)

    ###Plot_Half-life_comparison###
    gene_infor <- input_file[,1:InforColumn,with=F]
    
    half_life_column_1 <- comp_file_number[1]*(time_points + InforColumn + 23) #number
    half_1 <- as.numeric(input_file[[half_life_column_1]])
    
    half_life_column_2 <- comp_file_number[2]*(time_points + InforColumn + 23) #number
    half_2 <- as.numeric(input_file[[half_life_column_2]])
    
    half_data <- data.table(half1=half_1,half2=half_2)
    
    test <- lm(half2 ~ half1 + 0, data=half_data)
    coef <- as.numeric(test$coefficients)
    adjusted_half_2 <- half_2/coef
    adjusted_half_2_r <- NULL
    for(x in adjusted_half_2){
        if(is.na(x)){
        }else if(x >= 24){
            x <- 24
        }
        adjusted_half_2_r <- append(adjusted_half_2_r,x)
    }
    adjusted_half_2 <- adjusted_half_2_r # Vector
    
    adjusted_half_data <- data.table(Selected_HalfLife=adjusted_half_2)
    input_file <- input_file[,-half_life_column_2,with=F]
    output_data <- cbind(input_file,adjusted_half_data)
        
    write.table(output_data, file=OutputFile,sep="\t",row.names=FALSE)
}
