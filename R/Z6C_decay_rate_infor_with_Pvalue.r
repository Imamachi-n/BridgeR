#Title: Decay rate infor
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-10-09

###Decay_rate_Infor_function###
BridgeRHalfLifeComparisonWithPvalue <- function(filename = "BridgeR_7_HalfLife_Calibration.txt",
                                                group,
                                                hour,
                                                ComparisonFile,
                                                InforColumn = 4,
                                                LogScale=F,
                                                Calibration=F,
                                                OutputFig = "BridgeR_7_Halflife_comparison"){
    
    #ComparisonFile: The length of vector must be 2 => c("Control","Knockdown")
    ###Prepare_file_infor###
    time_points <- length(hour)
    group_number <- length(group)
    comp_file_number <- NULL
    for(a in 1:length(ComparisonFile)){
        comp_file_number <- append(comp_file_number, which(group == ComparisonFile[a]))
    }
    
    input_file <- fread(filename, header=T)
    figfile <- paste(OutputFig,"_",group[comp_file_number[1]],"_vs_",group[comp_file_number[2]],".png", sep="")
    png(filename=figfile,width = 1200, height = 1200)
    
    ###Plot_Half-life_comparison###
    half_life_column_1 <- time_points + InforColumn + 8 #number
    half_1 <- as.numeric(input_file[[half_life_column_1]])
    half_life_column_2 <- 2*(time_points + InforColumn) + 13 + 8 #number
    half_2 <- as.numeric(input_file[[half_life_column_2]])

    half_data <- data.table(half1=half_1,half2=half_2)
    
    if(Calibration == T){
        test <- lm(half2 ~ half1 + 0, data=half_data)
        coef <- as.numeric(test$coefficients)
        half_2 <- half_2/coef
        
        half_2_r <- NULL
        for(x in half_2){
            if(is.na(x)){
            }else if(x >= 24){
                x <- 24
            }
            half_2_r <- append(half_2_r,x)
        }
        half_2 <- half_2_r
    }

    gene_number <- length(half_1)
    half_1_fig <- NULL
    half_2_fig <- NULL
    factor_fig <- NULL
    for(x in 1:gene_number){
        Pvalue <- as.numeric(as.vector(as.matrix(input_file[x, 48, with=F])))
        div <- as.numeric(as.vector(as.matrix(input_file[x, 47, with=F])))
        
        if(is.na(half_1[x]) || is.na(half_2[x]) || Pvalue == "NaN"){
            next
        }
        half_1_fig <- append(half_1_fig, half_1[x])
        half_2_fig <- append(half_2_fig, half_2[x])
        
        if(div < 0 && Pvalue < 0.05){
            factor_fig <- append(factor_fig, 2) #Down-regulated
        }else if(div > 0 && Pvalue < 0.05){
            factor_fig <- append(factor_fig, 1) #up-regulated
        }else{
            factor_fig <- append(factor_fig, 0)
        }
    }
    if(LogScale){
        half_1_fig <- log2(half_1_fig)
        half_2_fig <- log2(half_2_fig)
    }
    up_genes <- length(which(factor_fig == 1))
    down_genes <- length(which(factor_fig == 2))
    
    plot_data <- data.frame(half_1_fig,half_2_fig,factor_fig)
    print_out <- paste("Plotted: ",length(plot_data[,1])," genes", sep="")
    print_out2 <- paste("At least 2-fold upregulated: ",up_genes," genes", sep="")
    print_out3 <- paste("At least 2-fold downregulated: ",down_genes," genes", sep="")
    print(print_out)
    print(print_out2)
    print(print_out3)
    
    p.scatter <- ggplot()
    p.scatter <- p.scatter + layer(data=plot_data, 
                                   mapping=aes(x=half_1_fig, y=half_2_fig,colour=factor(factor_fig)), 
                                   geom="point",
                                   #colour="black",
                                   size=2.5,
                                   alpha=0.5)
    
    p.scatter <- p.scatter + layer(data=plot_data, 
                                   mapping=aes(x=half_1_fig, y=half_2_fig),
                                   geom="smooth",
                                   geom_params=list(color = "blue", size=1.2),
                                   stat="smooth",
                                   stat_params=list(method="lm", se=F))
    
    p.scatter <- p.scatter + xlim(0,max(plot_data$half_1_fig)) + ylim(0,max(plot_data$half_2_fig))
    p.scatter <- p.scatter + ggtitle("Half-life comparison")
    name_xlab <- paste(group[comp_file_number[1]]," (Time)", sep="")
    name_ylab <- paste(group[comp_file_number[2]]," (Time)", sep="")
    p.scatter <- p.scatter + xlab(name_xlab)
    p.scatter <- p.scatter + ylab(name_ylab)
    p.scatter <- p.scatter + theme(legend.position="none") #Remove Guides
    p.scatter <- p.scatter + scale_colour_manual(values=c("black","red","blue")) #Change color
    plot(p.scatter)
        
    dev.off() #close_fig
    plot.new()
}
