#Title: Decay rate infor
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-10-09

###Decay_rate_Infor_function###
BridgeRHalfLifeComparison <- function(filename = "BridgeR_5_half-life_calculation.txt", group, hour, ComparisonFile, InforColumn = 4, LogScale=F, OutputFig = "BridgeR_6_Half-life_comparison"){
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
    half_life_column_1 <- comp_file_number[1]*(time_points + InforColumn + 8) #number
    half_1 <- input_file[[half_life_column_1]]
    half_life_column_2 <- comp_file_number[2]*(time_points + InforColumn + 8) #number
    half_2 <- input_file[[half_life_column_2]]
        
    gene_number <- length(half_1)
    half_1_fig <- NULL
    half_2_fig <- NULL
    factor_fig <- NULL
    for(x in 1:gene_number){
        if(is.na(half_1[x]) || is.na(half_2[x])){
            next
        }
        half_1_fig <- append(half_1_fig, half_1[x])
        half_2_fig <- append(half_2_fig, half_2[x])
        div <- half_1[x]/half_2[x]
        if(div <= 0.5){
            factor_fig <- append(factor_fig, 2)
        }else if(div >= 2){
            factor_fig <- append(factor_fig, 1)
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
    print_out2 <- paste("At least 2-fold upregulated: ",down_genes," genes", sep="")
    print_out3 <- paste("At least 2-fold downregulated: ",up_genes," genes", sep="")
    print(print_out)
    print(print_out2)
    print(print_out3)
    
    p.scatter <- ggplot()
    p.scatter <- p.scatter + layer(data=plot_data, 
                                   mapping=aes(x=half_1_fig, y=half_2_fig,colour=factor(factor_fig)), 
                                   geom="point",
                                   #colour="black",
                                   size=2.5,
                                   alpha=0.3)
    
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
    p.scatter <- p.scatter + scale_colour_manual(values=c("black","blue","red")) #Change color
    plot(p.scatter)
        
    dev.off() #close_fig
    plot.new()
}

BridgeRHalfLifeDistribution <- function(filename = "BridgeR_4_half-life_calculation.txt", group, hour, ComparisonFile, InforColumn = 4, OutputFig = "BridgeR_5_Half-life_distribution"){
    ###Prepare_file_infor###
    time_points <- length(hour)
    group_number <- length(group)
    sample_number <- length(ComparisonFile)
    comp_file_number <- NULL
    figfile <- OutputFig
    for(a in 1:sample_number){
        comp_file_number <- append(comp_file_number, which(group == ComparisonFile[a]))
        figfile <- paste(figfile,"_",ComparisonFile[a], sep="")
    }
    figfile <- paste(figfile,".png", sep="")
    
    input_file <- fread(filename, header=T)
    png(filename=figfile,width = 1200, height = 1200)
    
    ###Plot_Half-life_distribution###
    half_life_fig <- NULL
    for(x in 1:sample_number){
    #for(x in 1){
        data_file <- NULL
        half_life_column <- comp_file_number[x]*(time_points + InforColumn + 8)
        half_life_data <- input_file[[half_life_column]]
        half_life_data <- half_life_data[!is.na(half_life_data)]
        if(x == 1){
            Sample <-  as.factor(rep(ComparisonFile[x], length(half_life_data)))
            half_life_fig <- data.frame(half_life_data, Sample)
        }else{
            Sample <-  as.factor(rep(ComparisonFile[x], length(half_life_data)))
            half_life_fig <- rbind(half_life_fig, data.frame(half_life_data, Sample))
        }
    }

    p.scatter <- ggplot()
    p.scatter <- p.scatter + layer(data=half_life_fig, 
                                   mapping=aes(x=half_life_data, colour=Sample), 
                                   geom="freqpoly",
                                   binwidth=0.1,
                                   #geom="line",
                                   #stat="density",
                                   size=1.2,
                                   alpha=0.5)
    p.scatter <- p.scatter + xlim(0,25)
    p.scatter <- p.scatter + ggtitle("Half-life distribution")
    p.scatter <- p.scatter + xlab("half-life")
    p.scatter <- p.scatter + ylab("Transcripts #")
    #library(scales)
    #p.scatter <- p.scatter + scale_y_continuous(labels = percent)
    plot(p.scatter)
    
    dev.off() #close_fig
    plot.new()
}

BridgeRHalfLifeDifferenceHist <- function(filename = "BridgeR_4_half-life_calculation.txt", group, hour, ComparisonFile, InforColumn = 4, BinwidthFig = 0.01, OutputFig = "BridgeR_5_Half-life_difference_Histgram"){
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
    half_life_column_1 <- comp_file_number[1]*(time_points + InforColumn + 8) #number
    half_1 <- input_file[[half_life_column_1]]
    half_life_column_2 <- comp_file_number[2]*(time_points + InforColumn + 8) #number
    half_2 <- input_file[[half_life_column_2]]
    div_half <- log2(half_2/half_1)
    
    print(summary(half_1))
    print(summary(half_2))

    plot_data <- data.frame(div_half)

    p.scatter <- ggplot()
    p.scatter <- p.scatter + layer(data=plot_data, 
                                   mapping=aes(x=div_half), 
                                   geom="freqpoly",
                                   binwidth=BinwidthFig,
                                   #geom="line",
                                   #stat="density",
                                   size=1.2)
    p.scatter <- p.scatter + xlim(min(plot_data$div_half),max(plot_data$div_half))
    p.scatter <- p.scatter + ggtitle("Half-life difference")
    name_xlab <- paste("Relative half-life(",group[comp_file_number[1]],"_vs_",group[comp_file_number[2]],")",sep="")
    p.scatter <- p.scatter + xlab(name_xlab)
    p.scatter <- p.scatter + ylab("Transcripts #")
    plot(p.scatter)
    
    dev.off() #close_fig
    plot.new()
}

BridgeRHalfLifeDifferenceBox <- function(filename = "BridgeR_4_half-life_calculation.txt", group, hour, ComparisonFile, InforColumn = 4, OutputFig = "BridgeR_5_Half-life_difference_Boxplot"){
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
    fig_width <- 200*length(ComparisonFile)
    png(filename=figfile,width = fig_width, height = 1200)
    
    ###Plot_Half-life_comparison###
    half_life_column_1 <- comp_file_number[1]*(time_points + InforColumn + 8) #number
    half_1 <- input_file[[half_life_column_1]]
    half_1 <- half_1[!is.na(half_1)]
    half_1_number <- length(half_1)
    half_1_label <- rep(ComparisonFile[1],half_1_number)
    plot_data_1 <- data.frame(half_1,half_1_label)
    colnames(plot_data_1) <- c("half_data","label")
    
    half_life_column_2 <- comp_file_number[2]*(time_points + InforColumn + 8) #number
    half_2 <- input_file[[half_life_column_2]]
    half_2 <- half_2[!is.na(half_2)]
    half_2_number <- length(half_2)
    half_2_label <- rep(ComparisonFile[2],half_2_number)
    plot_data_2 <- data.frame(half_2,half_2_label)
    colnames(plot_data_2) <- c("half_data","label")
    
    print(summary(half_1))
    print(summary(half_2))
    
    plot_data <- rbind(plot_data_1,plot_data_2)

    p.boxplot <- ggplot()
    p.boxplot <- p.boxplot + layer(data=plot_data, 
                                   mapping=aes(x=label, y=half_data), 
                                   geom="boxplot")
    p.boxplot <- p.boxplot + ylim(0,24)
    p.boxplot <- p.boxplot + ggtitle("Half-life difference")
    p.boxplot <- p.boxplot + ylab("Half-life")
    plot(p.boxplot)
    
    dev.off() #close_fig
    plot.new()
}
