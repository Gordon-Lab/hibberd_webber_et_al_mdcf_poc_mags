library(ggplot2)
library(ggsignif)
library(tidyverse)
library('RColorBrewer')
library(ggpubr)
## Set WD
setwd("C:/Users/dmweb/Box/MDCF_POC_Microbiome_Paper/code for deposition/hibberd_webber_et_al_mdcf_poc_mags/glycan_plots")

# Load glycan data
(load("carbohydrate_quants.RData"))

####### Info about loaded objects:
# object, description, units
# -------------------------------------------------------------------------
# monosac, monosaccharide quantities, (ug monosaccharide)/(mg dry wt of ingredient or diet) 
# linkage, linkage quantities, linkage peak area normalized to input mass (reported as AUC  
# FITDOG, polysaccharide quantities, (ug polysaccharide)/(mg dry wt of ingredient or diet)

# Glycan linkages
###############################
  # 1. Calculate mean, SD, and perform T-test
  # Pivot data to tall 
  glycan_tall <- pivot_longer(linkage, col=2:ncol(linkage), names_to = "glycan", values_to = "AUC_glycan")
  # Calculate mean and SD
  glycan_mean_SD <- glycan_tall %>% dplyr::group_by(Sample_ID, glycan) %>% dplyr::summarise(mean = mean(AUC_glycan), SD=sd(AUC_glycan))
  # Narrow the glycan df down to MDCF2 and RUSF
  glycan_mean_SD_MDCF2_RUSF <- glycan_mean_SD[glycan_mean_SD$Sample_ID %in% c("MDCF_2","RUSF"),]
  glycan_MDCF2_vs_RUSF <-  glycan_tall[glycan_tall$Sample_ID %in% c("MDCF_2","RUSF"),]
  # Perform T-test
  glycan_MDCF2_vs_RUSF_stats <- compare_means(AUC_glycan ~ Sample_ID, data = glycan_MDCF2_vs_RUSF, group.by = "glycan", method = "t.test", paired = FALSE, p.adjust.method="fdr")
  #Pivot tall df with mean and SD to wide format so that it can be merged with stats tests
  glycan_mean_SD_wide <- pivot_wider(glycan_mean_SD, names_from = Sample_ID, values_from = c("mean","SD")) #Pivot from tall to wide format
  #Determine the fold change between MDCF2 and RUSF mono quants
  glycan_mean_SD_wide$log2MDCF2_RUSF <- log2(glycan_mean_SD_wide$mean_MDCF_2/glycan_mean_SD_wide$mean_RUSF)
  #Joint to mono stats df
  glycan_summary_wide <- left_join(glycan_mean_SD_wide, glycan_MDCF2_vs_RUSF_stats)
    #glycan_summary_wide <- left_join(glycan_summary_wide, LOD)
  #Get rid of unwanted columns
  glycan_summary_wide[,c(".y.","group1","group2")] <- NULL
  
  #Write glycan quant info to csv
  write.csv(glycan_summary_wide, "./out/glycan_summary_wide_210728.csv")

  # Convert AUC to millions
  glycan_mean_SD$mean_mil <-  glycan_mean_SD$mean /1000000
  glycan_mean_SD$SD_mil <- glycan_mean_SD$SD / 1000000

  # Narrow the glycan df down to MDCF2 vs RUSF, RUSF ingredients, and MDCF2 ingredients and set the order
  glycan_mean_SD_MDCF2_RUSF <- glycan_mean_SD[glycan_mean_SD$Sample_ID %in% c("RUSF","MDCF_2"),]
   glycan_mean_SD_MDCF2_RUSF$Sample_ID <- factor(glycan_mean_SD_MDCF2_RUSF$Sample_ID ,levels = c("RUSF","MDCF_2"))
  glycan_mean_SD_RUSF_ing <- glycan_mean_SD[glycan_mean_SD$Sample_ID %in% c("Lentil","Milk_powder","Rice"),]
   glycan_mean_SD_RUSF_ing$Sample_ID <- factor(glycan_mean_SD_RUSF_ing$Sample_ID ,levels = c("Lentil","Milk_powder","Rice"))
  glycan_mean_SD_MDCF2_ing <- glycan_mean_SD[glycan_mean_SD$Sample_ID %in% c("Chickpea","Green_banana","Peanut","Soybean"),]
   glycan_mean_SD_MDCF2_ing$Sample_ID <- factor(glycan_mean_SD_MDCF2_ing$Sample_ID, levels = c("Chickpea","Green_banana","Peanut","Soybean") )
  #######################################
  # Make a horizontal bar graph with separate x-axis scales for high and lower linkage quants 
  #######################################
    # A. Select and graph the linkages with larger values (upper portion of the plot)
    #####################################################
    #Set glycans to show in the figure
    glycan_above_cuttoff <- glycan_summary_wide[order(glycan_summary_wide$mean_MDCF_2,decreasing=TRUE)[1:20],c("glycan")]
    glycan_above_cuttoff$glycan <- rev(glycan_above_cuttoff$glycan) #Reverse the order
   
    linkage_filter_short <- c("T-Glucose","4-Glucose","T-Galactose")
    linkage_filter_long <- glycan_above_cuttoff[glycan_above_cuttoff$glycan %in% linkage_filter_short,][[1]]
    
    glycan_mean_SD_MDCF2_RUSF_top20 <- glycan_mean_SD_MDCF2_RUSF[glycan_mean_SD_MDCF2_RUSF$glycan %in% linkage_filter_long,] #Display the same set of 20 glycans in all 3 plots
    glycan_mean_SD_RUSF_ing_top20 <- glycan_mean_SD_RUSF_ing[glycan_mean_SD_RUSF_ing$glycan %in% linkage_filter_long,] #Display the same set of 20 glycans in all 3 plots
    glycan_mean_SD_MDCF2_ing_top20 <- glycan_mean_SD_MDCF2_ing[glycan_mean_SD_MDCF2_ing$glycan %in% linkage_filter_long,] #Display the same set of 20 glycans in all 3 plots
    
    #Set the X-axis order
    glycan_mean_SD_MDCF2_RUSF_top20$glycan <- factor(glycan_mean_SD_MDCF2_RUSF_top20$glycan, levels= glycan_above_cuttoff$glycan)
    glycan_mean_SD_RUSF_ing_top20$glycan <- factor(glycan_mean_SD_RUSF_ing_top20$glycan, levels= glycan_above_cuttoff$glycan)
    glycan_mean_SD_MDCF2_ing_top20$glycan <- factor(glycan_mean_SD_MDCF2_ing_top20$glycan, levels = glycan_above_cuttoff$glycan)
    
    # 3. Make glycan bar plots with normal axis
    # MDCF2 vs RUSF
    pdf("./out/glycan_quant_MDCF2_RUSF_upper.pdf", width = 2.5, height = 1.8)
    p <- ggplot(glycan_mean_SD_MDCF2_RUSF_top20, aes(x= glycan, y=mean_mil, fill=Sample_ID)) + 
      geom_bar(stat="identity", color="black", 
               position=position_dodge()) +
      geom_errorbar(aes(ymin=mean_mil-SD_mil, ymax=mean_mil+SD_mil), width=.4,
                    position=position_dodge(.9)) +
      theme_classic() +
      theme(axis.text.x = element_text(size = 6, vjust = 1, hjust=1)) + 
      theme(axis.text.y = element_text(size = 6)) + 
      theme(legend.position = "none") +
      ylim(NA,15) +
      scale_fill_manual(values=c("#FFFFFF","#CCCCCC"))
    p + coord_flip()
    dev.off()
    
    # RUSF
    pdf("./out/glycan_quant_RUSF_ing_upper.pdf", width = 2.5, height = 1.8)
    p <- ggplot(glycan_mean_SD_RUSF_ing_top20, aes(x= glycan, y=mean_mil, fill=Sample_ID)) + 
      geom_bar(stat="identity", color="black", 
               position=position_dodge()) +
      geom_errorbar(aes(ymin=mean_mil-SD_mil, ymax=mean_mil+SD_mil), width=.4,
                    position=position_dodge(.9)) +
      theme_classic() +
      theme(axis.text.x = element_text(size = 6, vjust = 1, hjust=1)) + 
      theme(axis.text.y = element_text(size = 6)) + 
      theme(legend.position = "none") +
      ylim(NA,15) +
      scale_fill_manual(values=c("#E41A1C","#377EB8","#666666"))
    p + coord_flip()
    dev.off()
    
    # MDCF2 
    pdf("./out/glycan_quant_MDCF2_ing_upper.pdf", width = 2.5, height = 1.8)
    p <- ggplot(glycan_mean_SD_MDCF2_ing_top20, aes(x= glycan, y=mean_mil, fill=Sample_ID)) + 
      geom_bar(stat="identity", color="black", 
               position=position_dodge()) +
      geom_errorbar(aes(ymin=mean_mil-SD_mil, ymax=mean_mil+SD_mil), width=.4,
                    position=position_dodge(.9)) +
      theme_classic() +
      theme(axis.text.x = element_text(size = 6, vjust = 1, hjust=1)) + 
      theme(axis.text.y = element_text(size = 6)) + 
      theme(legend.position = "none") +
      ylim(NA,15) +
      scale_fill_manual(values=c("#FFFF33","#4DAF4A","#A65628","#FF7F00"))
    p + coord_flip()
    dev.off()
  
    # B. Select and graph the linkages with smaller values (lower portion of the plot) 
    #####################################################
    #Set glycans to show in the figure
    glycan_above_cuttoff <- glycan_summary_wide[order(glycan_summary_wide$mean_MDCF_2,decreasing=TRUE)[1:20],c("glycan")]
    glycan_above_cuttoff$glycan <- rev(glycan_above_cuttoff$glycan) #Reverse the order
    
    linkage_filter_short <- c("T-Glucose","4-Glucose","T-Galactose")
    linkage_filter_long <- glycan_above_cuttoff[!glycan_above_cuttoff$glycan %in% linkage_filter_short,][[1]]
    
    glycan_mean_SD_MDCF2_RUSF_top20 <- glycan_mean_SD_MDCF2_RUSF[glycan_mean_SD_MDCF2_RUSF$glycan %in% linkage_filter_long,] #Display the same set of 20 glycans in all 3 plots
    glycan_mean_SD_RUSF_ing_top20 <- glycan_mean_SD_RUSF_ing[glycan_mean_SD_RUSF_ing$glycan %in% linkage_filter_long,] #Display the same set of 20 glycans in all 3 plots
    glycan_mean_SD_MDCF2_ing_top20 <- glycan_mean_SD_MDCF2_ing[glycan_mean_SD_MDCF2_ing$glycan %in% linkage_filter_long,] #Display the same set of 20 glycans in all 3 plots
    
    #Set the X-axis order
    glycan_mean_SD_MDCF2_RUSF_top20$glycan <- factor(glycan_mean_SD_MDCF2_RUSF_top20$glycan, levels= glycan_above_cuttoff$glycan)
    glycan_mean_SD_RUSF_ing_top20$glycan <- factor(glycan_mean_SD_RUSF_ing_top20$glycan, levels= glycan_above_cuttoff$glycan)
    glycan_mean_SD_MDCF2_ing_top20$glycan <- factor(glycan_mean_SD_MDCF2_ing_top20$glycan, levels = glycan_above_cuttoff$glycan)
    
    # 3. Make glycan bar plots with normal axis
    # MDCF2 vs RUSF
    pdf("./out/glycan_quant_MDCF2_RUSF_lower.pdf", width = 2.5, height = 5)
    p <- ggplot(glycan_mean_SD_MDCF2_RUSF_top20, aes(x= glycan, y=mean_mil, fill=Sample_ID)) + 
      geom_bar(stat="identity", color="black", 
               position=position_dodge()) +
      geom_errorbar(aes(ymin=mean_mil-SD_mil, ymax=mean_mil+SD_mil), width=.4,
                    position=position_dodge(.9)) +
      theme_classic() +
      theme(axis.text.x = element_text(size = 6, vjust = 1, hjust=1)) + 
      theme(axis.text.y = element_text(size = 6)) + 
      theme(legend.position = "none") +
      ylim(NA,2.5) +
      scale_fill_manual(values=c("#FFFFFF","#CCCCCC"))
    p + coord_flip()
    dev.off()
    
    # RUSF
    pdf("./out/glycan_quant_RUSF_ing_lower.pdf", width = 2.5, height = 5)
    p <- ggplot(glycan_mean_SD_RUSF_ing_top20, aes(x= glycan, y=mean_mil, fill=Sample_ID)) + 
      geom_bar(stat="identity", color="black", 
               position=position_dodge()) +
      geom_errorbar(aes(ymin=mean_mil-SD_mil, ymax=mean_mil+SD_mil), width=.4,
                    position=position_dodge(.9)) +
      theme_classic() +
      theme(axis.text.x = element_text(size = 6, vjust = 1, hjust=1)) + 
      theme(axis.text.y = element_text(size = 6)) + 
      theme(legend.position = "none") +
      ylim(NA,2.5) +
      scale_fill_manual(values=c("#E41A1C","#377EB8","#666666"))
    p + coord_flip()
    dev.off()
    
    # MDCF2 
    pdf("./out/glycan_quant_MDCF2_ing_lower.pdf", width = 2.5, height = 5)
    p <- ggplot(glycan_mean_SD_MDCF2_ing_top20, aes(x= glycan, y=mean_mil, fill=Sample_ID)) + 
      geom_bar(stat="identity", color="black", 
               position=position_dodge()) +
      geom_errorbar(aes(ymin=mean_mil-SD_mil, ymax=mean_mil+SD_mil), width=.4,
                    position=position_dodge(.9)) +
      theme_classic() +
      theme(axis.text.x = element_text(size = 6, vjust = 1, hjust=1)) + 
      theme(axis.text.y = element_text(size = 6)) + 
      theme(legend.position = "none") +
      ylim(NA,2.5) +
      scale_fill_manual(values=c("#FFFF33","#4DAF4A","#A65628","#FF7F00"))
    p + coord_flip()
    dev.off()
    
#FITDOG data:
#############################

  ###############################
  # 1. Calculate mean, SD, and perform T-test
  # Pivot data to tall 
  carb_tall <- pivot_longer(FITDOG, col=2:ncol(FITDOG), names_to = "carb", values_to = "ug_carb")
  #Narrow the df to the diets of interest
  carb_tall <- carb_tall[carb_tall$sample %in% c("chickpea","green_banana","lentil","MDCF2_human","milk_powder","peanut","rice","RUSF_human","soybean"),]   
  # Calculate mean and SD
  carb_mean_SD <- carb_tall %>% dplyr::group_by(sample, carb) %>% dplyr::summarise(mean = mean(ug_carb), SD=sd(ug_carb))
  carb_mean_SD_hold <- carb_mean_SD #Create a holder variable
  # Narrow the carb df down to MDCF2 and RUSF
  carb_mean_SD_MDCF2_RUSF <- carb_mean_SD[carb_mean_SD$sample %in% c("MDCF2_human","RUSF_human"),]
  carb_MDCF2_vs_RUSF <-  carb_tall[carb_tall$sample %in% c("MDCF2_human","RUSF_human"),]
  # Perform T-test
  carb_MDCF2_vs_RUSF_stats <- compare_means(ug_carb ~ sample, data = carb_MDCF2_vs_RUSF, group.by = "carb", method = "t.test", paired = FALSE, p.adjust.method= "fdr")
  #Pivot tall df with mean and SD to wide format so that it can be merged with stats tests
  carb_mean_SD_wide <- pivot_wider(carb_mean_SD, names_from = sample, values_from = c("mean","SD")) #Pivot from tall to wide format
  #Determine the fold change between MDCF2 and RUSF mono quants
  carb_mean_SD_wide$log2MDCF2_RUSF <- log2(carb_mean_SD_wide$mean_MDCF2_human/carb_mean_SD_wide$mean_RUSF_human)
  #Joint to mono stats df
  carb_summary_wide <- left_join(carb_mean_SD_wide, carb_MDCF2_vs_RUSF_stats)
  #carb_summary_wide <- left_join(carb_summary_wide, LOD)
  #Get rid of unwanted columns
  carb_summary_wide[,c(".y.","group1","group2")] <- NULL
  
  #Write glycan quant info to csv
  write.csv(carb_summary_wide, "./out/polysaccharide_summary_wide.csv")
  
  #######################################
  # Make a horizontal bar graph with separate x-axis scales for high and lower polysaccharide quants 
  #######################################
  # A. Select and graph non-starch FITDOG results (lower portion of the plot)
  #####################################################
    #Remove starch from the figure:
    carb_mean_SD <- carb_mean_SD_hold[carb_mean_SD_hold$carb %in% c("Arabinan","Cellulose","Galactan","Mannan","Xylan"),]

    # Narrow the carb df down to MDCF2 vs RUSF, RUSF ingredients, and MDCF2 ingredients and set the order
    carb_mean_SD_MDCF2_RUSF <- carb_mean_SD[carb_mean_SD$sample %in% c("RUSF_human","MDCF2_human"),]
    carb_mean_SD_MDCF2_RUSF$sample <- factor(carb_mean_SD_MDCF2_RUSF$sample ,levels = c("RUSF_human","MDCF2_human"))
    carb_mean_SD_RUSF_ing <- carb_mean_SD[carb_mean_SD$sample %in% c("lentil","milk_powder","rice"),]
    carb_mean_SD_RUSF_ing$sample <- factor(carb_mean_SD_RUSF_ing$sample ,levels = c("lentil","milk_powder","rice"))
    carb_mean_SD_MDCF2_ing <- carb_mean_SD[carb_mean_SD$sample %in% c("chickpea","green_banana","peanut","soybean"),]
    carb_mean_SD_MDCF2_ing$sample <- factor(carb_mean_SD_MDCF2_ing$sample, levels = c("chickpea","green_banana","peanut","soybean") )
  
    #Order the carbs:
    carb_mean_SD_MDCF2 <- carb_mean_SD[carb_mean_SD$sample %in% c("MDCF2_human"),]
    carb_order <- as.list(carb_mean_SD_MDCF2_RUSF[order(carb_mean_SD_MDCF2$mean, decreasing=FALSE),c("carb")])
    
    #Set the X-axis order
    carb_mean_SD_MDCF2_RUSF$carb <- factor(carb_mean_SD_MDCF2_RUSF$carb, levels= carb_order[[1]])
    carb_mean_SD_RUSF_ing$carb <- factor(carb_mean_SD_RUSF_ing$carb, levels= carb_order[[1]])
    carb_mean_SD_MDCF2_ing$carb <- factor(carb_mean_SD_MDCF2_ing$carb, levels = carb_order[[1]])
    
    # 3. Make polysaccharide bar plots for FITDOG data
    ############################
    
    #MDCF2 vs RUSF
    pdf("./out/polysaccharide_quant_MDCF2_RUSF_nonstarch.pdf", width = 2.1, height = 2)
    p <- ggplot(carb_mean_SD_MDCF2_RUSF, aes(x= carb, y=mean, fill=sample)) + 
      geom_bar(stat="identity", color="black", 
               position=position_dodge()) +
      geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD), width=.4,
                    position=position_dodge(.9)) +
      theme_classic() +
      theme(axis.text.x = element_text(size = 6, vjust = 1, hjust=1)) + 
      theme(axis.text.y = element_text(size = 6)) + 
      theme(legend.position = "none") +
      scale_fill_manual(values=c("#FFFFFF","#CCCCCC")) +
      ylim(0,15)
      p + coord_flip()
    dev.off()
    
    # RUSF ingredients
    pdf("./out/polysaccharide_quant_RUSF_ing_nonstarch.pdf", width = 2.1, height = 2)
    p <- ggplot(carb_mean_SD_RUSF_ing, aes(x= carb, y=mean, fill=sample)) + 
      geom_bar(stat="identity", color="black", 
               position=position_dodge()) +
      geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD), width=.4,
                    position=position_dodge(.9)) +
      theme_classic() +
      theme(axis.text.x = element_text(size = 6, vjust = 1, hjust=1)) + 
      theme(axis.text.y = element_text(size = 6)) + 
      theme(legend.position = "none") +
      scale_fill_manual(values=c("#E41A1C","#377EB8","#666666")) +
      ylim(0,15)
      p + coord_flip()
  
    dev.off()
    
    # MDCF2 ingredients
    pdf("./out/polysaccharide_quant_MDCF2_ing_nonstarch.pdf", width = 2.1, height = 2)
    p <- ggplot(carb_mean_SD_MDCF2_ing, aes(x= carb, y=mean, fill=sample)) + 
      geom_bar(stat="identity", color="black", 
               position=position_dodge()) +
      geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD), width=.4,
                    position=position_dodge(.9)) +
      theme_classic() +
      theme(axis.text.x = element_text(size = 6, vjust = 1, hjust=1)) + 
      theme(axis.text.y = element_text(size = 6)) + 
      theme(legend.position = "none") +
      scale_fill_manual(values=c("#FFFF33","#4DAF4A","#A65628","#FF7F00"))  +
      ylim(0,15)
      p + coord_flip()
    dev.off() 
    
    # B. Select and graph starch FITDOG results (upper portion of the plot)
    #####################################################
    #Remove starch from the figure:

    carb_mean_SD <- carb_mean_SD_hold[carb_mean_SD_hold$carb %in% c("Starch"),]
    
    # Narrow the carb df down to MDCF2 vs RUSF, RUSF ingredients, and MDCF2 ingredients and set the order
    carb_mean_SD_MDCF2_RUSF <- carb_mean_SD[carb_mean_SD$sample %in% c("RUSF_human","MDCF2_human"),]
    carb_mean_SD_MDCF2_RUSF$sample <- factor(carb_mean_SD_MDCF2_RUSF$sample ,levels = c("RUSF_human","MDCF2_human"))
    carb_mean_SD_RUSF_ing <- carb_mean_SD[carb_mean_SD$sample %in% c("lentil","milk_powder","rice"),]
    carb_mean_SD_RUSF_ing$sample <- factor(carb_mean_SD_RUSF_ing$sample ,levels = c("lentil","milk_powder","rice"))
    carb_mean_SD_MDCF2_ing <- carb_mean_SD[carb_mean_SD$sample %in% c("chickpea","green_banana","peanut","soybean"),]
    carb_mean_SD_MDCF2_ing$sample <- factor(carb_mean_SD_MDCF2_ing$sample, levels = c("chickpea","green_banana","peanut","soybean") )
    
    #Order the carbs:
    carb_mean_SD_MDCF2 <- carb_mean_SD[carb_mean_SD$sample %in% c("MDCF2_human"),]
    carb_order <- as.list(carb_mean_SD_MDCF2_RUSF[order(carb_mean_SD_MDCF2$mean, decreasing=FALSE),c("carb")])
    
    #Set the X-axis order
    carb_mean_SD_MDCF2_RUSF$carb <- factor(carb_mean_SD_MDCF2_RUSF$carb, levels= carb_order[[1]])
    carb_mean_SD_RUSF_ing$carb <- factor(carb_mean_SD_RUSF_ing$carb, levels= carb_order[[1]])
    carb_mean_SD_MDCF2_ing$carb <- factor(carb_mean_SD_MDCF2_ing$carb, levels = carb_order[[1]])
    
    # 3. Make polysaccharide bar plots for FITDOG data
    ############################
    
    #MDCF2 vs RUSF
    pdf("./out/polysaccharide_quant_MDCF2_RUSF_starch.pdf", width = 2.1, height = 1.1)
    p <- ggplot(carb_mean_SD_MDCF2_RUSF, aes(x= carb, y=mean, fill=sample)) + 
      geom_bar(stat="identity", color="black", 
               position=position_dodge()) +
      geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD), width=.4,
                    position=position_dodge(.9)) +
      theme_classic() +
      theme(axis.text.x = element_text(size = 6, vjust = 1, hjust=1)) + 
      theme(axis.text.y = element_text(size = 6)) + 
      theme(legend.position = "none") +
      scale_fill_manual(values=c("#FFFFFF","#CCCCCC")) +
      ylim(NA,500)
      p + coord_flip()
    dev.off()
    
    # RUSF ingredients
    pdf("./out/polysaccharide_quant_RUSF_ing_starch.pdf", width = 2.1, height = 1.1)
    p <- ggplot(carb_mean_SD_RUSF_ing, aes(x= carb, y=mean, fill=sample)) + 
      geom_bar(stat="identity", color="black", 
               position=position_dodge()) +
      geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD), width=.4,
                    position=position_dodge(.9)) +
      theme_classic() +
      theme(axis.text.x = element_text(size = 6, vjust = 1, hjust=1)) + 
      theme(axis.text.y = element_text(size = 6)) + 
      theme(legend.position = "none") +
      scale_fill_manual(values=c("#E41A1C","#377EB8","#666666")) +
      ylim(NA,500)
      p + coord_flip()
    dev.off()
    
    # MDCF2 ingredients
    pdf("./out/polysaccharide_quant_MDCF2_ing_starch.pdf", width = 2.1, height = 1.1)
    p <- ggplot(carb_mean_SD_MDCF2_ing, aes(x= carb, y=mean, fill=sample)) + 
      geom_bar(stat="identity", color="black", 
               position=position_dodge()) +
      geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD), width=.4,
                    position=position_dodge(.9)) +
      theme_classic() +
      theme(axis.text.x = element_text(size = 6, vjust = 1, hjust=1)) + 
      theme(axis.text.y = element_text(size = 6)) + 
      theme(legend.position = "none") +
      scale_fill_manual(values=c("#FFFF33","#4DAF4A","#A65628","#FF7F00"))  +
      ylim(NA,500)
      p + coord_flip()
    dev.off() 
    
############################## Monosaccharides ###############
  ############################################################
  #Extract monosaccharides LOD info from df
  LOD <- as.data.frame(t(monosac[nrow(monosac),2:ncol(monosac)]))
  LOD$monosaccharide <- rownames(LOD)
  colnames(LOD)[1] <- "LOD"
  rownames(LOD) <- NULL
  monosac <- monosac[1:nrow(monosac)-1,]
  #Monosaccharides
  ##############################
  # 1. Calculate mean, SD, and perform T-test
  # Pivot data to tall 
  monosac_tall <- pivot_longer(monosac, col=2:ncol(monosac), names_to = "monosaccharide", values_to = "ug_monosac")
  # Calculate mean and SD
  monosac_mean_SD <- monosac_tall %>% dplyr::group_by(Sample_ID, monosaccharide) %>% dplyr::summarise(mean = mean(ug_monosac), SD=sd(ug_monosac))
  
  # Narrow the monosac df down to MDCF2 and RUSF
  monosac_mean_SD_MDCF2_RUSF <- monosac_mean_SD[monosac_mean_SD$Sample_ID %in% c("MDCF_2","RUSF"),]
  monosac_MDCF2_vs_RUSF <-  monosac_tall[monosac_tall$Sample_ID %in% c("MDCF_2","RUSF"),]
  # Perform T-test
  monosac_MDCF2_vs_RUSF_stats <- compare_means(ug_monosac ~ Sample_ID, data = monosac_MDCF2_vs_RUSF, group.by = "monosaccharide", method = "t.test", paired = FALSE, p.adjust.method="fdr")
  #Pivot tall df with mean and SD to wide format so that it can be merged with stats tests
  monosac_mean_SD_wide <- pivot_wider(monosac_mean_SD, names_from = Sample_ID, values_from = c("mean","SD")) #Pivot from tall to wide format
  #Determine the fold change between MDCF2 and RUSF mono quants
  monosac_mean_SD_wide$log2MDCF2_RUSF <- log2(monosac_mean_SD_wide$mean_MDCF_2/monosac_mean_SD_wide$mean_RUSF)
  #Joint to mono stats df
  monosac_summary_wide <- left_join(monosac_mean_SD_wide, monosac_MDCF2_vs_RUSF_stats)
  monosac_summary_wide <- left_join(monosac_summary_wide, LOD)
  #Get rid of unwanted columns
  monosac_summary_wide[,c(".y.","group1","group2")] <- NULL
  
  #Write monosaccharide quant info to csv
  write.csv(monosac_summary_wide, "./out/monosac_summary_wide.csv")
  
  # 2. Filter data down to the set of monosacharides with the mean-SD of replicates > LOD in either MDCF2 and RUSF
  monosac_summary_wide$MDCF2_above_LOD <- FALSE #Initiate variable
  monosac_summary_wide[(monosac_summary_wide$mean_MDCF_2 - monosac_summary_wide$SD_MDCF_2) > (monosac_summary_wide$LOD),c("MDCF2_above_LOD")]  <- TRUE #Test to see if mean for each mono in MDCF2 is > LOD
  monosac_summary_wide$RUSF_above_LOD <- FALSE #Initiate variable
  monosac_summary_wide[(monosac_summary_wide$mean_RUSF - monosac_summary_wide$SD_RUSF) > (monosac_summary_wide$LOD),c("RUSF_above_LOD")]  <- TRUE #Test to see if mean for each mono in MDCF2 is > LOD
  monosac_summary_wide$MDCF2orRUSF_above_LOD <- FALSE
  monosac_summary_wide$MDCF2orRUSF_above_LOD <- monosac_summary_wide$RUSF_above_LOD == TRUE | monosac_summary_wide$MDCF2_above_LOD == TRUE
  monosac_above_LOD <- as.list(monosac_summary_wide[monosac_summary_wide$MDCF2orRUSF_above_LOD,c("monosaccharide")])

  # 3. Prepare df for bar plot: Narrow the glycan df down to MDCF2 vs RUSF, RUSF ingredients, and MDCF2 ingredients and set the order
  monosac_mean_SD_MDCF2_RUSF_flt <- monosac_mean_SD[monosac_mean_SD$Sample_ID %in% c("RUSF","MDCF_2"),]
  monosac_mean_SD_MDCF2_RUSF_flt$Sample_ID <- factor(monosac_mean_SD_MDCF2_RUSF_flt$Sample_ID ,levels = c("RUSF","MDCF_2"))
  monosac_mean_SD_RUSF_ing_flt <- monosac_mean_SD[monosac_mean_SD$Sample_ID %in% c("Lentil","Milk_powder","Rice"),]
  monosac_mean_SD_RUSF_ing_flt$Sample_ID <- factor(monosac_mean_SD_RUSF_ing_flt$Sample_ID ,levels = c("Lentil","Milk_powder","Rice"))
  monosac_mean_SD_MDCF2_ing_flt <- monosac_mean_SD[monosac_mean_SD$Sample_ID %in% c("Chickpea","Green_banana","Peanut","Soybean"),]
  monosac_mean_SD_MDCF2_ing_flt$Sample_ID <- factor(monosac_mean_SD_MDCF2_ing_flt$Sample_ID, levels = c("Chickpea","Green_banana","Peanut","Soybean") )
  
  #######################################
  # Make a horizontal bar graph with separate x-axis scales for high and lower monosaccharide quants 
  #######################################
  # A. Select and graph Glucose and Galactose results (upper portion of the plot)
  #####################################################
  
  #Set monosaccharide to show in the figure
  filtered <- c("Glucose","Galactose")
  monosaccharide_above_cuttoff <- monosac_summary_wide[monosac_summary_wide$monosaccharide %in% filtered,]
  monosaccharide_above_cuttoff <- monosaccharide_above_cuttoff[rev(order(monosaccharide_above_cuttoff$mean_RUSF, decreasing=TRUE)),c("monosaccharide")]
  monosac_mean_SD_MDCF2_RUSF <- monosac_mean_SD_MDCF2_RUSF_flt[monosac_mean_SD_MDCF2_RUSF_flt$monosaccharide %in% monosaccharide_above_cuttoff$monosaccharide,] #Display the same set of 10 mons above LOD in all 3 plots
  monosac_mean_SD_RUSF_ing <- monosac_mean_SD_RUSF_ing_flt[monosac_mean_SD_RUSF_ing_flt$monosaccharide %in% monosaccharide_above_cuttoff$monosaccharide,] #Display the same set of 10 mons above LOD in all 3 plots
  monosac_mean_SD_MDCF2_ing <- monosac_mean_SD_MDCF2_ing_flt[monosac_mean_SD_MDCF2_ing_flt$monosaccharide %in% monosaccharide_above_cuttoff$monosaccharide,] #Display the same set of 10 mons above LOD in all 3 plots
  
  #Set the X-axis order
  monosac_mean_SD_MDCF2_RUSF$monosaccharide <- factor(monosac_mean_SD_MDCF2_RUSF$monosaccharide, levels=monosaccharide_above_cuttoff$monosaccharide)
  monosac_mean_SD_RUSF_ing$monosaccharide <- factor(monosac_mean_SD_RUSF_ing$monosaccharide, levels= monosaccharide_above_cuttoff$monosaccharide)
  monosac_mean_SD_MDCF2_ing$monosaccharide <- factor(monosac_mean_SD_MDCF2_ing$monosaccharide, levels = monosaccharide_above_cuttoff$monosaccharide)
  
   # 3. Make monosaccharide bar plots with normal axis
  # MDCF2 vs RUSF
  pdf("./out/monosaccharide_quant_MDCF2_RUSF_GluGal.pdf", width = 2.5, height = 1.4)
  p <- ggplot(monosac_mean_SD_MDCF2_RUSF, aes(x= monosaccharide, y=mean, fill=Sample_ID)) + 
    geom_bar(stat="identity", color="black", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD), width=.4,
                  position=position_dodge(.9)) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 6, vjust = 1, hjust=1)) + 
    theme(axis.text.y = element_text(size = 6)) + 
    theme(legend.position = "none") +
    ylim(NA,1000) +
    scale_fill_manual(values=c("#FFFFFF","#CCCCCC"))
  p + coord_flip()
  dev.off()
  
  # RUSF
  pdf("./out/monosaccharide_quant_RUSF_ing_GluGal.pdf", width = 2.5, height = 1.4)
  p <- ggplot(monosac_mean_SD_RUSF_ing, aes(x= monosaccharide, y=mean, fill=Sample_ID)) + 
    geom_bar(stat="identity", color="black", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD), width=.4,
                  position=position_dodge(.9)) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 6, vjust = 1, hjust=1)) + 
    theme(axis.text.y = element_text(size = 6)) + 
    theme(legend.position = "none") +
    ylim(NA,1000) +
    scale_fill_manual(values=c("#E41A1C","#377EB8","#666666"))
  p + coord_flip()
  dev.off()
  
  # MDCF2 
  pdf("./out/monosaccharide_quant_MDCF2_ing_GluGal.pdf", width = 2.5, height = 1.4)
  p <- ggplot(monosac_mean_SD_MDCF2_ing, aes(x= monosaccharide, y=mean, fill=Sample_ID)) + 
    geom_bar(stat="identity", color="black", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD), width=.4,
                  position=position_dodge(.9)) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 6, vjust = 1, hjust=1)) + 
    theme(axis.text.y = element_text(size = 6)) + 
    theme(legend.position = "none") +
    ylim(NA,1000) +
    scale_fill_manual(values=c("#FFFF33","#4DAF4A","#A65628","#FF7F00"))
  p + coord_flip()
  dev.off()
  
  #######################################
  # B. Select and graph other monosaccharide results (lower portion of the plot)
  #####################################################
  
  #Set monosaccharide to show in the figure
  filtered <- c("Arabinose", "Fructose", "Fucose", "GalA", "Mannose", "Rhamnose", "Ribose", "Xylose") #Monos other than glucose and galactose (lower part of the panel)
  monosaccharide_above_cuttoff <- monosac_summary_wide[monosac_summary_wide$monosaccharide %in% filtered,]
  monosaccharide_above_cuttoff <- monosaccharide_above_cuttoff[rev(order(monosaccharide_above_cuttoff$mean_RUSF, decreasing=TRUE)),c("monosaccharide")]
  monosac_mean_SD_MDCF2_RUSF <- monosac_mean_SD_MDCF2_RUSF_flt[monosac_mean_SD_MDCF2_RUSF_flt$monosaccharide %in% monosaccharide_above_cuttoff$monosaccharide,] #Display the same set of 10 mons above LOD in all 3 plots
  monosac_mean_SD_RUSF_ing <- monosac_mean_SD_RUSF_ing_flt[monosac_mean_SD_RUSF_ing_flt$monosaccharide %in% monosaccharide_above_cuttoff$monosaccharide,] #Display the same set of 10 mons above LOD in all 3 plots
  monosac_mean_SD_MDCF2_ing <- monosac_mean_SD_MDCF2_ing_flt[monosac_mean_SD_MDCF2_ing_flt$monosaccharide %in% monosaccharide_above_cuttoff$monosaccharide,] #Display the same set of 10 mons above LOD in all 3 plots
  
  #Set the X-axis order
  monosac_mean_SD_MDCF2_RUSF$monosaccharide <- factor(monosac_mean_SD_MDCF2_RUSF$monosaccharide, levels=monosaccharide_above_cuttoff$monosaccharide)
  monosac_mean_SD_RUSF_ing$monosaccharide <- factor(monosac_mean_SD_RUSF_ing$monosaccharide, levels= monosaccharide_above_cuttoff$monosaccharide)
  monosac_mean_SD_MDCF2_ing$monosaccharide <- factor(monosac_mean_SD_MDCF2_ing$monosaccharide, levels = monosaccharide_above_cuttoff$monosaccharide)
  
  # 3. Make monosaccharide bar plots with normal axis
  # MDCF2 vs RUSF
  pdf("./out/monosaccharide_quant_MDCF2_RUSF_lower.pdf", width = 2.5, height = 3.3)
  p <- ggplot(monosac_mean_SD_MDCF2_RUSF, aes(x= monosaccharide, y=mean, fill=Sample_ID)) + 
    geom_bar(stat="identity", color="black", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD), width=.4,
                  position=position_dodge(.9)) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 6, vjust = 1, hjust=1)) + 
    theme(axis.text.y = element_text(size = 6)) + 
    theme(legend.position = "none") +
    ylim(NA,40) +
    scale_fill_manual(values=c("#FFFFFF","#CCCCCC"))
  p + coord_flip()
  dev.off()
  
  # RUSF
  pdf("./out/monosaccharide_quant_RUSF_ing_lower.pdf", width = 2.5, height = 3.3)
  p <- ggplot(monosac_mean_SD_RUSF_ing, aes(x= monosaccharide, y=mean, fill=Sample_ID)) + 
    geom_bar(stat="identity", color="black", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD), width=.4,
                  position=position_dodge(.9)) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 6, vjust = 1, hjust=1)) + 
    theme(axis.text.y = element_text(size = 6)) + 
    theme(legend.position = "none") +
    ylim(NA,40) +
    scale_fill_manual(values=c("#E41A1C","#377EB8","#666666"))
  p + coord_flip()
  dev.off()
  
  # MDCF2 
  pdf("./out/monosaccharide_quant_MDCF2_ing_lower.pdf", width = 2.5, height = 3.3)
  p <- ggplot(monosac_mean_SD_MDCF2_ing, aes(x= monosaccharide, y=mean, fill=Sample_ID)) + 
    geom_bar(stat="identity", color="black", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=mean-SD, ymax=mean+SD), width=.4,
                  position=position_dodge(.9)) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 6, vjust = 1, hjust=1)) + 
    theme(axis.text.y = element_text(size = 6)) + 
    theme(legend.position = "none") +
    ylim(NA,40) +
    scale_fill_manual(values=c("#FFFF33","#4DAF4A","#A65628","#FF7F00"))
  p + coord_flip()
  dev.off()