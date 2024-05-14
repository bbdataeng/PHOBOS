# Phobos study statistical analysis

# Libraries ---------------------------------------------------------
library(ComplexHeatmap)
library(tidyverse)
library(dplyr)
library(data.table)
library(varhandle)
library(writexl)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(gridExtra)
library(scales)
library(GGally)
library(reshape2)
library(survminer)
library(survival)
library(rms)
library(corrplot)
library(xlsx)
library(excel.link)
library(magrittr)
library(psych)
library(grid)
library(wesanderson)
library(reshape2)
library(scales)
library(Hmisc)
library(ggsci)
library(viridis)
library(table1)
library(summarytools)
library(pagedown)
library(tidytidbits)
library(ranger)
library(readxl)
library(pec)
library(R2HTML)
library(openxlsx)
library(gtable)
library(glmnet)
library(introdataviz)
library(cowplot)
library(ggpubr)
library(rstatix)
library(circlize)



setwd("path of your own folder")


# Load PHOBOS data and pre-processing -------------------------------

pw_message <- "Enter password for MS Excel sheet. 
               It is also important to note that the password must be entered for the tab that will open in the MS Excel file, otherwise the data will not be loaded."
Data <- xl.read.file("Syntetic_dataset.xlsx", password = rstudioapi::askForPassword(pw_message))

# change the name of the file when you receive it!!!!
# if no password is needed enter space in  the tab and press enter
# Proceed only in case you have entered password in both R and MS Excel tabs.



# Preprocess data into score

Data$pATM_N_p <- cut(Data$pATM_N_p, breaks = c(-Inf,0,10,50,80,+Inf), labels = c(0,1,2,3,4))
Data$pATM_C_p <- cut(Data$pATM_C_p, breaks = c(-Inf,0,10,50,80,+Inf), labels = c(0,1,2,3,4))
Data$pATR_N_p <- cut(Data$pATR_N_p, breaks = c(-Inf,0,10,50,80,+Inf), labels = c(0,1,2,3,4))
Data$pATR_C_p <- cut(Data$pATR_C_p, breaks = c(-Inf,0,10,50,80,+Inf), labels = c(0,1,2,3,4))
Data$pCHK1_N_p <- cut(Data$pCHK1_N_p, breaks = c(-Inf,0,10,50,80,+Inf),labels = c(0,1,2,3,4))
Data$pCHK1_C_p <- cut(Data$pCHK1_C_p, breaks = c(-Inf,0,10,50,80,+Inf), labels = c(0,1,2,3,4))
Data$pH2AX_N_p <- cut(Data$pH2AX_N_p, breaks = c(-Inf,0,10,50,80,+Inf), labels = c(0,1,2,3,4))
Data$pH2AX_C_p <- cut(Data$pH2AX_C_p, breaks = c(-Inf,0,10,50,80,+Inf), labels = c(0,1,2,3,4))
Data$pRPA32_N_p <- cut(Data$pRPA32_N_p, breaks = c(-Inf,0,10,50,80,+Inf), labels = c(0,1,2,3,4))
Data$pRPA32_C_p <- cut(Data$pRPA32_C_p, breaks = c(-Inf,0,10,50,80,+Inf), labels = c(0,1,2,3,4))
Data$pWEE_N_p <- cut(Data$pWEE_N_p, breaks = c(-Inf,0,10,50,80,+Inf), labels = c(0,1,2,3,4))
Data$pWEE_C_p <- cut(Data$pWEE_C_p, breaks = c(-Inf,0,10,50,80,+Inf), labels = c(0,1,2,3,4))

Data$pATM_N_p <- unfactor(Data$pATM_N_p)
Data$pATM_C_p <- unfactor(Data$pATM_C_p)
Data$pATR_N_p <- unfactor(Data$pATR_N_p)
Data$pATR_C_p <- unfactor(Data$pATR_C_p)
Data$pCHK1_N_p <- unfactor(Data$pCHK1_N_p)
Data$pCHK1_C_p <- unfactor(Data$pCHK1_C_p)
Data$pH2AX_N_p <- unfactor(Data$pH2AX_N_p)
Data$pH2AX_C_p <- unfactor(Data$pH2AX_C_p)
Data$pRPA32_N_p <- unfactor(Data$pRPA32_N_p)
Data$pRPA32_C_p <- unfactor(Data$pRPA32_C_p)
Data$pWEE_N_p <- unfactor(Data$pWEE_N_p)
Data$pWEE_C_p <- unfactor(Data$pWEE_C_p)

Data <- mutate(Data, ATM_N = pATM_N_p * pATM_N_s)
Data <- mutate(Data, ATM_C = pATM_C_p * pATM_C_s)
Data <- mutate(Data, ATR_N = pATR_N_p * pATR_N_s)
Data <- mutate(Data, ATR_C = pATR_C_p * pATR_C_s)
Data <- mutate(Data, CHK1_N = pCHK1_N_p * pCHK1_N_s)
Data <- mutate(Data, CHK1_C = pCHK1_C_p * pCHK1_C_s)
Data <- mutate(Data, H2AX_N = pH2AX_N_p * pH2AX_N_s)
Data <- mutate(Data, H2AX_C = pH2AX_C_p * pH2AX_C_s)
Data <- mutate(Data, RPA32_N = pRPA32_N_p * pRPA32_N_s)
Data <- mutate(Data, RPA32_C = pRPA32_C_p * pRPA32_C_s)
Data <- mutate(Data, WEE1_N = pWEE_N_p * pWEE_N_s)
Data <- mutate(Data, WEE1_C = pWEE_C_p * pWEE_C_s)



# Categorizing according to IRS Legend

Data$ATM_N_irs <- cut(Data$ATM_N, breaks = c(-Inf,1,3,8,+Inf),
                      labels = c("Negative","Positive, Weak","Positive, Intermediate","Positive, Strong"))
Data$ATM_C_irs <- cut(Data$ATM_C, breaks = c(-Inf,1,3,8,+Inf),
                      labels = c("Negative","Positive, Weak","Positive, Intermediate","Positive, Strong"))
Data$ATR_N_irs <- cut(Data$ATR_N, breaks = c(-Inf,1,3,8,+Inf),
                      labels = c("Negative","Positive, Weak","Positive, Intermediate","Positive, Strong"))
Data$ATR_C_irs <- cut(Data$ATR_C, breaks = c(-Inf,1,3,8,+Inf),
                      labels = c("Negative","Positive, Weak","Positive, Intermediate","Positive, Strong"))
Data$CHK1_N_irs <- cut(Data$CHK1_N, breaks = c(-Inf,1,3,8,+Inf),
                       labels = c("Negative","Positive, Weak","Positive, Intermediate","Positive, Strong"))
Data$CHK1_C_irs <- cut(Data$CHK1_C, breaks = c(-Inf,1,3,8,+Inf),
                       labels = c("Negative","Positive, Weak","Positive, Intermediate","Positive, Strong"))
Data$H2AX_N_irs <- cut(Data$H2AX_N, breaks = c(-Inf,1,3,8,+Inf),
                       labels = c("Negative","Positive, Weak","Positive, Intermediate","Positive, Strong"))
Data$H2AX_C_irs <- cut(Data$H2AX_C, breaks = c(-Inf,1,3,8,+Inf),
                       labels = c("Negative","Positive, Weak","Positive, Intermediate","Positive, Strong"))
Data$RPA32_N_irs <- cut(Data$RPA32_N, breaks = c(-Inf,1,3,8,+Inf),
                        labels = c("Negative","Positive, Weak","Positive, Intermediate","Positive, Strong"))
Data$RPA32_C_irs <- cut(Data$RPA32_C, breaks = c(-Inf,1,3,8,+Inf),
                        labels = c("Negative","Positive, Weak","Positive, Intermediate","Positive, Strong"))
Data$WEE1_N_irs <- cut(Data$WEE1_N, breaks = c(-Inf,1,3,8,+Inf),
                      labels = c("Negative","Positive, Weak","Positive, Intermediate","Positive, Strong"))
Data$WEE1_C_irs <- cut(Data$WEE1_C, breaks = c(-Inf,1,3,8,+Inf),
                      labels = c("Negative","Positive, Weak","Positive, Intermediate","Positive, Strong"))















# FIG_2 (All surv curves are not completed, go to end fot notes) ----



Data$Subtype <- ifelse(Data$Subtype == "LA_HR+/HER2-", "LA_HR+/HER2-",
                       ifelse(Data$Subtype == "LB_HR+/HER2-", "LB_HR+/HER2-",
                              ifelse(Data$Subtype == "LB_HR+/HER2+", "LB_HR+/HER2+",
                                     ifelse(Data$Subtype == "EN_HR-/HER2+", "EN_HR-/HER2+",
                                            ifelse(Data$Subtype == "TN_HR-/HER2-", "TN_HR-/HER2-", "Unknown")))))

Data$Her2_Level <- ifelse(Data$Her2 == "Her2_pos", "Her2_pos",
                          ifelse(Data$Her2 == "Her2_neg" & Data$Her2_ihc != 0, "Her2_low",
                                 ifelse(Data$Her2 == "Her2_neg" & Data$Her2_ihc == 0, "Her2_null", "Her2_unknown")))

Data$Her2_Level <- factor(Data$Her2_Level, levels = c("Her2_null", "Her2_low", "Her2_pos", "Her2_unknown"))

Data_Her_pos <- Data %>% filter(Her2_Level == "Her2_pos")




# Treatment DFS

DFS <- survfit(Surv(DFS_Time, DFS_State) ~ Treatment, data = Data)

cox_model_DFS <- coxph(Surv(Data$DFS_Time, Data$DFS_State) ~ Treatment, data = Data)
HR1 <-c(round(summary(cox_model_DFS)$conf.int[c(1,3,4)],2),round(summary(cox_model_DFS)$logtest[3],3))

DFS_surv <- ggsurvplot(DFS,
                       xlim = c(0, 336), 
                       xlab = "Months",
                       ylab = "Probability of Disease-Free Survival",
                       size = 0.8,
                       censor.size = 3.5,
                       palette = c("#3551DC", "#DB8623"),
                       break.time.by = 24,
                       conf.int = FALSE,
                       legend.title = "",
                       pval = FALSE,
                       pval.size = 4,
                       risk.table = TRUE,
                       risk.table.height = 0.18,
                       ggtheme = theme_classic(),
                       font.x = c(10, "bold"),
                       font.y = c(10, "bold"),
                       legend.labs = c("EC","D     EC"),
                       risk.table.y.text.col = TRUE,
                       risk.table.title = "N. patients at risk",
                       break.y.by = 0.1,
                       risk.table.y.text = FALSE,
                       risk.table.x.text = TRUE,
                       font.main = c(10,"plain"),
                       fontsize = 3.5,
                       font.legend = c(10,"bold"))

DFS_surv$plot <- DFS_surv$plot + 
        annotate("text", x=0, y=0.075, label=paste0("Hazard Ratio: ",HR1[1],"  (95% CI: ",HR1[2]," - ",HR1[3],")", "\np-value: ", HR1[4]), hjust = 0, fontface = "bold")

DFS_surv$table <- DFS_surv$table + labs(x="Months",y="", title = "N. patients at risk", caption = "DFS - all", fontsize = 10) +
                                   theme(plot.title = element_text(size = 10, face = "bold"),
                                         axis.title = element_text(size = 10, face = "bold"),
                                         plot.margin = margin(0.3, 0.25, 0.25, 0.05, "cm")) 

DFS_complete <- grid.arrange(DFS_surv$plot, DFS_surv$table, ncol = 1, heights = c(3,1))





# Treatment OS

OS <- survfit(Surv(OS_Time, OS_State) ~ Treatment, data = Data)

cox_model_OS <- coxph(Surv(Data$OS_Time, Data$OS_State) ~ Treatment, data = Data)
HR2 <-c(round(summary(cox_model_OS)$conf.int[c(1,3,4)],2),round(summary(cox_model_OS)$logtest[3],3))

OS_surv <- ggsurvplot(OS,
                      xlim = c(0, 336), 
                      xlab = "Months",
                      ylab = "Probability of Survival",
                      size = 0.8,
                      censor.size = 3.5,
                      palette = c("#3551DC", "#DB8623"),
                      break.time.by = 24,
                      conf.int = FALSE,
                      legend.title = "",
                      pval = FALSE,
                      pval.size = 4,
                      risk.table = TRUE,
                      risk.table.height = 0.18,
                      ggtheme = theme_classic(),
                      font.x = c(10, "bold"),
                      font.y = c(10, "bold"),
                      legend.labs = c("EC","D     EC"),
                      risk.table.y.text.col = TRUE,
                      risk.table.title = "N. patients at risk",
                      break.y.by = 0.1,
                      risk.table.y.text = FALSE,
                      risk.table.x.text = TRUE,
                      font.main = c(10,"plain"),
                      fontsize = 3.5,
                      font.legend = c(10,"bold"))

OS_surv$plot <- OS_surv$plot + 
      annotate("text", x=0, y=0.075, label=paste0("Hazard Ratio: ",HR2[1],"  (95% CI: ",HR2[2]," - ",HR2[3],")", "\np-value: ", HR2[4]), hjust = 0, fontface="bold")

OS_surv$table <- OS_surv$table + labs(x="Months",y="",title = "N. patients at risk", caption = "OS - all", fontsize = 10) +
                                 theme(plot.title = element_text(size = 10, face = "bold"),
                                       axis.title = element_text(size = 10, face = "bold"),
                                       plot.margin = margin(0.3, 0.25, 0.25, 0.05, "cm")) 

OS_complete <- grid.arrange(OS_surv$plot, OS_surv$table, ncol = 1, heights = c(3,1))






# Treatment DFS HER2-positive

DFS_pos <- survfit(Surv(DFS_Time, DFS_State) ~ Treatment, data = Data_Her_pos)

cox_model_DFS_pos <- coxph(Surv(Data_Her_pos$DFS_Time, Data_Her_pos$DFS_State) ~ Treatment, data = Data_Her_pos)
HR3 <- c(round(summary(cox_model_DFS_pos)$conf.int[c(1,3,4)],2), round(summary(cox_model_DFS_pos)$logtest[3],3))

DFS_surv_pos <- ggsurvplot(DFS_pos,
                           xlim = c(0, 336), 
                           xlab = "Months",
                           ylab = "Probability of Disease-Free Survival",
                           size = 0.8,
                           censor.size = 3.5,
                           palette = c("#3551DC", "#DB8623"),
                           break.time.by = 24,
                           conf.int = FALSE,
                           legend.title = "",
                           pval = FALSE,
                           pval.size = 4,
                           risk.table = TRUE,
                           risk.table.height = 0.18,
                           ggtheme = theme_classic(),
                           font.x = c(10, "bold"),
                           font.y = c(10, "bold"),
                           legend.labs = c("EC","D     EC"),
                           risk.table.y.text.col = TRUE,
                           risk.table.title = "N. patients at risk",
                           break.y.by = 0.1,
                           risk.table.y.text = FALSE,
                           risk.table.x.text = TRUE,
                           font.main = c(10,"plain"),
                           fontsize = 3.5,
                           font.legend = c(10,"bold"))

DFS_surv_pos$plot <- DFS_surv_pos$plot + 
  annotate("text", x=0, y=0.075, label=paste0("Hazard Ratio: ",HR3[1],"  (95% CI: ",HR3[2]," - ",HR3[3],")", "\np-value: ", HR3[4]), hjust = 0, fontface = "bold")

DFS_surv_pos$table <- DFS_surv_pos$table + labs(x="Months",y="",title = "N. patients at risk", caption = "DFS - Her2 pos", fontsize = 10) +
                                            theme(plot.title = element_text(size = 10, face = "bold"),
                                                  axis.title = element_text(size = 10, face = "bold"),
                                                  plot.margin = margin(0.3, 0.25, 0.25, 0.05, "cm")) 

DFS_complete_pos <- grid.arrange(DFS_surv_pos$plot, DFS_surv_pos$table, ncol = 1, heights = c(3,1))





# Treatment OS HER2-positive

OS_pos <- survfit(Surv(OS_Time, OS_State) ~ Treatment, data = Data_Her_pos)

cox_model_OS_pos <- coxph(Surv(Data_Her_pos$OS_Time, Data_Her_pos$OS_State) ~ Treatment, data = Data_Her_pos)
HR4 <- c(round(summary(cox_model_OS_pos)$conf.int[c(1,3,4)],2),round(summary(cox_model_OS_pos)$logtest[3],3))

OS_surv_pos <- ggsurvplot(OS_pos,
                          xlim = c(0, 336), 
                          xlab = "Months",
                          ylab = "Probability of Survival",
                          size = 0.8,
                          censor.size = 3.5,
                          palette = c("#3551DC", "#DB8623"),
                          break.time.by = 24,
                          conf.int = FALSE,
                          legend.title = "",
                          pval = FALSE,
                          pval.size = 4,
                          risk.table = TRUE,
                          risk.table.height = 0.18,
                          ggtheme = theme_classic(),
                          font.x = c(10, "bold"),
                          font.y = c(10, "bold"),
                          legend.labs = c("EC","D     EC"),
                          risk.table.y.text.col = TRUE,
                          risk.table.title = "N. patients at risk",
                          break.y.by = 0.1,
                          risk.table.y.text = FALSE,
                          risk.table.x.text = TRUE,
                          font.main = c(10,"plain"),
                          fontsize = 3.5,
                          font.legend = c(10,"bold"))

OS_surv_pos$plot <- OS_surv_pos$plot + 
  annotate("text", x=0, y=0.075, label=paste0("Hazard Ratio: ",HR4[1],"  (95% CI: ",HR4[2]," - ",HR4[3],")", "\np-value: ", HR4[4]), hjust = 0, fontface="bold")

OS_surv_pos$table <- OS_surv_pos$table + labs(x="Months",y="",title = "N. patients at risk", caption = "OS - Her2 pos", fontsize = 10) +
                                          theme(plot.title = element_text(size = 10, face = "bold"),
                                                axis.title = element_text(size = 10, face = "bold"),
                                                plot.margin = margin(0.3, 0.25, 0.25, 0.05, "cm")) 

OS_complete_pos <- grid.arrange(OS_surv_pos$plot, OS_surv_pos$table, ncol = 1, heights = c(3,1))





# Frequency plots and correlation plot

Data_sub <- Data[,65:76]
Dataset_c <- as.matrix(Data_sub, "numeric")
Dataset_c <- data.frame(na.omit(Dataset_c))
result <- cor(Dataset_c)
testRes <- cor.mtest(Dataset_c, conf.level = 0.95)
pal_0 <- viridis(8, direction = -1)



Data_sub1 <- Data[,77:88]
Dataset_d <- as.matrix(Data_sub1, "numeric")
Dataset_d <- data.frame(na.omit(Dataset_d))

data_N_irs <- Dataset_d %>% select(seq(1,12,2))
data_C_irs <- Dataset_d %>% select(seq(1,12,2)+1)
data_N <- Dataset_c %>% select(seq(1,12,2))
data_C <- Dataset_c %>% select(seq(1,12,2)+1)






# Plot Nuclear IRS category

data_N_irs_long <- pivot_longer(data_N_irs, cols = everything(), names_to = "Biomarker", values_to = "Value") %>% mutate(Value = as.factor(Value))
data_N_irs_long$Value <- factor(data_N_irs_long$Value, levels = c("Negative", "Positive, Weak", "Positive, Intermediate", "Positive, Strong"))
pal_1 <- viridis(length(levels(data_N_irs_long$Value))+1, direction = -1)
pal_1 <- pal_1[-1]

plot_N_irs_cat <-  ggplot(data_N_irs_long, 
                          aes(x = Biomarker, fill = Value, y = after_stat(count))) +
                          geom_bar() +
                          scale_fill_manual(values = pal_1) +  
                          labs(title = "Nuclear IRS category distribution", x = "Biomarker", y = "Absolute frequency") +
                          theme_bw() +  
                          theme(panel.border = element_blank(),  
                                panel.grid.major = element_blank(),  
                                panel.grid.minor = element_blank(),  
                                strip.background = element_blank(),  
                                axis.line = element_line(color = 'black', linewidth = 1.3), 
                                axis.title.x = element_text(color = "black", face = "bold", size = 13, margin = margin(t = 10)),
                                axis.title.y = element_text(color = "black", face = "bold", size = 13, margin = margin(r = 10)),
                                axis.ticks.x = element_line(color = "black", linewidth = 1),
                                axis.ticks.length.x = (unit(0.12, "cm")),
                                axis.ticks.length.y = (unit(0.12, "cm")),
                                plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  
                                text = element_text(size = 12, face = "bold"), 
                                legend.text = element_text(face = "bold", size = 10),
                                axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 12),  
                                axis.text.y = element_text(color = "black", face = "bold", size = 12),
                                legend.position = "right",
                                legend.direction = "vertical",
                                legend.box = "horizontal") + 
                          scale_x_discrete(labels = colnames(data_N))
plot_N_irs_cat




# Plot Nuclear IRS score

data_N_long <- pivot_longer(data_N, cols = everything(), names_to = "Biomarker", values_to = "Value") %>% mutate(Value = as.factor(Value))
pal_2 <- viridis(length(levels(data_N_long$Value))+1, direction = -1)
pal_2 <- pal_2[-1]

plot_N_irs <-  ggplot(data_N_long, 
                      aes(x = Biomarker, fill = Value, y = after_stat(count))) +
                      geom_bar() +
                      scale_fill_manual(values = pal_2) +  
                      labs(title = "Nuclear IRS score distribution", x = "Biomarker", y = "Absolute frequency") +
                      theme_bw() +  
                      theme(panel.border = element_blank(),  
                            panel.grid.major = element_blank(),  
                            panel.grid.minor = element_blank(),  
                            strip.background = element_blank(),  
                            axis.line = element_line(color = 'black', linewidth = 1.3), 
                            axis.title.x = element_text(color = "black", face = "bold", size = 13, margin = margin(t = 10)),
                            axis.title.y = element_text(color = "black", face = "bold", size = 13, margin = margin(r = 10)),
                            axis.ticks.x = element_line(color = "black", linewidth = 1),
                            axis.ticks.length.x = (unit(0.12, "cm")),
                            axis.ticks.length.y = (unit(0.12, "cm")),
                            plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  
                            text = element_text(size = 12, face = "bold"), 
                            legend.text = element_text(face = "bold", size = 10),
                            axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 12),  
                            axis.text.y = element_text(color = "black", face = "bold", size = 12),
                            legend.position = "right",
                            legend.direction = "vertical",
                            legend.box = "horizontal") + 
                      scale_x_discrete(labels = colnames(data_N))
plot_N_irs






# Plot Cytoplasmic IRS category

data_C_irs_long <- pivot_longer(data_C_irs, cols = everything(), names_to = "Biomarker", values_to = "Value") %>% mutate(Value = as.factor(Value))
data_C_irs_long$Value <- factor(data_C_irs_long$Value, levels = c("Negative", "Positive, Weak", "Positive, Intermediate", "Positive, Strong"))

plot_C_irs_cat <-  ggplot(data_C_irs_long, 
                          aes(x = Biomarker, fill = Value, y = after_stat(count))) +
                          geom_bar() +
                          scale_fill_manual(values = pal_1) +  
                          labs(title = "Cytoplasmic IRS category distribution", x = "Biomarker", y = "Absolute frequency") +
                          theme_bw() +  
                          theme(panel.border = element_blank(),  
                                panel.grid.major = element_blank(),  
                                panel.grid.minor = element_blank(),  
                                strip.background = element_blank(),  
                                axis.line = element_line(color = 'black', linewidth = 1.3), 
                                axis.title.x = element_text(color = "black", face = "bold", size = 13, margin = margin(t = 10)),
                                axis.title.y = element_text(color = "black", face = "bold", size = 13, margin = margin(r = 10)),
                                axis.ticks.x = element_line(color = "black", linewidth = 1),
                                axis.ticks.length.x = (unit(0.12, "cm")),
                                axis.ticks.length.y = (unit(0.12, "cm")),
                                plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  
                                text = element_text(size = 12, face = "bold"), 
                                legend.text = element_text(face = "bold", size = 10),
                                axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 12),  
                                axis.text.y = element_text(color = "black", face = "bold", size = 12),
                                legend.position = "right",
                                legend.direction = "vertical",
                                legend.box = "horizontal") + 
                          scale_x_discrete(labels = colnames(data_C))
plot_C_irs_cat




# Plot Cytoplasmic IRS score

data_C_long <- pivot_longer(data_C, cols = everything(), names_to = "Biomarker", values_to = "Value") %>% mutate(Value = as.factor(Value))
pal_3 <- viridis(length(levels(data_C_long$Value))+1, direction = -1)
pal_3 <- pal_3[-1]
plot_C_irs <-  ggplot(data_C_long, 
                      aes(x = Biomarker, fill = Value, y = after_stat(count))) +
                      geom_bar() +
                      scale_fill_manual(values = pal_3) +  
                      labs(title = "Cytoplasmic IRS score distribution", x = "Biomarker", y = "Absolute frequency") +
                      theme_bw() +  
                      theme(panel.border = element_blank(),  
                            panel.grid.major = element_blank(),  
                            panel.grid.minor = element_blank(),  
                            strip.background = element_blank(),  
                            axis.line = element_line(color = 'black', linewidth = 1.3), 
                            axis.title.x = element_text(color = "black", face = "bold", size = 13, margin = margin(t = 10)),
                            axis.title.y = element_text(color = "black", face = "bold", size = 13, margin = margin(r = 10)),
                            axis.ticks.x = element_line(color = "black", linewidth = 1),
                            axis.ticks.length.x = (unit(0.12, "cm")),
                            axis.ticks.length.y = (unit(0.12, "cm")),
                            plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  
                            text = element_text(size = 12, face = "bold"), 
                            legend.text = element_text(face = "bold", size = 11),
                            axis.text.x = element_text(angle = 45, hjust = 1, color = "black", face = "bold", size = 12),  
                            axis.text.y = element_text(color = "black", face = "bold", size = 12),
                            legend.position = "right",
                            legend.direction = "vertical",
                            legend.box = "horizontal") + 
                      scale_x_discrete(labels = colnames(data_C))
plot_C_irs


#Saving plots

svg(filename = "Fig_2I.svg")
cor.plot <- corrplot(result,
                  p.mat = testRes$p, 
                  method = 'color',
                  xlab = substitute(paste(bold('Etichetta X'))),
                  ylab = substitute(paste(bold('Etichetta Y'))),
                  diag = FALSE, 
                  type = 'lower',
                  sig.level = c(0.001, 0.01, 0.05), 
                  pch.cex = 1,
                  insig = 'label_sig', 
                  pch.col = 'black', 
                  order = 'hclust', 
                  hclust.method = 'average',
                  tl.srt = 90, 
                  tl.offset = 0.6, 
                  tl.col = 'black', 
                  tl.cex = 1.1,
                  tl.pos = 'ld', 
                  cl.cex = 1,
                  col = pal_0)     # MODIFICARE PRIMA DI LANCIARE COMANDO, VEDERE "colori.R"
dev.off()


ggsave(filename = "Fig_2E.svg", plot = plot_C_irs_cat, device = "svg", dpi = 600, height = 15, width = 17, units = "cm")
ggsave(filename = "Fig_2F.svg", plot = plot_C_irs, device = "svg", dpi = 600, height = 15, width = 14, units = "cm")
ggsave(filename = "Fig_2G.svg", plot = plot_N_irs_cat, device = "svg", dpi = 600, height = 15, width = 17, units = "cm")
ggsave(filename = "Fig_2H.svg", plot = plot_N_irs, device = "svg", dpi = 600, height = 15, width = 14, units = "cm")








# FIG_3  (Heatmap) --------------------------------------------------



# Pre-process data
Dataset_no_NA <- Data %>% tidyr::drop_na(ATM_N, ATM_C, ATR_N, ATR_C, CHK1_N, CHK1_C,
                                         H2AX_N, H2AX_C, RPA32_N, RPA32_C, WEE1_N, WEE1_C)

Dataset_no_NA$Treatment <- ifelse(Dataset_no_NA$Treatment == "EC", "EC", "D     EC")
Dataset_no_NA$Treatment <- factor(Dataset_no_NA$Treatment, levels = c("EC","D     EC"))
Dataset_no_NA$Stage <- factor(Dataset_no_NA$Stage, levels= c("2","3"))
Dataset_no_NA$Ki67_c <- ifelse(Dataset_no_NA$Ki67 > 20, "High", "Low")
Dataset_h <- Dataset_no_NA %>% select(c(65:76))
Dataset_h <- as.matrix(Dataset_h)


# Annotations
Age <- Dataset_no_NA$Age
Center <- Dataset_no_NA$Center
Menopause <- Dataset_no_NA$Menopause
Stage <- Dataset_no_NA$Stage
N_nr <- Dataset_no_NA$N_nr
Grade <- Dataset_no_NA$G
ER <- Dataset_no_NA$ER
PR <- Dataset_no_NA$PR
Her2 <- Dataset_no_NA$Her2
Her2_level <- Dataset_no_NA$Her2_Level
Ki67 <- Dataset_no_NA$Ki67
Ki67_c <- Dataset_no_NA$Ki67_c
Subtype <- Dataset_no_NA$Subtype
Treatment <- Dataset_no_NA$Treatment
DFS_State <- Dataset_no_NA$DFS_State
DFS_Time <- Dataset_no_NA$DFS_Time
OS_State <- Dataset_no_NA$OS_State
OS_Time <- Dataset_no_NA$OS_Time
Metastasis <- Dataset_no_NA$Metastasis

# Vectorize all future annotations
Age <- as.vector(Age)
Center <- as.vector(Center)
Menopause <- as.vector(Menopause)
Menopause <- as.character(Menopause)
Stage <- as.factor(Stage)
Stage <- as.character(Stage)
N_nr <- as.vector(N_nr)
Grade <- as.vector(Grade)
Grade <- as.character(Grade)
ER <- as.vector(ER)
PR <- as.vector(PR)
Her2 <- as.vector(Her2)
Her2_level <- as.vector(Her2_level)
Her2_level <- factor(Her2_level, levels = c("Her2_null", "Her2_low", "Her2_pos", "Her2_unk"))
Ki67 <- as.vector(Ki67)
Ki67_c <- as.vector(Ki67_c)
Subtype <- as.vector(Subtype)
Subtype <- factor(Subtype, levels = c("LA_HR+/HER2-", "LB_HR+/HER2-", "LB_HR+/HER2+", "EN_HR-/HER2+", "TN_HR-/HER2-", "Unknown"))
Treatment <- as.vector(Treatment)
DFS_State <- as.vector(DFS_State)
DFS_State <- as.character(DFS_State)
DFS_Time <- as.vector(DFS_Time)
OS_State <- as.vector(OS_State)
OS_State <- as.character(OS_State)
OS_Time <- as.vector(OS_Time)
Metastases <- as.vector(Metastasis)
Metastases <- as.character(Metastases)



# Constructing annotations

#An_a <- rowAnnotation(Age = anno_barplot(Age), annotation_name_gp = gpar(fontsize = 8, fontface = "bold"))
#An_b <- HeatmapAnnotation(Center = Center, which = "row", annotation_name_gp = gpar(fontsize = 8, fontface = "bold"))
An_c <- HeatmapAnnotation(Menopause = Menopause, which = "row", gap = unit(2, "cm"), annotation_name_gp = gpar(fontsize = 9, fontface = "bold"), annotation_name_offset = c("0.7cm"), annotation_name_rot = 60)
An_d <- HeatmapAnnotation(Stage = Stage, which = "row", annotation_name_gp = gpar(fontsize = 9, fontface = "bold"), annotation_name_offset = c("0.7cm"), annotation_name_rot = 60)
An_e <- rowAnnotation(N_nr = anno_barplot(N_nr), annotation_name_gp = gpar(fontsize = 9, fontface = "bold"), annotation_name_offset = c("0.7cm"), annotation_name_rot = 60)
An_f <- HeatmapAnnotation(Grade = Grade, which = "row", annotation_name_gp = gpar(fontsize = 9, fontface = "bold"), annotation_name_offset = c("0.7cm"), annotation_name_rot = 60)
An_g <- rowAnnotation(ER = anno_points(ER), annotation_name_gp = gpar(fontsize = 9, fontface = "bold"), annotation_name_rot = 60)
An_h <- rowAnnotation(PR = anno_points(PR), annotation_name_gp = gpar(fontsize = 9, fontface = "bold"), annotation_name_rot = 60)
An_i <- HeatmapAnnotation(Her2 = Her2, which = "row", annotation_name_gp = gpar(fontsize = 9, fontface = "bold"), annotation_name_rot = 60, annotation_name_offset = c("0.7cm"))
An_j <- HeatmapAnnotation(Her2_level = Her2_level, which = "row", annotation_name_gp = gpar(fontsize = 9, fontface = "bold"), annotation_name_offset = c("0.7cm"), annotation_name_rot = 60)
An_k <- rowAnnotation(Ki67 = anno_points(Ki67), annotation_name_gp = gpar(fontsize = 9, fontface = "bold"), annotation_name_rot = 60)
An_l <- HeatmapAnnotation(Ki67_c = Ki67_c, which = "row", annotation_name_gp = gpar(fontsize = 9, fontface = "bold"), annotation_name_offset = c("0.7cm"), annotation_name_rot = 60)
An_m <- HeatmapAnnotation(Subtype = Subtype, annotation_name_gp = gpar(fontsize = 9, fontface = "bold"), annotation_name_offset = c("0.7cm"), annotation_name_rot = 60,
                          col = list(Subtype = c("LA_HR+/HER2-" = "brown",
                                                 "LB_HR+/HER2-" = "orange", 
                                                 "LB_HR+/HER2+" = "yellow",
                                                 "EN_HR-/HER2+" = "green",
                                                 "TN_HR-/HER2-" = "blue",
                                                 "Unknown" = "lightgrey")),
                          which = "row")
An_n <- HeatmapAnnotation(Treatment = Treatment, which = "row", annotation_name_gp = gpar(fontsize = 9, fontface = "bold"), annotation_name_offset = c("0.7cm"), annotation_name_rot = 60)
An_o <- rowAnnotation(DFS_Time = anno_barplot(DFS_Time), annotation_name_gp = gpar(fontsize = 9, fontface = "bold"), annotation_name_rot = 60)
An_p <- HeatmapAnnotation(DFS_State = DFS_State, which = "row", annotation_name_gp = gpar(fontsize = 9, fontface = "bold"), annotation_name_offset = c("0.7cm"), annotation_name_rot = 60)
An_q <- rowAnnotation(OS_Time = anno_barplot(OS_Time), annotation_name_gp = gpar(fontsize = 9, fontface = "bold"), annotation_name_rot = 60)
An_r <- HeatmapAnnotation(OS_State = OS_State, which = "row", annotation_name_gp = gpar(fontsize = 9, fontface = "bold"), annotation_name_offset = c("0.7cm"), annotation_name_rot = 60)
An_s <- HeatmapAnnotation(Metastases = Metastases, which = "row", annotation_name_gp = gpar(fontsize = 9, fontface = "bold"), annotation_name_offset = c("0.7cm"), annotation_name_rot = 60)




# Heatmap and add annotation

pal <- wes_palette("Zissou1", 15, type = "continuous")

A <- ComplexHeatmap::Heatmap(Dataset_h, name = "Expression",
                             column_title = "Biomarkers", 
                             row_title = "Patients",
                             cluster_columns = TRUE,
                             clustering_distance_columns = "euclidean",
                             cluster_rows = TRUE,
                             clustering_distance_rows = "euclidean",
                             row_dend_width = unit(3, "cm"),
                             row_km = 2,
                             column_km = 3,
                             col = pal,
                             width = ncol(Dataset_h)*unit(4, "mm"), 
                             height = nrow(Dataset_h)*unit(1.2, "mm"),
                             column_names_gp = grid::gpar(fontsize = 8, fontface = "bold"),
                             column_names_rot = 60,
                             column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
                             row_title_gp = grid::gpar(fontsize = 12, fontface = "bold"))

A + An_c + An_d + An_e + An_f + An_g + An_h + An_i + An_j + An_k + An_m + An_n + An_o + An_p + An_q + An_r



#tiff("Fig_3.tiff", res = 1000, width = 28, height = 28, units = "cm")
#draw(A + An_c + An_d + An_e + An_f + An_g + An_h + An_i + An_j + An_k + An_m + An_n + An_o + An_p + An_q + An_r)
#dev.off()

svg("Fig_3.svg", width = 11, height = 13)
draw(A + An_c + An_d + An_e + An_f + An_g + An_h + An_i + An_j + An_k + An_m + An_n + An_o + An_p + An_q + An_r)
dev.off()






# Elastic-net with ridge penalization --------------------------------

ms_dataset <- data.frame("id"= seq(1:dim(Dataset_no_NA)[[1]]),Dataset_no_NA[,c(17,20:21,18:19,65:76,89,30:31,6)])
ms_dataset$RT_Adj <- as.factor(ms_dataset$RT_Adj)
ms_dataset$HT_Adj <- as.factor(ms_dataset$HT_Adj)
ms_dataset <- ms_dataset[-which(is.na(ms_dataset[,3])==TRUE),]




df <- ms_dataset[-193,7:18]                                                                                             # this observation was deleted due to the values reported in DFS_Time(value set to 0, this value does not allow to estimate the model)
new_df <- data.frame(matrix(nrow = nrow(df), ncol = 0))
for(i in 1:ncol(df)) {
  for(j in i:ncol(df)) {
    new_df[paste0(colnames(df)[i],"__",colnames(df)[j])] <- df[,i] * df[,j]
  }
}


new_df_2 <- cbind(df,new_df)
y_cox_en <- matrix(c(as.numeric(ms_dataset[-193,3]),as.numeric(ms_dataset[-193,4])),ncol=2)
colnames(y_cox_en) <- c("time","status")
x_cox_en <- as.matrix(scale(new_df_2))
x_cox_en <- replace(x_cox_en,which(is.na(x_cox_en)),0)



# Loop for elnet 

fit_cox_en <- glmnet(x_cox_en, y_cox_en, family = "cox", alpha = 0)
plot(fit_cox_en, xvar = "lambda")

set.seed(1234)
cvfit <- cv.glmnet(x_cox_en, y_cox_en, family = "cox",
                   type.measure = "C",
                   nfolds = 10,
                   alpha = 0)
plot(cvfit)
cvfit$lambda.min
cvfit$lambda.1se
cvfit 

fit_cox_en_1 <- glmnet(x_cox_en, y_cox_en, family = "cox", alpha = 0, lambda = cvfit$lambda.min)

elnet_coeff <- data.frame(as.matrix(round(coef(fit_cox_en_1),2)), var=rownames(coef(fit_cox_en_1)))
rownames(elnet_coeff) <- NULL

write.xlsx(elnet_coeff, paste0("elnet_ridge_coeff.xlsx"))

#test_elnet <- quantile(abs(elnet_coeff[,1]),seq(0,1,0.01))
#min_x_set_elnet <- which(abs(elnet_coeff[,1])>test_elnet[67])
#elnet_true <- elnet_coeff[min_x_set_elnet,]
#elnet_true
#x_cox_en <- x_cox_en[,elnet_true[,2]]
#fit_cox_en_2 <- glmnet(x_cox_en, y_cox_en, family = "cox", alpha = 0, lambda = cvfit$lambda.min)
#coef(fit_cox_en_2)

lp <- predict(fit_cox_en_1,
              newx=x_cox_en,
              type="link")

ms_data<-ms_dataset[-193,]
ms_data$prognosis <- ifelse(lp>0,"High Risk","Low Risk")
table(ms_data$prognosis)



#if (length(unique(ms_data$prognosis)) != 2) {
#  ms_data$prognosis  <- c(rep("High Risk",205),rep("Low Risk",6))
#}






# FIG_4 (Circular heatmap and box-plot in elnet) --------------------


# Preprocess data
ms_data_prova <- ms_data

ms_data_prova$Treatment <- as.factor(ifelse(ms_data_prova$Treatment == "EC", "EC", "D     EC"))
ms_data_prova$"Her2 Level" <- as.factor(ifelse(ms_data_prova$Her2_Level == "Her2_pos", "Her2 pos",
                                               ifelse(ms_data_prova$Her2_Level == "Her2_low", "Her2 low",
                                                      ifelse(ms_data_prova$Her2_Level == "Her2_null", "Her2 null", "Her2 unkn"))))
ms_data_prova$"Adjuvant RT" <- as.factor(ifelse(ms_data_prova$RT_Adj == 2, NA, ms_data_prova$RT_Adj))
ms_data_prova$"Adjuvant HT" <- as.factor(ifelse(ms_data_prova$HT_Adj == 2, NA, ms_data_prova$HT_Adj))
ms_data_prova$"Adjuvant RT" <- recode(ms_data_prova$"Adjuvant RT", "2"="No", "3"="Yes")
ms_data_prova$"Adjuvant HT" <- recode(ms_data_prova$"Adjuvant HT", "1"="No", "2"="Yes")
ms_data_prova$"DDR Risk Group" <- as.factor(ms_data_prova$prognosis)


ms_data_long <- ms_data_prova %>% pivot_longer(cols = -c(1:6,19:27),
                                               names_to = "variable",
                                               values_to = "value")

ms_data_long <- ms_data_long %>% mutate(BioType = case_when(str_detect(variable, "_N$") ~ "Nucl",
                                                            str_detect(variable, "_C$") ~ "Cyto",
                                                            TRUE ~ NA_character_))


ms_data_long$prognosis <- as.factor(ms_data_long$prognosis)
ms_data_long$numb <- seq(1,nrow(ms_data_long),1)
ms_data_long$col_risk <- recode(ms_data_long$prognosis, "High Risk"="#A6594B", "Low Risk"="#88B04B")


stat.test <- data.frame(ms_data_long) %>%
                      group_by(BioType,variable)%>% 
                      wilcox_test(value~prognosis) %>% 
                      adjust_pvalue(method = "fdr") %>%
                      add_significance("p.adj")



# Boxplots of biomarkers score in elnet risk groups

rainplot <- ggplot(ms_data_long,aes(x=prognosis,y=value)) +
                  geom_boxplot(aes(fill = prognosis), width = .6, alpha = 0.7, fatten = NULL, show.legend = TRUE, outlier.shape = NA) +
                  geom_point(aes(colour=prognosis,fill = prognosis), position = position_jitter(width = .2), pch = 21, size = .6)+
                  facet_wrap(~BioType+variable, nrow=2, scale="free_x", labeller = labeller(BioType = function(x) "", variable = identity))+
                  ylab('Biomarker score')+
                  xlab('')+
                  theme_cowplot()+
                  scale_y_continuous(breaks = seq(0, 12, by = 3),
                                     limits = c(0, 14.5)) +
                  guides(fill = "none") +
                  scale_fill_manual(values = c("#A6594B","#88B04B")) +
                  scale_color_manual(values = c("#A6594B","#88B04B")) +
                  theme(legend.position = "none",
                        text = element_text(face="bold"),
                        axis.title.x = element_text(size = 12,face="bold"),
                        axis.title.y = element_text(size = 12,face="bold"),
                        strip.text = element_text(size = 14,face="bold"),
                        strip.background = element_blank(),
                        panel.grid.major.y = element_line(linetype = "dashed", colour = "lightgrey"))+
                  geom_segment(aes(x = 1, xend = 2, y = 13.5, yend = 13.5), color = "black")+
                  geom_segment(aes(x = 1, xend = 1, y = 12.5, yend = 13.5), color = "black")+
                  geom_segment(aes(x = 2, xend = 2, y = 12.5, yend = 13.5), color = "black")+
                  geom_text(data = stat.test, aes(label = p.adj.signif, x = 1.5, y = 14.3), size = 5)
rainplot

ggsave("Fig_4B.tiff", plot = rainplot, device = "tiff", width = 32, height = 20, units = "cm", dpi = 600)
ggsave("Fig_4B.svg", plot = rainplot, device = "svg", width = 32, height = 20, units = "cm", dpi = 600)



# Circular hetmap of a specific biological interesting subset of elnet coefficients

mat0 <- matrix(elnet_coeff[,1],ncol=1)
subset_elnet <- which(elnet_coeff[,1]%in% c(-0.04,-0.05,0.04,0.05))
mat1<- as.matrix(mat0[subset_elnet[c(1,2,6,4,9,11,12,17)],])

testo_sostituito1 <- gsub("__", "*", elnet_coeff[subset_elnet[c(1,2,6,4,9,11,12,17)],2])

rownames(mat1) <- testo_sostituito1
colnames(mat1) <- c("Subset of elastic net coefficients")
mat1
col_funny <- colorRamp2(c(min(mat1), 0, max(mat1)), c("#88B04B", "white", "#A6594B"))


svg(filename = "Fig_4A.svg", height = 5, width = 5)
circplot <- circos.heatmap(mat1, col = col_funny, rownames.side = "outside")
circos.clear()
lgd = Legend(title = "Subset of elastic \nnet coefficients", col_fun = col_funny)
grid.draw(lgd)
dev.off()






# FIG_5 (Survival analysis after elnet) -----------------------------

# Survival plot effects on treatment based on elastic net

# DFS elnet - all
fit.surv.DFS <- survfit(Surv(DFS_Time, DFS_State) ~ prognosis, data = ms_data)

cox_model_DFS_en <- coxph(Surv(ms_data$DFS_Time, ms_data$DFS_State) ~ prognosis, data = ms_data)
HR5 <-c(round(summary(cox_model_DFS_en)$conf.int[c(1,3,4)],2),round(summary(cox_model_DFS_en)$logtest[3],3))


DFS_surv_AA <- ggsurvplot(fit.surv.DFS,
                          xlim = c(0, 336), 
                          xlab = "Months",
                          ylab = "Probability of Disease-Free Survival",
                          size = 0.8,
                          censor.size = 3.5,
                          palette = c("#A6594B","#88B04B"),
                          break.time.by = 24,
                          conf.int = FALSE,
                          legend.title = "",
                          pval = FALSE,
                          pval.size = 4,
                          risk.table = TRUE,
                          risk.table.height = 0.18,
                          ggtheme = theme_classic(),
                          font.x = c(10, "bold"),
                          font.y = c(10, "bold"),
                          legend.labs = c("High Risk","Low Risk"),
                          risk.table.y.text.col = TRUE,
                          risk.table.title = "N. patients at risk",
                          break.y.by = 0.1,
                          risk.table.y.text = FALSE,
                          risk.table.x.text = TRUE,
                          font.main = c(10,"plain"),
                          fontsize = 3.5,
                          font.legend = c(10,"bold"))

DFS_surv_AA$plot <- DFS_surv_AA$plot + 
  annotate("text", x=0, y=0.075, label=paste0("Hazard Ratio: ",HR5[1],"  (95% CI: ",HR5[2]," - ",HR5[3],")", "\np-value: ", "< 0.0001"), hjust = 0, fontface="bold")

DFS_surv_AA$table <- DFS_surv_AA$table + labs(x="Months",y="",title = "N. patients at risk",
                                              caption = "DFS risk groups based on elastic net", fontsize = 10) +
                                          theme(plot.title = element_text(size = 10, face = "bold"),
                                                axis.title = element_text(size = 10, face = "bold"),
                                                plot.margin = margin(0.3, 0.25, 0.25, 0.05, "cm")) 

DFS_surv_AA_complete <- grid.arrange(DFS_surv_AA$plot, DFS_surv_AA$table, ncol = 1, heights = c(3,1))





# OS elnet - all
fit.surv.OS <- survfit(Surv(OS_Time, OS_State) ~ prognosis, data = ms_data)

cox_model_OS_en <- coxph(Surv(ms_data$OS_Time, ms_data$OS_State) ~ prognosis, data = ms_data)
HR6 <-c(round(summary(cox_model_OS_en)$conf.int[c(1,3,4)],2),round(summary(cox_model_OS_en)$logtest[3],3))

OS_surv_AA <- ggsurvplot(fit.surv.OS,
                         xlim = c(0, 336), 
                         xlab = "Months",
                         ylab = "Probability of Survival",
                         size = 0.8,
                         censor.size = 3.5,
                         palette = c("#A6594B","#88B04B"),
                         break.time.by = 24,
                         conf.int = FALSE,
                         legend.title = "",
                         pval = FALSE,
                         pval.size = 4,
                         risk.table = TRUE,
                         risk.table.height = 0.18,
                         ggtheme = theme_classic(),
                         font.x = c(10, "bold"),
                         font.y = c(10, "bold"),
                         legend.labs = c("High Risk","Low Risk"),
                         risk.table.y.text.col = TRUE,
                         risk.table.title = "N. patients at risk",
                         break.y.by = 0.1,
                         risk.table.y.text = FALSE,
                         risk.table.x.text = TRUE,
                         font.main = c(10,"plain"),
                         fontsize = 3.5,
                         font.legend = c(10,"bold"))

OS_surv_AA$plot <- OS_surv_AA$plot + 
  annotate("text", x=0, y=0.075, label=paste0("Hazard Ratio: ",HR6[1],"  (95% CI: ",HR6[2]," - ",HR6[3],")", "\np-value: ", HR6[4]), hjust = 0, fontface="bold")

OS_surv_AA$table <- OS_surv_AA$table + labs(x="Months",y="",title = "N. patients at risk",
                                            caption = "OS risk groups based on elastic net", fontsize = 10) +
                                        theme(plot.title = element_text(size = 10, face = "bold"),
                                              axis.title = element_text(size = 10, face = "bold"),
                                              plot.margin = margin(0.3, 0.25, 0.25, 0.05, "cm")) 

OS_surv_AA_complete <- grid.arrange(OS_surv_AA$plot, OS_surv_AA$table, ncol = 1, heights = c(3,1))






# Survival plot effects on treatment for DDR risk - High Risk
ms_data_High <- ms_data %>% dplyr::filter(prognosis == "High Risk")

# DFS elnet - High
fit_DFS_High_elnet <- survfit(Surv(DFS_Time, DFS_State) ~ Treatment, data=ms_data_High)

cox_model_DFS_High_elnet <- coxph(Surv(ms_data_High$DFS_Time, ms_data_High$DFS_State) ~ Treatment, data = ms_data_High)
HR7 <-c(round(summary(cox_model_DFS_High_elnet)$conf.int[c(1,3,4)],2),round(summary(cox_model_DFS_High_elnet)$logtest[3],3))

DFS_High_elnet <- ggsurvplot(fit_DFS_High_elnet,
                             xlim = c(0, 336), 
                             xlab = "Months",
                             ylab = "Probability of Disease-Free Survival",
                             size = 0.8,
                             censor.size = 3.5,
                             palette = c("#3551DC", "#DB8623"),
                             break.time.by = 24,
                             conf.int = FALSE,
                             legend.title = "",
                             pval = FALSE,
                             pval.size = 4,
                             risk.table = TRUE,
                             risk.table.height = 0.18,
                             ggtheme = theme_classic(),
                             font.x = c(10, "bold"),
                             font.y = c(10, "bold"),
                             legend.labs = c("EC","D     EC"),
                             risk.table.y.text.col = TRUE,
                             risk.table.title = "N. patients at risk",
                             break.y.by = 0.1,
                             risk.table.y.text = FALSE,
                             risk.table.x.text = TRUE,
                             font.main = c(10,"plain"),
                             fontsize = 3.5,
                             font.legend = c(10,"bold"))

DFS_High_elnet$plot <- DFS_High_elnet$plot + 
  annotate("text", x=0, y=0.075, label=paste0("Hazard Ratio: ",HR7[1],"  (95% CI: ",HR7[2]," - ",HR7[3],")", "\np-value: ", HR7[4]), hjust = 0, fontface="bold")

DFS_High_elnet$table <- DFS_High_elnet$table + labs(x="Months",y="",title = "N. patients at risk",
                                                    caption = "DFS for treatment in high risk group based on elastic net", fontsize = 10) +
                                               theme(plot.title = element_text(size = 10, face = "bold"),
                                                     axis.title = element_text(size = 10, face = "bold"),
                                                     plot.margin = margin(0.3, 0.25, 0.25, 0.05, "cm")) 

DFS_complete_High_elnet <- grid.arrange(DFS_High_elnet$plot, DFS_High_elnet$table, ncol = 1, heights = c(3,1))





# OS elnet - High
fit_OS_High_elnet <- survfit(Surv(OS_Time, OS_State) ~ Treatment, data=ms_data_High)

cox_model_OS_High_elnet <- coxph(Surv(ms_data_High$OS_Time, ms_data_High$OS_State) ~ Treatment, data = ms_data_High)
HR8 <- c(round(summary(cox_model_OS_High_elnet)$conf.int[c(1,3,4)],2),round(summary(cox_model_OS_High_elnet)$logtest[3],3))

OS_High_elnet <- ggsurvplot(fit_OS_High_elnet,
                            xlim = c(0, 336), 
                            xlab = "Months",
                            ylab = "Probability of Survival",
                            size = 0.8,
                            censor.size = 3.5,
                            palette = c("#3551DC", "#DB8623"),
                            break.time.by = 24,
                            conf.int = FALSE,
                            legend.title = "",
                            pval = FALSE,
                            pval.size = 4,
                            risk.table = TRUE,
                            risk.table.height = 0.18,
                            ggtheme = theme_classic(),
                            font.x = c(10, "bold"),
                            font.y = c(10, "bold"),
                            legend.labs = c("EC","D     EC"),
                            risk.table.y.text.col = TRUE,
                            risk.table.title = "N. patients at risk",
                            break.y.by = 0.1,
                            risk.table.y.text = FALSE,
                            risk.table.x.text = TRUE,
                            font.main = c(10,"plain"),
                            fontsize = 3.5,
                            font.legend = c(10,"bold"))

OS_High_elnet$plot <- OS_High_elnet$plot + 
  annotate("text", x=0, y=0.075, label=paste0("Hazard Ratio: ",HR8[1],"  (95% CI: ",HR8[2]," - ",HR8[3],")", "\np-value: ", HR8[4]), hjust = 0, fontface="bold")

OS_High_elnet$table <- OS_High_elnet$table + labs(x="Months",y="",title = "N. patients at risk",
                                                  caption = "OS for treatment in high risk group based on elastic net", fontsize = 10) +
                                             theme(plot.title = element_text(size = 10, face = "bold"),
                                                   axis.title = element_text(size = 10, face = "bold"),
                                                   plot.margin = margin(0.3, 0.25, 0.25, 0.05, "cm")) 

OS_complete_High_elnet <- grid.arrange(OS_High_elnet$plot, OS_High_elnet$table, ncol = 1, heights = c(3,1))














# FIG_7 (Forest plot with Survival Risk after factor analysis) ------


# Multivariate COX estimates for DFS and OS 
fit_forest_DFS <- coxph(Surv(DFS_Time, DFS_State) ~ Treatment + Stage + `Her2 Level` + `DDR Risk Group` + `Adjuvant HT` + `Adjuvant RT`, data = ms_data_prova)
p_DFS <- ggforest(main= "Multivariate COX model - DFS (HR)", fit_forest_DFS, fontsize = 1.1, refLabel = "reference", data = ms_data_prova, cpositions = c(0.01, 0.18, 0.34))+ 
              theme(plot.margin = margin(0.5, 1.5, 0.5, 1, "cm"))


fit_forest_OS <- coxph(Surv(OS_Time, OS_State) ~ Treatment + Stage + `Her2 Level` + `DDR Risk Group` + `Adjuvant HT` + `Adjuvant RT`, data = ms_data_prova)
p_OS <- ggforest(main= "Multivariate COX model - OS (HR)", fit_forest_OS, fontsize = 1.1, refLabel = "reference", data = ms_data_prova, cpositions = c(0.01, 0.18, 0.34))+ 
  theme(plot.margin = margin(0.5, 0.5, 0.5, 1.5, "cm"), title = element_text(face = "bold")) 

p_all <- grid.arrange(p_DFS, p_OS, ncol = 2)

cox.zph(fit_forest_DFS)
cox.zph(fit_forest_OS)


# Saving forest plots
ggsave("Fig_7.tiff", plot = p_all, width = 52, height = 19, units = "cm", dpi = 600)
ggsave("Fig_7.svg", plot = p_all, width = 21, height = 9, device="svg")











# Supp_Fig_1 --------------------------------------------------------


# Survival plot effects on treatment for DDR risk - Her2 Risk
ms_data_Her2neg <- ms_data %>% dplyr::filter(Her2_Level %in% c("Her2_low","Her2_null"))
ms_data_Her2pos <- ms_data %>% dplyr::filter(Her2_Level == "Her2_pos")


# DFS - Her2 negative
fit_DFS_Her2neg_elnet <- survfit(Surv(DFS_Time, DFS_State) ~ Treatment, data=ms_data_Her2neg)

cox_model_DFS_Her2neg_elnet <- coxph(Surv(ms_data_Her2neg$DFS_Time, ms_data_Her2neg$DFS_State) ~ Treatment, data = ms_data_Her2neg)
HR11 <-c(round(summary(cox_model_DFS_Her2neg_elnet)$conf.int[c(1,3,4)],2), round(summary(cox_model_DFS_Her2neg_elnet)$logtest[3],3))

DFS_Her2neg_elnet <- ggsurvplot(fit_DFS_Her2neg_elnet,
                                xlim = c(0, 336), 
                                xlab = "Months",
                                ylab = "Probability of Disease-Free Survival",
                                size = 0.8,
                                censor.size = 3.5,
                                palette = c("#3551DC", "#DB8623"),
                                break.time.by = 24,
                                conf.int = FALSE,
                                legend.title = "",
                                pval = FALSE,
                                pval.size = 4,
                                risk.table = TRUE,
                                risk.table.height = 0.18,
                                ggtheme = theme_classic(),
                                font.x = c(10, "bold"),
                                font.y = c(10, "bold"),
                                legend.labs = c("EC","D     EC"),
                                risk.table.y.text.col = TRUE,
                                risk.table.title = "N. patients at risk",
                                break.y.by = 0.1,
                                risk.table.y.text = FALSE,
                                risk.table.x.text = TRUE,
                                font.main = c(10,"plain"),
                                fontsize = 3.5,
                                font.legend = c(10,"bold"))

DFS_Her2neg_elnet$plot <- DFS_Her2neg_elnet$plot + 
  annotate("text", x=0, y=0.075, label=paste0("Hazard Ratio: ",HR11[1],"  (95% CI: ",HR11[2]," - ",HR11[3],")", "\np-value: ", HR11[4]), hjust = 0, fontface="bold")

DFS_Her2neg_elnet$table <- DFS_Her2neg_elnet$table + labs(x="Months",y="",title = "N. patients at risk",
                                                          caption = "DFS for treatment in Her2 negative group", fontsize = 10) +
                                                      theme(plot.title = element_text(size = 10, face = "bold"),
                                                            axis.title = element_text(size = 10, face = "bold"),
                                                            plot.margin = margin(0.3, 0.25, 0.25, 0.05, "cm")) 

DFS_complete_Her2neg_elnet <- grid.arrange(DFS_Her2neg_elnet$plot, DFS_Her2neg_elnet$table, ncol = 1, heights = c(3,1))





# OS - Her2 negative
fit_OS_Her2neg_elnet <- survfit(Surv(OS_Time, OS_State) ~ Treatment, data=ms_data_Her2neg)

cox_model_OS_Her2neg_elnet <- coxph(Surv(ms_data_Her2neg$OS_Time, ms_data_Her2neg$OS_State) ~ Treatment, data = ms_data_Her2neg)
HR12 <-c(round(summary(cox_model_OS_Her2neg_elnet)$conf.int[c(1,3,4)],2),round(summary(cox_model_OS_Her2neg_elnet)$logtest[3],3))

OS_Her2neg_elnet <- ggsurvplot(fit_OS_Her2neg_elnet,
                               xlim = c(0, 336), 
                               xlab = "Months",
                               ylab = "Probability of Survival",
                               size = 0.8,
                               censor.size = 3.5,
                               palette = c("#3551DC", "#DB8623"),
                               break.time.by = 24,
                               conf.int = FALSE,
                               legend.title = "",
                               pval = FALSE,
                               pval.size = 4,
                               risk.table = TRUE,
                               risk.table.height = 0.18,
                               ggtheme = theme_classic(),
                               font.x = c(10, "bold"),
                               font.y = c(10, "bold"),
                               legend.labs = c("EC","D     EC"),
                               risk.table.y.text.col = TRUE,
                               risk.table.title = "N. patients at risk",
                               break.y.by = 0.1,
                               risk.table.y.text = FALSE,
                               risk.table.x.text = TRUE,
                               font.main = c(10,"plain"),
                               fontsize = 3.5,
                               font.legend = c(10,"bold"))

OS_Her2neg_elnet$plot <- OS_Her2neg_elnet$plot + 
  annotate("text", x=0, y=0.075, label=paste0("Hazard Ratio: ",HR12[1],"  (95% CI: ",HR12[2]," - ",HR12[3],")", "\np-value: ", HR12[4]), hjust = 0, fontface="bold")


OS_Her2neg_elnet$table <- OS_Her2neg_elnet$table + labs(x="Months",y="",title = "N. patients at risk",
                                                        caption = "OS for treatment in Her2 negative group", fontsize = 10) +
                                                    theme(plot.title = element_text(size = 10, face = "bold"),
                                                          axis.title = element_text(size = 10, face = "bold"),
                                                          plot.margin = margin(0.3, 0.25, 0.25, 0.05, "cm")) 

OS_complete_Her2neg_elnet <- grid.arrange(OS_Her2neg_elnet$plot, OS_Her2neg_elnet$table, ncol = 1, heights = c(3,1))








# Supp_Fig_2 --------------------------------------------------------


# Filter Data to exclude any rows with NA in 'Her2' and 'Treatment', and focus only on 'Her2_pos' and 'Her2_neg'
Data_filtered <- Data %>%
  filter(Her2 %in% c("Her2_pos", "Her2_neg", "Her2_null")) %>%
  drop_na(Her2, Treatment)  # removing rows with NA in both columns

# Ensure Treatment and Her2 columns are factors
Data_filtered$Treatment <- as.factor(Data_filtered$Treatment)
Data_filtered$Her2 <- as.factor(Data_filtered$Her2)

# Calculate percentages
Data_summary <- Data_filtered %>%
  group_by(Her2) %>%
  count(Treatment) %>%
  mutate(Percentage = n / sum(n))  # calculate the percentages

# Perform Fisher's exact test
contingency_table <- table(Data_filtered$Her2, Data_filtered$Treatment)
fisher_test_result <- fisher.test(contingency_table)
p_val_filt <- fisher_test_result$p.value

# Define specific color palettes for the treatments
# Select fewer colors to get darker tones. Here, we use only 5 darkest tones out of 9.
colors_EC <- brewer.pal(9, "Blues")[5:9]
colors_EC_D <- brewer.pal(9, "OrRd")[5:9]


# Calculate the formatted p-value string separately
formatted_p_value <- format.pval(p_val_filt, digits = 2, eps = 0.001)

annotation_text <- sprintf("Fisher's exact test\n p-value: %s", formatted_p_value)


# Create the plot
p_sup <- ggplot(Data_summary, aes(x = Her2, y = Percentage, fill = Treatment)) +
              geom_bar(stat = "identity", position = position_fill()) +
              geom_text(aes(label = scales::percent(Percentage)), fontface = "bold",
                        position = position_fill(vjust = 0.5), size = 4) +
              scale_fill_manual(values = c("#88B04B", "#A6594B"),
                                labels = c("EC", "D     EC")) + 
              labs(title = paste("Distribution of Treatments across Her2 Categories"),
                   x = "Her2 Category", 
                   y = "Percentage", 
                   fill = "Treatment") +
              scale_y_continuous(labels = scales::percent) +
              scale_x_discrete(labels = c("Her2 negative", "Her2 positive")) +
              theme_bw() +
              theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 0.6, unit = "cm")),
                    axis.title.x = element_text(size = 12, face = "bold", color = "black", margin = margin(t = 0.4, unit = "cm")),
                    axis.title.y = element_text(size = 12, face = "bold", color = "black"),
                    axis.text = element_text(size = 12, face = "bold"),
                    axis.text.x = element_text(hjust = 0.5, face = "bold", color = "black"),
                    axis.text.y = element_text(face = "bold", color = "black"),
                    legend.title = element_text(size = 12, face = "bold"),
                    legend.text = element_text(size = 12, face = "bold"),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank())+
              annotate("text", x = Inf, y = Inf, hjust = 1.3, vjust = 1.8, 
                       label = annotation_text,
                       size = 4, 
                       fontface = "bold.italic")

print(p_sup)


# Save Supplementary Figure 1
ggsave("Supp_Fig_2.tiff", plot = p_sup, dpi = 600, width = 15, height = 15, units = "cm")  
ggsave("Supp_Fig_2.svg", plot = p_sup, width = 6, height = 6)  









# Supp_Fig_3 --------------------------------------------------------


# Survival plot effects for DDR risk - Her2 Levels
ms_data_Her2_noUKN <- ms_data[-c(which(ms_data$Her2_Level=="Her2_unknown")),]

# DFS - Her2
fit_DFS_Her2_elnet <- survfit(Surv(DFS_Time, DFS_State) ~ Her2_Level, data=ms_data_Her2_noUKN)

cox_model_DFS_Her2_elnet <- coxph(Surv(ms_data_Her2_noUKN$DFS_Time, ms_data_Her2_noUKN$DFS_State) ~ Her2_Level, data = ms_data_Her2_noUKN)
HR15 <-c(round(summary(cox_model_DFS_Her2_elnet)$conf.int[c(1,3,4)],2),round(summary(cox_model_DFS_Her2_elnet)$logtest[3],3))
HER2groupDFS <- summary(cox_model_DFS_Her2_elnet)

DFS_Her2_elnet <- ggsurvplot(fit_DFS_Her2_elnet,
                             xlim = c(0, 336), 
                             xlab = "Months",
                             ylab = "Probability of Disease-Free Survival",
                             size = 0.8,
                             censor.size = 3.5,
                             palette = c("#8F0066", "#00ABAB", "#C2C208"),
                             break.time.by = 24,
                             conf.int = FALSE,
                             legend.title = "",
                             pval = FALSE,
                             pval.size = 4,
                             risk.table = TRUE,
                             
                             risk.table.height = 0.18,
                             ggtheme = theme_classic(),
                             font.x = c(10, "bold"),
                             font.y = c(10, "bold"),
                             legend.labs = c("Her2 null","Her2 low", "Her2 pos"),
                             risk.table.y.text.col = TRUE,
                             risk.table.title = "N. patients at risk",
                             break.y.by = 0.1,
                             risk.table.y.text = FALSE,
                             risk.table.x.text = TRUE,
                             font.main = c(10,"plain"),
                             fontsize = 3.5,
                             font.legend = c(10,"bold"))

HER2groupDFS_df <- round(data.frame("HR"=HER2groupDFS$conf.int[1:2,1],
                                    "lower.95"=HER2groupDFS$conf.int[1:2,3],
                                    "upper.95"=HER2groupDFS$conf.int[1:2,4],
                                    "p.value"= HER2groupDFS$coefficients[1:2,5]),2)
rownames(HER2groupDFS_df) <- c("HER2-low", "HER2-pos")
underTableDFS <- data.frame(rbind("HER2-null"= c("reference"),"p-value"=round(HER2groupDFS$logtest[3],3)))
colnames(underTableDFS) <- NULL

summary_table <- tableGrob(HER2groupDFS_df, theme = ttheme_minimal(core = list(fg_params = list(cex = 0.8)),
                                                                   rowhead = list(fg_params = list(cex = 0.8, fontface="bold")),
                                                                   colhead = list(fg_params = list(cex = 0.8))))
summary_underTableDFS <- tableGrob(underTableDFS, theme = ttheme_minimal(core = list(fg_params = list(cex = 0.8)),
                                                                         rowhead = list(fg_params = list(cex = 0.8, fontface="bold"))))
DFS_Her2_elnet$plot <- DFS_Her2_elnet$plot + 
  annotation_custom(grob = summary_table, xmin = 55, xmax = 156, ymin = 0.2, ymax = 0.31) + 
  annotation_custom(grob = summary_underTableDFS, xmin = 55, xmax = 40, ymin = 0, ymax = 0.1)

DFS_Her2_elnet$table <- DFS_Her2_elnet$table + labs(x="Months",y="",title = "N. patients at risk",
                                                    caption = "DFS for each Her2 level", fontsize = 10) +
                                                theme(plot.title = element_text(size = 10, face = "bold"),
                                                      axis.title = element_text(size = 10, face = "bold"),
                                                      plot.margin = margin(0.3, 0.25, 0.25, 0.05, "cm")) 

DFS_complete_Her2_elnet <- grid.arrange(DFS_Her2_elnet$plot, DFS_Her2_elnet$table, ncol = 1, heights = c(3,1))





# OS - Her2
fit_OS_Her2_elnet <- survfit(Surv(OS_Time, OS_State) ~ Her2_Level, data=ms_data_Her2_noUKN)

cox_model_OS_Her2_elnet <- coxph(Surv(ms_data_Her2_noUKN$OS_Time, ms_data_Her2_noUKN$OS_State) ~ Her2_Level, data = ms_data_Her2_noUKN)
HR16 <-c(round(summary(cox_model_OS_Her2_elnet)$conf.int[c(1,3,4)],2),round(summary(cox_model_OS_Her2_elnet)$logtest[3],3))
HER2group <- summary(cox_model_OS_Her2_elnet)

OS_Her2_elnet <- ggsurvplot(fit_OS_Her2_elnet,
                            xlim = c(0, 336), 
                            xlab = "Months",
                            ylab = "Probability of Survival",
                            size = 0.8,
                            censor.size = 3.5,
                            palette = c("#8F0066", "#00ABAB", "#C2C208"),
                            break.time.by = 24,
                            conf.int = FALSE,
                            legend.title = "",
                            pval = FALSE,
                            pval.size = 4,
                            risk.table = TRUE,
                            risk.table.height = 0.18,
                            ggtheme = theme_classic(),
                            font.x = c(10, "bold"),
                            font.y = c(10, "bold"),
                            legend.labs = c("Her2 null","Her2 low", "Her2 pos"),
                            risk.table.y.text.col = TRUE,
                            risk.table.title = "N. patients at risk",
                            break.y.by = 0.1,
                            risk.table.y.text = FALSE,
                            risk.table.x.text = TRUE,
                            font.main = c(10,"plain"),
                            fontsize = 3.5,
                            font.legend = c(10,"bold"))

HER2group_df <- round(data.frame("HR"=HER2group$conf.int[1:2,1],
                                 "lower.95"=HER2group$conf.int[1:2,3],
                                 "upper.95"=HER2group$conf.int[1:2,4],
                                 "p.value"= HER2group$coefficients[1:2,5]),2)
rownames(HER2group_df) <- c("HER2-low", "HER2-pos")
underTable <- data.frame(rbind("HER2-null"= c("reference"),"p-value"=round(HER2group$logtest[3],3)))
colnames(underTable) <- NULL

summary_table <- tableGrob(HER2group_df, theme = ttheme_minimal(core = list(fg_params = list(cex = 0.8)),
                                                                rowhead = list(fg_params = list(cex = 0.8, fontface="bold")),
                                                                colhead = list(fg_params = list(cex = 0.8))))
summary_underTable <- tableGrob(underTable, theme = ttheme_minimal(core = list(fg_params = list(cex = 0.8)),
                                                                               rowhead = list(fg_params = list(cex = 0.8, fontface="bold"))))
OS_Her2_elnet$plot <- OS_Her2_elnet$plot + 
                          annotation_custom(grob = summary_table, xmin = 55, xmax = 156, ymin = 0.2, ymax = 0.3) + 
                          annotation_custom(grob = summary_underTable, xmin = 55, xmax = 40, ymin = 0, ymax = 0.1)

OS_Her2_elnet$table <- OS_Her2_elnet$table + labs(x="Months",y="",title = "N. patients at risk",
                                                  caption = "OS for each Her2 level", fontsize = 10) +
                                              theme(plot.title = element_text(size = 10, face = "bold"),
                                                    axis.title = element_text(size = 10, face = "bold"),
                                                    plot.margin = margin(0.3, 0.25, 0.25, 0.05, "cm")) 

OS_complete_Her2_elnet <- grid.arrange(OS_Her2_elnet$plot, OS_Her2_elnet$table, ncol = 1, heights = c(3,1))




# Survival plot effects on treatment for DDR risk - Stage 

# DFS - Stage
fit_DFS_stage_elnet <- survfit(Surv(DFS_Time, DFS_State) ~ Stage, data=ms_data)

cox_model_DFS_stage_elnet <- coxph(Surv(ms_data$DFS_Time, ms_data$DFS_State) ~ Stage, data = ms_data)
HR17 <-c(round(summary(cox_model_DFS_stage_elnet)$conf.int[c(1,3,4)],2),round(summary(cox_model_DFS_stage_elnet)$logtest[3],3))

DFS_Stage_elnet <- ggsurvplot(fit_DFS_stage_elnet,
                                xlim = c(0, 336), 
                                xlab = "Months",
                                ylab = "Probability of Disease-Free Survival",
                                size = 0.8,
                                censor.size = 3.5,
                                palette = c("#4F4F4F", "#000000"),
                                break.time.by = 24,
                                conf.int = FALSE,
                                legend.title = "",
                                pval = FALSE,
                                pval.size = 4,
                                risk.table = TRUE,
                                risk.table.height = 0.18,
                                ggtheme = theme_classic(),
                                font.x = c(10, "bold"),
                                font.y = c(10, "bold"),
                                legend.labs = c("Stage 2","Stage 3"),
                                risk.table.y.text.col = TRUE,
                                risk.table.title = "N. patients at risk",
                                break.y.by = 0.1,
                                risk.table.y.text = FALSE,
                                risk.table.x.text = TRUE,
                                font.main = c(10,"plain"),
                                fontsize = 3.5,
                                font.legend = c(10,"bold"))

DFS_Stage_elnet$plot <- DFS_Stage_elnet$plot + 
  annotate("text", x=0, y=0.075, label=paste0("Hazard Ratio: ",HR17[1],"  (95% CI: ",HR17[2]," - ",HR17[3],")", "\np-value: ", HR17[4]), hjust = 0, fontface="bold")

DFS_Stage_elnet$table <- DFS_Stage_elnet$table + labs(x="Months",y="",title = "N. patients at risk",
                                                          caption = "DFS for Stage group", fontsize = 10) +
                                                  theme(plot.title = element_text(size = 10, face = "bold"),
                                                        axis.title = element_text(size = 10, face = "bold"),
                                                        plot.margin = margin(0.3, 0.25, 0.25, 0.05, "cm")) 

DFS_complete_Stage_elnet <- grid.arrange(DFS_Stage_elnet$plot, DFS_Stage_elnet$table, ncol = 1, heights = c(3,1))





# OS - Stage
fit_OS_stage_elnet <- survfit(Surv(OS_Time, OS_State) ~ Stage, data=ms_data)

cox_model_OS_stage_elnet <- coxph(Surv(ms_data$OS_Time, ms_data$OS_State) ~ Stage, data = ms_data)
HR18 <-c(round(summary(cox_model_OS_stage_elnet)$conf.int[c(1,3,4)],2), round(summary(cox_model_OS_stage_elnet)$logtest[3],3))

OS_Stage_elnet <- ggsurvplot(fit_OS_stage_elnet,
                               xlim = c(0, 336), 
                               xlab = "Months",
                               ylab = "Probability of Survival",
                               size = 0.8,
                               censor.size = 3.5,
                               palette = c("#4F4F4F", "#000000"),
                               break.time.by = 24,
                               conf.int = FALSE,
                               legend.title = "",
                               pval = FALSE,
                               pval.size = 4,
                               risk.table = TRUE,
                               risk.table.height = 0.18,
                               ggtheme = theme_classic(),
                               font.x = c(10, "bold"),
                               font.y = c(10, "bold"),
                               legend.labs = c("Stage 2","Stage 3"),
                               risk.table.y.text.col = TRUE,
                               risk.table.title = "N. patients at risk",
                               break.y.by = 0.1,
                               risk.table.y.text = FALSE,
                               risk.table.x.text = TRUE,
                               font.main = c(10,"plain"),
                               fontsize = 3.5,
                               font.legend = c(10,"bold"))

OS_Stage_elnet$plot <- OS_Stage_elnet$plot + 
  annotate("text", x=0, y=0.075, label=paste0("Hazard Ratio: ",HR18[1],"  (95% CI: ",HR18[2]," - ",HR18[3],")", "\np-value: ", "< 0.0001"), hjust = 0, fontface="bold")

OS_Stage_elnet$table <- OS_Stage_elnet$table + labs(x="Months",y="",title = "N. patients at risk",
                                                        caption = "OS for Stage group", fontsize = 10) +
                                                theme(plot.title = element_text(size = 10, face = "bold"),
                                                      axis.title = element_text(size = 10, face = "bold"),
                                                      plot.margin = margin(0.3, 0.25, 0.25, 0.05, "cm")) 

OS_complete_Stage_elnet <- grid.arrange(OS_Stage_elnet$plot, OS_Stage_elnet$table, ncol = 1, heights = c(3,1))









# Supp_Fig_4 --------------------------------------------------------

# Survival plot effects on treatment for DDR risk - Low Risk
ms_data_Low <- ms_data %>% dplyr::filter(prognosis == "Low Risk")

# DFS elnet - Low
fit_DFS_Low_elnet <- survfit(Surv(DFS_Time, DFS_State) ~ Treatment, data=ms_data_Low)

cox_model_DFS_Low_elnet <- coxph(Surv(ms_data_Low$DFS_Time, ms_data_Low$DFS_State) ~ Treatment, data = ms_data_Low)
HR9 <-c(round(summary(cox_model_DFS_Low_elnet)$conf.int[c(1,3,4)],2),round(summary(cox_model_DFS_Low_elnet)$logtest[3],3))
        
DFS_Low_elnet <- ggsurvplot(fit_DFS_Low_elnet,
                            xlim = c(0, 336), 
                            xlab = "Months",
                            ylab = "Probability of Disease-Free Survival",
                            size = 0.8,
                            censor.size = 3.5,
                            palette = c("#3551DC", "#DB8623"),
                            break.time.by = 24,
                            conf.int = FALSE,
                            legend.title = "",
                            pval = FALSE,
                            pval.size = 4,
                            risk.table = TRUE,
                            risk.table.height = 0.18,
                            ggtheme = theme_classic(),
                            font.x = c(10, "bold"),
                            font.y = c(10, "bold"),
                            legend.labs = c("EC","D     EC"),
                            risk.table.y.text.col = TRUE,
                            risk.table.title = "N. patients at risk",
                            break.y.by = 0.1,
                            risk.table.y.text = FALSE,
                            risk.table.x.text = TRUE,
                            font.main = c(10,"plain"),
                            fontsize = 3.5,
                            font.legend = c(10,"bold"))

DFS_Low_elnet$plot <- DFS_Low_elnet$plot + 
  annotate("text", x=0, y=0.075, label=paste0("Hazard Ratio: ",HR9[1],"  (95% CI: ",HR9[2]," - ",HR9[3],")", "\np-value: ", HR9[4]), hjust = 0, fontface="bold")

DFS_Low_elnet$table <- DFS_Low_elnet$table + labs(x="Months",y="",title = "N. patients at risk",
                                                  caption = "DFS for treatment in low risk group based on elastic net", fontsize = 10) +
                                              theme(plot.title = element_text(size = 10, face = "bold"),
                                                    axis.title = element_text(size = 10, face = "bold"),
                                                    plot.margin = margin(0.3, 0.25, 0.25, 0.05, "cm")) 

DFS_complete_Low_elnet <- grid.arrange(DFS_Low_elnet$plot, DFS_Low_elnet$table, ncol = 1, heights = c(3,1))





# OS elnet - Low
fit_OS_Low_elnet <- survfit(Surv(OS_Time, OS_State) ~ Treatment, data=ms_data_Low)

cox_model_OS_Low_elnet <- coxph(Surv(ms_data_Low$OS_Time, ms_data_Low$OS_State) ~ Treatment, data = ms_data_Low)
HR10 <-c(round(summary(cox_model_OS_Low_elnet)$conf.int[c(1,3,4)],2), round(summary(cox_model_OS_Low_elnet)$logtest[3],3))

OS_Low_elnet <- ggsurvplot(fit_OS_Low_elnet,
                           xlim = c(0, 336), 
                           xlab = "Months",
                           ylab = "Probability of Survival",
                           size = 0.8,
                           censor.size = 3.5,
                           palette = c("#3551DC", "#DB8623"),
                           break.time.by = 24,
                           conf.int = FALSE,
                           legend.title = "",
                           pval = FALSE,
                           pval.size = 4,
                           risk.table = TRUE,
                           risk.table.height = 0.18,
                           ggtheme = theme_classic(),
                           font.x = c(10, "bold"),
                           font.y = c(10, "bold"),
                           legend.labs = c("EC","D     EC"),
                           risk.table.y.text.col = TRUE,
                           risk.table.title = "N. patients at risk",
                           break.y.by = 0.1,
                           risk.table.y.text = FALSE,
                           risk.table.x.text = TRUE,
                           font.main = c(10,"plain"),
                           fontsize = 3.5,
                           font.legend = c(10,"bold"))

OS_Low_elnet$plot <- OS_Low_elnet$plot + 
  annotate("text", x=0, y=0.075, label=paste0("Hazard Ratio: ",HR10[1],"  (95% CI: ",HR10[2]," - ",HR10[3],")", "\np-value: ", HR10[4]), hjust = 0, fontface="bold")

OS_Low_elnet$table <- OS_Low_elnet$table + labs(x="Months",y="",title = "N. patients at risk",
                                                caption = "OS for treatment in low risk group based on elastic net", fontsize = 10) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 10, face = "bold"),
        plot.margin = margin(0.3, 0.25, 0.25, 0.05, "cm")) 

OS_complete_Low_elnet <- grid.arrange(OS_Low_elnet$plot, OS_Low_elnet$table, ncol = 1, heights = c(3,1))











# Supp_Fig_5 --------------------------------------------------------

#TCGA-RPPA


RPPA_ATM<- data.frame(read_excel("Prova_10_01_23.xlsx",1))
RPPA2<- data.frame(read_excel("Stage ALL.xlsx",4))
RPPA_STAGE<-merge(RPPA_ATM,RPPA2,by="Patient_ID")


# DFS median RPPA

expr_med<- median(RPPA_STAGE$ATM_ATM,na.rm = TRUE)
RPPA_STAGE$medianATM_ATM <- ifelse(RPPA_STAGE$ATM_ATM > expr_med, "High", "Low")
mod22 <- survfit(Surv(Disease_Free_Months,Disease_Free_Status)~ medianATM_ATM,data=RPPA_STAGE)

cox_model_mod22 <- coxph(Surv(RPPA_STAGE$Disease_Free_Months, RPPA_STAGE$Disease_Free_Status) ~ medianATM_ATM, data = RPPA_STAGE)
HR19 <-c(round(summary(cox_model_mod22)$conf.int[c(1,3,4)],2),round(summary(cox_model_mod22)$logtest[3],3))

DFS_RPPA_ATM <- ggsurvplot(mod22,
                           xlim = c(0, 336), 
                           xlab = "Months",
                           ylab = "Probability of Disease-Free Survival",
                           size = 0.8,
                           censor.size = 3.5,
                           palette = c("#88B04B", "#A6594B"),
                           break.time.by = 24,
                           conf.int = FALSE,
                           legend.title = "",
                           pval = FALSE,
                           pval.size = 4,
                           risk.table = TRUE,
                           risk.table.height = 0.18,
                           ggtheme = theme_classic(),
                           font.x = c(10, "bold"),
                           font.y = c(10, "bold"),
                           legend.labs = c("ATM-High","ATM-Low"),
                           risk.table.y.text.col = TRUE,
                           risk.table.title = "N. patients at risk",
                           break.y.by = 0.1,
                           risk.table.y.text = FALSE,
                           risk.table.x.text = TRUE,
                           font.main = c(10,"plain"),
                           fontsize = 3.5,
                           font.legend = c(10,"bold"))

DFS_RPPA_ATM$plot <- DFS_RPPA_ATM$plot + 
  annotate("text", x=0, y=0.075, label=paste0("Hazard Ratio: ",HR19[1],"  (95% CI: ",HR19[2]," - ",HR19[3],")", "\np-value: ", HR19[4]), hjust = 0, fontface="bold")

DFS_RPPA_ATM$table <- DFS_RPPA_ATM$table + labs(x="Months",y="",title = "N. patients at risk",
                                                caption = "RPPA - DFS based on ATM(median)", fontsize = 10) +
                                            theme(plot.title = element_text(size = 10, face = "bold"),
                                                  axis.title = element_text(size = 10, face = "bold"),
                                                  plot.margin = margin(0.3, 0.25, 0.25, 0.05, "cm")) 

DFS_RPPA_complete <- grid.arrange(DFS_RPPA_ATM$plot, DFS_RPPA_ATM$table, ncol = 1, heights = c(3,1))






# OS median RPPA


expr_med_OS<- median(RPPA_STAGE$ATM_ATM,na.rm = TRUE)
RPPA_STAGE$medianATM_ATM_OS <- ifelse(RPPA_STAGE$ATM_ATM > expr_med_OS, "High", "Low")
mod25 <- survfit(Surv(Overall_Survival_Months,Overall_Survival_Status)~ medianATM_ATM_OS,data=RPPA_STAGE)

cox_model_mod25 <- coxph(Surv(RPPA_STAGE$Overall_Survival_Months, RPPA_STAGE$Overall_Survival_Status) ~ medianATM_ATM_OS, data = RPPA_STAGE)
HR20 <-c(round(summary(cox_model_mod25)$conf.int[c(1,3,4)],2),round(summary(cox_model_mod25)$logtest[3],3))


OS_RPPA_ATM <- ggsurvplot(mod25,
                          xlim = c(0, 336), 
                          xlab = "Months",
                          ylab = "Probability of Survival",
                          size = 0.8,
                          censor.size = 3.5,
                          palette = c("#88B04B", "#A6594B"),
                          break.time.by = 24,
                          conf.int = FALSE,
                          legend.title = "",
                          pval = FALSE,
                          pval.size = 4,
                          risk.table = TRUE,
                          risk.table.height = 0.18,
                          ggtheme = theme_classic(),
                          font.x = c(10, "bold"),
                          font.y = c(10, "bold"),
                          legend.labs = c("ATM-High","ATM-Low"),
                          risk.table.y.text.col = TRUE,
                          risk.table.title = "N. patients at risk",
                          break.y.by = 0.1,
                          risk.table.y.text = FALSE,
                          risk.table.x.text = TRUE,
                          font.main = c(10,"plain"),
                          fontsize = 3.5,
                          font.legend = c(10,"bold"))

OS_RPPA_ATM$plot <- OS_RPPA_ATM$plot + 
  annotate("text", x=0, y=0.075, label=paste0("Hazard Ratio: ",HR20[1],"  (95% CI: ",HR20[2]," - ",HR20[3],")", "\np-value: ", HR20[4]), hjust = 0, fontface="bold")

OS_RPPA_ATM$table <- OS_RPPA_ATM$table + labs(x="Months",y="",title = "N. patients at risk",
                                              caption = "RPPA - OS based on ATM(median)", fontsize = 10) +
                                          theme(plot.title = element_text(size = 10, face = "bold"),
                                                axis.title = element_text(size = 10, face = "bold"),
                                                plot.margin = margin(0.3, 0.25, 0.25, 0.05, "cm")) 

OS_RPPA_complete <- grid.arrange(OS_RPPA_ATM$plot, OS_RPPA_ATM$table, ncol = 1, heights = c(3,1))








# Extras ------------------------------------------------------------

# Calculating Median Follow-Up time
Data_FollowUP_DFS <- Data %>% filter(!is.na(DFS_Time))
Data_FollowUP_DFS$DFS_State_rev <- ifelse(Data_FollowUP_DFS$DFS_State == 1, 0, 1)
surv_object_reverse_DFS <- survfit(Surv(DFS_Time, DFS_State_rev) ~ 1, data = Data_FollowUP_DFS)

Data_FollowUP_OS <- Data %>% filter(!is.na(OS_Time))
Data_FollowUP_OS$OS_State_rev <- ifelse(Data_FollowUP_OS$OS_State == 1, 0, 1)
surv_object_reverse_OS <- survfit(Surv(OS_Time, OS_State_rev) ~ 1, data = Data_FollowUP_OS)

# Print results
surv_object_reverse_DFS
summary(surv_object_reverse_DFS)[[2]]

surv_object_reverse_OS
summary(surv_object_reverse_OS)[[2]]

# Median Follow-Up time (both DFS-OS)
summary(sort(c(summary(surv_object_reverse_DFS)[[2]],summary(surv_object_reverse_OS)[[2]])))







# Adding extra annotations in survival curves -----------------------
library(ggfortify)



surv_info <- list(fortify(surv_object_reverse_DFS),
                  fortify(surv_object_reverse_OS),
                  fortify(DFS),
                  fortify(OS),                  
                  fortify(DFS_pos),             
                  fortify(OS_pos),              
                  fortify(fit.surv.DFS),         
                  fortify(fit.surv.OS),  
                  fortify(fit_DFS_High_elnet),      
                  fortify(fit_OS_High_elnet),      
                  fortify(fit_DFS_Her2neg_elnet),   
                  fortify(fit_OS_Her2neg_elnet), 
                  fortify(fit_DFS_Her2_elnet),      
                  fortify(fit_OS_Her2_elnet),       
                  fortify(fit_DFS_stage_elnet),     
                  fortify(fit_OS_stage_elnet),   
                  fortify(fit_DFS_Low_elnet),       
                  fortify(fit_OS_Low_elnet),        
                  fortify(mod22),            
                  fortify(mod25))



names(surv_info) <- c("data_FollowUp_DFS.xlsx",
                      "data_FollowUp_OS.xlsx",
                      "FIG_2A.xlsx",
                      "FIG_2B.xlsx",
                      "FIG_2C.xlsx",
                      "FIG_2D.xlsx",
                      "FIG_5A.xlsx",
                      "FIG_5B.xlsx",
                      "FIG_5C.xlsx",
                      "FIG_5D.xlsx",
                      "Supp_FIG_1A.xlsx",
                      "Supp_FIG_1B.xlsx",
                      "Supp_FIG_3A.xlsx",
                      "Supp_FIG_3B.xlsx",
                      "Supp_FIG_3C.xlsx",
                      "Supp_FIG_3D.xlsx",
                      "Supp_FIG_4A.xlsx",
                      "Supp_FIG_4B.xlsx",
                      "Supp_FIG_5A.xlsx",
                      "Supp_FIG_5B.xlsx")


surv_info2 <- list()

for (i in 3:length(surv_info)){
  aaa <- surv_info[[i]] %>% group_split(strata)
  names(aaa) <- levels(surv_info[[i]]$strata)
  surv_info2[[i-2]] <- aaa
}


surv_vals <- function(df, var, var2) {
  idx <- which.max(df[[var]][df[[var]] <= 260])
  return(df[[var2]][idx])
}

# Applica la funzione a ciascun data frame nella lista
surv_info_260<- lapply(surv_info2, function(aaa) lapply(aaa, surv_vals, var = "time", var2 = "surv"))


names(surv_info_260) <- c("DFS_complete",
                          "OS_complete",                  
                          "DFS_complete_pos",             
                          "OS_complete_pos",              
                          "DFS_surv_AA_complete",         
                          "OS_surv_AA_complete",  
                          "DFS_complete_High_elnet",      
                          "OS_complete_High_elnet",      
                          "DFS_complete_Her2neg_elnet",   
                          "OS_complete_Her2neg_elnet", 
                          "DFS_complete_Her2_elnet",      
                          "OS_complete_Her2_elnet",       
                          "DFS_complete_Stage_elnet",     
                          "OS_complete_Stage_elnet",   
                          "DFS_complete_Low_elnet",       
                          "OS_complete_Low_elnet",        
                          "DFS_RPPA_complete",            
                          "OS_RPPA_complete")

plot(OS_complete)



# DFS Survival curves

DFS_complete <-grid.arrange(DFS_surv$plot + 
                              geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[1]][1])), 
                                           linetype = "dashed", color = "#3551DC", size = 0.9)+
                              geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[1]][1]), xend = -15, yend = as.numeric(surv_info_260[[1]][1])),
                                           linetype = "dashed", color = "#3551DC", size = 0.9)+
                              geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[1]][2])), 
                                           linetype = "dashed", color = "#DB8623", size = 0.9)+
                              geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[1]][2]), xend = -15, yend = as.numeric(surv_info_260[[1]][2])),
                                           linetype = "dashed", color = "#DB8623", size = 0.9)+
                              geom_segment(aes(x=245, y=0.45, xend=260, yend=as.numeric(surv_info_260[[1]][1])), size=0.5)+ 
                              annotate(geom="text", x=243, y=0.45, label=paste0(round(as.numeric(surv_info_260[[1]][1])*100),"%"), color="#3551DC", hjust="right", size=4) +
                              geom_segment(aes(x=270, y=0.7, xend=260, yend=as.numeric(surv_info_260[[1]][2])), size=0.5)+ 
                              annotate(geom="text", x=272, y=0.7, label=paste0(round(as.numeric(surv_info_260[[1]][2])*100),"%"), color="#DB8623", hjust="left", size=4),
                            DFS_surv$table, ncol = 1, heights = c(3,1))


DFS_complete_pos <-grid.arrange(DFS_surv_pos$plot + 
                                  geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[3]][1])), 
                                               linetype = "dashed", color = "#3551DC", size = 0.9) +
                                  geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[3]][1]), xend = -15, yend = as.numeric(surv_info_260[[3]][1])),
                                               linetype = "dashed", color = "#3551DC", size = 0.9) +
                                  geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[3]][2])), 
                                               linetype = "dashed", color = "#DB8623", size = 0.9) +
                                  geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[3]][2]), xend = -15, yend = as.numeric(surv_info_260[[3]][2])),
                                               linetype = "dashed", color = "#DB8623", size = 0.9) +
                                  geom_segment(aes(x=270, y=0.4, xend=260, yend=as.numeric(surv_info_260[[3]][1])), size=0.5)+ 
                                  annotate(geom="text", x=272, y=0.4, label=paste0(round(as.numeric(surv_info_260[[3]][1])*100),"%"), color="#3551DC", hjust="left", size=4) +
                                  geom_segment(aes(x=270, y=0.7, xend=260, yend=as.numeric(surv_info_260[[3]][2])), size=0.5)+ 
                                  annotate(geom="text", x=272, y=0.7, label=paste0(round(as.numeric(surv_info_260[[3]][2])*100),"%"), color="#DB8623", hjust="left", size=4),
                                DFS_surv_pos$table, ncol = 1, heights = c(3,1))


DFS_surv_AA_complete <-grid.arrange(DFS_surv_AA$plot + 
                                  geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[5]][1])), 
                                               linetype = "dashed", color = "#A6594B", size = 0.9) +
                                  geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[5]][1]), xend = -15, yend = as.numeric(surv_info_260[[5]][1])),
                                               linetype = "dashed", color = "#A6594B", size = 0.9) +
                                  geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[5]][2])), 
                                               linetype = "dashed", color = "#88B04B", size = 0.9) +
                                  geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[5]][2]), xend = -15, yend = as.numeric(surv_info_260[[5]][2])),
                                               linetype = "dashed", color = "#88B04B", size = 0.9) +
                                  geom_segment(aes(x=280, y=0.4, xend=260, yend=as.numeric(surv_info_260[[5]][1])), size=0.5)+ 
                                  annotate(geom="text", x=282, y=0.4, label=paste0(round(as.numeric(surv_info_260[[5]][1])*100),"%"), color="#A6594B", hjust="left", size=4) +
                                  geom_segment(aes(x=280, y=0.8, xend=260, yend=as.numeric(surv_info_260[[5]][2])), size=0.5)+ 
                                  annotate(geom="text", x=282, y=0.8, label=paste0(round(as.numeric(surv_info_260[[5]][2])*100),"%"), color="#88B04B", hjust="left", size=4),
                                  DFS_surv_AA$table, ncol = 1, heights = c(3,1))


DFS_complete_High_elnet <-grid.arrange(DFS_High_elnet$plot + 
                                       geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[7]][1])), 
                                                    linetype = "dashed", color = "#3551DC", size = 0.9) +
                                       geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[7]][1]), xend = -15, yend = as.numeric(surv_info_260[[7]][1])),
                                                    linetype = "dashed", color = "#3551DC", size = 0.9) +
                                       geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[7]][2])), 
                                                    linetype = "dashed", color = "#DB8623", size = 0.9) +
                                       geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[7]][2]), xend = -15, yend = as.numeric(surv_info_260[[7]][2])),
                                                    linetype = "dashed", color = "#DB8623", size = 0.9)+
                                       geom_segment(aes(x=240, y=0.3, xend=260, yend=as.numeric(surv_info_260[[7]][1])), size=0.5)+ 
                                       annotate(geom="text", x=238, y=0.3, label=paste0(round(as.numeric(surv_info_260[[7]][1])*100),"%"), color="#3551DC", hjust="right", size=4) +
                                       geom_segment(aes(x=270, y=0.5, xend=260, yend=as.numeric(surv_info_260[[7]][2])), size=0.5)+ 
                                       annotate(geom="text", x=272, y=0.5, label=paste0(round(as.numeric(surv_info_260[[7]][2])*100),"%"), color="#DB8623", hjust="left", size=4),
                                       DFS_High_elnet$table, ncol = 1, heights = c(3,1))


DFS_complete_Her2neg_elnet <-grid.arrange(DFS_Her2neg_elnet$plot + 
                                   geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[9]][1])), 
                                                linetype = "dashed", color = "#3551DC", size = 0.9) +
                                   geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[9]][1]), xend = -15, yend = as.numeric(surv_info_260[[9]][1])),
                                                linetype = "dashed", color = "#3551DC", size = 0.9) +
                                   geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[9]][2])), 
                                                linetype = "dashed", color = "#DB8623", size = 0.9) +
                                   geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[9]][2]), xend = -15, yend = as.numeric(surv_info_260[[9]][2])),
                                                linetype = "dashed", color = "#DB8623", size = 0.9)+
                                   geom_segment(aes(x=280, y=0.75, xend=260, yend=as.numeric(surv_info_260[[9]][1])), size=0.5)+ 
                                   annotate(geom="text", x=282, y=0.75, label=paste0(round(as.numeric(surv_info_260[[9]][1])*100),"%"), color="#3551DC", hjust="left", size=4) +
                                   geom_segment(aes(x=240, y=0.45, xend=260, yend=as.numeric(surv_info_260[[9]][2])), size=0.5)+ 
                                   annotate(geom="text", x=240, y=0.45, label=paste0(round(as.numeric(surv_info_260[[9]][2])*100),"%"), color="#DB8623", hjust="right", size=4),
                                   DFS_Her2neg_elnet$table, ncol = 1, heights = c(3,1))


DFS_complete_Her2_elnet <-grid.arrange(DFS_Her2_elnet$plot + 
                                   geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[11]][1])), 
                                                linetype = "dashed", color = "#8F0066", size = 0.9) +
                                   geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[11]][1]), xend = -15, yend = as.numeric(surv_info_260[[11]][1])),
                                                linetype = "dashed", color = "#8F0066", size = 0.9) +
                                   geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[11]][2])), 
                                                linetype = "dashed", color = "#00ABAB", size = 0.9) +
                                   geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[11]][2]), xend = -15, yend = as.numeric(surv_info_260[[11]][2])),
                                                linetype = "dashed", color = "#00ABAB", size = 0.9) +
                                   geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[11]][3])), 
                                                linetype = "dashed", color = "#C2C208", size = 0.9) +
                                   geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[11]][3]), xend = -15, yend = as.numeric(surv_info_260[[11]][3])),
                                                linetype = "dashed", color = "#C2C208", size = 0.9) +
                                   geom_segment(aes(x=290, y=0.7, xend=260, yend=as.numeric(surv_info_260[[11]][1])), size=0.5)+ 
                                   annotate(geom="text", x=292, y=0.7, label=paste0(round(as.numeric(surv_info_260[[11]][1])*100),"%"), color="#8F0066", hjust="left", size=4) +
                                   geom_segment(aes(x=290, y=0.6, xend=260, yend=as.numeric(surv_info_260[[11]][2])), size=0.5)+ 
                                   annotate(geom="text", x=292, y=0.6, label=paste0(round(as.numeric(surv_info_260[[11]][2])*100),"%"), color="#00ABAB", hjust="left", size=4) +
                                   geom_segment(aes(x= 290, y=0.5, xend=260, yend=as.numeric(surv_info_260[[11]][3])), size=0.5)+ 
                                   annotate(geom="text", x=292, y=0.5, label=paste0(round(as.numeric(surv_info_260[[11]][3])*100),"%"), color="#C2C208", hjust="left", size=4),
                                   DFS_Her2_elnet$table, ncol = 1, heights = c(3,1))


DFS_complete_Stage_elnet <-grid.arrange(DFS_Stage_elnet$plot + 
                                  geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[13]][1])), 
                                               linetype = "dashed", color = "#4F4F4F", size = 0.9) +
                                  geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[13]][1]), xend = -15, yend = as.numeric(surv_info_260[[13]][1])),
                                               linetype = "dashed", color = "#4F4F4F", size = 0.9) +
                                  geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[13]][2])), 
                                               linetype = "dashed", color = "#000000", size = 0.9) +
                                  geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[13]][2]), xend = -15, yend = as.numeric(surv_info_260[[13]][2])),
                                               linetype = "dashed", color = "#000000", size = 0.9)+
                                  geom_segment(aes(x=275, y=0.75, xend=260, yend=as.numeric(surv_info_260[[13]][1])), size=0.5)+ 
                                  annotate(geom="text", x=277, y=0.75, label=paste0(round(as.numeric(surv_info_260[[13]][1])*100),"%"), color="#4F4F4F", hjust="left", size=4) +
                                  geom_segment(aes(x=275, y=0.55, xend=260, yend=as.numeric(surv_info_260[[13]][2])), size=0.5)+ 
                                  annotate(geom="text", x=277, y=0.55, label=paste0(round(as.numeric(surv_info_260[[13]][2])*100),"%"), color="#000000", hjust="left", size=4),
                                  DFS_Stage_elnet$table, ncol = 1, heights = c(3,1))


DFS_complete_Low_elnet <-grid.arrange(DFS_Low_elnet$plot + 
                                          geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[15]][1])), 
                                                       linetype = "dashed", color = "#3551DC", size = 0.9) +
                                          geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[15]][1]), xend = -15, yend = as.numeric(surv_info_260[[15]][1])),
                                                       linetype = "dashed", color = "#3551DC", size = 0.9) +
                                          geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[15]][2])), 
                                                       linetype = "dashed", color = "#DB8623", size = 0.9) +
                                          geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[15]][2]), xend = -15, yend = as.numeric(surv_info_260[[15]][2])),
                                                       linetype = "dashed", color = "#DB8623", size = 0.9)+
                                          geom_segment(aes(x=275, y=0.8, xend=260, yend=as.numeric(surv_info_260[[15]][1])), size=0.5)+ 
                                          annotate(geom="text", x=277, y=0.8, label=paste0(round(as.numeric(surv_info_260[[15]][1])*100),"%"), color="#3551DC", hjust="left", size=4) +
                                          geom_segment(aes(x=235, y=0.6, xend=260, yend=as.numeric(surv_info_260[[15]][2])), size=0.5)+ 
                                          annotate(geom="text", x=233, y=0.6, label=paste0(round(as.numeric(surv_info_260[[15]][2])*100),"%"), color="#DB8623", hjust="right", size=4),
                                      DFS_Low_elnet$table, ncol = 1, heights = c(3,1))


# OS Survival curves

OS_complete <-grid.arrange(OS_surv$plot + 
                             geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[2]][1])), 
                                          linetype = "dashed", color = "#3551DC", size = 0.9)+
                             geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[2]][1]), xend = -15, yend = as.numeric(surv_info_260[[2]][1])),
                                          linetype = "dashed", color = "#3551DC", size = 0.9)+
                             geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[2]][2])), 
                                          linetype = "dashed", color = "#DB8623", size = 0.9)+
                             geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[2]][2]), xend = -15, yend = as.numeric(surv_info_260[[2]][2])),
                                          linetype = "dashed", color = "#DB8623", size = 0.9)+
                             geom_segment(aes(x=286, y=0.80, xend=260, yend=as.numeric(surv_info_260[[2]][1])), size=0.5)+ 
                             annotate(geom="text", x=288, y=0.8, label=paste0(round(as.numeric(surv_info_260[[2]][1])*100),"%"), color="#3551DC", hjust="left", size=4) +
                             geom_segment(aes(x=242, y=0.5, xend=260, yend=as.numeric(surv_info_260[[2]][2])), size=0.5)+ 
                             annotate(geom="text", x=240, y=0.5, label=paste0(round(as.numeric(surv_info_260[[2]][2])*100),"%"), color="#DB8623", hjust="right", size=4),
                           OS_surv$table, ncol = 1, heights = c(3,1))


OS_complete_pos <-grid.arrange(OS_surv_pos$plot + 
                                 geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[4]][1])), 
                                              linetype = "dashed", color = "#3551DC", size = 0.9) +
                                 geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[4]][1]), xend = -15, yend = as.numeric(surv_info_260[[4]][1])),
                                              linetype = "dashed", color = "#3551DC", size = 0.9) +
                                 geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[4]][2])), 
                                              linetype = "dashed", color = "#DB8623", size = 0.9) +
                                 geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[4]][2]), xend = -15, yend = as.numeric(surv_info_260[[4]][2])),
                                              linetype = "dashed", color = "#DB8623", size = 0.9) +
                                 geom_segment(aes(x=270, y=0.4, xend=260, yend=as.numeric(surv_info_260[[4]][1])), size=0.5)+ 
                                 annotate(geom="text", x=272, y=0.4, label=paste0(round(as.numeric(surv_info_260[[4]][1])*100),"%"), color="#3551DC", hjust="left", size=4) +
                                 geom_segment(aes(x=270, y=0.7, xend=260, yend=as.numeric(surv_info_260[[4]][2])), size=0.5)+ 
                                 annotate(geom="text", x=272, y=0.7, label=paste0(round(as.numeric(surv_info_260[[4]][2])*100),"%"), color="#DB8623", hjust="left", size=4),
                               OS_surv_pos$table, ncol = 1, heights = c(3,1))


OS_surv_AA_complete <-grid.arrange(OS_surv_AA$plot + 
                                     geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[6]][1])), 
                                                  linetype = "dashed", color = "#A6594B", size = 0.9) +
                                     geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[6]][1]), xend = -15, yend = as.numeric(surv_info_260[[6]][1])),
                                                  linetype = "dashed", color = "#A6594B", size = 0.9) +
                                     geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[6]][2])), 
                                                  linetype = "dashed", color = "#88B04B", size = 0.9) +
                                     geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[6]][2]), xend = -15, yend = as.numeric(surv_info_260[[6]][2])),
                                                  linetype = "dashed", color = "#88B04B", size = 0.9) +
                                     geom_segment(aes(x=280, y=0.4, xend=260, yend=as.numeric(surv_info_260[[6]][1])), size=0.5)+ 
                                     annotate(geom="text", x=282, y=0.4, label=paste0(round(as.numeric(surv_info_260[[6]][1])*100),"%"), color="#A6594B", hjust="left", size=4) +
                                     geom_segment(aes(x=280, y=0.8, xend=260, yend=as.numeric(surv_info_260[[6]][2])), size=0.5)+ 
                                     annotate(geom="text", x=282, y=0.8, label=paste0(round(as.numeric(surv_info_260[[6]][2])*100),"%"), color="#88B04B", hjust="left", size=4),
                                   OS_surv_AA$table, ncol = 1, heights = c(3,1))


OS_complete_High_elnet <-grid.arrange(OS_High_elnet$plot + 
                                        geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[8]][1])), 
                                                     linetype = "dashed", color = "#3551DC", size = 0.9) +
                                        geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[8]][1]), xend = -15, yend = as.numeric(surv_info_260[[8]][1])),
                                                     linetype = "dashed", color = "#3551DC", size = 0.9) +
                                        geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[8]][2])), 
                                                     linetype = "dashed", color = "#DB8623", size = 0.9) +
                                        geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[8]][2]), xend = -15, yend = as.numeric(surv_info_260[[8]][2])),
                                                     linetype = "dashed", color = "#DB8623", size = 0.9)+
                                        geom_segment(aes(x=270, y=0.75, xend=260, yend=as.numeric(surv_info_260[[8]][1])), size=0.5)+ 
                                        annotate(geom="text", x=272, y=0.75, label=paste0(round(as.numeric(surv_info_260[[8]][1])*100),"%"), color="#3551DC", hjust="left", size=4) +
                                        geom_segment(aes(x=270, y=0.4, xend=260, yend=as.numeric(surv_info_260[[8]][2])), size=0.5)+ 
                                        annotate(geom="text", x=272, y=0.4, label=paste0(round(as.numeric(surv_info_260[[8]][2])*100),"%"), color="#DB8623", hjust="left", size=4),
                                      OS_High_elnet$table, ncol = 1, heights = c(3,1))


OS_complete_Her2neg_elnet <-grid.arrange(OS_Her2neg_elnet$plot + 
                                           geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[10]][1])), 
                                                        linetype = "dashed", color = "#3551DC", size = 0.9) +
                                           geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[10]][1]), xend = -15, yend = as.numeric(surv_info_260[[10]][1])),
                                                        linetype = "dashed", color = "#3551DC", size = 0.9) +
                                           geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[10]][2])), 
                                                        linetype = "dashed", color = "#DB8623", size = 0.9) +
                                           geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[10]][2]), xend = -15, yend = as.numeric(surv_info_260[[10]][2])),
                                                        linetype = "dashed", color = "#DB8623", size = 0.9)+
                                           geom_segment(aes(x=280, y=0.85, xend=260, yend=as.numeric(surv_info_260[[10]][1])), size=0.5)+ 
                                           annotate(geom="text", x=282, y=0.85, label=paste0(round(as.numeric(surv_info_260[[10]][1])*100),"%"), color="#3551DC", hjust="left", size=4) +
                                           geom_segment(aes(x=280, y=0.68, xend=260, yend=as.numeric(surv_info_260[[10]][2])), size=0.5)+ 
                                           annotate(geom="text", x=282, y=0.68, label=paste0(round(as.numeric(surv_info_260[[10]][2])*100),"%"), color="#DB8623", hjust="left", size=4),
                                         OS_Her2neg_elnet$table, ncol = 1, heights = c(3,1))


OS_complete_Her2_elnet <-grid.arrange(OS_Her2_elnet$plot + 
                                        geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[12]][1])), 
                                                     linetype = "dashed", color = "#8F0066", size = 0.9) +
                                        geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[12]][1]), xend = -15, yend = as.numeric(surv_info_260[[12]][1])),
                                                     linetype = "dashed", color = "#8F0066", size = 0.9) +
                                        geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[12]][2])), 
                                                     linetype = "dashed", color = "#00ABAB", size = 0.9) +
                                        geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[12]][2]), xend = -15, yend = as.numeric(surv_info_260[[12]][2])),
                                                     linetype = "dashed", color = "#00ABAB", size = 0.9) +
                                        geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[12]][3])), 
                                                     linetype = "dashed", color = "#C2C208", size = 0.9) +
                                        geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[12]][3]), xend = -15, yend = as.numeric(surv_info_260[[12]][3])),
                                                     linetype = "dashed", color = "#C2C208", size = 0.9) +
                                        geom_segment(aes(x=290, y=0.82, xend=260, yend=as.numeric(surv_info_260[[12]][1])), size=0.5)+ 
                                        annotate(geom="text", x=292, y=0.82, label=paste0(round(as.numeric(surv_info_260[[12]][1])*100),"%"), color="#8F0066", hjust="left", size=4) +
                                        geom_segment(aes(x=290, y=0.6, xend=260, yend=as.numeric(surv_info_260[[12]][2])), size=0.5)+ 
                                        annotate(geom="text", x=292, y=0.6, label=paste0(round(as.numeric(surv_info_260[[12]][2])*100),"%"), color="#00ABAB", hjust="left", size=4) +
                                        geom_segment(aes(x= 290, y=0.48, xend=260, yend=as.numeric(surv_info_260[[12]][3])), size=0.5)+ 
                                        annotate(geom="text", x=292, y=0.48, label=paste0(round(as.numeric(surv_info_260[[12]][3])*100),"%"), color="#C2C208", hjust="left", size=4),
                                      OS_Her2_elnet$table, ncol = 1, heights = c(3,1))


OS_complete_Stage_elnet <-grid.arrange(OS_Stage_elnet$plot + 
                                         geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[14]][1])), 
                                                      linetype = "dashed", color = "#4F4F4F", size = 0.9) +
                                         geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[14]][1]), xend = -15, yend = as.numeric(surv_info_260[[14]][1])),
                                                      linetype = "dashed", color = "#4F4F4F", size = 0.9) +
                                         geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[14]][2])), 
                                                      linetype = "dashed", color = "#000000", size = 0.9) +
                                         geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[14]][2]), xend = -15, yend = as.numeric(surv_info_260[[14]][2])),
                                                      linetype = "dashed", color = "#000000", size = 0.9)+
                                         geom_segment(aes(x=275, y=0.85, xend=260, yend=as.numeric(surv_info_260[[14]][1])), size=0.5)+ 
                                         annotate(geom="text", x=277, y=0.85, label=paste0(round(as.numeric(surv_info_260[[14]][1])*100),"%"), color="#4F4F4F", hjust="left", size=4) +
                                         geom_segment(aes(x=275, y=0.6, xend=260, yend=as.numeric(surv_info_260[[14]][2])), size=0.5)+ 
                                         annotate(geom="text", x=277, y=0.6, label=paste0(round(as.numeric(surv_info_260[[14]][2])*100),"%"), color="#000000", hjust="left", size=4),
                                       OS_Stage_elnet$table, ncol = 1, heights = c(3,1))


OS_complete_Low_elnet <-grid.arrange(OS_Low_elnet$plot + 
                                       geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[16]][1])), 
                                                    linetype = "dashed", color = "#3551DC", size = 0.9) +
                                       geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[16]][1]), xend = -15, yend = as.numeric(surv_info_260[[16]][1])),
                                                    linetype = "dashed", color = "#3551DC", size = 0.9) +
                                       geom_segment(aes(x = 260, y = 0, xend = 260, yend = as.numeric(surv_info_260[[16]][2])), 
                                                    linetype = "dashed", color = "#DB8623", size = 0.9) +
                                       geom_segment(aes(x = 260, y = as.numeric(surv_info_260[[16]][2]), xend = -15, yend = as.numeric(surv_info_260[[16]][2])),
                                                    linetype = "dashed", color = "#DB8623", size = 0.9)+
                                       geom_segment(aes(x=275, y=0.82, xend=260, yend=as.numeric(surv_info_260[[16]][1])), size=0.5)+ 
                                       annotate(geom="text", x=277, y=0.82, label=paste0(round(as.numeric(surv_info_260[[16]][1])*100),"%"), color="#3551DC", hjust="left", size=4) +
                                       geom_segment(aes(x=275, y=0.6, xend=260, yend=as.numeric(surv_info_260[[16]][2])), size=0.5)+ 
                                       annotate(geom="text", x=277, y=0.6, label=paste0(round(as.numeric(surv_info_260[[16]][2])*100),"%"), color="#DB8623", hjust="left", size=4),
                                     OS_Low_elnet$table, ncol = 1, heights = c(3,1))



ggsave(filename = "Fig_2A.svg" , plot = DFS_complete, device = "svg", dpi = 1000, height = 6, width = 6)
ggsave(filename = "Fig_2B.tiff" , plot = OS_complete, device = "tiff", dpi = 600, height = 6, width = 6)
ggsave(filename = "Fig_2C.tiff" , plot = DFS_complete_pos, device = "tiff", dpi = 600, height = 6, width = 6)
ggsave(filename = "Fig_2D.tiff" , plot = OS_complete_pos, device = "tiff", dpi = 600, height = 6, width = 6)
ggsave(filename = "Fig_5A.tiff" , plot = DFS_surv_AA_complete, device = "tiff", dpi = 600, height = 6, width = 6)
ggsave(filename = "Fig_5B.tiff" , plot = OS_surv_AA_complete, device = "tiff", dpi = 600, height = 6, width = 6)
ggsave(filename = "Fig_5C.tiff" , plot = DFS_complete_High_elnet, device = "tiff", dpi = 600, height = 6, width = 6)
ggsave(filename = "Fig_5D.tiff" , plot = OS_complete_High_elnet, device = "tiff", dpi = 600, height = 6, width = 6)
ggsave(filename = "Supp_Fig_1A.tiff" , plot = DFS_complete_Her2neg_elnet, device = "tiff", dpi = 600, height = 6, width = 6)
ggsave(filename = "Supp_Fig_1B.tiff" , plot = OS_complete_Her2neg_elnet, device = "tiff", dpi = 600, height = 6, width = 6)
ggsave(filename = "Supp_Fig_3A.tiff" , plot = DFS_complete_Her2_elnet, device = "tiff", dpi = 600, height = 6, width = 6)
ggsave(filename = "Supp_Fig_3B.tiff" , plot = OS_complete_Her2_elnet, device = "tiff", dpi = 600, height = 6, width = 6)
ggsave(filename = "Supp_Fig_3C.tiff" , plot = DFS_complete_Stage_elnet, device = "tiff", dpi = 600, height = 6, width = 6)
ggsave(filename = "Supp_Fig_3D.tiff" , plot = OS_complete_Stage_elnet, device = "tiff", dpi = 600, height = 6, width = 6)
ggsave(filename = "Supp_Fig_4A.tiff" , plot = DFS_complete_Low_elnet, device = "tiff", dpi = 600, height = 6, width = 6)
ggsave(filename = "Supp_Fig_4B.tiff" , plot = OS_complete_Low_elnet, device = "tiff", dpi = 600, height = 6, width = 6)
ggsave(filename = "Supp_Fig_5A.tiff" , plot = DFS_RPPA_complete, device = "tiff", dpi = 600, height = 6, width = 6)
ggsave(filename = "Supp_Fig_5B.tiff" , plot = OS_RPPA_complete, device = "tiff", dpi = 600, height = 6, width = 6)


ggsave(filename = "Fig_2A.svg" , plot = DFS_complete, device = "svg", dpi = 1000, height = 6, width = 6)
ggsave(filename = "Fig_2B.svg" , plot = OS_complete, device = "svg", dpi = 600, height = 6, width = 6)
ggsave(filename = "Fig_2C.svg" , plot = DFS_complete_pos, device = "svg", dpi = 600, height = 6, width = 6)
ggsave(filename = "Fig_2D.svg" , plot = OS_complete_pos, device = "svg", dpi = 600, height = 6, width = 6)
ggsave(filename = "Fig_5A.svg" , plot = DFS_surv_AA_complete, device = "svg", dpi = 600, height = 6, width = 6)
ggsave(filename = "Fig_5B.svg" , plot = OS_surv_AA_complete, device = "svg", dpi = 600, height = 6, width = 6)
ggsave(filename = "Fig_5C.svg" , plot = DFS_complete_High_elnet, device = "svg", dpi = 600, height = 6, width = 6)
ggsave(filename = "Fig_5D.svg" , plot = OS_complete_High_elnet, device = "svg", dpi = 600, height = 6, width = 6)
ggsave(filename = "Supp_Fig_1A.svg" , plot = DFS_complete_Her2neg_elnet, device = "svg", dpi = 600, height = 6, width = 6)
ggsave(filename = "Supp_Fig_1B.svg" , plot = OS_complete_Her2neg_elnet, device = "svg", dpi = 600, height = 6, width = 6)
ggsave(filename = "Supp_Fig_3A.svg" , plot = DFS_complete_Her2_elnet, device = "svg", dpi = 600, height = 6, width = 6)
ggsave(filename = "Supp_Fig_3B.svg" , plot = OS_complete_Her2_elnet, device = "svg", dpi = 600, height = 6, width = 6)
ggsave(filename = "Supp_Fig_3C.svg" , plot = DFS_complete_Stage_elnet, device = "svg", dpi = 600, height = 6, width = 6)
ggsave(filename = "Supp_Fig_3D.svg" , plot = OS_complete_Stage_elnet, device = "svg", dpi = 600, height = 6, width = 6)
ggsave(filename = "Supp_Fig_4A.svg" , plot = DFS_complete_Low_elnet, device = "svg", dpi = 600, height = 6, width = 6)
ggsave(filename = "Supp_Fig_4B.svg" , plot = OS_complete_Low_elnet, device = "svg", dpi = 600, height = 6, width = 6)
ggsave(filename = "Supp_Fig_5A.svg" , plot = DFS_RPPA_complete, device = "svg", dpi = 600, height = 6, width = 6)
ggsave(filename = "Supp_Fig_5B.svg" , plot = OS_RPPA_complete, device = "svg", dpi = 600, height = 6, width = 6)
