setwd("~/Documents/QuarantaLab/PC9_cFP_JR/")
# checking to see if diprate package is installed, and install if not detected
# devtools::install_github("QuLab-VU/dipDRC", subdir="diprate", 
#                          user="coreyhayford", auth_token="fdf1ef1048bde5db83cb7c5e4d1d28a15e26068e", 
#                          dependencies=TRUE)
library(lubridate)
library(gplots)
library(diprate)
library(tidyxl)
library(readxl)
require(ggplot2)
require(Hmisc)
library(plyr)
library(reshape2)
library(stringr)
library(ggpubr)
source("SumSE.R") # Summarize function to get the summary statistics;
source("functionsDRC.R")
setwd("PC9_cFP/")
# ============================================================================================================
# ============================================================================================================

#All of these samples were seeded at 1 cell per well, grown for 1 week in RPMI 1640, 
#and 3um Erlotinib was added. The plates were imaged once a day for 1 week.
files = list.files(pattern = "M.csv")
CFPdata = data.frame()
CFPNL2 = data.frame()

#Make dataframe with all 384 well entries and samples from files
for(i in 1:length(files)){
  df = cellCountCV(read.csv(files[i], header=T))
  f = strsplit(files[i],"_")
  df$Sample = as.character(f[[1]][2])
  lwell = unique(df$Well)
  for(k in 1:length(lwell)){
    df[df$Well == lwell[k],"Seq"] = k 
  }
  #Get the log2 normalized counts for EACH well in a sample
  # dfNL2 = compNL2(df, ntimepoint = 3)
  dfNL2 = compNL2(df, ntimepoint = 4)
  #Remove outer wells
  dfNL2 = subset(dfNL2,!(dfNL2$Row%in%c("R01","R16")|dfNL2$Column%in%c("C01","C24")))
  #INF filter
  d <- dfNL2
  dnew <- do.call(data.frame,lapply(d, function(x) replace(x, is.infinite(x),NA)))
  dfNL2 <- dnew[complete.cases(dnew),]
  #Filter out entries that had less than 50 cells initially
  ind <- which(dfNL2[dfNL2$Time==0,]$Count > 50)
  d_norm_i <- dfNL2[dfNL2$Time==0,]
  Well_keep <- d_norm_i[ind,]$Well
  dfNL2_n <- dfNL2[dfNL2$Well %in% Well_keep,]
  #Combine into a master DF
  #CFPdata = rbind(CFPdata,df)
  CFPNL2 = rbind(CFPNL2,dfNL2_n)
  rm(f,df,dfNL2,lwell)
}

#Pruning data for the times that were after the normalization
# adCFPNL2 = subset(CFPNL2, CFPNL2$Time > 40)
adCFPNL2 = subset(CFPNL2, CFPNL2$Time > 65)
adCFPNL2$Time = as.numeric(adCFPNL2$Time)
adCFPNL2$nl2 = as.numeric(adCFPNL2$nl2)
adCFPNL2$Seq = as.numeric(adCFPNL2$Seq)

adCFPNL2_sub = subset(adCFPNL2, !(adCFPNL2$nl2 > 4 & adCFPNL2$Sample == "PC9-DS7"))
adCFPNL2_sub = adCFPNL2_sub[,c("Count","Time","Well","Sample",
                               "Seq", "l2", "nl2")]
adCFPNL2_sub$Sample <- replace(as.character(adCFPNL2_sub$Sample), 
                               adCFPNL2_sub$Sample == "PC9-parental", "PC9-VU")
adCFPNL2_sub$Sample <- replace(as.character(adCFPNL2_sub$Sample), 
                               adCFPNL2_sub$Sample == "PC9-DS1", "DS1")
adCFPNL2_sub$Sample <- replace(as.character(adCFPNL2_sub$Sample), 
                               adCFPNL2_sub$Sample == "PC9-DS3", "DS3")
adCFPNL2_sub$Sample <- replace(as.character(adCFPNL2_sub$Sample), 
                               adCFPNL2_sub$Sample == "PC9-DS4", "DS4")
adCFPNL2_sub$Sample <- replace(as.character(adCFPNL2_sub$Sample), 
                               adCFPNL2_sub$Sample == "PC9-DS6", "DS6")
adCFPNL2_sub$Sample <- replace(as.character(adCFPNL2_sub$Sample), 
                               adCFPNL2_sub$Sample == "PC9-DS7", "DS7")
adCFPNL2_sub$Sample <- replace(as.character(adCFPNL2_sub$Sample), 
                               adCFPNL2_sub$Sample == "PC9-DS9", "DS9")

DS1_times <- round(unique(subset(adCFPNL2_sub, Sample == "DS1")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "DS1")$Time))[1]
DS3_times <- round(unique(subset(adCFPNL2_sub, Sample == "DS3")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "DS3")$Time))[1]
DS4_times <- round(unique(subset(adCFPNL2_sub, Sample == "DS4")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "DS4")$Time))[1]
DS6_times <- round(unique(subset(adCFPNL2_sub, Sample == "DS6")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "DS6")$Time))[1]
DS7_times <- round(unique(subset(adCFPNL2_sub, Sample == "DS7")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "DS7")$Time))[1]
DS9_times <- round(unique(subset(adCFPNL2_sub, Sample == "DS9")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "DS9")$Time))[1]

sublines <- c("DS1", "DS3", "DS4", "DS6", "DS7", "DS9")
cFP_sublines <- subset(adCFPNL2_sub, Sample %in% sublines)

DS1_exp <- subset(cFP_sublines, Sample == "DS1")
DS1_exp$Time <- rep(DS1_times, rep = length(unique(DS1_exp$Time)))
DS1_exp$variable = as.factor(rep(seq(length(unique(DS1_exp$Well))), each = length(unique(DS1_exp$Time))))
colnames(DS1_exp)[4] <- "Subline"
DS1_exp$Type <- "Experimental"
DS1_exp <- DS1_exp[,c("Time", "variable", "nl2", "Subline", "Type")]

DS3_exp <- subset(cFP_sublines, Sample == "DS3")
DS3_exp$Time <- rep(DS3_times, rep = length(unique(DS3_exp$Time)))
DS3_exp$variable = as.factor(rep(seq(length(unique(DS3_exp$Well))), each = length(unique(DS3_exp$Time))))
colnames(DS3_exp)[4] <- "Subline"
DS3_exp$Type <- "Experimental"
DS3_exp <- DS3_exp[,c("Time", "variable", "nl2", "Subline", "Type")]

DS4_exp <- subset(cFP_sublines, Sample == "DS4")
DS4_exp$Time <- rep(DS4_times, rep = length(unique(DS4_exp$Time)))
DS4_exp$variable = as.factor(rep(seq(length(unique(DS4_exp$Well))), each = length(unique(DS4_exp$Time))))
colnames(DS4_exp)[4] <- "Subline"
DS4_exp$Type <- "Experimental"
DS4_exp <- DS4_exp[,c("Time", "variable", "nl2", "Subline", "Type")]

DS6_exp <- subset(cFP_sublines, Sample == "DS6")
DS6_exp$Time <- rep(DS6_times, rep = length(unique(DS6_exp$Time)))
DS6_exp$variable = as.factor(rep(seq(length(unique(DS6_exp$Well))), each = length(unique(DS6_exp$Time))))
colnames(DS6_exp)[4] <- "Subline"
DS6_exp$Type <- "Experimental"
DS6_exp <- DS6_exp[,c("Time", "variable", "nl2", "Subline", "Type")]

DS7_exp <- subset(cFP_sublines, Sample == "DS7")
DS7_exp$Time <- rep(DS7_times, rep = length(unique(DS7_exp$Time)))
DS7_exp$variable = as.factor(rep(seq(length(unique(DS7_exp$Well))), each = length(unique(DS7_exp$Time))))
colnames(DS7_exp)[4] <- "Subline"
DS7_exp$Type <- "Experimental"
DS7_exp <- DS7_exp[,c("Time", "variable", "nl2", "Subline", "Type")]

DS9_exp <- subset(cFP_sublines, Sample == "DS9")
DS9_exp$Time <- rep(DS9_times, rep = length(unique(DS9_exp$Time)))
DS9_exp$variable = as.factor(rep(seq(length(unique(DS9_exp$Well))), each = length(unique(DS9_exp$Time))))
colnames(DS9_exp)[4] <- "Subline"
DS9_exp$Type <- "Experimental"
DS9_exp <- DS9_exp[,c("Time", "variable", "nl2", "Subline", "Type")]


####

setwd("~/git/ThreeStateMelanoma/")
DS1_trajs <- read.csv('trajectories_DS1_G50.csv', row.names = 1)
DS1_trajs <- DS1_trajs[,colSums(is.na(DS1_trajs))<nrow(DS1_trajs)]
names(DS1_trajs) <- seq(ncol(DS1_trajs))
DS1_trajs$Time <- seq(nrow(DS1_trajs)) -1
DS1_trajs <- melt(DS1_trajs, id = c("Time"), value.name = "nl2")
DS1_trajs$Subline <- "DS1"
DS1_trajs$Type <- "Simulation"
DS1_trajs <- subset(DS1_trajs, Time %in% DS1_times)

DS3_trajs <- read.csv('trajectories_DS3_G50.csv', row.names = 1)
DS3_trajs <- DS3_trajs[,colSums(is.na(DS3_trajs))<nrow(DS3_trajs)]
names(DS3_trajs) <- seq(ncol(DS3_trajs))
DS3_trajs$Time <- seq(nrow(DS3_trajs)) -1
DS3_trajs <- melt(DS3_trajs, id = c("Time"), value.name = "nl2")
DS3_trajs$Subline <- "DS3"
DS3_trajs$Type <- "Simulation"
DS3_trajs <- subset(DS3_trajs, Time %in% DS3_times)

DS4_trajs <- read.csv('trajectories_DS4_G50.csv', row.names = 1)
DS4_trajs <- DS4_trajs[,colSums(is.na(DS4_trajs))<nrow(DS4_trajs)]
names(DS4_trajs) <- seq(ncol(DS4_trajs))
DS4_trajs$Time <- seq(nrow(DS4_trajs)) -1
DS4_trajs <- melt(DS4_trajs, id = c("Time"), value.name = "nl2")
DS4_trajs$Subline <- "DS4"
DS4_trajs$Type <- "Simulation"
DS4_trajs <- subset(DS4_trajs, Time %in% DS4_times)


DS6_trajs <- read.csv('trajectories_DS6_G50.csv', row.names = 1)
DS6_trajs <- DS6_trajs[,colSums(is.na(DS6_trajs))<nrow(DS6_trajs)]
names(DS6_trajs) <- seq(ncol(DS6_trajs))
DS6_trajs$Time <- seq(nrow(DS6_trajs)) -1
DS6_trajs <- melt(DS6_trajs, id = c("Time"), value.name = "nl2")
DS6_trajs$Subline <- "DS6"
DS6_trajs$Type <- "Simulation"
DS6_trajs <- subset(DS6_trajs, Time %in% DS6_times)


DS7_trajs <- read.csv('trajectories_DS7_G50.csv', row.names = 1)
DS7_trajs <- DS7_trajs[,colSums(is.na(DS7_trajs))<nrow(DS7_trajs)]
names(DS7_trajs) <- seq(ncol(DS7_trajs))
DS7_trajs$Time <- seq(nrow(DS7_trajs)) -1
DS7_trajs <- melt(DS7_trajs, id = c("Time"), value.name = "nl2")
DS7_trajs$Subline <- "DS7"
DS7_trajs$Type <- "Simulation"
DS7_trajs <- subset(DS7_trajs, Time %in% DS7_times)


DS9_trajs <- read.csv('trajectories_DS9_G50.csv', row.names = 1)
DS9_trajs <- DS9_trajs[,colSums(is.na(DS9_trajs))<nrow(DS9_trajs)]
names(DS9_trajs) <- seq(ncol(DS9_trajs))
DS9_trajs$Time <- seq(nrow(DS9_trajs)) -1
DS9_trajs <- melt(DS9_trajs, id = c("Time"), value.name = "nl2")
DS9_trajs$Subline <- "DS9"
DS9_trajs$Type <- "Simulation"
DS9_trajs <- subset(DS9_trajs, Time %in% DS9_times)

DS3_4_data <- rbind(DS3_trajs, DS4_trajs, DS3_exp, DS4_exp)
DS1_6_7_9_data <- rbind(DS1_trajs, DS6_trajs, DS7_trajs, DS9_trajs,
                         DS1_exp, DS6_exp, DS7_exp, DS9_exp)

#CFP Population Curves for each PC9 cell line
cols_DS3_4 <- c("DS3" = "brown", "DS4" = "deepskyblue")
labs_DS3_4 <- c("DS3", "DS4")

DS3_4_sim <- rbind(DS3_trajs, DS4_trajs)
DS3_4_exp <- rbind(DS3_exp, DS4_exp)

n_DF <- data.frame(Subline = c("DS3", "DS4"),
                   label=c(paste("n =", as.character(length(unique(DS3_exp$variable)))),
                           paste("n =", as.character(length(unique(DS4_exp$variable))))))


ggplot(data = DS3_4_exp, aes(x=Time, y=nl2, group = variable, color = Subline)) +
  geom_line(size = 0.5, alpha = 0.25) +
  labs(x = "Time (hours) Post Drug Penetrance", y = "Normalized Log2 Cell Count") +
  facet_wrap(.~Subline, ncol = 1) +
  geom_text(data = n_DF,
            aes(x = 100, y = 2, label = label),
            inherit.aes=FALSE, parse = FALSE, size = 5) +
  scale_color_manual(values = c("grey", "grey")) +
  ylim(-2.5,2.5) +
  theme_bw() +
  theme(
    # axis.line.x = element_line(colour = "black", size = 1),
    # axis.line.y = element_line(colour = "black", size = 1),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    # legend.position = c(0.95, -0.15),
    # legend.justification = c(1, 0),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
    legend.title = element_text(size=14), axis.title=element_text(size=14),
    strip.text = element_text(size = 14)) +
  ggsave("PC9_cFPs_experimental_grey_DS3_4.pdf", width = 5, height = 8)

kdiv_DF <- data.frame(Subline = c("DS3", "DS4"),
                      label=c(paste("k[division] == 0.030"), #"=", as.character(0.030)),
                              paste("k[division] == 0.035"))) #"=", as.character(0.035))))
kdth_DF <- data.frame(Subline = c("DS3", "DS4"),
                      label=c(paste("k[death] == 0.03075"),# "=", as.character(0.03075)),
                              paste("k[death] == 0.03175")))# "=", as.character(0.03175))))

ggplot(data = DS3_4_sim, aes(x=Time, y=nl2, group = variable, color = Subline)) +
  # geom_line(data = DS3_4_exp, aes(x=Time, y=nl2, group = variable), color = "grey", 
  #           size = 0.5, alpha = 0.25) +
  geom_line(size = 0.5, alpha = 0.25) +
  labs(x = "Time (hours) Post Drug Penetrance", y = "Normalized Log2 Cell Count") +
  geom_text(data= kdiv_DF,
            aes(x = 5,y = 2.3, label = label), 
            inherit.aes=FALSE, parse = TRUE, size = 5,
            hjust = 0) +
  geom_text(data= kdth_DF,
            aes(x = 5,y = 1.8, label = label), 
            inherit.aes=FALSE, parse = TRUE, size = 5,
            hjust = 0) +
  scale_color_manual(values = c("brown", "deepskyblue")) +
  facet_wrap(.~Subline, ncol = 1) +
  ylim(-2.5,2.5) +
  theme_bw() +
  theme(
    # axis.line.x = element_line(colour = "black", size = 1),
    # axis.line.y = element_line(colour = "black", size = 1),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    # legend.position = c(0.95, -0.15),
    # legend.justification = c(1, 0),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
    legend.title = element_text(size=14), axis.title=element_text(size=14),
    strip.text = element_text(size = 14)) +
  ggsave("PC9_cFPs_simulated_color_DS3_4_withRates_G50.pdf", width = 5, height = 8)

DS3_4_distribution <- read.csv("distributions_G50.csv", row.names = 1)[,c(2,3)]
DS3_4_distribution <- melt(DS3_4_distribution)
names(DS3_4_distribution) <- c("Subline", "DIP")
DS3_4_distribution$Type = "Simulated"

DS3_4_distribution_exp <- read.csv("cFP_rates_VUlines.csv", row.names = 1)
DS3_4_distribution_exp <- subset(DS3_4_distribution_exp, Cell_Line %in% c("PC9-DS3", "PC9-DS4"))
DS3_4_distribution_exp$Subline <- str_split_fixed(DS3_4_distribution_exp$Cell_Line, "-", 2)[,2]
DS3_4_distribution_exp <- DS3_4_distribution_exp[,c("DIP_Rate", "Subline")]
names(DS3_4_distribution_exp) <- c("DIP", "Subline")
DS3_4_distribution_exp$Type = "Experimental"

DS3_4_distribution_data <- rbind(DS3_4_distribution, DS3_4_distribution_exp)

pvalues <- read.csv("pvalues_G50.csv", row.names = 1)
p_DF <- data.frame(Subline = c("DS3", "DS4"),
                   label=c(paste("p =", as.character(round(pvalues["DS3"],3))),
                           paste("p =", as.character(round(pvalues["DS4"],3)))))


ggplot(DS3_4_distribution, aes(x=DIP, fill=Subline, color=Subline)) +
  geom_density(alpha=.25, size = 0.5) +  
  geom_density(data = DS3_4_distribution_exp, aes(x=DIP), 
               fill = "grey90", color = "grey10",
               alpha=0.75, size = 0.5) +
  theme_bw() +
  scale_color_manual(values = c("brown", "deepskyblue")) +
  scale_fill_manual(values = c("brown", "deepskyblue")) +
  geom_vline(xintercept = 0, size = 0.5, colour = "black",
             linetype = "dashed") +
  labs(x = "DIP Rate", y = "Density")+ #xlim(-0.03, 0.05) +
  # facet_grid(.~Sample) +
  facet_wrap(.~Subline, ncol = 1) + xlim(-0.025, 0.020) +
  geom_text(data= p_DF,
            aes(x = -0.015,y = 125, label = label), 
            inherit.aes=FALSE, parse = FALSE, size = 5) +
  theme(
    # axis.line.x = element_line(colour = "black", size = 1),
    # axis.line.y = element_line(colour = "black", size = 1),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    # legend.position = c(0.95, -0.15),
    # legend.justification = c(1, 0),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
    legend.title = element_text(size=14), axis.title=element_text(size=14),
    strip.text = element_text(size = 14)) +
  ggsave("DS3_4_experimental_simulated_NEW.pdf", width = 5, height = 8)

# p <- ggplot(DS3_4_distribution_data, aes(x=DIP, fill=Type, color=Type)) +
#   geom_density(alpha=.25, size = 0.5) + theme_bw() +
#   scale_color_manual(values = c("grey10", "deepskyblue")) +
#   scale_fill_manual(values = c("grey90", "deepskyblue")) +
#   theme(
#     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#     legend.position = "right",
#     # legend.position = c(0.95, -0.15),
#     # legend.justification = c(1, 0),
#     plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
#     legend.title = element_text(size=14), axis.title=element_text(size=14),
#     strip.text = element_text(size = 14))
# leg <- get_legend(p)
# as_ggplot(leg) + ggsave("legend_expSim_DS4.pdf", width = 1, height = 2)


#########################################################################
#########################################################################


DS1_6_7_9_data <- rbind(DS1_trajs, DS6_trajs, DS7_trajs, DS9_trajs,
                        DS1_exp, DS6_exp, DS7_exp, DS9_exp)

cols_DS1_6_7_9 <- c("DS1" = "coral", "DS6" = "deeppink",
                "DS7" = "darkorchid", "DS9" = "gold")
labs_DS1_6_7_9 <- c("DS1", "DS6", "DS7", 
                "DS9")

DS1_6_7_9_sim <- rbind(DS1_trajs, DS6_trajs, DS7_trajs, DS9_trajs)
DS1_6_7_9_exp <- rbind(DS1_exp, DS6_exp, DS7_exp, DS9_exp)

n_DF_others <- data.frame(Subline = c("DS1", "DS6", "DS7", "DS9"),
                   label=c(paste("n =", as.character(length(unique(DS1_exp$variable)))),
                           paste("n =", as.character(length(unique(DS6_exp$variable)))),
                           paste("n =", as.character(length(unique(DS7_exp$variable)))),
                           paste("n =", as.character(length(unique(DS9_exp$variable))))))


ggplot(data = DS1_6_7_9_exp, aes(x=Time, y=nl2, group = variable, color = Subline)) +
  geom_line(size = 0.5, alpha = 0.25) +
  labs(x = "Time (hours) Post Drug Penetrance", y = "Normalized Log2 Cell Count") +
  facet_wrap(.~Subline, ncol = 2) +
  geom_text(data = n_DF_others,
            aes(x = 100, y = 2, label = label),
            inherit.aes=FALSE, parse = FALSE, size = 5) +
  scale_color_manual(values = c("grey", "grey", "grey", "grey")) +
  ylim(-2.5,2.5) +
  theme_bw() +
  theme(
    # axis.line.x = element_line(colour = "black", size = 1),
    # axis.line.y = element_line(colour = "black", size = 1),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
    legend.title = element_text(size=14), axis.title=element_text(size=14),
    strip.text = element_text(size = 14)) +
  ggsave("PC9_cFPs_experimental_grey_DS1_6_7_9.pdf", width = 5, height = 16)

kdiv_DF <- data.frame(Subline = c("DS1", "DS6", "DS7", "DS9"),
                      label=c(paste("k[division] == 0.032"),
                              paste("k[division] == 0.030"),
                              paste("k[division] == 0.028"),
                              paste("k[division] == 0.025"))) 
kdth_DF <- data.frame(Subline = c("DS1", "DS6", "DS7", "DS9"),
                      label=c(paste("k[death] == 0.03050"),
                              paste("k[death] == 0.02940"),
                              paste("k[death] == 0.02625"),
                              paste("k[death] == 0.02375")))

ggplot(data = DS1_6_7_9_sim, aes(x=Time, y=nl2, group = variable, color = Subline)) +
  geom_line(size = 0.5, alpha = 0.25) +
  labs(x = "Time (hours) Post Drug Penetrance", y = "Normalized Log2 Cell Count") +
  geom_text(data= kdiv_DF,
            aes(x = 5,y = 2.3, label = label), 
            inherit.aes=FALSE, parse = TRUE, size = 5,
            hjust = 0) +
  geom_text(data= kdth_DF,
            aes(x = 5,y = 1.8, label = label), 
            inherit.aes=FALSE, parse = TRUE, size = 5,
            hjust = 0) +
  scale_color_manual(values = c("coral", "deeppink", "darkorchid", "gold")) +
  facet_wrap(.~Subline, ncol = 1) +
  ylim(-2.5,2.5) +
  theme_bw() +
  theme(
    # axis.line.x = element_line(colour = "black", size = 1),
    # axis.line.y = element_line(colour = "black", size = 1),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    # legend.position = c(0.95, -0.15),
    # legend.justification = c(1, 0),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
    legend.title = element_text(size=14), axis.title=element_text(size=14),
    strip.text = element_text(size = 14)) +
  ggsave("PC9_cFPs_simulated_color_DS1_6_7_9_withRates_G50.pdf", width = 5, height = 16)


DS1_6_7_9_distribution <- read.csv("distributions_G50.csv", row.names = 1)[,c(1,4,5,6)]
DS1_6_7_9_distribution <- melt(DS1_6_7_9_distribution)
names(DS1_6_7_9_distribution) <- c("Subline", "DIP")
DS1_6_7_9_distribution$Type = "Simulated"

DS1_6_7_9_distribution_exp <- read.csv("cFP_rates_VUlines.csv", row.names = 1)
DS1_6_7_9_distribution_exp <- subset(DS1_6_7_9_distribution_exp, 
                                     Cell_Line %in% c("PC9-DS1", "PC9-DS6",
                                                      "PC9-DS7", "PC9-DS9"))
DS1_6_7_9_distribution_exp$Subline <- str_split_fixed(DS1_6_7_9_distribution_exp$Cell_Line, "-", 2)[,2]
DS1_6_7_9_distribution_exp <- DS1_6_7_9_distribution_exp[,c("DIP_Rate", "Subline")]
names(DS1_6_7_9_distribution_exp) <- c("DIP", "Subline")
DS1_6_7_9_distribution_exp$Type = "Experimental"

DS1_6_7_9_distribution_data <- rbind(DS1_6_7_9_distribution, DS1_6_7_9_distribution_exp)

pvalues <- read.csv("pvalues_G50.csv", row.names = 1)
p_DF <- data.frame(Subline = c("DS1", "DS6", "DS7", "DS9"),
                   label=c(paste("p =", as.character(round(pvalues["DS1"],3))),
                           paste("p =", as.character(round(pvalues["DS6"],3))),
                           paste("p =", as.character(round(pvalues["DS7"],3))),
                           paste("p =", as.character(round(pvalues["DS9"],3)))))


ggplot(DS1_6_7_9_distribution, aes(x=DIP, fill=Subline, color=Subline)) +
  geom_density(alpha=.25, size = 0.5) +  
  geom_density(data = DS1_6_7_9_distribution_exp, aes(x=DIP), 
               fill = "grey90", color = "grey10",
               alpha=0.75, size = 0.5) +
  theme_bw() +
  scale_color_manual(values = c("coral", "deeppink", "darkorchid", "gold")) +
  scale_fill_manual(values = c("coral", "deeppink", "darkorchid", "gold")) +
  geom_vline(xintercept = 0, size = 0.5, colour = "black",
             linetype = "dashed") +
  labs(x = "DIP Rate", y = "Density")+ xlim(-0.025, 0.020) +
  # facet_grid(.~Sample) +
  facet_wrap(.~Subline, ncol = 1) +
  geom_text(data= p_DF,
            aes(x = -0.015,y = 125, label = label), 
            inherit.aes=FALSE, parse = FALSE, size = 5) +
  theme(
    # axis.line.x = element_line(colour = "black", size = 1),
    # axis.line.y = element_line(colour = "black", size = 1),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    # legend.position = c(0.95, -0.15),
    # legend.justification = c(1, 0),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
    legend.title = element_text(size=14), axis.title=element_text(size=14),
    strip.text = element_text(size = 14)) +
  ggsave("DS1_6_7_9_experimental_simulated_NEW.pdf", width = 5, height = 16)

# p <- ggplot(DS1_6_7_9_distribution_data, aes(x=DIP, fill=Type, color=Type)) +
#   geom_density(alpha=.25, size = 0.5) + theme_bw() +
#   scale_color_manual(values = c("grey10", "gold")) +
#   scale_fill_manual(values = c("grey90", "gold")) +
#   theme(
#     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#     legend.position = "right",
#     # legend.position = c(0.95, -0.15),
#     # legend.justification = c(1, 0),
#     plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
#     legend.title = element_text(size=14), axis.title=element_text(size=14),
#     strip.text = element_text(size = 14))
# leg <- get_legend(p)
# as_ggplot(leg) + ggsave("legend_expSim_DS9.pdf", width = 1, height = 2)

##### 
## Recreate distributions.csv
#####
# d1 <- DS3_4_distribution
# d1$Type = NULL
# d2 <- DS1_6_7_9_distribution
# d2$Type = NULL
# 
# d1$id <- c(1:220, 1:220)
# d3 <- dcast(data = d1,formula = id~Subline,fun.aggregate = sum,value.var = "DIP")
# 
# d2$id <- c(1:220, 1:220, 1:220, 1:220)
# d4 <- dcast(data = d2,formula = id~Subline,fun.aggregate = sum,value.var = "DIP")
# 
# d5 <- cbind(d4[,c(1,2)],d3[,c(2,3)], d4[,c(3,4,5)])
# d5$id <- NULL
# 
# write.csv(d5, "distributions.csv")
#################
# DS8
################

adCFPNL2_sub$Sample <- replace(as.character(adCFPNL2_sub$Sample), 
                               adCFPNL2_sub$Sample == "PC9-DS8", "DS8")

DS8_times <- round(unique(subset(adCFPNL2_sub, Sample == "DS8")$Time)) - 
  round(unique(subset(adCFPNL2_sub, Sample == "DS8")$Time))[1]

DS8_exp <- subset(adCFPNL2_sub, Sample == "DS8")
DS8_exp$Time <- rep(DS8_times, rep = length(unique(DS8_exp$Time)))
DS8_exp$variable = as.factor(rep(seq(length(unique(DS8_exp$Well))), each = length(unique(DS8_exp$Time))))
colnames(DS8_exp)[4] <- "Subline"
DS8_exp$Type <- "Experimental"
DS8_exp <- DS8_exp[,c("Time", "variable", "nl2", "Subline", "Type")]

setwd("~/git/ThreeStateMelanoma/")
DS8_trajs <- read.csv('trajectories_DS8_G50.csv', row.names = 1)
DS8_trajs <- DS8_trajs[,colSums(is.na(DS8_trajs))<nrow(DS8_trajs)]
names(DS8_trajs) <- seq(ncol(DS8_trajs))
DS8_trajs$Time <- seq(nrow(DS8_trajs)) -1
DS8_trajs <- melt(DS8_trajs, id = c("Time"), value.name = "nl2")
DS8_trajs$Subline <- "DS8"
DS8_trajs$Type <- "Simulation"
DS8_trajs <- subset(DS8_trajs, Time %in% DS8_times)

n_DF_DS8 <- data.frame(Subline = c("DS8"),
                          label=paste("n =", as.character(length(unique(DS8_exp$variable)))))


ggplot(data = DS8_exp, aes(x=Time, y=nl2, group = variable, color = Subline)) +
  geom_line(size = 0.5, alpha = 0.25) +
  labs(x = "Time (hours) Post Drug Penetrance", y = "Normalized Log2 Cell Count") +
  facet_wrap(.~Subline, ncol = 1) +
  geom_text(data = n_DF_DS8,
            aes(x = 100, y = 2, label = label),
            inherit.aes=FALSE, parse = FALSE, size = 5) +
  scale_color_manual(values = c("grey")) +
  ylim(-2.5,2.5) +
  theme_bw() +
  theme(
    # axis.line.x = element_line(colour = "black", size = 1),
    # axis.line.y = element_line(colour = "black", size = 1),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    # legend.position = c(0.95, -0.15),
    # legend.justification = c(1, 0),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
    legend.title = element_text(size=14), axis.title=element_text(size=14),
    strip.text = element_text(size = 14)) +
  ggsave("PC9_cFPs_experimental_grey_DS8.pdf", width = 5, height = 4)

kdiv1_DF <- data.frame(Subline = c("DS8"),
                      label=c(paste("k[division] == 0.032"))) 
kdiv2_DF <- data.frame(Subline = c("DS8"),
                       label=c(paste("k[division] == 0.033"))) 
kdth1_DF <- data.frame(Subline = c("DS8"),
                       label=c(paste("k[death] == 0.0311"))) 
kdth2_DF <- data.frame(Subline = c("DS8"),
                       label=c(paste("k[death] == 0.0265"))) 


ggplot(data = DS8_trajs, aes(x=Time, y=nl2, group = variable, color = Subline)) +
  geom_line(size = 0.5, alpha = 0.25) +
  labs(x = "Time (hours) Post Drug Penetrance", y = "Normalized Log2 Cell Count") +
  geom_text(data= kdiv1_DF,
            aes(x = 5,y = 2.3, label = label), 
            inherit.aes=FALSE, parse = TRUE, size = 5,
            hjust = 0) +
  geom_text(data= kdth1_DF,
            aes(x = 5,y = 1.8, label = label), 
            inherit.aes=FALSE, parse = TRUE, size = 5,
            hjust = 0) +
  geom_text(data= kdiv2_DF,
            aes(x = 5,y = -1.8, label = label), 
            inherit.aes=FALSE, parse = TRUE, size = 5,
            hjust = 0) +
  geom_text(data= kdth2_DF,
            aes(x = 5,y = -2.3, label = label), 
            inherit.aes=FALSE, parse = TRUE, size = 5,
            hjust = 0) +
  scale_color_manual(values = c("seagreen")) +
  facet_wrap(.~Subline, ncol = 1) +
  ylim(-2.5,2.5) +
  theme_bw() +
  theme(
    # axis.line.x = element_line(colour = "black", size = 1),
    # axis.line.y = element_line(colour = "black", size = 1),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    # legend.position = c(0.95, -0.15),
    # legend.justification = c(1, 0),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
    legend.title = element_text(size=14), axis.title=element_text(size=14),
    strip.text = element_text(size = 14)) +
  ggsave("PC9_cFPs_simulated_color_DS8_withRates_G50.pdf", width = 5, height = 4)


DS8_distribution <- read.csv("distributions_DS8_G50.csv", row.names = 1)[,c(1)]
DS8_distribution <- as.data.frame(DS8_distribution)
DS8_distribution$Subline <- "DS8"
DS8_distribution$Type = "Simulated"
names(DS8_distribution) <- c("DIP", "Subline", "Type")

DS8_distribution_exp <- read.csv("cFP_rates_VUlines.csv", row.names = 1)
DS8_distribution_exp <- subset(DS8_distribution_exp, 
                                     Cell_Line %in% c("PC9-DS8"))
DS8_distribution_exp$Subline <- str_split_fixed(DS8_distribution_exp$Cell_Line, "-", 2)[,2]
DS8_distribution_exp <- DS8_distribution_exp[,c("DIP_Rate", "Subline")]
names(DS8_distribution_exp) <- c("DIP", "Subline")
DS8_distribution_exp$Type = "Experimental"

DS8_distribution_data <- rbind(DS8_distribution, DS8_distribution_exp)

# pvalues <- read.csv("pvalues.csv", row.names = 1)
p_DF <- data.frame(Subline = c("DS8"),
                   label=c(paste("p =", as.character(round(0.3098197140751014,3)))))


ggplot(DS8_distribution, aes(x=DIP, fill=Subline, color=Subline)) +
  geom_density(alpha=.25, size = 0.5) +  
  geom_density(data = DS8_distribution_exp, aes(x=DIP), 
               fill = "grey90", color = "grey10",
               alpha=0.75, size = 0.5) +
  theme_bw() +
  scale_color_manual(values = c("seagreen")) +
  scale_fill_manual(values = c("seagreen")) +
  geom_vline(xintercept = 0, size = 0.5, colour = "black",
             linetype = "dashed") +
  labs(x = "DIP Rate", y = "Density")+ xlim(-0.025, 0.020) + ylim(0,210) +
  # facet_grid(.~Sample) +
  facet_wrap(.~Subline, ncol = 1) +
  geom_text(data= p_DF,
            aes(x = -0.015,y = 125, label = label), 
            inherit.aes=FALSE, parse = FALSE, size = 5) +
  theme(
    # axis.line.x = element_line(colour = "black", size = 1),
    # axis.line.y = element_line(colour = "black", size = 1),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "none",
    # legend.position = c(0.95, -0.15),
    # legend.justification = c(1, 0),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
    legend.title = element_text(size=14), axis.title=element_text(size=14),
    strip.text = element_text(size = 14)) +
  ggsave("DS8_experimental_simulated_NEW.pdf", width = 5, height = 4)

# pval = 0.3098197140751014

p <- ggplot(DS8_distribution_data, aes(x=DIP, fill=Type, color=Type)) +
  geom_density(alpha=.25, size = 0.5) + theme_bw() +
  scale_color_manual(values = c("grey10", "seagreen")) +
  scale_fill_manual(values = c("grey90", "seagreen")) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.position = "right",
    # legend.position = c(0.95, -0.15),
    # legend.justification = c(1, 0),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"), axis.text=element_text(size=14),
    legend.title = element_text(size=14), axis.title=element_text(size=14),
    strip.text = element_text(size = 14))
leg <- get_legend(p)
as_ggplot(leg) + ggsave("legend_expSim_DS8.pdf", width = 1, height = 2)


#### Median test all sublines
DSLines_DIP <- read.csv("cFP_rates_VUlines.csv", row.names = 1)
DSLines_DIP <- subset(DSLines_DIP, Cell_Line %in% c("PC9-DS1", "PC9-DS3", 
                                                 "PC9-DS4", "PC9-DS6", 
                                                 "PC9-DS7", "PC9-DS8", 
                                                 "PC9-DS9"))
library(RVAideMemoire)

mood.medtest(DIP_Rate ~ Cell_Line,
             data  = DSLines_DIP,
             exact = FALSE)
summary(aov(DIP_Rate ~ Cell_Line, data = DSLines_DIP))