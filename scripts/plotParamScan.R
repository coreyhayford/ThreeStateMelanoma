library(ggplot2)
setwd('~/git/ThreeStateMelanoma/')
tile <- read.csv('all_cellLine_tile.csv')
tile_new <- subset(tile, !cell.line == "PC9.DS8") # DS8 cFP did not match model prediction
cols <- c("not.assigned" = "grey90", "PC9.DS1" = "coral", "PC9.DS3" = "brown",
          "PC9.DS4" = "deepskyblue", "PC9.DS6" = "deeppink", "PC9.DS7" = "darkorchid",
          "PC9.DS9" = "gold")
pops <- c("not.assigned", "PC9.DS1", "PC9.DS3", "PC9.DS4", 
          "PC9.DS6","PC9.DS7", "PC9.DS9")
labs <- c("Not Assigned", "DS1", "DS3", "DS4", "DS6", "DS7", "DS9")
ggplot(tile_new) +
  geom_tile(data = subset(tile_new, cell.line.new == 'PC9.DS1'), 
            aes(DIP.rate, division.rate, fill = cell.line.new, alpha = p.value), color = "coral") +
  geom_tile(data = subset(tile_new, cell.line.new == 'PC9.DS3'), 
            aes(DIP.rate, division.rate, fill = cell.line.new, alpha = p.value), color = "brown") +
  geom_tile(data = subset(tile_new, cell.line.new == 'PC9.DS4'), 
            aes(DIP.rate, division.rate, fill = cell.line.new, alpha = p.value), color = "deepskyblue") +
  geom_tile(data = subset(tile_new, cell.line.new == 'PC9.DS6'), 
            aes(DIP.rate, division.rate, fill = cell.line.new, alpha = p.value), color = "deeppink") +
  geom_tile(data = subset(tile_new, cell.line.new == 'PC9.DS7'), 
            aes(DIP.rate, division.rate, fill = cell.line.new, alpha = p.value), color = "darkorchid") +
  geom_tile(data = subset(tile_new, cell.line.new == 'PC9.DS9'), 
            aes(DIP.rate, division.rate, fill = cell.line.new, alpha = p.value), color = "gold") +
  theme_bw() + 
  scale_fill_manual(values = cols,
                    breaks = pops,
                    labels = labs,
                    name = "Subline") +
  scale_alpha_continuous(name = "p-value") +
  xlab("DIP Rate") + ylab("Division Rate") +
  theme(legend.text = element_text(size = 12), legend.position = "right", 
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), axis.text=element_text(size=12),
        legend.title = element_text(size=12), axis.title=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # ggtitle("Subline Parameter Sweep") +
  ggsave("parameterScan_Tile.pdf", width = 10, height = 5)

tile_forMainFig <- subset(tile, cell.line == c("PC9.DS3", "PC9.DS4")) 
cols <- c("PC9.DS3" = "brown", "PC9.DS4" = "deepskyblue")
pops <- c("PC9.DS3", "PC9.DS4")
labs <- c("DS3", "DS4")
ggplot(tile_forMainFig) +
  geom_tile(data = subset(tile_new, cell.line.new == 'PC9.DS3'),
            aes(DIP.rate, division.rate, fill = cell.line.new, alpha = p.value),
            color = "brown") + #, fill = "brown") +
  geom_tile(data = subset(tile_new, cell.line.new == 'PC9.DS4'),
            aes(DIP.rate, division.rate, fill = cell.line.new, alpha = p.value), 
            color = "deepskyblue") + #, fill = "deepskyblue") +
  theme_bw() + 
  scale_fill_manual(values = cols,
                    breaks = pops,
                    labels = labs,
                    name = "Subline") +
  scale_alpha_continuous(name = "p-value", range = c(0,1), limits = c(0,1)) +
  xlab("Divison Rate - Death Rate") + ylab("Division Rate") +
  theme(legend.text = element_text(size = 12), legend.position = "right", 
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), axis.text=element_text(size=12),
        legend.title = element_text(size=12), axis.title=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # ggtitle("Subline Parameter Sweep") +
  ggsave("parameterScan_Tile_DS3_4_forMainFigs.pdf", width = 10, height = 5)

tile_forSuppFig <- subset(tile, cell.line %in% c("PC9.DS1", "PC9.DS6", "PC9.DS7", "PC9.DS9")) 
cols <- c("PC9.DS1" = "coral", "PC9.DS6" = "deeppink", 
          "PC9.DS7" = "darkorchid", "PC9.DS9" = "gold")
pops <- c("PC9.DS1", "PC9.DS6", "PC9.DS7", "PC9.DS9")
labs <- c("DS1", "DS6", "DS7", "DS9")
ggplot(tile_forSuppFig) +
  geom_tile(data = subset(tile_new, cell.line.new == 'PC9.DS1'),
            aes(DIP.rate, division.rate, fill = cell.line.new, alpha = p.value), color = "coral") + #, fill = "coral") +
  geom_tile(data = subset(tile_new, cell.line.new == 'PC9.DS6'),
            aes(DIP.rate, division.rate, fill = cell.line.new, alpha = p.value), color = "deeppink") + #, fill = "deeppink") +
  geom_tile(data = subset(tile_new, cell.line.new == 'PC9.DS7'),
            aes(DIP.rate, division.rate, fill = cell.line.new, alpha = p.value), color = "darkorchid") + #, fill = "darkorchid") +
  geom_tile(data = subset(tile_new, cell.line.new == 'PC9.DS9'),
            aes(DIP.rate, division.rate, fill = cell.line.new, alpha = p.value), color = "gold") + #, fill = "gold") +
  theme_bw() +
  scale_fill_manual(values = cols,
                    breaks = pops,
                    labels = labs,
                    name = "Subline") +
  scale_alpha_continuous(name = "p-value", range = c(0,1), limits = c(0,1)) +
  xlim(-0.0002, 0.0022) +
  xlab("Division Rate - Death Rate") + ylab("Division Rate") +
  theme(legend.text = element_text(size = 12), legend.position = "right", 
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), axis.text=element_text(size=12),
        legend.title = element_text(size=12), axis.title=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # ggtitle("Subline Parameter Sweep") +
  ggsave("parameterScan_Tile_DS1-6-7-9_forSuppFigs.pdf", width = 7.5, height = 4.5)
  # ggsave("parameterScan_Tile_DS9_colorAlpha.pdf", width = 7.5, height = 4.5)


#### DS8 two state

# tile_DS8 <- read.csv('DS8_twoState_tile_expandRange.csv')
# tile_DS8_largerDiv <- read.csv('DS8_twoState_tile_largerDiv.csv')
# tile_DS8_largerDIP <- read.csv('DS8_twoState_tile_largerDIPs.csv')
# tile_DS8_wideRange <- read.csv('DS8_twoState_tile_wideRangeTest.csv')

# tile_DS8_wholeRange_lowDivs <- read.csv('DS8_twoState_lowDivs.csv')
# tile_DS8_wholeRange_lowMedDivs <- read.csv('DS8_twoState_lowMedDivs.csv')
# tile_DS8_wholeRange_medHighDivs <- read.csv('DS8_twoState_medHighDivs.csv')
# tile_DS8_wholeRange_highDivs <- read.csv('DS8_twoState_highDivs.csv')
# tile_DS8_wholeRange <- rbind(tile_DS8_wholeRange_lowDivs,
#                              tile_DS8_wholeRange_lowMedDivs,
#                              tile_DS8_wholeRange_medHighDivs,
#                              tile_DS8_wholeRange_highDivs)
# write.csv(tile_DS8_wholeRange, "DS8_wholeRange.csv")
tile_DS8_wholeRange <- read.csv('DS8_twoState_tile_wholeRange.csv')



# DS8_state1 <- subset(tile_DS8, Cell.Line == 'PC9-DS8.1')
# DS8_state2 <- subset(tile_DS8, Cell.Line == 'PC9-DS8.2')
# DS8_largerDiv_state1 <- subset(tile_DS8_largerDiv, Cell.Line == 'PC9-DS8.1')
# DS8_largerDiv_state2 <- subset(tile_DS8_largerDiv, Cell.Line == 'PC9-DS8.2')
# DS8_largerDIP_state1 <- subset(tile_DS8_largerDIP, Cell.Line == 'PC9-DS8.1')
# DS8_largerDIP_state2 <- subset(tile_DS8_largerDIP, Cell.Line == 'PC9-DS8.2')
# DS8_WRT_state1 <- subset(tile_DS8_wideRange, Cell.Line == 'PC9-DS8.1')
# DS8_WRT_state2 <- subset(tile_DS8_wideRange, Cell.Line == 'PC9-DS8.2')
DS8_WR_state1 <- subset(tile_DS8_wholeRange, Cell.Line == 'PC9-DS8.1')
DS8_WR_state2 <- subset(tile_DS8_wholeRange, Cell.Line == 'PC9-DS8.2')

# DS8_state1$ParamPair <- as.character(seq(1:nrow(DS8_state1)))
# DS8_state2$ParamPair <- as.character(seq(1:nrow(DS8_state2)))
# DS8_largerDiv_state1$ParamPair <- as.character(seq(1:nrow(DS8_largerDiv_state1)))
# DS8_largerDiv_state2$ParamPair <- as.character(seq(1:nrow(DS8_largerDiv_state2)))
# DS8_largerDIP_state1$ParamPair <- as.character(seq(1:nrow(DS8_largerDIP_state1)))
# DS8_largerDIP_state2$ParamPair <- as.character(seq(1:nrow(DS8_largerDIP_state2)))
# DS8_WRT_state1$ParamPair <- as.character(seq(1:nrow(DS8_WRT_state1)))
# DS8_WRT_state2$ParamPair <- as.character(seq(1:nrow(DS8_WRT_state2)))
DS8_WR_state1$ParamPair <- as.character(seq(1:nrow(DS8_WR_state1)))
DS8_WR_state2$ParamPair <- as.character(seq(1:nrow(DS8_WR_state2)))

# DS8 <- rbind(DS8_state1, DS8_state2, DS8_largerDiv_state1, DS8_largerDiv_state2)
# DS8_WRT <- rbind(DS8_WRT_state1, DS8_WRT_state2)
DS8_WR <- rbind(DS8_WR_state1, DS8_WR_state2)


DS8_WR_df <- data.frame(
  p.value = DS8_WR$p.value,
  Cell.Line = DS8_WR$Cell.Line,
  DIP.Rate = jitter(DS8_WR$DIP.Rate,2),
  Division.Rate = jitter(DS8_WR$Division.Rate,2),
  Death.Rate = jitter(DS8_WR$Death.Rate,2),
  ParamPair = as.numeric(DS8_WR$ParamPair)
)

# DS8_1 <- rbind(DS8_largerDiv_state1, DS8_state1)
# DS8_2 <- rbind(DS8_largerDiv_state2, DS8_state2)
# 
# DS8_largerDIP <- rbind(DS8_largerDIP_state1, DS8_largerDIP_state2)


cols_DS8 <- c("PC9-DS8.1" = "seagreen", "PC9-DS8.2" = "seagreen")
pops_DS8 <- c("PC9-DS8.1", "PC9-DS8.2")
labs_DS8 <- c("DS8 State 1", "DS8 State 2")

plt <- ggplot(DS8_WR_df, aes(DIP.Rate, Division.Rate, color = Cell.Line, 
                             alpha = p.value, fill = Cell.Line)) +#, fill = ParamPair)) +
  geom_point(aes(alpha = p.value), color = "seagreen", fill = "seagreen") +
  # geom_density_2d() +
  # geom_tile(aes(width = 0.0002, height = 0.001), size = 1.1) +
  theme_bw()  + ylim(0.01, 0.045) + xlim(0, 0.008) + 
  scale_color_manual(values = cols_DS8,
                    breaks = pops_DS8,
                    labels = labs_DS8,
                    name = "Subline") +
  # scale_fill_gradientn(colors = rainbow(nrow(DS8_df)/2)) +
  scale_fill_manual(values = cols_DS8,
                    breaks = pops_DS8,
                    labels = labs_DS8,
                    name = "Subline") +
  scale_alpha_continuous(name = "p-value", range = c(0,1), limits = c(0,1)) +
  xlab("Division Rate - Death Rate") + ylab("Division Rate") +  #facet_grid(.~Cell.Line, scales = "free_x") +#facet_wrap(~ParamPair, ncol = 5) + 
  theme(legend.text = element_text(size = 12), legend.position = "none", 
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), axis.text=element_text(size=12),
        legend.title = element_text(size=12), axis.title=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.spacing = unit(2, "lines"))

  # ggtitle("DS8 Parameter Sweep") +
plt + ggsave("DS8_twoState_wholeRange_colorCellLine_pointJitterAlphaPvalue_legendRight.pdf", width = 5, height = 3)
leg_plt <- get_legend(plt)
as_ggplot(leg_plt) + ggsave("DS8_PS_legend.svg", width = 3, height = 6)






ggplot(tile_DS8) +
  # geom_point(data = DS8_state1, position = position_jitterdodge(jitter.width = 0.002, jitter.height = 0.002),
  #           aes(DIP.Rate, Division.Rate, color = Cell.Line)) +
  # geom_point(data = DS8_state2, position = position_jitterdodge(jitter.width = 0.002, jitter.height = 0.002),
  #           aes(DIP.Rate, Division.Rate, color = Cell.Line)) +
  # # geom_density_2d(data = subset(tile_DS8, Cell.Line == 'PC9-DS8.1'), aes(DIP.Rate, Division.Rate, color = Cell.Line)) +
  # # geom_density_2d(data = subset(tile_DS8, Cell.Line == 'PC9-DS8.2'), aes(DIP.Rate, Division.Rate, colour = Cell.Line)) +
  geom_tile(data = subset(tile_DS8, Cell.Line == 'PC9-DS8.1'), 
            aes(DIP.Rate, Division.Rate, fill = Cell.Line, alpha = p.value, width = 0.001, height = 0.001)) +
  geom_tile(data = subset(tile_DS8, Cell.Line == 'PC9-DS8.2'), 
            aes(DIP.Rate, Division.Rate, fill = Cell.Line, alpha = p.value, width = 0.001, height = 0.001)) +
  # stat_density_2d(data = DS8_state1, aes(DIP.Rate, Division.Rate, fill = ..density..), geom = "raster", contour = FALSE) +
  # stat_density_2d(data = DS8_state2, aes(DIP.Rate, Division.Rate, fill = ..density..), geom = "raster", contour = FALSE) +
  # 
  theme_bw() + xlim(0, 0.008) + ylim(0.01, 0.045) + 
  scale_fill_manual(values = cols_DS8,
                    breaks = pops_DS8,
                    labels = labs_DS8,
                    name = "Subline") +
  scale_alpha_continuous(name = "p-value", range = c(0,1), limits = c(0,1)) +
  # scale_color_manual(values = cols_DS8,
  #                   breaks = pops_DS8,
  #                   labels = labs_DS8,
  #                   name = "Subline") +
  scale_alpha_continuous(name = "p-value") +
  xlab("DIP Rate") + ylab("Division Rate") +
  theme(legend.text = element_text(size = 12), legend.position = "right", 
        plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), axis.text=element_text(size=12),
        legend.title = element_text(size=12), axis.title=element_text(size=12),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("DS8 Parameter Sweep") 

+
  ggsave("parameterScan_Tile_DS8_redBlue.pdf", width = 10, height = 4)
