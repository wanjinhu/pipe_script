library(ggplot2)
library(tidyr)
library(ggpubr)
library(dplyr)
setwd("F:/03-Dipro项目/11-张瑞-宏基因组-犬类/analysis/03_diff_compare/")
colDef <- c("#dc4c43","#4b8ac0","#e6811d","#228a58","#8e48d8","#6d6d6d")
# 差异细菌双向箱形图
data <- read.delim("diff_sp_g_s.xls", 
                   sep="\t", 
                   header = T, 
                   check.names = F)
names(data)[1:2] <- c("Genus","Species")
group <- read.table("../group1.txt", 
                    sep="\t", 
                    header = T)
names(group)[1] <- "sampleID"
data <- data %>% melt() 
colnames(data)[3] <- "sampleID"
ndf <- left_join(data,group,by="sampleID")
ndf$count <- ifelse(ndf$value == 0, ndf$value, ndf$value)
# ndf$count <- log2(ndf$value)

ndf$count[ndf$group=="pet_dog"] <- ndf$count[ndf$group=="pet_dog"]*-1
order <- data.frame(Genus=ndf$Genus,Species=ndf$Species)
order <- order %>% filter(!duplicated(Species))
ndf$Species <- factor(ndf$Species,levels=rev(order$Species))
legend_order <- order %>% filter(!duplicated(Genus))
ndf$Genus <- factor(ndf$Genus, levels = legend_order$Genus)

p <- 
  ggplot(ndf)+
  geom_boxplot(aes(x=Species,y=count,fill=Genus))+
  facet_wrap(~group,scale="free_x")+
  coord_flip()+
  theme_bw()+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10, face = "plain"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, face = "plain"),
        #panel.border = element_blank(), axis.line = element_line(colour = "black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.placement = "outside")+
  guides(fill = guide_legend( ncol = 1, byrow = TRUE))

ggsave("diff_species.pdf",p)
