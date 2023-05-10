library(tidyverse)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(reshape2)

setwd("F:/03-Dipro项目/11-张瑞-宏基因组-犬类/analysis/04_diff_KEGG/")
pathway <- read.delim("../pathway_samples.xls", 
                      header=T, 
                      check.names = F,
                      sep = "\t")
group <- read.table("../group1.txt", 
                    header=T,
                    check.names = F,
                    sep = "\t")
name <- c()
t_result <- c()
# wilcox_result <- c()
for (i in 1:nrow(pathway)){
  target <- pathway[i,c(4:ncol(pathway))]
  target_name <- as.character(target["pathway"][1])
  test <- as.data.frame(t(target))
  test$sample <- rownames(test)
  test <- as.data.frame(test[-1,])
  names(test) <- c(target_name,"sample")
  df_merge <- merge(test,group,by="sample")
  df_merge[,2] <- as.numeric(df_merge[,2])
  df_merge[,3] <- as.factor(df_merge[,3])
  t_stat <- round((t.test(df_merge[,2]~df_merge[,3]))$p.value,4)
  # wilcox_stat <- round((wilcox.test(df_merge[,2]~df_merge[,3]))$p.value,4)
  name <- c(name,target_name)
  t_result <- c(t_result, t_stat)
  # wilcox_result <- c(wilcox_result,wilcox_stat)
}

result <- data.frame(id=name,  
                     t_value=t_result,
                     # wilcox_value=wilcox_result,
                     t_determine=ifelse(t_result <= 0.05, "Diff", "Non-diff"))

# 筛选t检验显著差异的id以及数据表
diff <- result[result$t_determine=="Diff", "id"]
diff_df <- pathway %>% filter(pathway %in% diff)
diff_sig <- merge(result,diff_df,by.x = "id",by.y = "pathway")
names(diff_sig)[1] <- "pathway"
write.table(diff_sig,
            file = "diff_pathway_sig.xls",
            sep = "\t",
            col.names = NA,
            quote = F)
data <- diff_df %>% melt() 
names(data)[5] <- "sample"
ndf <- left_join(data,group,by="sample")
ndf$value_log10 <- ifelse(ndf$value == 0, ndf$value, round(log10(ndf$value),4))
ndf$value_log10[ndf$group=="pet_dog"] <- ndf$value_log10[ndf$group=="pet_dog"]*-1
order_data <- data.frame(Level1=ndf$level1,
                    Pathway=ndf$pathway)
order_data <- order_data %>% filter(!duplicated(Pathway))
order_data <- order_data[order(order_data$Level1),]
ndf$pathway <- factor(ndf$pathway,levels=rev(order_data$Pathway))
legend_order <- order_data %>% filter(!duplicated(Level1))
ndf$level1 <- factor(ndf$level1, levels = legend_order$Level1)

p <- 
  ggplot(ndf)+
  geom_boxplot(aes(x=pathway,y=value_log10,fill=level1))+
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

ggsave("diff_pathway.pdf",p,width = 15,height = 13)
