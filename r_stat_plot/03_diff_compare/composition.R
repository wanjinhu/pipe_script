library(tidyverse)
library(ggsci)
library(ggplot2)
setwd("F:/03-Dipro项目/11-张瑞-宏基因组-犬类/analysis/03_diff_compare/")
df <- read.table("../06_metaphlan_species.txt",
                 header = T,
                 sep = "\t",
                 check.names = F)
split_data <- strsplit(df$clade_name, ";")
split_data <- do.call("rbind",split_data)
cols <- df[, c(2:ncol(df))]
df_merge <- do.call("cbind", list(split_data,cols))
names(df_merge)[1:7] <- c("kindom","phylum","class","order","family","genus","species")
group <- read.table("../group.txt",
                    header = T,
                    sep = "\t",
                    check.names = F)
spe <- df_merge[-c(1:7)]
tax <- df_merge[1:7]

# 2.1.1 物种组成数据按照门进行汇总
## spe和tax数据表中物种排序一致
phy <- spe %>%
  group_by(tax$phylum) %>% # 使用tax中的门水平进行分类
  summarise_all(sum) %>%
  rename(phylum = `tax$phylum`) %>%
  gather(key="samples",value = "abun",-phylum) %>% # 数据形式转换：“宽”转“长”
  left_join(group,by=c("samples"="sample"))

# 2.1.2 颜色
col=pal_d3("category20")(20)
col2 = pal_d3("category20",alpha = 0.5)(20)
mypal=c(col,col2)

# 2.1.3  物种组成堆叠柱形图[这里的数据已经是百分比值了]
p1 <-
  ggplot(phy,aes(samples,abun,fill=phylum)) + 
  geom_bar(stat = "identity",position = "fill") +
  facet_grid(. ~ group1, scales = "free", space = "free_x") +
  labs(x="",y="Percent")+
  scale_fill_manual(values = mypal)+labs(fill="")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(expand=c(0,0))
pdf("phylum_percent_2groups.pdf",width = 12,height = 10,family="Times")
p1
dev.off()

p2 <-
  ggplot(phy,aes(samples,abun,fill=phylum)) + 
  geom_bar(stat = "identity",position = "fill") +
  facet_grid(. ~ group2, scales = "free", space = "free_x") +
  labs(x="",y="Percent")+
  scale_fill_manual(values = mypal)+labs(fill="")+
  theme_bw()+
  theme(legend.title=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(expand=c(0,0))
pdf("phylum_percent_4groups.pdf",width = 12,height = 10,family="Times")
p2
dev.off()
