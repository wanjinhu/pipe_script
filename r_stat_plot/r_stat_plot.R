## 目的: 整理早就应该整理的数据处理-绘图的一些代码
## 方式: 定义具体的函数,处理具体问题,函数命名规则如下:
# 主要功能为数据处理: stat_xx_1
# 主要功能为绘图: plot_xx_1

#定义函数
library(picante)       #picante 包加载时默认同时加载 vegan

alpha <- function(x, tree = NULL, base = exp(1)) {
  Shannon <- diversity(x, index = 'shannon', base = base)
  Simpson <- diversity(x, index = 'simpson')
  result <- data.frame(Shannon, Simpson)
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
  }
  result
}

#现在直接使用定义好的命令 alpha()，一步得到多种 Alpha 多样性指数
#加载 OTU 丰度表和进化树文件
setwd("F:/03-Dipro项目/11-张瑞-宏基因组-犬类/analysis/00_alpha/")
otu <- read.delim('../06_metaphlan_species.txt', 
                  row.names = 1, 
                  sep = '\t', 
                  stringsAsFactors = FALSE, 
                  check.names = FALSE)
otu <- t(otu)
tree <- read.tree('otu_tree.tre')

#不包含谱系多样性，无需指定进化树；Shannon 公式的 log 底数我们使用 2
alpha_all <- alpha(otu, base = 2)
#包含谱系多样性时，指定进化树文件；Shannon 公式的 log 底数我们使用 2
alpha_all <- alpha(otu, tree, base = 2)

#输出保存在本地
write.csv(alpha_all, 'alpha.csv', quote = FALSE)

# --------------------------
# alpha_vilon作图
# --------------------------
library(ggplot2)
library(tidyr)
library(ggpubr)
setwd("F:/03-Dipro项目/11-张瑞-宏基因组-犬类/analysis/00_alpha/")
data <- read.csv("alpha.csv", header = T, check.names = F)
group <- read.table("../group2.txt", sep="\t", header = T)
dataCom<-merge(data,group,by.x='sample', by.y='sample')
colDef=c("#dc4c43","#4b8ac0","#e6811d","#228a58","#8e48d8","#6d6d6d")

# shannon指数
p1 <- 
  ggplot(dataCom, aes(as.factor(group), Shannon, fill=group)) + 
  # geom_violin(trim = FALSE, color="white") + 
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_boxplot(aes(fill = group), width = 0.2, show.legend = F) + 
  scale_fill_manual(values=colDef) + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, face = "plain"),
        axis.title.y = element_text(size = 16, face = "plain"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, face = "plain"),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(colour = "black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  labs(y = "Shannon", col = "group") + 
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c('Chinese Pastoral Dog','German shepherds'),
                                        c('Chinese Pastoral Dog','Golden Retriever'),
                                        c('Chinese Pastoral Dog','Kunming Dog'),
                                        c('German shepherds','Golden Retriever'),
                                        c('German shepherds','Kunming Dog'),
                                        c('Golden Retriever','Kunming Dog')),
                     label = "p.signif") +
  stat_compare_means(label.y = 8)
# Simpson指数
p2 <- 
  ggplot(dataCom, aes(as.factor(group), Simpson, fill=group)) + 
  # geom_violin(trim = FALSE, color="white") + 
  stat_boxplot(geom = "errorbar", width = 0.1) +
  geom_boxplot(aes(fill = group), width = 0.2, show.legend = F) + 
  scale_fill_manual(values=colDef) + 
  theme_bw() +
  theme(axis.text.y = element_text(size = 12, face = "plain"),
        axis.title.y = element_text(size = 16, face = "plain"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, face = "plain"),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(), axis.line = element_line(colour = "black", size = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  labs(y = "Simpson", col = "group") + 
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c('Chinese Pastoral Dog','German shepherds'),
                                        c('Chinese Pastoral Dog','Golden Retriever'),
                                        c('Chinese Pastoral Dog','Kunming Dog'),
                                        c('German shepherds','Golden Retriever'),
                                        c('German shepherds','Kunming Dog'),
                                        c('Golden Retriever','Kunming Dog')),
                     label = "p.signif") +
  stat_compare_means(label.y = 1.4)

p <- ggarrange(p1, p2, ncol = 2, nrow = 1)
pdf("alpha_species.pdf",width = 20)
p
dev.off()
