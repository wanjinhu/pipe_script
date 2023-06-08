#!/usr/bin/Rscript
# -*- encoding: utf-8 -*-

#@File    :   alpha.R
#@Time    :   2023/06/08 17:45:22
#@Author  :   Wanjin.Hu 
#@Version :   1.0
#@Contact :   wanjin.hu@diprobio.com
#@Description :

stat_alpha <- function(x, tree = NULL, base = 2) {
  "
  @功能: alpha多样性计算
  @参数x: 物种组成表
  @参数tree: 物种序列进化树[如果有的话]
  @return: alpha结果
  "
  Shannon <- diversity(x, index = 'shannon', base = base)
  Simpson <- diversity(x, index = 'simpson')
  result <- data.frame(Shannon, Simpson)
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
  }
  return(result)
}

plot_bar <- function(data,group_file,target) {
  "
  @功能: 柱状图
  @参数data: factor组成表,第一列为样本,后面几列是factor数据
  @参数group_file: 分组文件,只支持1个分组方式,即第1列为样本,第2列为分组方案
  @参数target: 对哪个factor作图,即参数data的列名
  @return: 柱状图结果
  "
  names(data)[1] <- c("sample")
  names(group_file) <- c("sample","group")
  dataCom<-merge(data,group_file,by.x='sample', by.y='sample')
  group_info <- unique(dataCom$group)
  dataCom$group <- factor(dataCom$group, 
                          levels = group_info)
  # 生成两两比较的list
  num <- length(group_info)
  compare_list <- vector("list", 0)
  n <- 0
  for (i in 1:(num-1)) {
    for (j in (i+1):num) {
      n <- n + 1
      compare_list[[n]] <- c(group_info[[i]],
                             group_info[[j]])
    }
  }
  colDef=c("#dc4c43","#4b8ac0","#e6811d","#228a58","#8e48d8","#6d6d6d")
  p1 <- 
    ggplot(dataCom, aes(as.factor(group), dataCom[[target]], fill=group)) + 
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
    labs(y = target, col = "group") + 
    stat_compare_means(method = "wilcox.test", 
                       comparisons = compare_list,
                       label = "p.signif")
  return(p1)
}

library(getopt)
command <- matrix(c(
  "help", "h", 0, "logical", "Help docs",
  "data", "d", 1, "character", "data file",
  "group", "g", 1, "character", "group file"),
  byrow = TRUE, ncol = 5
)
args <- getopt(command)
if (! is.null(args$help) || 
    is.null(args$data) || 
    is.null(args$group)) {
  cat(paste(getopt(command, usage = TRUE), "\n"))
  q(status = 1)
}

library(vegan)
library(picante)
library(ggplot2)
library(tidyr)
library(ggpubr)
otu <- read.delim(args$data, 
                  header = T,
                  row.names = 1, 
                  sep = '\t', 
                  stringsAsFactors = FALSE, 
                  check.names = FALSE)
group <- read.delim(args$group, 
                    sep="\t", 
                    header = T)
otu <- t(otu)
alpha_all <- stat_alpha(otu, base = 2)
write.csv(alpha_all,"alpha.csv", quote = FALSE)
message("Finish alpha stat!")
alpha_all <- read.csv("alpha.csv", header = T, check.names = F)
simpson_plot <- plot_bar(data = alpha_all, 
                         group_file = group,
                         target = "Simpson")
ggsave("Simpson_boxplot.pdf",simpson_plot,width = 12,height = 10)
shannon_plot <- plot_bar(data = alpha_all, 
                         group_file = group,
                         target = "Shannon")
ggsave("Shannon_boxplot.pdf",shannon_plot,width = 12,height = 10)
message("Finish alpha plot!")

## 本地测试
# setwd("F:/Wanjin-Project/pipe_script/general/data/")
# otu <- read.delim("06_metaphlan_species.txt", 
#                   header = T,
#                   row.names = 1, 
#                   sep = '\t', 
#                   stringsAsFactors = FALSE, 
#                   check.names = FALSE)
# group <- read.delim("group.txt", sep="\t", header = T)
# names(group) <- c("sample","group")
# otu <- t(otu)
# alpha_all <- stat_alpha(otu, base = 2)
# write.csv(alpha_all, paste0("alpha.csv"), quote = FALSE)
# alpha_all <- read.csv("alpha.csv", header = T, check.names = F)
# simpson_plot <- plot_bar(data = alpha_all, 
#                          group_file = group,
#                          target = "Simpson")
# ggsave("Simpson_boxplot.pdf",simpson_plot,width = 12,height = 10)
# shannon_plot <- plot_bar(data = alpha_all, 
#                          group_file = group,
#                          target = "Shannon")
# ggsave("Shannon_boxplot.pdf",shannon_plot,width = 12,height = 10)
