# #####################################################
# ################# 函数的demo测试 ####################
# #####################################################
## -------------------------------------------------------------------------------
## 所有的demo测试数据都在此文件夹中
setwd("F:/Wanjin-Project/pipe_script/r_stat_plot/demo_data/")
## 如果函数保留导出结果功能，则保存在自己设定的目录中,否则保存在该目录中
demo_result <- "F:/Wanjin-Project/pipe_script/r_stat_plot/demo_data/demo_result/"
source("F:/Wanjin-Project/pipe_script/r_stat_plot/r_stat_plot.R")

## -------------------------------------------------------------------------------
## 1. alpha多样性计算
otu <- read.delim("06_metaphlan_species.txt", 
                  row.names = 1, 
                  sep = '\t', 
                  stringsAsFactors = FALSE, 
                  check.names = FALSE)
otu <- t(otu)
alpha_all <- stat_alpha_1(otu, base = 2)
write.csv(alpha_all, paste0(demo_result,"alpha.csv"), quote = FALSE)
## -------------------------------------------------------------------------------
## 2. beta_pcoa绘图
data <- read.delim("06_metaphlan_species.txt",
                   sep = "\t",
                   header = T,
                   check.names = FALSE,
                   row.names = 1)
# 读取分组文件
group <- read.delim("group.txt",
                    header = T,
                    sep = "\t",
                    check.names = F)
# 选择哪个分组方案[sample group],分析前做好选择,可以是其他方法
group_target <- data.frame(group[,c(1:2)])
plot_pcoa_1(data = data, sd = group_target, prefix = "test")
## -------------------------------------------------------------------------------
## 3. 柱状图[显著性标注,针对1里面计算得到的alpha数据表]
data <- read.csv("demo_result/alpha.csv", header = T, check.names = F)
group <- read.table("group.txt", sep="\t", header = T)
# 选择哪个分组方案[sample group],分析前做好选择,可以是其他方法
group_target <- data.frame(group[,c(1,3)])
box_plot <- plot_bar_1(data = data, 
                       group_file = group_target,
                       target = "Shannon",
                       ymax = 8)
ggsave(paste0(demo_result,"Shannon_boxplot.pdf"),
       box_plot,width = 12,height = 10)
## -------------------------------------------------------------------------------
