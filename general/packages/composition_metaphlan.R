#!/usr/bin/Rscript
# -*- encoding: utf-8 -*-

#@File    :   composition_metaphlan.R
#@Time    :   2023/06/09 14:30:50
#@Author  :   Wanjin.Hu 
#@Version :   1.0
#@Contact :   wanjin.hu@diprobio.com
#@Description :

composition <- function(tax_lev) {
  "
  @功能: 对格式为<06_metaphlan_species.txt>的物种表做柱状组成图
  @参数tax_lev: 在哪个水平做组成分析
  "
  split_data <- strsplit(df$clade_name, ";")
  split_data <- do.call("rbind",split_data)
  cols <- df[, c(2:ncol(df))]
  df_merge <- do.call("cbind", list(split_data,cols))
  names(df_merge)[1:7] <- c("kindom","phylum","class","order","family","genus","species")
  spe <- df_merge[-c(1:7)]
  tax <- df_merge[1:7]
  phy <- spe %>%
    group_by(tax[[tax_lev]]) %>% # 使用tax中的X水平进行分类
    summarise_all(sum) %>%
    rename_at(1, ~ tax_lev) %>%
    gather(key="samples",value = "abun",- tax_lev) %>% # 数据形式转换：“宽”转“长”
    left_join(group,by=c("samples"="sample"))
  p1 <-
    ggplot(phy,aes(samples,abun,fill=.data[[tax_lev]])) + 
    geom_bar(stat = "identity",position = "fill") +
    facet_grid(. ~ group, scales = "free", space = "free_x") +
    labs(x="",y="Percent")+
    scale_fill_manual(values = mypal)+labs(fill="")+
    theme_bw()+
    theme(legend.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, color = "black"))+
    scale_y_continuous(expand=c(0,0))
  return(p1)
}

library(getopt)
command <- matrix(c(
  "help", "h", 0, "logical", "Help docs",
  "data", "d", 1, "character", "data file",
  "group", "g", 1, "character", "group file",
  "tax", "t", 1, "character", "which lev you graph"),
  byrow = TRUE, ncol = 5
)
args <- getopt(command)
if (! is.null(args$help) || 
    is.null(args$data) || 
    is.null(args$group) ||
    is.null(args$tax)) {
  cat(paste(getopt(command, usage = TRUE), "\n"))
  q(status = 1)
}

library(tidyverse)
library(ggsci)
library(ggplot2)
df <- read.table(args$data,
                 header = T,
                 sep = "\t",
                 check.names = F)
group <- read.table(args$group,
                    header = T,
                    sep = "\t",
                    check.names = F)
tax <- args$tax
# 颜色设置
colors_lancet <- pal_lancet("lanonc",alpha = 1)(9)
colors_npg <- pal_npg("nrc",alpha = 1)(10)
colors_d3_20 <- pal_d3("category20",alpha = 1)(20)
colors_igv <- pal_igv("default",alpha = 1)(51)
mypal=c(colors_lancet,colors_npg,colors_d3_20,colors_igv)
p1 <- composition(tax_lev = tax)
pdf(paste0(tax,"_percent.pdf"),width = 12,height = 10,family="Times")
p1
dev.off()
message("Finish composition plot!")

## 本地测试
# df <- read.table("06_metaphlan_species.txt",
#                  header = T,
#                  sep = "\t",
#                  check.names = F)
# group <- read.table("group.txt",
#                     header = T,
#                     sep = "\t",
#                     check.names = F)
# col=pal_d3("category20")(20)
# col2 = pal_d3("category20",alpha = 0.5)(20)
# mypal=c(col,col2)
# p1 <- composition(tax_lev = "phylum")

# pdf("phylum_percent.pdf",width = 12,height = 10,family="Times")
# p1
# dev.off()
