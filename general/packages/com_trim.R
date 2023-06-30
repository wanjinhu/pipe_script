#!/usr/bin/Rscript
# -*- encoding: utf-8 -*-

#@File    :   com_trim.R
#@Time    :   2023/06/20 14:30:50
#@Author  :   Wanjin.Hu 
#@Version :   1.0
#@Contact :   wanjin.hu@diprobio.com
#@Description :

com_plot <- function(tax_lev) {
  "
  @功能: 对格式为<06_metaphlan_species.txt>的物种表做柱状组成图,SGB信息归到Others
  @参数tax_lev: 在哪个水平做组成分析
  "
  df_1 <- df %>% select(clade_name,group$sample)
  split_data <- strsplit(df_1$clade_name, ";")
  split_data <- do.call("rbind",split_data)
  cols <- df_1[, c(2:ncol(df_1))]
  df_merge <- do.call("cbind", list(split_data,cols))
  names(df_merge)[1:7] <- c("kindom",
                            "phylum",
                            "class",
                            "order",
                            "family",
                            "genus",
                            "species")
  spe <- df_merge[-c(1:7)]
  tax <- df_merge[1:7]
  if (tax_lev == "species") {
    mer_lev = "_SGB"
  } else if (tax_lev == "genus") {
    mer_lev = "_GGB"
  } else if (tax_lev == "family") {
    mer_lev = "_FGB"
  } else if (tax_lev == "order") {
    mer_lev = "_OFGB"
  } else if (tax_lev == "class") {
    mer_lev = "_CFGB"
  } else  {
    mer_lev = "-"
  }
  # phylum不需要处理，展示全部的
  if (tax_lev == "phylum") {
    phy <- spe %>%
      group_by(tax[[tax_lev]]) %>% 
      summarise_all(sum) %>%
      rename_at(1, ~ "class")
    phy <- phy %>% arrange(desc(rowSums(phy[,2:ncol(phy)])))
  } else {
    # 举例如果是class水平，合并包含"_CFGB"为Others
    phy1 <- spe %>%
      group_by(tax[[tax_lev]]) %>%
      summarise_all(sum) %>%
      rename_at(1, ~ "class") %>%
      filter(str_detect(class, mer_lev)) %>%
      summarise(across(-1, sum)) %>%
      add_column(class = "Others", .before = 1)
    # 举例如果是class水平，除了"_CFGB"之外的
    phy2 <- spe %>%
      group_by(tax[[tax_lev]]) %>%
      summarise_all(sum) %>%
      rename_at(1, ~ "class") %>%
      filter(!str_detect(class, mer_lev))
    phy2 <- phy2 %>% arrange(desc(rowSums(phy2[,2:ncol(phy2)])))
    phy <- rbind(phy1,phy2)
  }
  write.table(phy,
              file = paste0("trim_",tax_lev,".xls"),
              sep = "\t",
              col.names = NA,
              quote = F)
  phy_plot <- phy %>% 
    gather(key="samples", value = "abun", - "class") %>% # 数据形式转换：“宽”转“长”
    left_join(group,by=c("samples"="sample"))
  phy_plot$class <- factor(phy_plot$class,levels = phy$class)
  p1 <-
    ggplot(phy_plot,aes(samples,abun,fill=.data[["class"]])) + 
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
  "data", "d", 1, "character", "data file[06_metaphlan_species.txt]",
  "group", "g", 1, "character", "group file[group.txt]",
  "tax", "t", 1, "character", "which lev[phylum,class,order,family,genus,species]"),
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
df <- read.table(args$data, # 06_metaphlan_species.txt
                 header = T,
                 sep = "\t",
                 check.names = F)
group <- read.table(args$group,
                    header = T,
                    sep = "\t",
                    check.names = F)
names(group) <- c("sample","group")
tax <- args$tax

# 颜色设置
colors_lancet <- pal_lancet("lanonc",alpha = 1)(9)
colors_npg <- pal_npg("nrc",alpha = 1)(10)
colors_d3_20 <- pal_d3("category20",alpha = 1)(20)
colors_igv <- pal_igv("default",alpha = 1)(51)
if (tax == "phylum") {
  mypal=c(colors_npg[1:8],colors_d3_20[1:10],colors_igv)
} else {
  mypal=c("#C7C7C7FF",colors_npg[1:8],colors_d3_20[1:10],colors_igv)
}

p1 <- com_plot(tax_lev = tax)
pdf(paste0(tax,"_percent.pdf"),width = 12,height = 10,family="Times")
p1
dev.off()
message(paste0("Finish ",tax," analysis!"))

