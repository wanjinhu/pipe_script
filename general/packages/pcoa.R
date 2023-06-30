#!/usr/bin/Rscript
# -*- encoding: utf-8 -*-

#@File    :   pcoa.R
#@Time    :   2023/06/08 18:02:48
#@Author  :   Wanjin.Hu 
#@Version :   1.0
#@Contact :   wanjin.hu@diprobio.com
#@Description :

plot_pcoa <- function(data,sd) {
  "
  @功能: beta多样性计算/pcoa绘图
  @参数data: 读入的组成表
  @参数sd: 读入的分组文件表
  @return: list(pcoa_no_label,pcoa_with_label,pcoa_permanova)
  "
  dataT <- t(data)
  dist <- vegdist(dataT, method="bray")
  dist <- as.matrix(dist)
  adist<-as.dist(dist)
  rownames(sd) <- as.character(sd[,1])
  sd[,2] <- as.character(sd[,2])
  dist <- dist[as.character(sd[,1]),][,as.character(sd[,1])]
  pc_num <- as.numeric(unlist(strsplit("1-2","-")))
  pc_x <- pc_num[1]
  pc_y <- pc_num[2]
  pca <- cmdscale(dist, k=3, eig=TRUE)
  pc12 <- pca$points[,pc_num]
  pc <- round(pca$eig/sum(pca$eig)*100,digits = 2)
  pc12 <- as.data.frame(pc12)
  colnames(pc12) <- c("pc_x","pc_y")
  pc12['sample'] <- rownames(pc12)
  colnames(sd)[1:2] <- c("sample","group")
  sd$group<-factor(sd$group,levels=sd$group[!duplicated(sd$group)])
  pc12 <- merge(pc12,sd,by="sample")
  mex <- 0.2*abs(max(pc12$pc_x)-min(pc12$pc_x))
  mey <- 0.2*abs(max(pc12$pc_y)-min(pc12$pc_y))
  mytheme <- theme(panel.border = element_rect(color = "black",linewidth = 1.0,fill = NA),
                   # axis.text.x = element_blank(),
                   # axis.text.y = element_blank(),
                   #axis.ticks = element_blank(),
                   text = element_text(size=18))
  # plot.margin = unit(c(1, 0.5, 1, 1), "lines"))
  pc12$group<-factor(pc12$group,levels=levels(sd$group))
  # p1出的结果是没有label的
  p1 <- ggscatter(pc12, x = "pc_x", y = "pc_y",
                  color = "group", shape = "group", palette = "npg", size=3)+
    ylab(paste0("PCoA",pc_y,"(",round(pc[pc_y],2),"%",")"))+
    xlab(paste0("PCoA",pc_x,"(",round(pc[pc_x],2),"%",")"))+
    geom_hline(yintercept = 0, color = '#B3B3B3', linetype = "longdash")+
    geom_vline(xintercept = 0, color = '#B3B3B3', linetype = "longdash")+
    theme(legend.position = "right",legend.title = element_blank())
  p1 <- p1 + stat_ellipse(aes(x = pc_x, y = pc_y, color = group), 
                          linetype = 1, level = 0.95) + mytheme
  # p2出的结果是有label的
  p2 <- ggscatter(pc12, x = "pc_x", y = "pc_y",
                  color = "group", shape = "group", palette = "npg", size=2)+
    ylab(paste0("PCoA",pc_y,"(",round(pc[pc_y],2),"%",")"))+
    xlab(paste0("PCoA",pc_x,"(",round(pc[pc_x],2),"%",")"))+
    geom_hline(yintercept = 0, color = '#B3B3B3', linetype = "longdash")+
    geom_vline(xintercept = 0, color = '#B3B3B3', linetype = "longdash")+
    geom_text_repel(aes(pc_x, pc_y, label = sample),max.overlaps = 100)+
    theme(legend.position = "right",legend.title = element_blank())
  p2 <- p2 +  stat_ellipse(aes(x = pc_x, y = pc_y, color = group), 
                           linetype = 1, level = 0.95) + mytheme
  # 利用adonis函数进行PERMANOVA分析
  ADONIS<-adonis(dist~sd$group)
  TEST<-ADONIS$aov.tab$`Pr(>F)`[1]
  R2adonis<-round(ADONIS$aov.tab$R2[1],digits = 3)
  # sink(paste0(prefix,'_adonis.txt'))
  # print(ADONIS)
  # sink()
  # 将PERMANOVA计算结果分别加入到p1和p2的结果图中
  xpos<-ggplot_build(p1)$layout$panel_scales_x[[1]]$range$range
  ypos<-ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range
  p1 <- p1 + geom_text(aes(x=xpos[1],y=ypos[2]*1.1),
                       label=paste("PERMANOVA, P","=",TEST,sep = ''),
                       size=6,hjust=0)
  xpos<-ggplot_build(p2)$layout$panel_scales_x[[1]]$range$range
  ypos<-ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range
  p2 <- p2 + geom_text(aes(x=xpos[1],y=ypos[2]*1.1),
                       label=paste("PERMANOVA, P","=",TEST,sep = ''),
                       size=6,hjust=0)
  return(list(pcoa_no_label=p1,
              pcoa_with_label=p2,
              pcoa_permanova=ADONIS))
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
library(ggpubr)
library(reshape2)
library(patchwork)
library(ggrepel)
data <- read.delim(args$data,
                   sep = "\t",
                   header = T,
                   check.names = FALSE,
                   row.names = 1)
group <- read.delim(args$group,
                    header = T,
                    sep = "\t",
                    check.names = F)
names(group) <- c("sample","group")
result <- plot_pcoa(data=data, sd=group)
ggsave("Pcoa_withlabel.pdf",result$pcoa_with_label,width = 10,height = 10)
ggsave("Pcoa_nolabel.pdf",result$pcoa_no_label,width = 10,height = 10)
sink("Pcoa_permanova.txt")
print(result$pcoa_permanova)
sink()
message("Finish beta pcoa plot!")

## 本地测试
# data <- read.delim("06_metaphlan_species.txt",
#                    sep = "\t",
#                    header = T,
#                    check.names = FALSE,
#                    row.names = 1)
# group <- read.delim("group.txt",
#                     header = T,
#                     sep = "\t",
#                     check.names = F)
# names(group) <- c("sample","group")
# result <- plot_pcoa(data=data, sd=group)
# ggsave("Pcoa_withlabel.pdf",result$pcoa_with_label,width = 10,height = 10)
# ggsave("Pcoa_nolabel.pdf",result$pcoa_no_label,width = 10,height = 10)
# sink("Pcoa_permanova.txt")
# print(result$pcoa_permanova)
# sink()

