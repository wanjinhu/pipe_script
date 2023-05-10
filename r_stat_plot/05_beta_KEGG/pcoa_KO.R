
library(vegan)
library(ggpubr)
library(reshape2)
library(patchwork)
library(ggrepel)

#读取样本文件并且生成距离矩阵
setwd("F:/03-Dipro项目/11-张瑞-宏基因组-犬类/analysis/05_beta_KEGG/")
data <- read.table("../KO_samples.xls",
                   sep="\t",
                   header = T,
                   row.names = 3,
                   check.names = FALSE,
                   quote = "")

data <- data[,-c(1:2)]
dataT <- t(data)
dist <- vegdist(dataT, method="bray")
dist <- as.matrix(dist)
adist<-as.dist(dist)

#读取分组文件
options(stringsAsFactors=F)
sd <- read.table("../group2.txt",
                 head=T,
                 sep="\t",
                 comment.char = "",
                 check.names = FALSE)

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

mytheme <- theme(panel.border = element_rect(color = "black",size = 1.0,fill = NA),
                 # axis.text.x = element_blank(),
                 # axis.text.y = element_blank(),
                 #axis.ticks = element_blank(),
                 text = element_text(size=18))
#plot.margin = unit(c(1, 0.5, 1, 1), "lines"))
pc12$group<-factor(pc12$group,levels=levels(sd$group))

#p1出的结果是没有label的
p1 <- ggscatter(pc12, x = "pc_x", y = "pc_y",
                color = "group", shape = "group", palette = "npg", size=3)+
  ylab(paste0("PCoA",pc_y,"(",round(pc[pc_y],2),"%",")"))+
  xlab(paste0("PCoA",pc_x,"(",round(pc[pc_x],2),"%",")"))+
  geom_hline(yintercept = 0, color = '#B3B3B3', linetype = "longdash")+
  geom_vline(xintercept = 0, color = '#B3B3B3', linetype = "longdash")+
  theme(legend.position = "bottom",legend.title = element_blank())

p1 <- p1 +  stat_ellipse(aes(x = pc_x, y = pc_y, color = group), linetype = 1, level = 0.95)
p1 <- p1 + mytheme
p1

#p2出的结果是有label的
p2 <- ggscatter(pc12, x = "pc_x", y = "pc_y",
                color = "group", shape = "group", palette = "npg", size=2)+
  ylab(paste0("PCoA",pc_y,"(",round(pc[pc_y],2),"%",")"))+
  xlab(paste0("PCoA",pc_x,"(",round(pc[pc_x],2),"%",")"))+
  geom_hline(yintercept = 0, color = '#B3B3B3', linetype = "longdash")+
  geom_vline(xintercept = 0, color = '#B3B3B3', linetype = "longdash")+
  geom_text_repel(aes(pc_x, pc_y, label = sample),max.overlaps = 100)+
  theme(legend.position = "bottom",legend.title = element_blank())

p2 <- p2 +  stat_ellipse(aes(x = pc_x, y = pc_y, color = group), linetype = 1, level = 0.95)
p2 <- p2 + mytheme
p2

#利用adonis函数进行PERMANOVA分析，并保存计算文件
ADONIS<-adonis(dist~sd$group)
TEST<-ADONIS$aov.tab$`Pr(>F)`[1]
R2adonis<-round(ADONIS$aov.tab$R2[1],digits = 3)
sink('adonis.txt')
print(ADONIS)
sink()

#将PERMANOVA计算结果分别加入到p1和p2的结果图中，保存为pdf格式的文件
xpos<-ggplot_build(p1)$layout$panel_scales_x[[1]]$range$range
ypos<-ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range
pdf("pcoaNotWithLabel.pdf",width = 10)
p1+geom_text(aes(x=xpos[1],y=ypos[2]*1.1),label=paste("PERMANOVA, P","=",TEST,sep = ''),size=6,hjust=0)
dev.off()

xpos<-ggplot_build(p2)$layout$panel_scales_x[[1]]$range$range
ypos<-ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range
pdf("pcoaWithLabel.pdf",width = 10)
p2+geom_text(aes(x=xpos[1],y=ypos[2]*1.1),label=paste("PERMANOVA, P","=",TEST,sep = ''),size=6,hjust=0)
dev.off()
