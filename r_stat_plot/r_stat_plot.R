## 目的: 整理早就应该整理的数据处理-绘图的一些代码
## 方式: 定义具体的函数,处理具体问题,函数命名规则如下:
# 主要功能为数据处理: stat_xx_1
# 主要功能为绘图: plot_xx_1
# 函数都设计为返回值的形式,避免直接生成文件

# #####################################################
# ################# 主题和颜色 ########################
# #####################################################



# #####################################################
# ################# 定义的函数 ########################
# #####################################################
## 1. alpha多样性计算
stat_alpha_1 <- function(x, tree = NULL, base = 2) {
  "
  @功能: alpha多样性计算
  @参数x: 物种组成表
  @参数tree: 物种序列进化树[如果有的话]
  "
  library(vegan)
  library(picante)
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

## 2. beta_pcoa绘图
plot_pcoa_1 <- function(data,sd) {
  "
  @功能: beta多样性计算/pcoa绘图
  @参数data: 读入的组成表
  @参数sd: 读入的分组文件表
  "
  library(vegan)
  library(ggpubr)
  library(reshape2)
  library(patchwork)
  library(ggrepel)
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
    theme(legend.position = "right",legend.title = element_blank())
  p1 <- p1 + stat_ellipse(aes(x = pc_x, y = pc_y, color = group), 
                           linetype = 1, level = 0.95) + mytheme
  
  #p2出的结果是有label的
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
  #利用adonis函数进行PERMANOVA分析，并保存计算文件
  ADONIS<-adonis(dist~sd$group)
  TEST<-ADONIS$aov.tab$`Pr(>F)`[1]
  R2adonis<-round(ADONIS$aov.tab$R2[1],digits = 3)
  sink(paste0(prefix,'_adonis.txt'))
  print(ADONIS)
  sink()
  #将PERMANOVA计算结果分别加入到p1和p2的结果图中，保存为pdf格式的文件
  xpos<-ggplot_build(p1)$layout$panel_scales_x[[1]]$range$range
  ypos<-ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range
  p1 <- p1 + geom_text(aes(x=xpos[1],y=ypos[2]*1.1),
                       label=paste("PERMANOVA, P","=",TEST,sep = ''),
                       size=6,hjust=0)
  ggsave(paste0(prefix,"_pcoaNotWithLabel.pdf"),p1,width = 12,height = 12)
  
  xpos<-ggplot_build(p2)$layout$panel_scales_x[[1]]$range$range
  ypos<-ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range
  p2 <- p2 + geom_text(aes(x=xpos[1],y=ypos[2]*1.1),
                       label=paste("PERMANOVA, P","=",TEST,sep = ''),
                       size=6,hjust=0)
  ggsave(paste0(prefix,"_pcoaWithLabel.pdf"),p2,width = 12,height = 12)
}

## 3. 柱状图[显著性标注,针对1里面计算得到的alpha数据表]
plot_bar_1 <- function(data,group_file,target,ymax = 10) {
  "
  @功能: 柱状图
  @参数data: factor组成表,第一列为样本,后面几列是factor数据
  @参数group_file: 分组文件,只支持1个分组方式,即第1列为样本,第2列为分组方案
  @参数target: 对哪个factor作图,即参数data的列名
  @参数ymax: 绘制多组时,有kw检验结果,放置图上垂直的位置
  "
  library(ggplot2)
  library(tidyr)
  library(ggpubr)
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
  # 确定多组kw检验的p值在图上的摆放位置
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
  
  if (num == 2) {
    return(p1)
  } else {
    p1 <- p1 + stat_compare_means(label.y = ymax)
    return(p1)
  }
}

