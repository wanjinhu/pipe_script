## 函数集合
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
    ggplot(dataCom, aes(as.factor(group), .data[[target]], fill=group)) + 
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
          axis.line = element_line(colour = "black", linewidth = 1),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + 
    labs(y = target, col = "group") + 
    stat_compare_means(method = "wilcox.test", 
                       comparisons = compare_list,
                       label = "p.signif")
  return(p1)
}

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

com_plot <- function(tax_lev) {
  "
  @功能: 对格式为<06_metaphlan_species.txt>的物种表做柱状组成图,SGB信息归到Others
  @参数tax_lev: 在哪个水平做组成分析
  "
  # 颜色设置
  colors_lancet <- pal_lancet("lanonc",alpha = 1)(9)
  colors_npg <- pal_npg("nrc",alpha = 1)(10)
  colors_d3_20 <- pal_d3("category20",alpha = 1)(20)
  colors_igv <- pal_igv("default",alpha = 1)(51)
  if (tax_lev == "phylum") {
    mypal=c(colors_npg[1:8],colors_d3_20[1:10],colors_igv)
  } else {
    mypal=c("#C7C7C7FF",colors_npg[1:8],colors_d3_20[1:10],colors_igv)
  }
  split_data <- strsplit(df$clade_name, ";")
  split_data <- do.call("rbind",split_data)
  cols <- df[, c(2:ncol(df))]
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
  # group中不一定包含全部样本，上一步会为不存在的样本增加NA的group，这里删除NA行
  phy_plot_noNA <- phy_plot[!is.na(phy_plot$group),]
  phy_plot_noNA$class <- factor(phy_plot_noNA$class,levels = phy$class)
  p1 <-
    ggplot(phy_plot_noNA,aes(samples,abun,fill=.data[["class"]])) + 
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

