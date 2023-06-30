library(microeco)
library(magrittr)
library(ggplot2)
library(ggh4x)
library(ggalluvial)
library(ggtern)
library(ggnested)
library(ggradar)
library(FSA)
library(rcompanion)
library(agricolae)
library(ggdendro)
library(file2meco)

# load example data
data(sample_info_16S)
data(otu_table_16S)
data(taxonomy_table_16S)
data(phylo_tree_16S)
data(env_data_16S)
set.seed(123)
theme_set(theme_bw())
# clean taxonomy data
taxonomy_table_16S %<>% tidy_taxonomy
# Let's create a microtable object with more information
dataset <- microtable$new(sample_table = sample_info_16S, 
                          otu_table = otu_table_16S, 
                          tax_table = taxonomy_table_16S, 
                          phylo_tree = phylo_tree_16S)
dataset
## -----------------------------------------------
## composition class
## -----------------------------------------------
# abundance
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 8)
t1$plot_bar(others_color = "grey70", facet = "Group", xtext_keep = FALSE, legend_text_italic = FALSE)
# return a ggplot2 object
# require package ggh4x, first run install.packages("ggh4x") if not installed
t1$plot_bar(others_color = "grey70", facet = c("Group", "Type"), xtext_keep = FALSE, legend_text_italic = FALSE, barwidth = 1)

t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 8)
# require ggalluvial package
# use_alluvium = TRUE make the alluvial plot, clustering =TRUE can be used to reorder the samples by clustering
# bar_type = "notfull" can discard 'others'; select another color palette
p <- t1$plot_bar(bar_type = "notfull", use_alluvium = TRUE, clustering = TRUE, xtext_angle = 30, xtext_size = 3, color_values = RColorBrewer::brewer.pal(8, "Set2"))
# The groupmean parameter can be used to obtain the group-mean barplot.
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 10, groupmean = "Group")
g1 <- t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
g1 + theme_classic() + theme(axis.title.y = element_text(size = 18))
# show 15 taxa at Class level
t1 <- trans_abund$new(dataset = dataset, taxrank = "Class", ntaxa = 15)
t1$plot_box(group = "Group", xtext_angle = 30)
# show 40 taxa at Genus level
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 40)
t1$plot_heatmap(facet = "Group", xtext_keep = FALSE, withmargin = FALSE)
# Line chart
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 5)
t1$plot_line()
t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 5, groupmean = "Type")
t1$plot_line(position = position_dodge(0.3), xtext_angle = 0)
# pie chart
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 6, groupmean = "Group")
# all pie chart in one row
t1$plot_pie(facet_nrow = 1)
t1$plot_pie(facet_nrow = 1, add_label = TRUE)
# donut chart
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 8, groupmean = "Group")
t1$plot_donut(label = FALSE)
t1$plot_donut(label = TRUE)
# radar charts
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 8, groupmean = "Group")
t1$plot_radar(values.radar = c("0%", "25%", "50%"), grid.min = 0, grid.mid = 0.25, grid.max = 0.5)
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 8, groupmean = "Type")
t1$plot_radar(values.radar = c("0%", "25%", "50%"), grid.min = 0, grid.mid = 0.25, grid.max = 0.5)
# ternary plot
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 8, groupmean = "Group")
t1$plot_tern()
# require ggnested package; see https://chiliubio.github.io/microeco_tutorial/intro.html#dependence
test1 <- trans_abund$new(dataset = dataset, taxrank = "Class", ntaxa = 10, high_level = "Phylum", prefix = "\\|")
test1$plot_bar(ggnested = TRUE, facet = c("Group", "Type"), xtext_angle = 30)
# fixed number in each phylum
test1 <- trans_abund$new(dataset = dataset, taxrank = "Class", ntaxa = 30, show = 0, high_level = "Phylum", high_level_fix_nsub = 4)
test1$plot_bar(ggnested = TRUE, xtext_angle = 30, facet = c("Group", "Type"))
# sum others in each phylum
test1 <- trans_abund$new(dataset = dataset, taxrank = "Class", ntaxa = 20, show = 0, high_level = "Phylum", high_level_fix_nsub = 3, prefix = "\\|")
test1$plot_bar(ggnested = TRUE, high_level_add_other = TRUE, xtext_angle = 30, facet = c("Group", "Type"))
# clustering plot in the bar plot
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 10, groupmean = "Group")
g1 <- t1$plot_bar(coord_flip = TRUE)
g1 <- g1 + theme_classic() + theme(axis.title.x = element_text(size = 16), axis.ticks.y = element_blank(), axis.line.y = element_blank())
g1
g1 <- t1$plot_bar(clustering_plot = TRUE)
# In this case, g1 (aplot object) is the combination of different ggplot objects
# to adjust the main plot, please select g1[[1]]
g1[[1]] <- g1[[1]] + theme_classic() + theme(axis.title.x = element_text(size = 16), axis.ticks.y = element_blank(), axis.line.y = element_blank())
g1
## venn
# merge samples as one community for each group
dataset1 <- dataset$merge_samples(use_group = "Group")
# dataset1 is a new microtable object
# create trans_venn object
t1 <- trans_venn$new(dataset1, ratio = NULL)
t1$plot_venn()
# create venn plot with more information
t1 <- trans_venn$new(dataset1, ratio = "seqratio")
t1$plot_venn()
# The integer is OTU number
# The percentage data is the sequence number/total sequence number
# use "Type" column in sample_table
dataset1 <- dataset$merge_samples(use_group = "Type")
t1 <- trans_venn$new(dataset1)
t1$plot_venn(petal_plot = TRUE, petal_color = RColorBrewer::brewer.pal(8, "Dark2"))
t1$plot_venn(petal_plot = TRUE, petal_center_size = 50, petal_r = 1.5, petal_a = 3, petal_move_xy = 3.8, petal_color_center = "#BEBADA")
## UpSet
tmp <- dataset$merge_samples(use_group = "Type")
tmp
t1 <- trans_venn$new(dataset = tmp)
# only show some sets with large intersection numbers
t1$data_summary %<>% .[.[, 1] > 20, ]
g1 <- t1$plot_bar(left_plot = TRUE, bottom_height = 0.5, left_width = 0.15, up_bar_fill = "grey50", left_bar_fill = "grey50", bottom_point_color = "black")
g1
# g1 is aplot class and can be saved with ggplot2::ggsave, aplot::ggsave or cowplot::save_plot function
# as g1 is comprised of several sub-plots, please adjust the details for each sub-plot
g1[[1]]
g1[[2]]
dataset1 <- dataset$merge_samples(use_group = "Group")
t1 <- trans_venn$new(dataset1)
# transform venn results to the sample-species table, here do not consider abundance, only use presence/absence.
t2 <- t1$trans_comm(use_frequency = TRUE)
# t2 is a new microtable class, each part is considered a sample
class(t2)
# calculate taxa abundance, that is, the frequency
t2$cal_abund()
# transform and plot
t3 <- trans_abund$new(dataset = t2, taxrank = "Genus", ntaxa = 8)
t3$plot_bar(bar_type = "part", legend_text_italic = T, xtext_angle = 30, color_values = RColorBrewer::brewer.pal(8, "Set2"),
            order_x = c("IW", "CW", "TW", "IW&CW", "IW&TW", "CW&TW", "IW&CW&TW")) + ylab("Frequency (%)")
t3 <- trans_abund$new(dataset = t2, taxrank = "Phylum", ntaxa = 8)
t3$data_abund$Sample %<>% factor(., levels = unique(.))
t3$plot_pie(facet_nrow = 3, color_values = c(RColorBrewer::brewer.pal(8, "Dark2"), "grey50"))
# 
test <- dataset$merge_samples(use_group = "Type")
test$sample_table %<>% .[c("YML", "NE", "NW", "NC", "QTP", "SC"), , drop = FALSE]
test$tidy_dataset()
# The columns of otu_table can also be reordered according to the sample_table after running tidy_dataset function
t1 <- trans_venn$new(test)
t1$plot_bar(sort_samples = FALSE)
test <- dataset$merge_samples(use_group = "Type")
t1 <- trans_venn$new(test)
# remove left bar in the UpSet plot
t1$plot_bar(left_plot = FALSE)
# sort samples in the axis according to the number of features
t1$plot_bar(sort_samples = TRUE)
# original orders in test$sample_table
t1$plot_bar(sort_samples = FALSE)

## -----------------------------------------------
## diversity class
## -----------------------------------------------
## alpha
t1 <- trans_alpha$new(dataset = dataset, group = "Group")
# return t1$data_stat
head(t1$data_stat)
# alpha的差异检验
t1$cal_diff(method = "KW")
# return t1$res_diff
head(t1$res_diff)
t1$cal_diff(method = "KW_dunn")
# return t1$res_diff
head(t1$res_diff)
# more options
t1$cal_diff(method = "KW_dunn", KW_dunn_letter = FALSE)
head(t1$res_diff)
t1$cal_diff(method = "wilcox")
head(t1$res_diff)
t1$cal_diff(method = "t.test")
head(t1$res_diff)
t1$cal_diff(method = "anova")
# return t1$res_diff
head(t1$res_diff)
# two-way anova
t1 <- trans_alpha$new(dataset = dataset, group = "Group")
t1$cal_diff(method = "anova", formula = "Group+Type")
head(t1$res_diff)
# plot_alpha
t1$cal_diff(method = "anova")
t1$plot_alpha(measure = "Chao1")
t1$plot_alpha(measure = "Chao1", boxplot_add = "jitter", order_x_mean = TRUE, add_sig_text_size = 6)
t1$cal_diff(method = "wilcox")
t1$plot_alpha(measure = "Chao1", shape = "Group")
# remove the 'ns' in the label
t1$res_diff %<>% base::subset(Significance != "ns")
t1$plot_alpha(measure = "Chao1", boxplot_add = "dotplot", xtext_size = 15)
# group => by_group
t1 <- trans_alpha$new(dataset = dataset, group = "Type", by_group = "Group")
t1$cal_diff(method = "wilcox")
t1$plot_alpha(measure = "Shannon")
## beta
# pcoa
# create an trans_beta object
# measure parameter must be one of names(dataset$beta_diversity)
data(dataset)
t1 <- trans_beta$new(dataset = dataset, measure = "bray", group = "Group")
# PCoA, PCA and NMDS are available
t1$cal_ordination(ordination = "PCoA")
# t1$res_ordination is the ordination result list
class(t1$res_ordination)
# plot the PCoA result with confidence ellipse
t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "ellipse"))
# pcoa相关的调整
t1$plot_ordination(plot_color = "Type", plot_type = "point")
t1$plot_ordination(plot_color = "Group", point_size = 5, point_alpha = .2, plot_type = c("point", "ellipse"), ellipse_chull_fill = FALSE)
t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "centroid"))
t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "ellipse", "centroid"))
t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "chull"))
t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "chull", "centroid"))
t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("chull", "centroid"))
t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = c("point", "chull", "centroid"), add_sample_label = "SampleID")
t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = "centroid")
t1$plot_ordination(plot_color = "Group", plot_shape = "Group", plot_type = "centroid", centroid_segment_alpha = 0.9, centroid_segment_size = 1, centroid_segment_linetype = 1)
t1$plot_ordination(plot_type = c("point", "centroid"), plot_color = "Type", centroid_segment_linetype = 1)
t1$plot_ordination(plot_color = "Saline", point_size = 5, point_alpha = .2, plot_type = c("point", "chull"), ellipse_chull_fill = FALSE, ellipse_chull_alpha = 0.1)
t1$plot_ordination(plot_color = "Group") + theme(panel.grid = element_blank()) + geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2)
# beta多样性柱状图
# calculate and plot sample distances within groups
t1$cal_group_distance(within_group = TRUE)
# return t1$res_group_distance
# perform Wilcoxon Rank Sum and Signed Rank Tests
t1$cal_group_distance_diff(method = "wilcox")
# plot_group_order parameter can be used to adjust orders in x axis
t1$plot_group_distance(boxplot_add = "mean")
# calculate and plot sample distances between groups
t1$cal_group_distance(within_group = FALSE)
t1$cal_group_distance_diff(method = "wilcox")
t1$plot_group_distance(boxplot_add = "mean")
# Clustering plot 
# extract a part of data
d1 <- clone(dataset)
d1$sample_table %<>% subset(Group %in% c("CW", "TW"))
d1$tidy_dataset()
t1 <- trans_beta$new(dataset = d1, group = "Group")
# use replace_name to set the label name, group parameter used to set the color
t1$plot_clustering(group = "Type", replace_name = c("Type"))
# PerMANOVA
# manova for all groups when manova_all = TRUE
t1$cal_manova(manova_all = TRUE)
t1$res_manova
# manova for each paired groups
t1$cal_manova(manova_all = FALSE)
t1$res_manova
# manova for specified group set: such as "Group + Type"
t1$cal_manova(manova_set = "Group + Type")
t1$res_manova
# for the whole comparison and for each paired groups
t1$cal_betadisper()

## -----------------------------------------------
## diff class
## -----------------------------------------------
## lefse
t1 <- trans_diff$new(dataset = dataset,
                     taxa_level = "all",
                     method = "lefse", 
                     group = "Group", 
                     alpha = 0.01, 
                     lefse_subgroup = NULL)
# see t1$res_diff for the result
# From v0.8.0, threshold is used for the LDA score selection.
t1$plot_diff_bar(threshold = 4)
# we show 20 taxa with the highest LDA (log10)
t1$plot_diff_bar(use_number = 1:30, width = 0.8, group_order = c("CW", "IW", "TW"))
# show part of the table
t1$res_diff[1:5, c(1, 3, 4, 6)]
# plot the abundance of biomarkers detected by LEfSe
t1$plot_diff_abund(use_number = 1:30, group_order = c("CW", "IW", "TW"))
# clade_label_level 5 represent phylum level in this analysis
# require ggtree package
t1$plot_diff_cladogram(use_taxa_num = 200, use_feature_num = 50, clade_label_level = 5, group_order = c("CW", "IW", "TW"))
# choose some taxa according to the positions in the previous picture; those taxa labels have minimum overlap
use_labels <- c("c__Deltaproteobacteria", "c__Actinobacteria", "o__Rhizobiales", "p__Proteobacteria", "p__Bacteroidetes", 
                "o__Micrococcales", "p__Acidobacteria", "p__Verrucomicrobia", "p__Firmicutes", 
                "p__Chloroflexi", "c__Acidobacteria", "c__Gammaproteobacteria", "c__Betaproteobacteria", "c__KD4-96",
                "c__Bacilli", "o__Gemmatimonadales", "f__Gemmatimonadaceae", "o__Bacillales", "o__Rhodobacterales")
# then use parameter select_show_labels to show
t1$plot_diff_cladogram(use_taxa_num = 200, use_feature_num = 50, select_show_labels = use_labels)
## random forest
# use Genus level for parameter taxa_level, if you want to use all taxa, change to "all"
# nresam = 1 and boots = 1 represent no bootstrapping and use all samples directly
t1 <- trans_diff$new(dataset = dataset, method = "rf", group = "Group", taxa_level = "Genus")
# plot the MeanDecreaseGini bar
# group_order is designed to sort the groups
g1 <- t1$plot_diff_bar(use_number = 1:20, group_order = c("TW", "CW", "IW"))
# plot the abundance using same taxa in g1
g2 <- t1$plot_diff_abund(group_order = c("TW", "CW", "IW"), select_taxa = t1$plot_diff_bar_taxa)
# now the y axis in g1 and g2 is same, so we can merge them
# remove g1 legend; remove g2 y axis text and ticks
g1 <- g1 + theme(legend.position = "none")
g2 <- g2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
gridExtra::grid.arrange(g1, g2, ncol = 2, nrow = 1, widths = c(2, 1.7))
## wilcox
t1 <- trans_diff$new(dataset = dataset, method = "wilcox", group = "Group", taxa_level = "Genus", filter_thres = 0.001)
# filter something not needed to show
t1$res_diff %<>% subset(Significance %in% "***")
t1$plot_diff_abund(use_number = 1:10, add_sig = T, add_sig_label = "Significance")
## anova
t1 <- trans_diff$new(dataset = dataset, method = "anova", group = "Group", taxa_level = "Genus", filter_thres = 0.001)
t1$plot_diff_abund(use_number = 1:10, add_sig = T, coord_flip = F)
## metastat
# metastat analysis at Genus level
t1 <- trans_diff$new(dataset = dataset, method = "metastat", group = "Group", taxa_level = "Genus")
# t1$res_diff is the differential test result
# t1$res_abund is the group abundance
# select_group should be one of groups in t1$res_diff$Comparison
t1$plot_diff_abund(use_number = 1:20, select_group = "CW - TW", coord_flip = F)
# examples for the methods ‘KW’, ‘KW_dunn’, ‘wilcox’, ‘t.test’ and ‘anova’.
# Kruskal-Wallis Rank Sum Test for all groups (>= 2)
t1 <- trans_diff$new(dataset = dataset, method = "KW", group = "Group", taxa_level = "all", filter_thres = 0.001)
t1$plot_diff_abund(use_number = 1:20)
# Dunn's Kruskal-Wallis Multiple Comparisons when group number > 2; require FSA package
t1 <- trans_diff$new(dataset = dataset, method = "KW_dunn", group = "Group", taxa_level = "Genus", filter_thres = 0.0001)
t1$plot_diff_abund(use_number = 1:10, add_sig = T, coord_flip = F)
# Wilcoxon Rank Sum and Signed Rank Tests for all paired groups
t1 <- trans_diff$new(dataset = dataset, method = "wilcox", group = "Group", taxa_level = "Genus", filter_thres = 0.001)
t1$plot_diff_bar(use_number = 1:20, select_group = "CW - TW")
t1$plot_diff_abund(use_number = 1:20, select_group = "CW - TW", group_order = c("TW", "CW"))
# t.test
t1 <- trans_diff$new(dataset = dataset, method = "t.test", group = "Group", taxa_level = "all", filter_thres = 0.001)
# anova
t1 <- trans_diff$new(dataset = dataset, method = "anova", group = "Group", taxa_level = "Phylum", filter_thres = 0.001)
head(t1$res_diff)

## --------------------
## 其他
## --------------------
library(microeco)
library(file2meco)
library(magrittr)
?humann2meco
sample_file_path <- system.file("extdata", "example_metagenome_sample_info.tsv", package="file2meco")
match_file_path <- system.file("extdata", "example_metagenome_match_table.tsv", package="file2meco")
# use KEGG pathway based HUMAnN result
abund_file_path <- system.file("extdata", "example_HUMAnN_KEGG_abund.tsv", package="file2meco")
# match_table parameter can be used to replace sample names
test <- humann2meco(abund_file_path, db = "KEGG", sample_table = sample_file_path, match_table = match_file_path)
# remove the unclassified pathway in Level.1
test$tax_table %<>% subset(Level.1 != "unclassified")
test$tidy_dataset()
# rel = FALSE donot use relative abundance, use the raw RPK
test$cal_abund(select_cols = 1:3, rel = FALSE)
# use_percentage = FALSE disable percentage for relative abundance
test1 <- trans_abund$new(test, taxrank = "Level.2", ntaxa = 10, use_percentage = FALSE)
test1$plot_bar(facet = "Group", bar_type = "notfull", xtext_angle = 30) + ggplot2::ylab("Abundance (RPK)")

