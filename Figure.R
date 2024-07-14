library(microeco)
#?microtable # 展示所有的功能和细节描述
library(magrittr)
library(ggplot2)
library(scales)
library(randomcoloR)
library(circlize)#引用包
library(reshape2)
library(tidyr)
library(tidyverse)
theme_set(theme_bw())

#import data
otu_table_16S_1<-read.csv("16s2.csv",row.names = 1)
#filter group
#otu_table_16S_1 <- otu_table_16S_1[,c(1:7,26:42)]
sample_info_16S_1<-read.csv("SampleIDZ25.csv",row.names = 1)
#filter group
#sample_info_16S_1 <- sample_info_16S_1[c(1:7,26:42),]

taxonomy_table_16S_1<-read.csv("tax.csv",row.names = 1)

dataset <- microtable$new(sample_table = sample_info_16S_1, 
                          otu_table = otu_table_16S_1,
                          tax_table = taxonomy_table_16S_1)

#group
#dataset$sample_table$Group <-  factor(dataset$sample_table$Group, levels = c("LFD","HFD","1711","17_3","2016_49_7","C1_A31","Z25","ZS2058"))

dataset$sample_table$Group <-  factor(dataset$sample_table$Group, levels = c("ZS2058","C1_A31","2016_49_7","Z25","1711","17_3","HFD","LFD"))
#touch data
dataset$tidy_dataset()
dataset$cal_abund()
class(dataset$taxa_abund)
#创建丰度表格
#dir.create("taxa_abund")
#dataset$save_abund(dirpath = "taxa_abund")


#color
#color_compaired <- c("#09f9f5", "#d561dd", "#c93f00", "#ddd53e","#4aef7b", "#e86502", "#931635", "#373bbf", "#9ed84e", "#8249aa", "#99db27", "#e07233", "#ff523f","#ce2523", "#f7aa5d", "#cebb10", "#03827f", "#a1ce4c", "#ef3bb6", "#d66551","#1a918f", "#ff66fc", "#2927c4", "#7149af" ,"#57e559" ,"#8e3af4" ,"#f9a270" ,"#22547f", "#db5e92","#edd05e", "#6f25e8", "#0dbc21", "#280f7a", "#6373ed", "#5b910f" ,"#7b34c1" ,"#0cf29a" ,"#d80fc1","#dd27ce", "#07a301", "#167275", "#391c82", "#2baeb5","#925bea", "#63ff4f")

#color_compaired <- distinctColorPalette(30)
#alpha_dicersity
dataset$tidy_dataset()
print(dataset)
#dataset$rarefy_samples(sample.size = 10000)
#dataset$sample_sums() %>% range
#dataset$cal_alphadiv(PD = FALSE)
#dir.create("alpha_diversity")
#dataset$save_alphadiv(dirpath = "alpha_diversity")

#kingdom
t1 <- trans_abund$new(dataset = dataset, taxrank = "Kingdom", ntaxa = 5)

t1$plot_pie(color_values = c(color_compaired),legend_text_italic = FALSE)

t1$plot_bar(others_color = "grey70", color_values = c(color_compaired),legend_text_italic = FALSE)


#Phylum
t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 20, groupmean = "Group")

t1$plot_pie(color_values = c(color_compaired),legend_text_italic = FALSE)

#color_compaired <- distinctColorPalette(10)
t1$plot_bar(others_color = "grey70", color_values = c(color_compaired),legend_text_italic = FALSE)

t1$plot_bar(others_color = "grey70", color_values = c(color_compaired),legend_text_italic = FALSE, use_alluvium = T)


#Phylum circlize

data <- trans_abund$new(dataset = dataset, taxrank = "Phylum",groupmean = "Group")
data$data_abund
data2 <- as.data.frame(data$data_abund)
data3 <- data2[,-c(3,5,6,7)]

# transform the data into circlize input format
data5 <- dcast(data3,Taxonomy~Sample) #data6 <- spread(data3,key = Sample,value = Abundance)
data5 <- as.data.frame(data5)
rownames(data5) <- data5$Taxonomy
data5 <- data5[,-1]

##cirzlize P1  
phylum2 <- as.matrix(data5)
chordDiagram(phylum2)

colgroup <- distinctColorPalette(length(colnames(data5)))
collink <- distinctColorPalette(length(rownames(data5)))
grid.col <- c(colgroup,collink)
group.col <- collink
#grid.col <- distinctColorPalette(sum(length(rownames(data5)),length(colnames(data5))))
#group.col <- c(grid.col[(length(colnames(data5)) + 1):length(grid.col)])
chordDiagram(phylum2, grid.col = grid.col, 
             column.col = group.col, directional = -1)##directional=1或-1 可在vector与links之间再加一层扇形

##cirzlize P2
phylum2 <- as.matrix(data5)

df <- data.frame(from = rep(rownames(phylum2), ncol(phylum2)),
                 to = rep(colnames(phylum2), each = nrow(phylum2)),
                 value = as.vector(phylum2))

##颜色设定
color <- NULL
color[colnames(phylum2)] <- distinctColorPalette(length(colnames(data5)))
color[rownames(phylum2)] <-  distinctColorPalette(length(rownames(data5)))

chordDiagram(df, 
             grid.col =color,#颜色设置
             grid.border=NULL,#边框颜色设置，设置为NULL则默认与填充色一致
             transparency = 0.2,#连接颜色透明度
             link.lwd = 0.01,#线条宽度
             link.lty = 1,    # 线路类型
             link.border = 0,#边框颜色
             directional = -1,#表示线条的方向，0代表没有方向，1代表正向，-1代表反向，2代表双向
             diffHeight = mm_h(2),#外圈和中间连线的间隔
             direction.type = c("diffHeight","arrows"), #线条是否带有箭头
             link.arr.type = "big.arrow")#箭头类型

legend(x=1.2,
       y=1,###位置 
       pch=20,
       title = "Phylum",
       legend=rownames(phylum2),
       col=color[rownames(phylum2)],
       cex=1,
       pt.cex=3,
       border="black",
       xpd = T) 




#Species
color_compaired <- distinctColorPalette(100)
show_col(color_compaired)
t1 <- trans_abund$new(dataset = dataset, taxrank = "Species", ntaxa = 99, groupmean = "Group")
t1$plot_bar(others_color = "grey70", color_values = c(color_compaired),legend_text_italic = FALSE)

#######species circlize#########
data <- trans_abund$new(dataset = dataset, taxrank = "Species",groupmean = "Group", ntaxa = 100)
data$data_abund
data2 <- as.data.frame(data$data_abund)
data3 <- data2[,-c(3,5,6,7)]

# transform the data into circlize input format
data4 <- dcast(data3,Taxonomy~Sample) #data6 <- spread(data3,key = Sample,value = Abundance)
data4 <- as.data.frame(data4)
rownames(data4) <- data4$Taxonomy
data4 <- data4[,-1]

#filter top 100 in species level
data4$totalabundance <- rowSums(data4)
data5 <- top_n(data4,10,abs(totalabundance))
data5 <- data5[,-4]
##cirzlize P1  
species <- as.matrix(data5)
chordDiagram(species)

colgroup <- distinctColorPalette(length(colnames(species)))
collink <- distinctColorPalette(length(rownames(species)))
grid.col <- c(colgroup,collink)
group.col <- collink
#grid.col <- distinctColorPalette(sum(length(rownames(data5)),length(colnames(data5))))
#group.col <- c(grid.col[(length(colnames(data5)) + 1):length(grid.col)])
chordDiagram(species, grid.col = grid.col, 
             column.col = group.col, directional = -1)##directional=1或-1 可在vector与links之间再加一层扇形

##cirzlize P2
species <- as.matrix(data5)

df <- data.frame(from = rep(rownames(species), ncol(species)),
                 to = rep(colnames(species), each = nrow(species)),
                 value = as.vector(species))

##颜色设定
color <- NULL
color[colnames(species)] <- distinctColorPalette(length(colnames(data5)))
color[rownames(species)] <-  distinctColorPalette(length(rownames(data5)))

chordDiagram(df, 
             grid.col =color,#颜色设置
             grid.border=NULL,#边框颜色设置，设置为NULL则默认与填充色一致
             transparency = 0.2,#连接颜色透明度
             link.lwd = 0.01,#线条宽度
             link.lty = 1,    # 线路类型
             link.border = 0,#边框颜色
             directional = -1,#表示线条的方向，0代表没有方向，1代表正向，-1代表反向，2代表双向
             diffHeight = mm_h(2),#外圈和中间连线的间隔
             direction.type = c("diffHeight","arrows"), #线条是否带有箭头
             link.arr.type = "big.arrow")#箭头类型

legend(x=1.2,
       y=1,###位置 
       pch=20,
       legend=rownames(species),
       col=color[rownames(species)],
       cex=1,
       pt.cex=3,
       border="black",
       xpd = T) 

######class#####
color_compaired3 <- distinctColorPalette(10)
t1 <- trans_abund$new(dataset = dataset, taxrank = "Class", ntaxa = 30)
t1$plot_box(group = "Group",color_values = c(color_compaired3))


######genus#####
color_compaired3 <- distinctColorPalette(60)

t1 <- trans_abund$new(dataset = dataset, taxrank = "Genus", ntaxa = 100,groupmean = "Group")
t1$plot_bar(legend_text_italic = FALSE, color_values = c(color_compaired3),xtext_size = 6, others_color = "grey70", use_alluvium = F)



#######venn#####
dataset1 <- dataset$merge_samples(use_group = "Group")
t1 <- trans_venn$new(dataset1, ratio = "seqratio")
t1$plot_venn()
data <- print(t1[["data_details"]][["NC31"]])
data
write.csv(data,"NC31_venn.csv") 


#Lefse
t1 <- trans_diff$new(dataset = dataset, method = "lefse", group = "Group", alpha = 0.01, lefse_subgroup = NULL, p_adjust_method = "none", taxa_level = "Genus")
# t1$res_lefse is the LEfSe result
# t1$res_abund is the abundance information
t1$res_lefse[1:5, ]
t1$plot_diff_abund(use_number = 1:30)
library(ggtree)
t1$plot_lefse_cladogram(use_taxa_num = 200, use_feature_num = 50, clade_label_level = 5)


#RandomForest
t1 <- trans_diff$new(dataset = dataset, method = "rf", group = "Group", taxa_level = "Genus",alpha = 0.01,p_adjust_method = "none")
g1 <- t1$plot_diff_bar(use_number = 1:30)
g2 <- t1$plot_diff_abund(select_taxa = t1$plot_diff_bar_taxa)
g1 <- g1 + theme(legend.position = "none")
g2 <- g2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
p1 <- gridExtra::grid.arrange(g1, g2, ncol = 2, nrow = 1, widths = c(2, 1.7))


#Multi-group compaired anova
colorgroup <- c("blue3","red2","green2","purple2","orange2")
t1 <- trans_diff$new(dataset = dataset, method = "anova", group = "Group", taxa_level = "Species", filter_thres = 0.001)
t1$plot_diff_abund(use_number = 1:15, add_sig = T, coord_flip = F, color_values = c(colorgroup))


#group-group compaired
#plot for supprotinformation
t1 <- trans_diff$new(dataset = dataset, method = "metastat", group = "Group", taxa_level = "Species")
t1$plot_diff_abund(use_number = 1:20, select_group = "HF - NC31", coord_flip = F, color_values = c("blue3","orange2"))
t1$plot_diff_abund(use_number = 1:20, select_group = "HF - NY15", coord_flip = F, color_values = c("blue3","red2"))
t1$plot_diff_abund(use_number = 1:20, select_group = "HF - HeNa1610", coord_flip = F, color_values = c("blue3","green2"))
t1$plot_diff_abund(use_number = 1:20, select_group = "HF - FXJWS23M45", coord_flip = F, color_values = c("blue3","purple2"))
