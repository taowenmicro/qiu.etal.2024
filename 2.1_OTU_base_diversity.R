# #清空内存#######
# rm(list=ls()) 

# source("./coding/OTU_parameter.R")

# #-分析共四个部分，这是第二个部分，代号2，第一部分是扩增子原始序列处理


#- 0.1 基本表格保存#----------
tabpath = paste(otupath,"/report_table/",sep = "")
dir.create(tabpath)
#--raw otu tab
otu = as.data.frame(t(ggClusterNet::vegan_otu(ps)))
head(otu)
FileName <- paste(tabpath,"/otutab.raw.csv", sep = "")
write.csv(otu,FileName,sep = "")
# tax table
tax = as.data.frame((ggClusterNet::vegan_tax(ps)))
head(tax)``
FileName <- paste(tabpath,"/tax.raw.csv", sep = "")
write.csv(otu,FileName,sep = "")

ps_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) );ps_rela 
#--norm otu tab
otu_norm = as.data.frame(t(ggClusterNet::vegan_otu(ps_rela)))
FileName <- paste(tabpath,"/otutab.norm.rela.csv", sep = "")
write.csv(otu_norm,FileName,sep = "")

otutax <- cbind(as.data.frame(t(ggClusterNet::vegan_otu(ps_rela))),as.data.frame((ggClusterNet::vegan_tax(ps_rela))))
FileName <- paste(tabpath,"/otu.tax.norm.rela.csv", sep = "")
write.csv(otutax,FileName,sep = "")


for (i in 2: length(phyloseq::rank_names(ps))) {
  psi  <- ggClusterNet::tax_glom_wt(ps = ps,ranks = phyloseq::rank_names(ps)[i])
  #--raw otu tab
  otu = as.data.frame(t(ggClusterNet::vegan_otu(psi)))
  FileName <- paste(tabpath,"/otutab",phyloseq::rank_names(ps)[i],".csv", sep = "")
  write.csv(otu,FileName,sep = "")
  # tax table
  tax = as.data.frame((ggClusterNet::vegan_tax(ps)))
  FileName <- paste(tabpath,"/tax",phyloseq::rank_names(ps)[i],".csv", sep = "")
  write.csv(otu,FileName,sep = "")
  
  psi_rela  = phyloseq::transform_sample_counts(psi, function(x) x / sum(x) );psi_rela 
  #--norm otu tab
  otu_norm = as.data.frame(t(ggClusterNet::vegan_otu(psi_rela)))
  FileName <- paste(tabpath,"/otutab.norm.rela",phyloseq::rank_names(psi)[i],".csv", sep = "")
  write.csv(otu_norm,FileName,sep = "")
  
  otutax <- cbind(as.data.frame(t(ggClusterNet::vegan_otu(psi_rela))),as.data.frame((ggClusterNet::vegan_tax(psi_rela))))
  FileName <- paste(tabpath,"/otutax.norm.rela",phyloseq::rank_names(ps)[i],".csv", sep = "")
  write.csv(otutax,FileName,sep = "")
}

# out.ps.data(ps,
#             path = tabpath,
#             mark = "16S"
#             )


#如果导入数据是测序原始count，可以选择后面这几种标准化方法# sampling, log,TMM,RLE,upperquartile
# 下面是运行案例，直接降otu表格导出即可。
# ps1  = ps %>%
#   ggClusterNet::scale_micro(method = "TMM" )
# otu = as.data.frame(t(ggClusterNet::vegan_otu(ps1)))
# head(otu)


#--1.1 alpha多样性#---------
alppath = paste(otupath,"/alpha/",sep = "")
dir.create(alppath)


# source("../micro/alpha-diversity.R")

#---多种指标alpha多样性分析加出图-标记显著性
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/alpha-diversity.R")
index = c("Shannon","Inv_Simpson","Pielou_evenness","Simpson_evenness" ,"Richness" ,"Chao1","ACE" )

#--多种组合alpha分析和差异分析出图
alp = alpha(ps = ps,inde="Shannon",group = "Group",Plot = TRUE )
index= alp
head(index)

#--提取三个代表指标作图
all.alpha = TRUE

if (all.alpha) {
  sel = c(match("Inv_Simpson",colnames(index)),
          match("Pielou_evenness",colnames(index)),
          match("Simpson_evenness",colnames(index)),
          match("Richness",colnames(index)),
          match("Chao1",colnames(index)),
          match("ACE",colnames(index)),
          match("Shannon",colnames(index))
          
  )
  h = 3
} else{
  sel = c(match("Shannon",colnames(index)),match("Richness",colnames(index)),
          match("Pielou_evenness",colnames(index)))
  h = 1
}

n = length(sel) + 3


data = cbind(data.frame(ID = 1:length(index$Group),group = index$Group),index[sel])
head(data)


result = EasyStat::MuiKwWlx2(data = data,num = c(3:(n -1)))

FileName <- paste(alppath,"/alpha_diversity_different_label.csv", sep = "")
write.csv(result,FileName,sep = "")
FileName <- paste(alppath,"/alpha_diversity_index.csv", sep = "")
write.csv(index,FileName,sep = "")


result1 = EasyStat::FacetMuiPlotresultBox(data = data,num = c(3:(n -1)),
                                result = result,
                                sig_show ="abc",ncol = 1 )
p1_1 = result1[[1]] + 
  ggplot2::scale_x_discrete(limits = axis_order) + 
  mytheme2 +
  ggplot2::guides(fill = guide_legend(title = NULL)) +
  ggplot2::scale_fill_manual(values = colset1)
p1_1

#如何升级展示-提取数据用小提琴图展示
p1_0 = result1[[2]] %>% ggplot(aes(x=group , y=dd )) + 
  geom_violin(alpha=1, aes(fill=group)) +
  geom_jitter( aes(color = group),position=position_jitter(0.17), size=3, alpha=0.5)+
  labs(x="", y="")+
  facet_wrap(.~name,scales="free_y",ncol  = 3) +
  # theme_classic()+
  geom_text(aes(x=group , y=y ,label=stat)) +
  ggplot2::scale_x_discrete(limits = axis_order) + 
  mytheme1 +
  guides(color=guide_legend(title = NULL),
         shape=guide_legend(title = NULL),
         fill = guide_legend(title = NULL)
         ) +
  ggplot2::scale_fill_manual(values = colset1)
p1_0


res = EasyStat::FacetMuiPlotresultBar(data = data,num = c(3:(n -1)),result = result,sig_show ="abc",ncol = 3)
p1_2 = res[[1]]+ scale_x_discrete(limits = axis_order) + guides(color = FALSE) +
  mytheme1+ 
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p1_2

res = EasyStat::FacetMuiPlotReBoxBar(data = data,num = c(3:(n -1)),result = result,sig_show ="abc",ncol = 3)
p1_3 = res[[1]]+ scale_x_discrete(limits = axis_order) + 
  mytheme1 + 
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p1_3


FileName <- paste(alppath,"Alpha_Facet_vl", ".pdf", sep = "")
ggsave(FileName, p1_0, width = ((1 + gnum) * 3), height =4 *h,limitsize = FALSE)
FileName <- paste(alppath,"Alpha_Facet_vl", ".png", sep = "")
ggsave(FileName, p1_0, width = ((1 + gnum) * 3), height =4 *h,limitsize = FALSE)


FileName <- paste(alppath,"Alpha_Facet_box", ".pdf", sep = "")
ggsave(FileName, p1_1, width = 20, height =80,limitsize = FALSE)

FileName <- paste(alppath,"Alpha_Facet_bar", ".pdf", sep = "")
ggsave(FileName, p1_2, width = ((1 + gnum) * 3), height = 4 *h,limitsize = FALSE)

FileName <- paste(alppath,"Alpha_Facet_boxbar", ".pdf", sep = "")
ggsave(FileName, p1_3, width = ((1 + gnum) * 3), height = 4 *h,limitsize = FALSE)

FileName <- paste(alppath,"Alpha_Facet_box", ".jpg", sep = "")
ggsave(FileName, p1_1, width = ((1 + gnum) * 3), height =4*h,limitsize = FALSE)

FileName <- paste(alppath,"Alpha_Facet_bar", ".jpg", sep = "")
ggsave(FileName, p1_2, width = ((1 + gnum) * 3), height = 4*h,limitsize = FALSE)

FileName <- paste(alppath,"Alpha_Facet_boxbar", ".jpg", sep = "")
ggsave(FileName, p1_3, width = ((1 + gnum) * 3), height = 4*h,limitsize = FALSE)


#--总体差异检测alpha多样性
krusk1 = ggpubr::compare_means( Shannon ~ group, data=data, method = "kruskal.test")
krusk2 = ggpubr::compare_means( Richness ~ group, data=data, method = "kruskal.test")
krusk3 = ggpubr::compare_means( Pielou_evenness ~ group, data=data, method = "kruskal.test")

dat = rbind(krusk1,krusk2,krusk3) %>% as.data.frame()
FileName <- paste(alppath,"/alpha_diversity_index_all_p_Kruskal-Wallis.csv", sep = "")
write_csv(dat,FileName)

#--1.2 系统发育多样性#------
alppath = paste(otupath,"/alpha/",sep = "")
dir.create(alppath)
com_2020 <- ps %>% vegan_otu() %>% 
  as.data.frame()
library(ape)
library(picante) 

## 加载系统发育分析包
rooted <- phy_tree(ps)  
## 读取系统发育树文件
# rooted <- root(phylo,outgroup = "Amborella_trichopoda",resolve.root = T)  ## 转化为有根树
## 转化为距离矩阵

cover2020.pd<-pd(com_2020,rooted,include.root=F) 

map = sample_data(ps)
head(map)

data = cbind(map[,c("ID","Group")],pd = cover2020.pd[,1])
head(data)
colnames(data)[2] = "group"
data$group = as.factor(data$group)
result = EasyStat::MuiKwWlx2(data = data,num = c(3))

FileName <- paste(alppath,"/pd_diversity_different_label.csv", sep = "")
write_csv(result,FileName)
FileName <- paste(alppath,"/pd_diversity_index.csv", sep = "")
write_csv(data,FileName)

result1 = EasyStat::FacetMuiPlotresultBox(data = data,num = c(3),
                                          result = result,
                                          sig_show ="abc",ncol = 1 )
p1_1 = result1[[1]] + 
  ggplot2::scale_x_discrete(limits = axis_order) + 
  mytheme2 +
  ggplot2::guides(fill = guide_legend(title = NULL)) +
  ggplot2::scale_fill_manual(values = colset1)
p1_1

#如何升级展示-提取数据用小提琴图展示
p1_1 = result1[[2]] %>% ggplot(aes(x=group , y=dd )) + 
  geom_violin(alpha=1, aes(fill=group)) +
  geom_jitter( aes(color = group),position=position_jitter(0.17), size=3, alpha=0.5)+
  labs(x="", y="")+
  facet_wrap(.~name,scales="free_y",ncol  = 3) +
  # theme_classic()+
  geom_text(aes(x=group , y=y ,label=stat)) +
  ggplot2::scale_x_discrete(limits = axis_order) + 
  mytheme2 +
  guides(color=guide_legend(title = NULL),
         shape=guide_legend(title = NULL),
         fill = guide_legend(title = NULL)
  ) +
  ggplot2::scale_fill_manual(values = colset1)
p1_1

res = EasyStat::FacetMuiPlotresultBar(data = data,num = c(3),result = result,sig_show ="abc",ncol = 1)
p1_2 = res[[1]]+ scale_x_discrete(limits = axis_order) + guides(color = FALSE) +
  mytheme2+ 
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p1_2

res = EasyStat::FacetMuiPlotReBoxBar(data = data,num = c(3),result = result,sig_show ="abc",ncol =1)
p1_3 = res[[1]]+ scale_x_discrete(limits = axis_order) + 
  mytheme2 + 
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p1_3


FileName <- paste(alppath,"pd_Facet_box", ".pdf", sep = "")
ggsave(FileName, p1_1, width = ((1 + gnum) ), height =4,limitsize = FALSE)

FileName <- paste(alppath,"pd_Facet_bar", ".pdf", sep = "")
ggsave(FileName, p1_2, width = ((1 + gnum) ), height = 4,limitsize = FALSE)

FileName <- paste(alppath,"pd_Facet_boxbar", ".pdf", sep = "")
ggsave(FileName, p1_3, width = ((1 + gnum) ), height = 4,limitsize = FALSE)

FileName <- paste(alppath,"pd_Facet_box", ".jpg", sep = "")
ggsave(FileName, p1_1, width = ((1 + gnum) ), height =4,limitsize = FALSE)

FileName <- paste(alppath,"pd_Facet_bar", ".jpg", sep = "")
ggsave(FileName, p1_2, width = ((1 + gnum) ), height = 4,limitsize = FALSE)

FileName <- paste(alppath,"pd_Facet_boxbar", ".jpg", sep = "")
ggsave(FileName, p1_3, width = ((1 + gnum) ), height = 4,limitsize = FALSE)


#--1.3 alpha稀释曲线绘制#-----

rare <- mean(phyloseq::sample_sums(ps))/10

source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\alpha_rare_all.R",encoding = "utf-8")
result = alpha_rare_all(ps = ps, group = "Group", method = "Richness", start = 100, step = rare)
#--提供单个样本溪稀释曲线的绘制
p2_1 <- result[[1]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)

## 提供数据表格，方便输出
raretab <- result[[2]]
head(raretab)


#--按照分组展示稀释曲线
p2_2 <- result[[3]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
#--按照分组绘制标准差稀释曲线
p2_3 <- result[[4]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)

# p2_4 <- result[[5]] +
#   mytheme1 +
#   guides(fill = guide_legend(title = NULL))+
#   scale_fill_manual(values = colset1)


FileName <- paste(alppath,"Alpha_rare_sample", ".pdf", sep = "")
ggsave(FileName, p2_1, width = 8, height =6)

FileName <- paste(alppath,"Alpha_rare_group", ".pdf", sep = "")
ggsave(FileName, p2_2, width = 8, height =6)
FileName <- paste(alppath,"Alpha_rare_groupwithSD", ".pdf", sep = "")
ggsave(FileName, p2_3, width = 8, height =6)
FileName <- paste(alppath,"Alpha_rare_sample", ".jpg", sep = "")
ggsave(FileName, p2_1, width = 8, height =6)
FileName <- paste(alppath,"Alpha_rare_group", ".jpg", sep = "")
ggsave(FileName, p2_2, width = 8, height =6)
FileName <- paste(alppath,"Alpha_rare_groupwithSD", ".jpg", sep = "")
ggsave(FileName, p2_3, width = 8, height =6)


FileName <- paste(alppath,"/Alpha_rare_data.csv", sep = "")
write.csv(raretab,FileName,sep = "")

#---2.1 排序分析beta-diversity#----
betapath = paste(otupath,"/beta/",sep = "")
dir.create(betapath)

# "unifrac" "wunifrac" "dpcoa" "jsd" "manhattan" "euclidean"   "canberra" "bray" "kulczynski" 
# "jaccard" "gower" "altGower" "morisita" "horn" "mountford"  "raup" "binomial" 
# "chao"  "cao" "w"  "-1"  "c" "wb"  "r"   "I"  "e" "t" "me"   "j"  "sor"  "m"   "-2"  "co"
# DCA, CCA, RDA, NMDS, MDS, PCoA, PCA, LDA tsne 

source("E:\\Shared_Folder\\Function_local\\R_function\\micro/BetaDiv.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/MicroTest.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/pairMicroTest.R")

methodlist = c("NMDS","PCoA", "PCA")

# methodlist = c("LDA")
for (method in methodlist) {
  result = BetaDiv(ps = ps, group = "Group", dist = "bray",
                   method = method, Micromet = "anosim", pvalue.cutoff = 0.05,
                   pair = F)
  p3_1 = result[[1]] + 
    scale_fill_manual(values = colset1)+
    scale_color_manual(values = colset1,guide = F) +
    mytheme1 
    # theme(legend.position = c(0.2,0.2))
  p3_1
  #带标签图形出图
  p3_2 = result[[3]] +
    scale_fill_manual(values = colset1)+
    scale_color_manual(values = colset1,guide = F) + 
    mytheme1 
    # theme(legend.position = c(0.2,0.2))
  p3_2
  
  FileName <- paste(betapath,"/a2_",method,"bray.pdf", sep = "")
  ggsave(FileName, p3_1, width = 8, height = 7)
  FileName1 <- paste(betapath,"/a2_",method,"",method,"bray.jpg", sep = "")
  ggsave(FileName1 , p3_1, width = 12, height = 11)
  
  FileName <- paste(betapath,"/a2_",method,"bray_label.pdf", sep = "")
  ggsave(FileName, p3_2, width = 12, height = 12)
  FileName1 <- paste(betapath,"/a2_",method,"bray_label.jpg", sep = "")
  ggsave(FileName1 , p3_2, width = 12, height = 11)
  
  # 提取出图数据
  plotdata = result[[2]]
  FileName <-  paste(betapath,"/a2_",method,"bray.csv", sep = "")
  write.csv(plotdata,FileName)
  #---------排序-精修图
  plotdata =result[[2]]
  head(plotdata)
  # 求均值
  cent <- aggregate(cbind(x,y) ~Group, data = plotdata, FUN = mean)
  cent
  # 合并到样本坐标数据中
  segs <- merge(plotdata, setNames(cent, c('Group','oNMDS1','oNMDS2')),
                by = 'Group', sort = FALSE)
  
  # p2$layers[[2]] = NULL
  # library(ggcor)
  library(ggsci)
  p3_3 = p3_1 +geom_segment(data = segs,
                            mapping = aes(xend = oNMDS1, yend = oNMDS2,color = Group),show.legend=F) + # spiders
    geom_point(mapping = aes(x = x, y = y),data = cent, size = 5,pch = 24,color = "black",fill = "yellow") +
    scale_fill_manual(values = colset1)+
    scale_color_manual(values = colset1,guide = F) + 
    mytheme1 
    # theme(legend.position = c(0.2,0.2))
  p3_3
  
  FileName <- paste(betapath,"/a2_",method,"bray_star.pdf", sep = "")
  ggsave(FileName, p3_3, width = 8, height = 7)
  FileName1 <- paste(betapath,"/a2_",method,"bray_star.jpg", sep = "")
  ggsave(FileName1 , p3_3, width = 8, height = 7)
  
}

map

#提取总体比较
TResult =result[[5]]
head(TResult)

# 提取两两检测结果
pair = result[[4]]
pair
FileName <- paste(betapath,"Pair_anosim.csv", sep = "")
write.csv(pair,FileName)
FileName <- paste(betapath,"Total_anosim.csv", sep = "")
write.csv(TResult,FileName)

#--换用adonis差异分析
title1 = MicroTest(ps = ps, Micromet = "adonis", dist = "bray")
title1
FileName <- paste(betapath,"Total_adonis.csv", sep = "")
write.csv(title1,FileName)
# pairResult = pairMicroTest(ps = ps, Micromet = "adonis")
# FileName <- paste(betapath,"Pair_anosim.csv", sep = "")
# write.csv(pair,FileName)

##多因素多元置换检验#--------
# ps1_rela  = phyloseq::transform_sample_counts(ps, function(x) x / sum(x) );ps1_rela
# # library(vegan)
# map = as.data.frame(phyloseq::sample_data(ps1_rela))
# unif = phyloseq::distance(ps, method= "bray")
# ado =  vegan::adonis2(unif ~ map$Group *map$species,permutations = 999)
# sink(paste(betapath,"/Mui.adonis.txt",sep = ""))
# print(ado)
# sink()


#---2.2 普氏分析#------
map = phyloseq::sample_data(ps)
samegroup = map$Group %>% table() %>% as.data.frame() %>% .$Freq %>% unique() %>% length() == 1

if (samegroup) {
  source("E:\\Shared_Folder\\Function_local\\R_function\\micro/matel_pro_plot.R")
  #--门特尔检验-普氏分析
  if (length(unique(phyloseq::sample_data(ps)$Group)) < 6) {
    # library(vegan)
    size = combn(unique(phyloseq::sample_data(ps)$Group),2) %>% dim()
    size
    result <- mantal.micro(ps = ps,
                           method =  "spearman",
                           group = "Group",
                           ncol = size[2],
                           nrow = 1
    )
    data <- result[[1]]
    
    p3_7 <- result[[2]] +  mytheme1 
    p3_7
    
    FileName <- paste(betapath,"mantel_pro.csv", sep = "")
    write.csv(data,FileName)
    FileName1 <- paste(betapath,"/a2_","Mantel_Pro.pdf", sep = "")
    ggsave(FileName1 , p3_7, width = size[2] *6, height = 6,limitsize = FALSE)
    FileName1 <- paste(betapath,"/a2_","Mantel_Pro.jpg", sep = "")
    ggsave(FileName1 , p3_7, width = size[2] *6, height =6,limitsize = FALSE)
  }
}


#--3.1 物种分类树分析#----
# detach("package:treeio")
# detach("package:vegan")
library(ggstar)
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/phy_tree_micro.R")

barpath = paste(otupath,"/phy_tree_micro/",sep = "");print(barpath)
dir.create(barpath)
library(ggClusterNet)
# rank.names(ps)
# Top_micro = 150

result <- phy_tree_micro(ps = ps,Top = 150)

p0 = result[[1]]
p1 = result[[2]]
p2 = result[[3]]
p3 = result[[4]]
p4 = result[[5]]
p5 = result[[6]]
p6 = result[[7]]
p7 = result[[8]]

detach("package:ggstar")

FileName <- paste(barpath,Top_micro,"phy_tree_micro1", ".pdf", sep = "")
ggsave(FileName, p0, width = 4, height = 4)
FileName <- paste(barpath,Top_micro,"phy_tree_micro2", ".pdf", sep = "")
ggsave(FileName, p1, width = 6, height =4)
FileName <- paste(barpath,Top_micro,"phy_tree_micro3", ".pdf", sep = "")
ggsave(FileName, p2, width = 6, height = 4)

FileName <- paste(barpath,Top_micro,"phy_tree_micro4", ".pdf", sep = "")
ggsave(FileName, p3, width = 12, height = 10)
FileName <- paste(barpath,Top_micro,"phy_tree_micro5", ".pdf", sep = "")
ggsave(FileName, p4, width = 15, height = 13)
FileName <- paste(barpath,Top_micro,"phy_tree_micro5_2", ".pdf", sep = "")
ggsave(FileName, p5, width = 18, height = 15)
FileName <- paste(barpath,Top_micro,"phy_tree_micro6", ".pdf", sep = "")
ggsave(FileName, p6, width = 5, height = 4)
FileName <- paste(barpath,Top_micro,"phy_tree_micro7", ".pdf", sep = "")
ggsave(FileName, p7, width = 6, height = 5)

# library(cowplot)
# save_plot(FileName, p2, base_height = 7, base_width =7)

FileName <- paste(barpath,Top_micro,"phy_tree_micro1", ".png", sep = "")
ggsave(FileName, p0, width = 4, height = 4)
FileName <- paste(barpath,Top_micro,"phy_tree_micro2", ".png", sep = "")
ggsave(FileName, p1, width = 6, height =4)
FileName <- paste(barpath,Top_micro,"phy_tree_micro3", ".png", sep = "")
ggsave(FileName, p2, width = 7, height = 7,dpi = 72)
FileName <- paste(barpath,Top_micro,"phy_tree_micro4", ".png", sep = "")
ggsave(FileName, p3, width = 12, height = 12,dpi = 72)
FileName <- paste(barpath,Top_micro,"phy_tree_micro5", ".png", sep = "")
ggsave(FileName, p4, width = 15, height = 15,dpi = 72)
FileName <- paste(barpath,Top_micro,"phy_tree_micro5_2", ".png", sep = "")
ggsave(FileName, p5, width = 18, height = 15,dpi = 72)
FileName <- paste(barpath,Top_micro,"phy_tree_micro6", ".png", sep = "")
ggsave(FileName, p6, width = 5, height = 4,dpi = 72)
FileName <- paste(barpath,Top_micro,"phy_tree_micro7", ".png", sep = "")
ggsave(FileName, p7,width = 6, height = 5,dpi = 72)


#--3.2 和弦图#-----
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/cir_plot.R")
# otupath = paste(res1path,"/OTU/",sep = "");otupath
# dir.create(otupath)

barpath =  paste(otupath,"/circle_plot/",sep = "");otupath
dir.create(barpath)

for (i in 2:7) {
  cir_plot(ps  = ps,Top = 12,rank = i,
           path = barpath)
}



#--3.3 物种组成展示#---------
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/barMainplot.R")
barpath = paste(otupath,"/Microbial_composition2/",sep = "")
dir.create(barpath)

phyloseq::rank_names(ps)
pst = ps %>% subset_taxa.wt("Species","Unassigned",TRUE)
pst = pst %>% subset_taxa.wt("Genus","Unassigned",TRUE)

for (j in c("Phylum" , "Class" ,  "Order"  , "Family" , "Genus","Species")) {
  result = barMainplot(ps = pst,
                       j = j,
                       # axis_ord = axis_order,
                       label = FALSE,
                       sd = FALSE,
                       Top = 35)
  p4_1 <- result[[1]] + 
    # scale_fill_brewer(palette = "Paired") + 
    scale_fill_manual(values = colset2) +
    scale_x_discrete(limits = axis_order) +
    mytheme1
  p4_1

  p4_2  <- result[[3]] + 
    # scale_fill_brewer(palette = "Paired") + 
    scale_fill_manual(values = colset2) +
    scale_x_discrete(limits = axis_order) + 
    mytheme1
  p4_2
  
  databar <- result[[2]] %>% group_by(Group,aa) %>%
    dplyr::summarise(sum(Abundance)) %>% as.data.frame()
  head(databar)
  colnames(databar) = c("Group",j,"Abundance(%)")
  
  
  FileName1 <- paste(barpath,"/a2_",j,"_barflow",".pdf", sep = "")
  ggsave(FileName1, p4_2, width = (5+ gnum), height =18,limitsize = FALSE)
  FileName2 <- paste(barpath,"/a2_",j,"_barflow",".jpg", sep = "")
  ggsave(FileName2, p4_2, width = (5+ gnum), height =18,limitsize = FALSE)
  
  FileName1 <- paste(barpath,"/a2_",j,"_bar",".pdf", sep = "")
  ggsave(FileName1, p4_1, width = (5+ gnum), height =18 ,limitsize = FALSE)
  FileName2 <- paste(barpath,"/a2_",j,"_bar",".jpg", sep = "")
  ggsave(FileName2, p4_1, width = (5+ gnum), height =18,limitsize = FALSE)
  
  FileName <- paste(barpath,"/a2_",j,"_bar_data",".csv", sep = "")
  write.csv(databar,FileName,quote = F)
}

detach("package:ggalluvial")



#--3.4 微生物桑基图#--------
library(ggClusterNet)
snapath =  paste(otupath,"/sankeyNetwork/",sep = "");otupath
dir.create(snapath)
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\sankey.m.Group.R",encoding = "utf-8")



result = sankey.m.Group(
  ps = ps  %>% subset_taxa.wt("Species","Unassigned",TRUE),
  rank = 6,
  Top = 50

)

p = result[[1]]
p
dat = result[[2]]

FileName <-paste(snapath,"/sankey_Group.csv", sep = "")
write.csv(dat,FileName,sep = "")

saveNetwork(p,paste(snapath,"/sankey_Group.html", sep = ""))
library(webshot)
# webshot::install_phantomjs()
# webshot(paste(snapath,"/sankey1.html", sep = "") ,paste(snapath,"/sankey1.png", sep = ""))
webshot(paste(snapath,"/sankey_Group.html", sep = "") , paste(snapath,"/sankey_Group.pdf", sep = ""))


source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\sankey.micro.R",encoding = "utf-8")

sankey.micro (ps = ps,
              rank = 6,# 参数目前不可修改
              Top = 50,
              snapath = snapath
              
              )


#--3.5 环状物种堆叠柱状图#----
barpath = paste(otupath,"/circle_Micro_strack_bar/",sep = "");print(barpath)
dir.create(barpath)

source("E:\\Shared_Folder\\Function_local\\R_function\\micro/circle_stark_bar.R")
library(ggtree)
# j = "Phylum"
p2 = circle_starc_bar(
  ps = ps,
  Top = 15,
  dist = "bray",
  cuttree = 3,
  hcluter_method = "complete")

FileName2 <- paste(barpath,"/a2_","_bar",".jpg", sep = "")
ggsave(FileName2, p2, width = 10, height =8 )

FileName2 <- paste(barpath,"/a2_","_bar",".pdf", sep = "")
ggsave(FileName2, p2, width = 10, height =8 )

#--3.6 距离和丰度合并#-------
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/cluMicro.bar.R")


for (j in c("Phylum" , "Class" ,  "Order"  , "Family" , "Genus")) {
  
  result <-  cluMicro.bar (dist = "bray",
                           ps = ps,
                           j = j,
                           Top = Top, # 提取丰度前十的物种注释
                           tran = TRUE, # 转化为相对丰度值
                           hcluter_method = "complete",
                           Group = "Group",
                           cuttree = length(unique(phyloseq::sample_data(ps)$Group))
  )
  
  p5_1 <- result[[1]]
  p5_2 <- result[[2]]
  p5_3 <- result[[3]]
  p5_4 <- result[[4]]
  clubardata <- result[[5]]
  
  
  FileName1 <- paste(barpath,"/a2_",j,"_cluster_sample",".pdf", sep = "")
  ggsave(FileName1, p5_1, width = 6, height = dim(phyloseq::sample_data(ps))[1]/4,limitsize = FALSE)
  FileName1 <- paste(barpath,"/a2_",j,"_cluster_bar_sample",".pdf", sep = "")
  ggsave(FileName1, p5_2, width = 12, height = dim(phyloseq::sample_data(ps))[1]/4 ,limitsize = FALSE)
  FileName1 <- paste(barpath,"/a2_",j,"_cluster_sample",".jpg", sep = "")
  ggsave(FileName1, p5_1, width = 6, height =dim(phyloseq::sample_data(ps))[1]/4 ,limitsize = FALSE)
  FileName1 <- paste(barpath,"/a2_",j,"_cluster_bar_sample",".jpg", sep = "")
  ggsave(FileName1, p5_2, width = 12, height =dim(phyloseq::sample_data(ps))[1]/4 ,limitsize = FALSE)
  
  FileName1 <- paste(barpath,"/a2_",j,"_cluster_Group",".pdf", sep = "")
  ggsave(FileName1, p5_3, width = 6, height = gnum,limitsize = FALSE)
  FileName1 <- paste(barpath,"/a2_",j,"_cluster_bar_Group",".pdf", sep = "")
  ggsave(FileName1, p5_4, width = 12, height = gnum ,limitsize = FALSE)
  FileName1 <- paste(barpath,"/a2_",j,"_cluster_Group",".jpg", sep = "")
  ggsave(FileName1, p5_3, width = 6, height = gnum ,limitsize = FALSE)
  FileName1 <- paste(barpath,"/a2_",j,"_cluster_bar_Group",".jpg", sep = "")
  ggsave(FileName1, p5_4, width = 12, height = gnum ,limitsize = FALSE)
  
  FileName <- paste(barpath,"/a2_",j,"_cluster_bar_data",".csv", sep = "")
  write.csv(clubardata,FileName)
  
}



#--3.7 三元图#--------
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/Micro_tern.R")
library(ggtern)
# 设定大于三组以后才可以使用
if (gnum  >= 3) {
  ternpath = paste(otupath,"/ggtern/",sep = "")
  dir.create(ternpath)
  
  ps1 = ps %>% filter_OTU_ps(500)
  ps1
  Micro_tern(ps = ps1,ternpath = ternpath )
  detach("package:ggtern")
  
}



#---共有微生物特有微生物
#--3.8 flower plot#-----
flowpath = paste(otupath,"/flowplot/",sep = "")
dir.create(flowpath)

source("E:\\Shared_Folder\\Function_local\\R_function\\micro/ggflowerplot.R")
p0_1 <- ggflower(ps = ps ,
                 # rep = 1,
                 group = "ID",
                 start = 1, # 风车效果
                 m1 = 1, # 花瓣形状，方形到圆形到棱形，数值逐渐减少。
                 a = 0.2, # 花瓣胖瘦
                 b = 1, # 花瓣距离花心的距离
                 lab.leaf = 1, # 花瓣标签到圆心的距离
                 col.cir = "yellow",
                 N = 0.5
) 

p0_1 

# p + scale_fill_brewer(palette = "Paired")
FileName1 <- paste(flowpath,"ggflowerID.pdf", sep = "")
ggsave(FileName1, p0_1, width = 8, height = 8)
FileName2 <- paste(flowpath,"ggflowerID.jpg", sep = "")
ggsave(FileName2, p0_1, width = 8, height = 8 )



p0_2 <- ggflower(ps = ps,
                 # rep = 1,
                 group = "Group",
                 start = 1, # 风车效果
                 m1 = 1.8, # 花瓣形状，方形到圆形到棱形，数值逐渐减少
                 a = 0.3, # 花瓣胖瘦
                 b = 1, # 花瓣距离花心的距离
                 lab.leaf = 1, # 花瓣标签到圆心的距离
                 col.cir = "yellow",
                 N = 0.1
) + scale_fill_manual(values = colset1) 
p0_2

FileName1 <- paste(flowpath,"ggflowerGroup.pdf", sep = "")
ggsave(FileName1, p0_2, width = 14, height = 14)
FileName2 <- paste(flowpath,"ggflowerGroup.jpg", sep = "")
ggsave(FileName2, p0_2, width = 14, height = 14 )


#3.9 ggplot升级版本韦恩图和Upset#-------

source("E:\\Shared_Folder\\Function_local\\R_function\\micro/Ven.Upset.gg.R")

if (gnum < 6) {
  Venpath = paste(otupath,"/Ven_Upset_super/",sep = "")
  dir.create(Venpath)
  
  library(ggVennDiagram)
  res = Ven.Upset(ps =  ps,
                  group = "Group",
                  N = 0.5,
                  size = 3)
  
  p1 = res[[1]]
  p2 = res[[2]]
  
  filename3 <- paste(Venpath,"Ven_gg.pdf", sep = "")
  ggsave(filename3, p1, width = 8, height = 8)
  filename3 <- paste(Venpath,"Upset_gg.pdf", sep = "")
  ggsave(filename3, p2, width = 8, height = 8)
}


#-3.10 Ven-Upset#----------
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/VenSeper.R")
source("E:/Shared_Folder/Function_local/R_function/micro/barMainplot.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/Ven-Upset.R")
# j = "Genus"
group = "Group"
ps_Ven = ps
# BiocManager::install("VennDiagram")
# otutab = as.data.frame(otu_table(ps))

map = as.data.frame(phyloseq::sample_data(ps_Ven))
gnumven <- map[,group] %>% unique() %>% dim()

if (gnumven[1] < 6) {
  

  Venpath = paste(otupath,"/Ven_Upset_super/",sep = "")
  dir.create(Venpath)

  
  result = VenUpset(ps = ps_Ven,
                    group = group,
                    path = Venpath
  )
  
  #---每个部分
  result = VenSeper(ps = ps_Ven,
                   path = Venpath,
                   group = group,
                   j = j,
                   Top = 10
                   
  )
  # 提取韦恩图中全部部分的otu极其丰度做门类柱状图
  p7_1 <- result[[1]]
  #每个部分序列的数量占比，并作差异
  p8 <- result[[2]]
  # 每部分的otu门类冲积图
  p7_2 <- result[[3]]
  
  
  FileName <- paste(Venpath,j,"count_Facet_ven", ".pdf", sep = "")
  ggsave(FileName, p7_1, width = 15, height = 12)
  
  FileName <- paste(Venpath,j,"diff_count_box", ".pdf", sep = "")
  ggsave(FileName, p8, width = 15, height = 12)
  
  FileName <- paste(Venpath,j,"count_Facet_ven_flow", ".pdf", sep = "")
  ggsave(FileName, p7_2, width = 15, height = 12)
  
  FileName <- paste(Venpath,j,"count_Facet_ven", ".jpg", sep = "")
  ggsave(FileName, p7_1, width = 15, height = 12)
  
  FileName <- paste(Venpath,j,"diff_count_box", ".jpg", sep = "")
  ggsave(FileName, p8, width = 15, height = 12)
  
  FileName <- paste(Venpath,j,"count_Facet_ven_flow", ".jpg", sep = "")
  ggsave(FileName, p7_2, width = 15, height = 12)
  
}




#-3.11 维恩网络#-------
library(ggClusterNet)
library(phyloseq)
library(ggrepel)

biospath = paste(otupath,"/biospr_network_Ven/",sep = "")
dir.create(biospath)
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/ven.network.R")
result = ven.network(
  ps = ps,
  N = 0.5,
  fill = "Phylum"
)

p  = result[[1]]

data = result[[2]]

filename = paste(biospath,"/","biostr_Ven_network.species.several.pdf",sep = "")
ggsave(filename,p,width = (15),height = (12))
filename = paste(biospath,"/","biostr_Ven_network.jpg",sep = "")
ggsave(filename,p,width = (15),height = (12))

filename = paste(biospath,"Ven.network.all.csv",sep = "")
write.csv(data,filename)


detach("package:ggClusterNet")
detach("package:phyloseq")


#-4.1 差异分析edgeR#----
diffpath = paste(otupath,"/diff_tax/",sep = "")
dir.create(diffpath)

# 准备脚本
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\EdgerSuper.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\DESep2_micro.R")
# source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\Plot.CompareWithCK.R",encoding = "utf-8")

diffpath.2 = paste(diffpath,"/EDgeR/",sep = "")
dir.create(diffpath.2)

res = EdgerSuper(ps = ps,group  = "Group",artGroup = NULL,
                 j = "OTU",
                 path = diffpath.2
)
head(res)
filename = paste(diffpath.2,"/","_","OTU","_","edger_all.csv",sep = "")
write.csv(res,filename)

for (j in c(2:6)) {
  res = EdgerSuper(ps = ps,group  = "Group",artGroup = NULL,
                   j = j,
                   path = diffpath.2
  )
  head(res)
  filename = paste(diffpath.2,"/","_",j,"_","edger_all.csv",sep = "")
  write.csv(res,filename)
}

#-4.1 差异分析DEsep2#----

diffpath.1 = paste(diffpath,"/DEsep2/",sep = "")
dir.create(diffpath.1)
res = DESep2_Meta2(ps = ps,
                   group  = "Group",
                   artGroup =NULL,
                   j = "OTU",
                   path = diffpath.1
)

head(res)
tail(res)
filename = paste(diffpath.1,"/","_","OTU","_","DESep2_all.csv",sep = "")
write.csv(res,filename,quote = F)

for (j in 2:6) {
  res = DESep2_Meta2(ps = ps,
                     group  = "Group",
                     artGroup =NULL,
                     j = j,
                     path = diffpath.1
  )
  head(res)
  
  filename = paste(diffpath.1,"/","_",j,"_","DESep2_all.csv",sep = "")
  write.csv(res,filename,quote = F)
}



#-4.2 edger--曼哈顿图绘制#-------
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\edge_Manhattan.R")

diffpath = paste(otupath,"/diff_Manhattan/",sep = "")
dir.create(diffpath)

edge_Manhattan(
  ps = ps,
  pvalue = 0.05,
  lfc = 0,
  diffpath = diffpath 
)



# 4.3 -相对对比柱状图-更新完成#-----


source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\Plot.CompareWithCK.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\EdgerSuper.R")

library(ggClusterNet)
library(phyloseq)
library(parallel)
library(tidyverse)


j = "OTU" 
res = Plot.CompareWithCK(ps = ps,
                         CK = "WT",
                         j = j,
                         abun = 0.001,
                         cpu = 6)

p = res[[1]]
p
data = res[[2]]

filename = paste(diffpath,"/","edger_001_diff_bio_plot.pdf",sep = "")
ggsave(filename,p,width = 10,height = dim(data)[1]/8)
filename = paste(diffpath,"/","edger_001_diff_bio_plot.jpg",sep = "")
ggsave(filename,p,width = 10,height = dim(data)[1]/8)

filename = paste(diffpath,"/","edger_Top_2000_plotdata.csv",sep = "")
write.csv(data,filename)


j = "Genus"
res = Plot.CompareWithCK(ps = ps,
                         CK = "WT",
                         j = j,
                         abun = 0,
                         cpu = 6)

p = res[[1]]
p
data = res[[2]]

filename = paste(diffpath,"/",j,"edger_diff_bio_plot.pdf",sep = "")
ggsave(filename,p,width = 10,height = dim(data)[1]/8)
filename = paste(diffpath,"/",j,"edger_001_diff_bio_plot.jpg",sep = "")
ggsave(filename,p,width = 10,height = dim(data)[1]/8)

filename = paste(diffpath,"/",j,"edger_Top_2000_plotdata.csv",sep = "")
write.csv(data,filename)


#-4.4 stemp_差异分析#-------
library(ggClusterNet)
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/stemp_diff.R")

# library(phyloseq)
# sample_data(ps)$Group

diffpath = paste(otupath,"/stemp_diff/",sep = "")
dir.create(diffpath)
# https://mp.weixin.qq.com/s/DTOz37JgH80kuLNi6Ae6-g
#---分组两两提取
map = phyloseq::sample_data(ps)
map$ID = row.names(map)

# map$Group = as.factor(map$Group)
sample_data(ps) = map

allgroup <- combn(unique(map$Group),2)
for (i in 1:dim(allgroup)[2]) {
  ps_sub <- phyloseq::subset_samples(ps,Group %in% allgroup[,i]);ps_sub
  
  for (j in 2:7) {
    p <- stemp_diff(ps = ps_sub,Top = 20,ranks = j)
    p
    # filename = paste(diffpath,"/",paste(allgroup[,i][1],allgroup[,i][2],sep = "_"),"stemp_P_plot.csv",sep = "")
    # write.csv(diff.mean,filename)
    filename = paste(diffpath,"/",paste(allgroup[,i][1],
                                        allgroup[,i][2],sep = "_"),phyloseq::rank.names(ps)[j],"stemp_plot.pdf",sep = "")
    ggsave(filename,p,width = 14,height = 6)
    
    filename = paste(diffpath,"/",paste(allgroup[,i][1],
                                        allgroup[,i][2],sep = "_"),phyloseq::rank.names(ps)[j],"stemp_plot.jpg",sep = "")
    ggsave(filename,p,width = 14,height = 6)
  }
}

detach("package:patchwork")

library(ggClusterNet)
library(tidyverse)
library(phyloseq)


source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\EdgerSuper.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\EdgerSuper2.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\Mui.cluster-group.volcano.R")

# 4.5 聚类差异火山图指定分组#-------
# source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\Plot.CompareWithCK.R",encoding = "utf-8")
id = sample_data(ps)$Group %>% unique()
# group1 = c("Group1","Group2")
# b= data.frame(group1)
diffpath.1 = paste(otupath,"/Mui.cluster.v/",sep = "")
dir.create(diffpath.1)

i= 1
for (i in 1:length(id)) {
  aaa = combn(id,2)
  diffpathv = paste(diffpath.1,"/Mui.cluster.v",paste(aaa[,i][1],aaa[,i][2],sep = "_"),sep = "")
  dir.create(diffpathv)

  group2 = c(aaa[1,i],aaa[2,i])
  b= data.frame(group2)

  res = EdgerSuper(ps = ps,
                   group  = "Group",
                   artGroup = b,
                   j = "OTU",
                   path = diffpathv
  )

  head(res)
  

  result = Mui.cluster.volcano(res = res,rs.k = 6)
  p = result[[1]]
  # p

  p1 = result[[2]]
  # p
  filename = paste(diffpathv,"/","Mui.cluster.volcano.label.pdf",sep = "")
  ggsave(filename,p,width = 8,height = 4)
  filename = paste(diffpathv,"/","Mui.cluster.volcano.pdf",sep = "")
  ggsave(filename,p1,width = 8,height = 4)
}


#-4.6 多组差异分析火山图#------
library(ggrepel)
diffpath.1 = paste(otupath,"/Mui.Group.v/",sep = "")
dir.create(diffpath.1)
res = EdgerSuper2 (ps = ps,group  = "Group",artGroup =NULL,
                   j = "OTU",
                   path = diffpath.1
)

head(res)


result = Mui.Group.volcano (res = res)
p = result[[2]]
p
filename = paste(diffpath.1,"/","Mui.group.volcano.pdf",sep = "")
ggsave(filename,p,width = 12,height = 6,limitsize = FALSE)


#--3.12 热图和气泡图展示微生物丰度--------

heatpath = paste(otupath,"/heapmap_boplot/",sep = "")
dir.create(heatpath)

#--注意map文件中一定要有ID列
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\Microheatmap.R",encoding = "utf-8")

#提取丰度最高的前20个OTU做展示
# ps_rela  = ps %>% scale_micro(method = "TMM")
# j = 2
# rank.names(ps)[j]
heatnum = 30
map = phyloseq::sample_data(ps)
map$ID = row.names(map)
phyloseq::sample_data(ps) = map
for (j in 2:6) {
  ps_tem = ps %>% 
    ggClusterNet::scale_micro(method = "TMM") %>%
    ggClusterNet::tax_glom_wt(ranks = j) 
  # rowSD = function(x){
  #   apply(x,1, sd)
  # }
  # 
  # rowCV = function(x){
  #   rowSD(x)/rowMeans(x)
  # }
  
  id <- ps %>% 
    ggClusterNet::scale_micro(method = "TMM") %>%
    ggClusterNet::tax_glom_wt(ranks = j) %>%
    ggClusterNet::filter_OTU_ps(100) %>%
    ggClusterNet::vegan_otu() %>%
    t() %>% as.data.frame() %>%rowCV %>%
    sort(decreasing = TRUE) %>%
    head(heatnum) %>%
    names()
  
  result <- Microheatmap(ps_rela = ps_tem,id = id ,col_cluster = FALSE,row_cluster = FALSE)
  
  p1 <- result[[1]] 
  p1
  # p1 +  scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Set3"))(60))
  p2 <- result[[2]]
  p2
  filename = paste(heatpath,"/",phyloseq::rank.names(ps)[j],"Topggheatmap.pdf",sep = "")
  ggsave(filename,p1,width = 14,height = (6 + heatnum/10))
  
  filename = paste(heatpath,phyloseq::rank.names(ps)[j],"Topggbubble.pdf",sep = "")
  ggsave(filename,p2,width = 14,height = (6 + heatnum/10))
  
  filename = paste(heatpath,"/",phyloseq::rank.names(ps)[j],"Topggheatmap.png",sep = "")
  ggsave(filename,p1,width = 14,height = (6 + heatnum/10))
  
  filename = paste(heatpath,phyloseq::rank.names(ps)[j],"Topggbubble.png",sep = "")
  ggsave(filename,p2,width = 14,height = (6 + heatnum/10))
  
  # filename = paste(heatpath,"/",rank.names(ps)[j],"/","Topggheatmap.jpg",sep = "")
  # ggsave(filename,p1,width = 14,height = (6 + heatnum/10))
  # filename = paste(heatpath,"/",rank.names(ps)[j],"/","Topggbubble.jpg",sep = "")
  # ggsave(filename,p2,width = 14,height = (6 + heatnum/10))
}


#-5.1 lefse--py-在win子系统中可以运行#--------
if (T) {
  source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\lefse_py_pre.R",encoding = "utf-8")
  lefpath = paste(otupath,"/lefse_py/",sep = "")
  dir.create(lefpath)
  
  # library(phyloseq)
  # library(EasyMicrobiome)
  # library("tidyverse")
  
  tablefse <- lefse_py_pre(ps = ps,taxGlomRank = "Genus",filter = 250)
  dim(tablefse)
  filename = paste(lefpath,"/LEFSE_to_run_G_level.txt",sep = "")
  write.table(tablefse,filename,append = F, quote = F,col.names= F,sep = "\t")
  
  # #文件预处理
  # conda activate qiime1
  # lefse-format_input.py LEFSE_to_run_G_level.txt pri_lefse.in -c 1 -u 2 -o 1000000
  # run_lefse.py pri_lefse.in pri_lefse_2.res  -l 2
  # lefse-plot_res.py pri_lefse_2.res lefse_barplot.pdf --format pdf
  # lefse-plot_cladogram.py pri_lefse_2.res lefse_tree.pdf --format pdf
  # mkdir biomarkers_raw_images
  # lefse-plot_features.py pri_lefse.in pri_lefse_2.res biomarkers_raw_images/
  
}



lefsepath = paste(otupath,"/lefse_R_plot/",sep = "")
dir.create(lefsepath)

# library(ggpubr)
# library(patchwork)
# library(MicrobiotaProcess)
library(ggtree)
source("E:/Shared_Folder/Function_local/R_function/micro/R_lefse_SAV.R",encoding = "utf-8")
# source("../micro/R_lefse_SAV.R",encoding = "utf-8")


p1 <- p_base(ps,Top = 100)
p1$data

tablda = LDA_Micro(ps = ps,
                   Top = 100,
                   p.lvl = 0.05,
                   lda.lvl = 1,
                   seed = 11, 
                   adjust.p = F)

p <- lefse_bar(taxtree = tablda[[2]])

tem = tablda[[2]]


FileName <- paste(lefsepath,"bar_lefse", ".pdf", sep = "")
ggsave(FileName, p, width = 25, height =dim(tem)/4,limitsize = FALSE)
FileName <- paste(lefsepath,"bar_lefse", ".png", sep = "")
ggsave(FileName, p, width = 60, height =dim(tem)/4,dpi = 72,limitsize = FALSE)

res = tablda[[2]]
FileName <- paste(lefsepath,"tree_lefse_data", ".csv", sep = "")
write.csv(res,FileName,quote = F)

# 注释树
library(patchwork)
p2 <- clade.anno_wt(p1, tablda[[1]], alpha=0.3, anno.depth = 7)
FileName <- paste(lefsepath,"tree_lefse", ".pdf", sep = "")
ggsave(FileName,p2,width = 15,height = 10)

FileName <- paste(lefsepath,"tree_lefse", ".png", sep = "")
ggsave(FileName,p2,width = 15,height = 10,dpi = 72)


source("E:/Shared_Folder/Function_local/R_function/micro/R_lefse_allRank.R",encoding = "utf-8")

lefsepath = paste(otupath,"/lefse_R_plot/",sep = "")
dir.create(lefsepath)


for (j in 2:6) {
  
  p1 <- p_base(ps,Top = 200,ranks =j)
  p1
  
  tablda = LDA_Micro(ps = ps,
                     Top = 200,
                     ranks = j,
                     p.lvl = 0.05,
                     lda.lvl = 2,
                     seed = 11,
                     adjust.p = F)
  
  p2 <- clade.anno_wt(p1, tablda[[1]], alpha=0.3,anno.depth = 2)
  p2
  FileName <- paste(lefsepath,j,"_tree_lefse", ".pdf", sep = "")
  ggsave(FileName,p2,width = 15,height = 10)
  FileName <- paste(lefsepath,j,"_tree_lefse", ".png", sep = "")
  ggsave(FileName,p2,width = 15,height = 10)
  p <- lefse_bar(taxtree = tablda[[2]])
  FileName <- paste(lefsepath,j,"_bar_lefse", ".pdf", sep = "")
  ggsave(FileName, p, width = 15, height =9)
  
  FileName <- paste(lefsepath,j,"_bar_lefse", ".png", sep = "")
  ggsave(FileName, p, width = 15, height =9)
  
  res = tablda[[2]]
  FileName <- paste(lefsepath,j,"_tree_lefse_data", ".csv", sep = "")
  write.csv(res,FileName,quote = F)
}




#-5.2 机器学习--两两分类建模 #-----


source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\MicroMachine_learning.R")

library(randomForest)
library(caret)
library(ROCR) ##用于计算ROC
library(e1071)

id = sample_data(ps)$Group %>% unique()
aaa = combn(id,2)
i= 1
group = c(aaa[1,i],aaa[2,i])
b= data.frame(group)

ps = ps %>% subset_taxa.wt("Family","Unassigned",T)
ps = ps %>% subset_taxa.wt("Order","Unassigned",T)
ps = ps %>% subset_taxa.wt("Genus","Unassigned",T)
ps = ps %>% subset_taxa.wt("Phylum","Unassigned",T)
ps = ps %>% subset_taxa.wt("class","Unassigned",T)


for (i in 1:length(aaa[1,])) {
  matpath = paste(otupath,"/Machine_learing.",paste(aaa[,i][1],aaa[,i][2],sep = "."),sep = "")
  dir.create(matpath )
  
  pst = ps %>% subset_samples.wt("Group",group) %>%
    filter_taxa(function(x) sum(x ) > 10, TRUE)
  randomforest.wt(
    pst = pst,
    ROC = F,
    rfcv = F,
    optimal = 50,
    matpath=matpath)
}





#-6.1 分组之间距离比较和可视化#-------
#--一下调整必须在map文件中包含两列，
# 一列是时间梯度，使用数字表示，不能带有特殊表述，一列是treat，就是分组信息，这些分组信息必须是有时间梯度的。
if (6 >gnum & gnum > 2) {
  
  group = "Group"
  map = as.data.frame(sample_data(ps))
  alppath = paste(otupath,"/distance/",sep = "")
  dir.create(alppath)
  gro = map[,group] %>% unique()
  colnames(gro) = "group"
  conbgroup = combn(gro$group,2)
  # 计算包括终点均值的所有样品bray距离
  bray_curtis = vegan::vegdist(vegan_otu(ps), method = "bray")
  bray_curtis = as.matrix(bray_curtis)
  
  for (i in 1:dim(conbgroup)[2]) {
    a = conbgroup[,i]
    map = as.data.frame(sample_data(ps))
    head(map)
    
    chose1 = map[as.matrix(map[,group]) %>% as.vector() == a[1],] %>% row.names()
    chose2 = map[as.matrix(map[,group]) %>% as.vector() == a[2],] %>% row.names()
    
    dat = data.frame(group = paste(a[1],a[2],sep = "_VS_"), Distance =bray_curtis[chose1,chose2] %>% as.dist() %>% as.vector() )
    head(dat)
    
    if (i == 1) {
      table = dat
    }
    
    if (i != 1) {
      table = rbind(table,dat)
    }
  }
  
  head(table)
  table$id = 1:dim(table)[1]
  data <- table %>% dplyr::select(id,everything())
  
  
  result = EasyStat::MuiKwWlx(data = data,num = c(3))
  FileName <- paste(alppath,"/distance_label.csv", sep = "")
  write.csv(result,FileName,sep = "")
  FileName <- paste(alppath,"/distance_index.csv", sep = "")
  write.csv(data,FileName,sep = "")
  
  result1 = EasyStat::FacetMuiPlotresultBox(data = data,num = c(3),result = result,sig_show ="abc",ncol = 1 )
  p1_1 = result1[[1]] + 
    # scale_x_discrete(limits = axis_order) + 
    mytheme2 +
    guides(fill = guide_legend(title = NULL)) +
    scale_fill_manual(values = colset4)
  p1_1
  res = EasyStat::FacetMuiPlotresultBar(data = data,num = c(3),result = result,sig_show ="abc",ncol = 1)
  p1_2 = res[[1]]+ 
    # scale_x_discrete(limits = axis_order) + 
    guides(color = FALSE) +
    mytheme2+ 
    guides(fill = guide_legend(title = NULL))+
    scale_fill_manual(values = colset4)
  p1_2
  res = EasyStat::FacetMuiPlotReBoxBar(data = data,num = c(3),result = result,sig_show ="abc",ncol = 1)
  p1_3 = res[[1]]+ 
    # scale_x_discrete(limits = axis_order) + 
    mytheme2 + 
    guides(fill = guide_legend(title = NULL))+
    scale_fill_manual(values = colset4)
  p1_3
  FileName <- paste(alppath,"distance_box", ".pdf", sep = "")
  ggsave(FileName, p1_1, width = ((4+ gnum) ), height =8,limitsize = FALSE)
  
  FileName <- paste(alppath,"distance_bar", ".pdf", sep = "")
  ggsave(FileName, p1_2, width = ((4+gnum) ), height = 8,limitsize = FALSE)
  
  FileName <- paste(alppath,"distance_boxbar", ".pdf", sep = "")
  ggsave(FileName, p1_3, width = ((4+gnum) ), height = 8,limitsize = FALSE)
  
  FileName <- paste(alppath,"distance_box", ".jpg", sep = "")
  ggsave(FileName, p1_1, width = (( 4+gnum) ), height =8,limitsize = FALSE)
  
  FileName <- paste(alppath,"distance_bar", ".jpg", sep = "")
  ggsave(FileName, p1_2, width = ((4+ gnum) ), height = 8,limitsize = FALSE)
  
  FileName <- paste(alppath,"distance_boxbar", ".jpg", sep = "")
  ggsave(FileName, p1_3, width = ((4+ gnum) ), height = 8,limitsize = FALSE)
  
}


#-3.13 maptree#-------
maptpath = paste(otupath,"/maptree/",sep = "")
dir.create(maptpath)
#--OTU 水平的maptree
library(ggClusterNet)
library(phyloseq)
library(tidyverse)
library(ggraph)
library(data.tree)
library(igraph)

source("E:\\Shared_Folder\\Function_local\\R_function\\micro/maptree_micro.R")

tax = ps %>% vegan_tax() %>% 
  as.data.frame()
head(tax)
#--如果数据中存在两个king，则无法运行，所以，这里去掉古菌OTU。
ps.bac <- ps %>% subset_taxa.wt("Kingdom", "Bacteria") 
ps.bac

p = micro_maptree(ps = ps.bac,
                  Top = 100,
                  labtab =  NULL,
                  seed = 11)
p

FileName <- paste(maptpath,"maptree", ".pdf", sep = "")
ggsave(FileName, p,width = 18,height = 16)


FileName <- paste(maptpath,"maptree", ".png", sep = "")
ggsave(FileName, p,width = 18,height = 16,dpi = 72)



#--6.2 FEAST溯源分析#-----
betapath = paste(otupath,"/Feast/",sep = "")
dir.create(betapath)


# source("E:\\Shared_Folder\\Function_local\\R_function\\Liu_project\\amplicon-master\\R\\开发花絮\\FEAST-master\\FEAST_src\\src.R")
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\FEAST.R",encoding = "UTF-8")


result = FEAST.micro(ps = ps,
               group = "Group",
               sinkG = "Group1",
               sourceG = c("Group2","Group3"),
               path = "E:/Shared_Folder/Function_local/R_function/micro/" # 注意按照自己设定的路径进行修改
               )

# result
p <- Plot_FEAST(data = result)
p
p2 = MuiPlot_FEAST(data = result)

FileName <- paste(betapath,"Feast_data.csv", sep = "")
write.csv(result,FileName, quote = F)

FileName <- paste(betapath,"Feast_Group.pdf", sep = "")
ggsave(FileName, p, width = 5, height = 5)

p2 = MuiPlot_FEAST(data = result)
p2
FileName <- paste(betapath,"Feast_sample.pdf", sep = "")
ggsave(FileName, p2, width = 15, height = 15)


#-7.2 FAPROTAX#----

Bca = T


#----FAPROTAX#---
Bca = TRUE
if (Bca ==T ) {
  # 加载R包
  library(microeco)
  # 加载ggplot2绘图包并设置样式
  library(ggplot2)
  library("WGCNA")
  library(tidyverse)
  library(ggtree)
  library("SpiecEasi")
  library(ggClusterNet)
  library(phyloseq)
  library(magrittr)
  p_list = c("ggplot2", "BiocManager", "devtools","picante", "GUniFrac", 
             "ggalluvial", "rgexf")
  # for(p in p_list){if (!requireNamespace(p)){install.packages(p)}
  #   library(p, character.only = T, quietly = T, warn.conflicts = F)}
  
  # ps = readRDS("./data/dataNEW/ps_16s.rds")
  
  
  
  otu = ps %>% vegan_otu() %>%
    t() %>%
    as.data.frame()
  
  tax = ps %>% vegan_tax() %>%
    as.data.frame()
  # 构建分析对象
  dataset = microtable$new(sample_table = sample_data(ps), otu_table = otu, tax_table = tax)
  
  t2 = trans_func$new(dataset)
  t2$cal_spe_func()
  t2$res_spe_func[1:5, 1:6]
  
  data = t2$res_spe_func
  
  data = data[rowSums(data)> 0,]
  
  betapath = paste(otupath,"/FAPROTAX/",sep = "")
  dir.create(betapath)
  
  
  # dir.create("./result_and_plot/Base_diversity_16s//OTU_220921//FAPROTAX")
  
  write.csv(data,paste(betapath,"/FAPROTAX.csv",sep = ""))
  # 查看功能 分组列表
  t2$func_group_list
  # 查看某一类
  t2$show_prok_func("methanotrophy")
  
}

#  元素循环注释结果之鞥呢判断这个微生物是否有这个能力，而不能判断能力强弱
head(data)
source("E:\\Shared_Folder\\Function_local\\R_function\\micro/barMainplot.R")
barpath = paste(otupath,"/FAPROTAX.plot/",sep = "")
dir.create(barpath)


ids = colnames(data)
for (i in 1:length(ids)) {
  id = ids[i]
  dat2  = data %>% select(id) %>% rownames_to_column("ID") %>%
    filter(!!sym(id) > 0)
  head(dat2)
  
  if (dim(dat2)[1] != 0) {
    
    pst =ps %>% subset_taxa.wt("OTU",dat2$ID)
    
    phyloseq::rank_names(ps)
    pst = ps %>% subset_taxa.wt("Species","Unassigned",TRUE)
    
    
    result = barMainplot(ps = pst,
                         j = "Genus",
                         # axis_ord = axis_order,
                         label = FALSE,
                         sd = FALSE,
                         tran = F,
                         Top = 12)
    p4_1 <- result[[1]] + 
      # scale_fill_brewer(palette = "Paired") + 
      scale_fill_manual(values = colset2) +
      scale_x_discrete(limits = axis_order) +
      labs(title = id)+
      mytheme1
    p4_1
    
    p4_2  <- result[[3]] + 
      # scale_fill_brewer(palette = "Paired") + 
      scale_fill_manual(values = colset2) +
      scale_x_discrete(limits = axis_order) + 
      labs(title = id)+
      mytheme1
    p4_2
    
    databar <- result[[2]] %>% group_by(Group,aa) %>%
      summarise(sum(Abundance)) %>% as.data.frame()
    head(databar)
    colnames(databar) = c("Group",j,"Abundance(%)")
    
    
    FileName1 <- paste(barpath,"/a2_",id,"_barflow",".pdf", sep = "")
    ggsave(FileName1, p4_2, width = (5+ gnum), height =8,limitsize = FALSE)
    FileName2 <- paste(barpath,"/a2_",id,"_barflow",".jpg", sep = "")
    ggsave(FileName2, p4_2, width = (5+ gnum), height =8,limitsize = FALSE)
    
    FileName1 <- paste(barpath,"/a2_",id,"_bar",".pdf", sep = "")
    ggsave(FileName1, p4_1, width = (5+ gnum), height =8 ,limitsize = FALSE)
    FileName2 <- paste(barpath,"/a2_",id,"_bar",".jpg", sep = "")
    ggsave(FileName2, p4_1, width = (5+ gnum), height =8,limitsize = FALSE)
    
    FileName <- paste(barpath,"/a2_",j,"_bar_data",".csv", sep = "")
    write.csv(databar,FileName,quote = F)
  }
  
}


head(data)

dat3 = data.frame(ID = colnames(data),class = 'Function')
head(dat3)


dat4 = t2$func_group_list %>% unlist() %>% as.data.frame()
head(dat4)

A = c()

for (i in 1:length(dat3$ID)) {
  A[i] = (t2$show_prok_func(dat3$ID[i]))
}
dat3$class = A

FileName <- paste(barpath,"/class_data",".csv", sep = "")
write_csv(dat3,FileName)






