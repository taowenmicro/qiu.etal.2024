#清空内存#######
rm(list=ls()) 
# 参数文件; Tao Wen, 2024.1.31#------
library(tidyverse)
library(phyloseq)
library(ggClusterNet)

#--输入和设定文件
# ps01 = base::readRDS("./data/dataNEW/ps.rds")
ps = ps01 = base::readRDS("./data/dataNEW//ps_16s.rds")%>%
  filter_taxa(function(x) sum(x ) > 0 , TRUE)
path.id = "16S_micro"

# ps = ps02 = base::readRDS("./data/dataNEW//ps_ITS.rds")%>% 
#   filter_taxa(function(x) sum(x ) > 0 , TRUE)
# path.id = "ITS_micro"


#--最终确定的phyloseq对象定义为ps


#--提取分组因子数量
gnum = phyloseq::sample_data(ps)$Group %>% unique() %>% length()
gnum
#--设定排序顺序1：按照ps对象中map文件顺序进行
axis_order =  phyloseq::sample_data(ps)$Group %>%unique();axis_order
# axis_order = c("KO","WT","OE")
#设定排序顺序2：设定顺序按照已有分组文件的顺序，读入map文件
#--物种分类树展示的OTU数量
Top_micro = 150
#--堆叠柱状图展示前Top的微生物,j 展示的微生物分类等级
Top = 12
#jj = j = "Phylum"
#-主题--颜色等
source("E:\\Shared_Folder\\Function_local\\R_function\\micro\\total_amplicon.R")
res = theme_my()
mytheme1 = res[[1]]
mytheme2 = res[[2]]; 
colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]
#  主题颜色更新-可改method可选anhei,出来暗黑风格的图片，Top就是堆叠柱状图展示的数量
#  升级颜色挑选方式
set.seed(122)
res = theme.col(ps = ps,method = "no",Top =c(Top + 1),gnum = gnum)
mytheme1 = res$mytheme[sample(c(1:4),size=1)]#可选1-4
mytheme2 = res$mytheme[sample(c(5:8),size=1)]; #可选5-8
colset1 = res$col.group[[sample(c(1:5),size=1)]];#1-5可选
colset2 = res$col.bar[[sample(c(1:3),size=1)]];#1-3可选
colset3 = res$col.time[[sample(c(1:6),size=1)]];#1-6可选
colset4 = colset3

# 构建保存结果文件夹
result<- dir.amp(smart = F)#--如果出现错误，设定smart = F；是由于注释信息不符合本代码标准导致的
res1path = result[[1]];res1path
#-构建子文件夹保存:设定当前日期保存
a = Sys.Date() %>% as.character()
a = gsub("-",".",a)
otupath = paste(res1path,"/",path.id,a,"/",sep = "");otupath
dir.create(otupath)
print(otupath)




































#  附件一些其他函数#--------


## R安装问题：修改R包安装路径#-------
# ps01 = base::readRDS("./Error/221121/ps_its.rds")
# ps01
# .libPaths()
# .libPaths(new="C:/Program Files/R/R-4.1.1/library")


# 检查七个等级注释名字
# change.rank.name(ps01)


## 用于筛选微生物#-----
# ps01 %>% scale_micro() %>%
#   filter_taxa(function(x) sum(x ) > 0.01 , TRUE)
# #低丰度微生物挑选
# ps01 %>% scale_micro() %>%
# filter_taxa(function(x) sum(x ) < 0.0001 , TRUE)
# #--如何筛选样本:去除sample1
# ps_sub <- subset_samples.wt(ps01,"ID",c("sample1"),T);ps_sub

# #--如何筛选微生物
# ps01 <- ps01 %>% subset_taxa.wt("Kingdom", id) ;ps01
# # 是否需要仅仅注释为细菌或者真菌的做进一步分析
# ps01 <- ps01 %>% subset_taxa.wt("Kingdom", id) ;ps01


# 更多的phylsoeq对象操作内容参考我的个人笔记#-------
# 【有道云笔记】ggClusterNet开发一系列小工具润滑微生物组数据分析.md
# 点击这个链接：https://note.youdao.com/s/Sv6dChHV



#  一些其他参数#------

#韦恩网络设置过滤阈值
ps_biost = ggClusterNet::filter_OTU_ps(ps = ps,Top = 500)
#--差异分析设定两两比对
# group1 = c("Gro1","Gro2")
# group2 = c("Gro1","Gro2")
# b= data.frame(group1,group2)
# b
b = NULL# 如果每两个组之间都做差异，那就指定b为NULL
# 热图展示的OTU数量
heatnum　=　30
#--R语言做lefse的过滤数量
ps_Rlefse = ggClusterNet::filter_OTU_ps(ps = ps,Top = 400)
#--机器学习部分
ROC = FALSE# ROC是三种机器学习的ROC曲线，但是只能跑两个分组，如果两组，可以选择为T。
rfcv = FALSE# 是否做交叉检验
optimal = 40 # 选择多少个重要变量


#-----选择性功能
#设置CK，用于双向柱状图绘制-目前不绘制
CK = unique(phyloseq::sample_data(ps)$Group)[1]
# 用python做lefse
lefse.py = T
if (lefse.py) {
  lefsenum = 0
  ps_lefse <- ps %>%
    phyloseq::subset_taxa(
      # Kingdom == "Fungi"
      Kingdom == id
      # Genus  == "Genus1"
      # Species %in%c("species1") 
      # row.names(tax_table(ps01))%in%c("OTU1")
    )
  
  ps_lefse = ggClusterNet::filter_OTU_ps(ps = ps_lefse,Top = 400)
  
}
