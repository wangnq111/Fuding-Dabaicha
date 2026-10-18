library(magrittr)
library(LorMe)
library(dplyr)
library(FactoMineR)
library(factoextra)
library(vegan)
library(tidyverse)
library(agricolae)
library(patchwork)

my_col=color_scheme("Plan2") %>% .[c(1,3,4)]
names(my_col)=c("CK","CF","BP")
soil=read.table("./input/soil.txt",header = T,sep="\t")
tea=read.table("./input/tea.txt",header = T,sep="\t")
metafile=read.table("./input/metafile.txt",header = T,sep="\t")


##soil PCA####
res.pca_soil <- PCA(soil[,3:12], 
               scale.unit = TRUE,# 默认对数据进行scale标准化，最后响应变量间的均值为0，sd=1。
               ##相同尺度的标准化避免了一些变量仅仅因为其较大的度量单位而成为主导变量。使变量具有可比性。
               ncp = 11, # 结果中默认保留5，nrow(env[-ind.sup,])-1(样本数)和ncol(env[2:14][-c(quanti.sup,quali.sup)])(响应变量数)中的最小值个主成分
               #ind.sup = 1, # 附加样本，1表示数据表的第一行，不是必须的数据，默认为NULL。附加样本的坐标将根据PCA分析结果进行预测。
               quanti.sup = NULL, # 附加的定量响应变量数据，不是必须的，默认为NULL。基于PCA分析结果预测坐标。
               quali.sup = 1:2, # 附加定性变量，分类信息最好不要用数字表示，在行数区间的数字，运行会报错。
               ##1:2表示定性变量在env[2:14]中的列数。不是必须的，默认为NULL。基于PCA分析结果预测坐标，可用于对样本进行分组着色。
               row.w = NULL, col.w = NULL, # 分别给样本和响应变量设定权重，默认所有样本或所有响应变量的权重都是一样的。默认相当于对变量数据进行了中心化，变量均值接近0。可以设置为1个数字，也可以设置为分别与行数或列数相等的矢量。
               graph = FALSE)

sam_soil=data.frame(res.pca_soil$ind$coord) %>% .[1:2] %>% data.frame(.,metafile[,c(1:2)])
get_eigenvalue(res.pca_soil) 

cent = aggregate(sam_soil[, 1:2], by = list(sam_soil$Treatment),FUN = mean)
colnames(cent) = c("Treatment", "cent1", "cent2")
PCAframe_soil = left_join(data.frame(sam_soil, joint = sam_soil$Treatment), 
                     data.frame(cent[, 2:3], joint = cent$Treatment)) %>% suppressMessages()
PCAframe_soil = PCAframe_soil[, -which(colnames(PCAframe_soil) == "joint")]
set.seed(0)
adonis2((soil[,3:12])~Treatment+Block,data=metafile,method = "euclidean",by = "margin")
set.seed(0)
vegdist(soil[,3:12], method = "euclidean")  %>% 
  betadisper(
    group = obj_plan1$groupfile$Treatment,
    type = "median"
  ) %>% permutest(permutations = 999)
F1A=ggplot(PCAframe_soil,aes(x=Dim.1,y=Dim.2,fill=factor(Treatment)))+
  geom_hline(yintercept=0,linetype=4,color="grey") +
  geom_vline(xintercept=0,linetype=4,color="grey") +
  geom_segment(aes(xend = cent1, yend = cent2,color=Treatment),show.legend = FALSE, size = 0.3, alpha = 0.6)+
  geom_point(size=2.5,pch=21,alpha=.8)+
  scale_color_manual(values = my_col)+
  scale_fill_manual(values = my_col)+
  labs(x="PC1: 76.59%",y="PC2: 22.82%",fill="")+
  # stat_ellipse(level=0.95)+
  theme_zg()+guides(fill="none")
## soil heat####
var_names <- setdiff(colnames(soil), c("Sample", "Treatment"))

get_sign_letter <- function(v) {
  fit <- aov(as.formula(paste(v, "~ Treatment")), data = soil)
  lsd <- LSD.test(fit, "Treatment", p.adj = "bonferroni") 
  df_let <- as.data.frame(lsd$groups)
  df_let$Treatment <- rownames(df_let)
  df_let$Index <- v
  df_let[, c("Treatment", "Index", "groups")] %>% rename(letter = groups)
}
letter_df_soil <- map_dfr(var_names, get_sign_letter)

heat_mean_soil <- soil %>%
  pivot_longer(-c(Sample, Treatment), names_to = "Index", values_to = "Value") %>%
  group_by(Index) %>%
  mutate(z_score = scale(Value)[, 1]) %>%
  group_by(Treatment, Index) %>%
  summarise(z_mean = mean(z_score), .groups = "drop") %>%
  left_join(letter_df_soil, by = c("Treatment", "Index"))

index_label_soil <- c(
  pH = "Soil pH",
  OM = "Organic matter",
  TN = "Total N",
  NO3N = "Nitrate N",
  NH4N = "Ammonium N",
  AP = "Available P",
  AK = "Available K",
  URE = "Urease",
  ACP = "Acid phosphatase",
  SUC = "Sucrase"
)

heat_mean_soil$Index=factor(heat_mean_soil$Index,levels = colnames(soil)[-c(1:2)])
F2A=ggplot(heat_mean_soil, aes(x = Treatment, y = Index)) +
  geom_tile(aes(fill = z_mean), color = "white", linewidth = 1) +
  geom_text(aes(label = letter), size = 4, fontface = "bold") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  scale_y_discrete(labels = index_label_soil) +
  #labs(x = "Treatments", y = "Soil indicators", fill = "Standardized Z-score") +
  coord_flip()+
  theme_zg()+
  theme(
    axis.ticks = element_blank(),
    panel.background= element_blank(),
    axis.title = element_blank()
  )+guides(fill="none")

##tea PCA####
res.pca_tea <- PCA(tea[,3:15], 
                   scale.unit = TRUE,# 默认对数据进行scale标准化，最后响应变量间的均值为0，sd=1。
                   ##相同尺度的标准化避免了一些变量仅仅因为其较大的度量单位而成为主导变量。使变量具有可比性。
                   ncp = 11, # 结果中默认保留5，nrow(env[-ind.sup,])-1(样本数)和ncol(env[2:14][-c(quanti.sup,quali.sup)])(响应变量数)中的最小值个主成分
                   #ind.sup = 1, # 附加样本，1表示数据表的第一行，不是必须的数据，默认为NULL。附加样本的坐标将根据PCA分析结果进行预测。
                   quanti.sup = NULL, # 附加的定量响应变量数据，不是必须的，默认为NULL。基于PCA分析结果预测坐标。
                   quali.sup = 1:2, # 附加定性变量，分类信息最好不要用数字表示，在行数区间的数字，运行会报错。
                   ##1:2表示定性变量在env[2:14]中的列数。不是必须的，默认为NULL。基于PCA分析结果预测坐标，可用于对样本进行分组着色。
                   row.w = NULL, col.w = NULL, # 分别给样本和响应变量设定权重，默认所有样本或所有响应变量的权重都是一样的。默认相当于对变量数据进行了中心化，变量均值接近0。可以设置为1个数字，也可以设置为分别与行数或列数相等的矢量。
                   graph = FALSE)

sam_tea=data.frame(res.pca_tea$ind$coord) %>% .[1:2] %>% data.frame(.,metafile[,c(1:2)])
get_eigenvalue(res.pca_tea) 

cent = aggregate(sam_tea[, 1:2], by = list(sam_tea$Treatment),FUN = mean)
colnames(cent) = c("Treatment", "cent1", "cent2")
PCAframe_tea = left_join(data.frame(sam_tea, joint = sam_tea$Treatment), 
                         data.frame(cent[, 2:3], joint = cent$Treatment)) %>% suppressMessages()
PCAframe_tea = PCAframe_tea[, -which(colnames(PCAframe_tea) == "joint")]

adonis2((tea[,3:15])~Treatment,data=metafile,method = "euclidean")
F3A=ggplot(PCAframe_tea,aes(x=Dim.1,y=Dim.2,fill=factor(Treatment)))+
  geom_hline(yintercept=0,linetype=4,color="grey") +
  geom_vline(xintercept=0,linetype=4,color="grey") +
  geom_segment(aes(xend = cent1, yend = cent2,color=Treatment),show.legend = FALSE, size = 0.3, alpha = 0.6)+
  geom_point(size=2.5,pch=21,alpha=.8)+
  scale_color_manual(values = my_col)+
  scale_fill_manual(values = my_col)+
  labs(x="PC1: 74.42%",y="PC2: 19.32%",fill="")+
  # stat_ellipse(level=0.95)+
  theme_zg()+guides(fill="none")
## tea heat####
var_names <- setdiff(colnames(tea), c("SampleID", "Treatment"))

get_sign_letter <- function(v) {
  fit <- aov(as.formula(paste(v, "~ Treatment")), data = tea)
  lsd <- LSD.test(fit, "Treatment", p.adj = "bonferroni")
  df_let <- as.data.frame(lsd$groups)
  df_let$Treatment <- rownames(df_let)
  df_let$Index <- v
  df_let[, c("Treatment", "Index", "groups")] %>% rename(letter = groups)
}
letter_df_tea <- map_dfr(var_names, get_sign_letter)

heat_mean_tea <- tea %>%
  pivot_longer(-c(SampleID, Treatment), names_to = "Index", values_to = "Value") %>%
  group_by(Index) %>%
  mutate(z_score = scale(Value)[, 1]) %>%
  group_by(Treatment, Index) %>%
  summarise(z_mean = mean(z_score), .groups = "drop") %>%
  left_join(letter_df_tea, by = c("Treatment", "Index"))

index_label_tea <- c(
  TN = "Leaf total nitrogen",
  TP = "Leaf total phosphorus",
  TK = "Leaf total potassium",
  TPP = "Tea polyphenols",
  FAA = "Free amino acids",
  WE = "Water extracts",
  SS = "Soluble sugars",
  CAF = "Caffeine",
  EC = "(-)-Epicatechin",
  EGC = "(-)-Epigallocatechin",
  EGCG = "(-)-Epigallocatechin gallate",
  C = "(+)-Catechin",
  ECG = "(-)-Epicatechin gallate"
)

heat_mean_tea$Index=factor(heat_mean_tea$Index,levels = colnames(tea)[-c(1:2)])
F4A=ggplot(heat_mean_tea, aes(x = Treatment, y = Index)) +
  geom_tile(aes(fill = z_mean), color = "white", linewidth = 1) +
  geom_text(aes(label = letter), size = 4, fontface = "bold") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  scale_y_discrete(labels = index_label_tea) +
  labs(x = "Treatments", y = "tea indicators", fill = "Standardized Z-score") +
  coord_flip()+
  theme_zg()+
  theme(
    axis.ticks = element_blank(),
    panel.background= element_blank(),
    axis.title = element_blank()
  )+guides(fill="none")

((F1A|F2A)+plot_layout(widths = c(1,3)))/
  ((F3A|F4A)+plot_layout(widths = c(1,3)))


ggsave("./Figure1/Fig1.pdf",width = 8,height = 4.5)
save(soil,tea,file="./outdata/pheodata.rda")


agdata=read.csv("./input/agdata.csv",header = T)
inputdata<- agdata 
facet_compare<-data.frame(compare=0,Letters=0,type=0,Mean=0,std=0,letterp=0,facetlabel=0)[0,]
for(i in unique(inputdata[,2])){
  sub_facet<-inputdata[which((inputdata[,2])==i),] 
  results<-auto_signif_test(data =sub_facet,treatment_col =1,value_col =3,prior = T,comparison_method = "bonferroni")
  facet_compare<-rbind(facet_compare,data.frame(results$comparison_letters[,1:5],letterp=max(1.3*(results$comparison_letters %>% .[,'Mean'])+1.3*(results$comparison_letters %>% .[,'std'])),facetlabel=i))
}

colnames(facet_compare)[7]=colnames(inputdata)[2]
mean_frame<- aggregate(inputdata[,3],by=list(inputdata[,1],inputdata[,2]),FUN=mean)
Sd<- aggregate(inputdata[,3],by=list(inputdata[,1],inputdata[,2]),FUN=sd)%>% .[,'x']
Treatment_Name<- mean_frame$Group.1
N<- table(inputdata[,1]) %>% as.numeric()
Mean<- mean_frame[,'x']
SEM<- Sd/(N^0.5)
input_mean_frame=data.frame(Treatment_Name,N,Mean,Sd,SEM,mean_frame[,'Group.2'])
colnames(input_mean_frame)[6]=colnames(inputdata)[2]
ggplot(input_mean_frame,aes(x=as.factor(Treatment_Name),y=Mean))+
  theme_zg()+
  geom_bar(stat = 'identity',size=0.2,width=0.5,aes(fill=as.factor(Treatment_Name)),color='#000000',alpha=0.8)+
  geom_errorbar(aes(ymin=Mean-Sd,ymax=Mean+Sd),size=0.2,width=0.2)+
  scale_y_continuous(expand = c(0,0.01))+
  scale_color_manual(values=my_col)+
  scale_fill_manual(values=my_col)+
  labs(x='',fill='',y=' ')+
  facet_wrap(~get(colnames(inputdata)[2]),scales ='free_y',strip.position ='top',as.table =T)+
  geom_text(data=facet_compare,aes(x=as.factor(compare),y=letterp,label=Letters),size=6)+
  theme(legend.position = 'right',
  )+
  guides(color='none',fill='none')
ggsave("./Supp/Sup_fig2.pdf",width = 6,height = 3)

