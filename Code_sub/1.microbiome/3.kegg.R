library(LorMe)
library(dplyr)
library(stringr)
library(tidyverse)
library(agricolae)
library(patchwork)
library(vegan)
library(magrittr)
metafile=read.table("./input/metafile.txt",header = T,sep="\t")
my_col=color_scheme("Plan2") %>% .[c(1,3,4)]
names(my_col)=c("CK","CF","BP")
#kegg profile####
load("./input/kegg_cato.rda")
kegg=read.csv("./input/kegg.profile.L3.csv",header = T)
kegg_anno=left_join(data.frame(Category3=kegg$Category),kegg_cato)
kegg_anno[is.na(kegg_anno)]="Other"
kegg_anno=data.frame(KEGGID=paste0("Kegg",1:nrow(kegg_anno)),anno=paste(kegg_anno$Category1,kegg_anno$Category2,kegg_anno$Category3,kegg_anno$Category3,kegg_anno$Category3,kegg_anno$Category3,kegg_anno$Category3,sep = "_"))
kegg_anno$anno=str_replace(kegg_anno$anno,": ","-")

kegg_obj=tax_summary(groupfile = metafile,inputtable = kegg[,2:10],reads = T,taxonomytable = kegg_anno,into = "standard",outputtax = "standard",sep = "_")
kegg_obj_plan1=object_config(taxobj = kegg_obj,treat_location = 2,treat_col = my_col,treat_order = c("CK","CF","BP")) 


beta_kegg=structure_plot(taxobj = kegg_obj_plan1,taxlevel = "Class",diagram = "stick")
F3A=beta_kegg$PCoA_Plot
beta_kegg$PERMANOVA_statistics
set.seed(999)
vegdist(t(kegg_obj_plan1$Genus[,-1]), method = "bray")  %>% 
  betadisper(
    group = kegg_obj_plan1$groupfile$Treatment,
    type = "median"
  ) %>% permutest(permutations = 999)

n=15
L3_pct=kegg_obj_plan1$Class_percent
L3_anno=kegg_obj_plan1$Class_taxonomy

L3_all=data.frame(L3_pct,L3_anno)
L3_metab=L3_all[L3_all$Domain=="Metabolism",]
L3_metab_top=L3_metab[order(rowMeans(L3_metab[,2:10]),decreasing = T) %>% .[1:n],]
L3_metab_top_trans=t(L3_metab_top[,2:10]) %>% data.frame(metafile[,1:2])
colnames(L3_metab_top_trans)[1:n]=L3_metab_top$Class


var_names <- L3_metab_top$Class 

get_sign_letter <- function(v) {
  fit <- aov( L3_metab_top_trans[,v]~Treatment, data = L3_metab_top_trans)
  lsd <- LSD.test(fit, "Treatment")
  df_let <- as.data.frame(lsd$groups)
  df_let$Treatment <- rownames(df_let)
  df_let$Index <- v
  df_let[, c("Treatment", "Index", "groups")] %>% rename(letter = groups)
}
letter_df_kegg <- map_dfr(var_names, get_sign_letter)

heat_mean_kegg <- L3_metab_top_trans %>%
  pivot_longer(-c(SampleID, Treatment), names_to = "Index", values_to = "Value") %>%
  group_by(Index) %>%
  mutate(z_score = scale(Value)[, 1]) %>%
  group_by(Treatment, Index) %>%
  summarise(z_mean = mean(z_score), .groups = "drop") %>%
  left_join(letter_df_kegg, by = c("Treatment", "Index"))

heat_mean_kegg$Index=factor(heat_mean_kegg$Index,levels = var_names %>% rev())
heat_mean_kegg$Treatment=factor(heat_mean_kegg$Treatment,levels = c("CK","CF","BP"))
F3B=ggplot(heat_mean_kegg, aes(x = Treatment, y = Index)) +
  geom_tile(aes(fill = z_mean), color = "white", linewidth = 1) +
  geom_text(aes(label = letter), size = 4, fontface = "bold") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  theme_zg()+
  theme(
    axis.ticks = element_blank(),
    panel.background= element_blank(),
    axis.title = element_blank()
  )+guides(fill="none")

L3_metab_top_mean=data.frame(KEGGL3=L3_metab_top$Class,mean=rowMeans(L3_metab_top[,2:10]))
L3_metab_top_mean$KEGGL3=factor(L3_metab_top_mean$KEGGL3,levels = var_names %>% rev())

F3C=ggplot(L3_metab_top_mean,aes(x=mean,y=KEGGL3))+
  geom_vline(xintercept = c(0.05,0.1,0.15,0.2),size=.2,lty=2)+
  geom_col()+
  scale_x_continuous(expand = c(0,0),labels = scales::label_percent())+
  theme_zg()+
  theme(axis.ticks.y =element_blank(),
        axis.text.y=element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(color = "black",size=0.2),
        axis.title=element_blank())

##overall LSD####
n=47
L3_pct=kegg_obj_plan1$Class_percent
L3_anno=kegg_obj_plan1$Class_taxonomy

L3_all=data.frame(L3_pct,L3_anno)
L3_metab=L3_all[L3_all$Domain=="Metabolism",]
L3_metab_top=L3_metab[order(rowMeans(L3_metab[,2:10]),decreasing = T) %>% .[1:n],]
L3_metab_top_trans=t(L3_metab_top[,2:10]) %>% data.frame(metafile[,1:2])
colnames(L3_metab_top_trans)[1:n]=L3_metab_top$Class


var_names <- L3_metab_top$Class 

get_sign_letter <- function(v) {
  fit <- aov( L3_metab_top_trans[,v]~Treatment, data = L3_metab_top_trans)
  lsd <- LSD.test(fit, "Treatment")
  df_let <- as.data.frame(lsd$groups)
  df_let$Treatment <- rownames(df_let)
  df_let$Index <- v
  df_let[, c("Treatment", "Index", "groups")] %>% rename(letter = groups)
}
letter_df_kegg <- map_dfr(var_names, get_sign_letter)

kept_list_kegg=as.character()
for(i in unique(letter_df_kegg$Index) ){
  subdata=letter_df_kegg[letter_df_kegg$Index==i,]
  if(subdata[subdata$Treatment=="BP",3]=="a"&subdata[subdata$Treatment=="CF",3]!="a"&subdata[subdata$Treatment=="CK",3]!="a"){
    kept_list_kegg=c(kept_list_kegg,i)
  }
}


L3_metab_top_long=L3_metab_top[L3_metab_top$Class %in% kept_list_kegg,1:10] %>% gather("Type","Value",-c(Class))
L3_metab_top_long$Class=factor(L3_metab_top_long$Class,levels = kept_list_kegg)
L3_metab_top_long$Treatment=c(rep("CK",30),rep("CF",30),rep("BP",30))
L3_metab_top_long$Treatment=factor(L3_metab_top_long$Treatment,levels =c("CK","CF","BP") )
#N-cycling profile####
ncyc=read.csv("./input/Ncyc.csv",header = T)
ncyc_anno=data.frame(ID=ncyc$GeneID,anno=paste(ncyc$Type,ncyc$SubType,ncyc$GeneID,ncyc$GeneID,ncyc$GeneID,ncyc$GeneID,ncyc$GeneID,sep = ";"))

ncyc_obj=tax_summary(groupfile = metafile,inputtable = ncyc[,2:10],reads = T,taxonomytable = ncyc_anno,into = "standard",outputtax = "standard",sep = ";")
ncyc_obj_plan1=object_config(taxobj = ncyc_obj,treat_location = 2,treat_col = my_col,treat_order = c("CK","CF","BP")) 

beta_ncyc=structure_plot(taxobj = ncyc_obj_plan1,taxlevel = "Class",diagram = "stick")
F3D=beta_ncyc$PCoA_Plot+guides(fill="none")
beta_ncyc$PERMANOVA_statistics

com_ncyc=community_plot(taxobj =ncyc_obj_plan1,taxlevel = "Domain",n =7,palette = "Paired")

get_sign_letter <- function(v) {
  fit <- aov( as.numeric(com_ncyc$Top7Domain[com_ncyc$Top7Domain$Domain==v,-1])~Treatment, data = metafile)
  lsd <- LSD.test(fit, "Treatment")
  df_let <- as.data.frame(lsd$groups)
  df_let$Treatment <- rownames(df_let)
  df_let$Index <- v
  df_let[, c("Treatment", "Index", "groups")] %>% rename(letter = groups)
}
letter_df_kegg <- map_dfr(com_ncyc$Top7Domain$Domain, get_sign_letter)
letter_df_kegg


F3E=com_ncyc$alluvialplot+guides(fill="none")
ncyc_type_long=com_ncyc$Grouped_Top7Domain

ncyc_subtype=ncyc_obj_plan1$Phylum_percent

subtype_stat=data.frame()
for(i in ncyc_subtype$Phylum){
  subdata=ncyc_subtype[ncyc_subtype$Phylum==i,2:10] %>% as.numeric()
  sub_stat=aov(subdata~Treatment,data=metafile) %>% summary()
  pval=sub_stat[[1]] %>% .[1,5]
  meanval=mean(subdata)
  subtype_stat=rbind(subtype_stat,data.frame(Subtype=i,pvalue=pval,Mean_abundance=meanval))
}


###
ncyc_subtype_trans=t(ncyc_subtype[,-1]) %>% data.frame(metafile[,1:2])
colnames(ncyc_subtype_trans)[1:55]=ncyc_subtype$Phylum
var_names <- ncyc_subtype$Phylum 

get_sign_letter <- function(v) {
  fit <- aov( ncyc_subtype_trans[,v]~Treatment, data = ncyc_subtype_trans)
  lsd <- LSD.test(fit, "Treatment")
  df_let <- as.data.frame(lsd$groups)
  df_let$Treatment <- rownames(df_let)
  df_let$Index <- v
  df_let[, c("Treatment", "Index", "groups")] %>% rename(letter = groups)
}
letter_df_kegg <- map_dfr(var_names, get_sign_letter)
kept_list=as.character()
for(i in ncyc_subtype$Phylum ){
  subdata=letter_df_kegg[letter_df_kegg$Index==i,]
  if(subdata[subdata$Treatment=="BP",3]=="a"&subdata[subdata$Treatment=="CF",3]!="a"&subdata[subdata$Treatment=="CK",3]!="a"){
    kept_list=c(kept_list,i)
  }
}

kept_list=kept_list[which(kept_list %in%  subtype_stat$Subtype[subtype_stat$pvalue<0.05])]

kept_ncyc_long=ncyc_obj_plan1$Phylum_percent[ncyc_obj_plan1$Phylum_percent$Phylum %in% kept_list,] %>% 
  gather("Tax","Rel",-c(Phylum))
kept_ncyc_long$Treatment=c(rep("CK",12),rep("CF",12),rep("BP",12))

inputdata<- kept_ncyc_long 
facet_compare<-data.frame(compare=0,Letters=0,type=0,Mean=0,std=0,letterp=0,facetlabel=0)[0,]
for(i in unique(inputdata[,1])){
  sub_facet<-inputdata[which((inputdata[,1])==i),] 
  results<-auto_signif_test(data =sub_facet,treatment_col =4,value_col =3,prior = T)
  facet_compare<-rbind(facet_compare,
                       data.frame(results$comparison_letters[,1:5],letterp=max(1.2*(results$comparison_letters %>% .[,'Mean'])+1.3*(results$comparison_letters %>% .[,'std'])),facetlabel=i))
}

colnames(facet_compare)[7]=colnames(inputdata)[1]
mean_frame<- aggregate(inputdata[,3],by=list(inputdata[,4],inputdata[,1]),FUN=mean)
Sd<- aggregate(inputdata[,3],by=list(inputdata[,4],inputdata[,1]),FUN=sd)%>% .[,'x']
Treatment_Name<- mean_frame$Group.1
N<- table(inputdata[,4]) %>% as.numeric()
Mean<- mean_frame[,'x']
SEM<- Sd/(N^0.5)
input_mean_frame=data.frame(Treatment_Name,N,Mean,Sd,SEM,mean_frame[,'Group.2'])
colnames(input_mean_frame)[6]=colnames(inputdata)[1]
input_mean_frame$Treatment_Name=factor(input_mean_frame$Treatment_Name,levels = c("CK","CF","BP"))
facet_compare$compare=factor(facet_compare$compare,levels = c("CK","CF","BP"))
F3F=ggplot(input_mean_frame,aes(x=as.factor(Treatment_Name),y=Mean))+
  theme_zg()+
  geom_bar(stat = 'identity',size=0.2,width=0.5,aes(fill=as.factor(Treatment_Name)),color='#000000',alpha=0.8)+
  geom_errorbar(aes(ymin=Mean-Sd,ymax=Mean+Sd),size=0.2,width=0.2)+
  scale_y_continuous(expand = c(0,0))+
  scale_color_manual(values=my_col)+
  scale_fill_manual(values=my_col)+
  labs(x='',fill='',y='')+
  facet_wrap(~Phylum,scales ='free_y',strip.position ='top',as.table =T,nrow = 2)+
  geom_text(data=facet_compare,aes(x=as.factor(compare),y=letterp,label=Letters),size=6)+
  theme(legend.position = 'right')+guides(fill="none")


F3up=(F3A|plot_spacer()|F3C)*guides(fill="none")+plot_layout(widths = c(0.8,1.5,1.3))
F3down=(F3D|F3E|F3F)+plot_layout(widths = c(1,1,2))
F3up/
  F3down

ggsave("./Figure3/Figure3.pdf",width = 12,height = 6.5)


ncyc_select=ncyc_obj_plan1$Phylum_percent[ncyc_obj_plan1$Phylum_percent$Phylum %in% kept_list,]
save(kegg_obj_plan1,ncyc_select,kept_list_kegg,file="./outdata/kegg.rda")
