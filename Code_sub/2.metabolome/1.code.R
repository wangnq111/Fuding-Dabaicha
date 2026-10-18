library(ropls)
library(tidyverse)
library(LorMe)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(vegan)
library(tidyverse)
library(magrittr)
library(agricolae)

metabolite=read.csv("./input/metabolome.csv",header = T,encoding = "UTF-8",as.is = TRUE)
metafile=read.table("./input/metafile.txt",header = T,sep="\t")
my_col=color_scheme("Plan2") %>% .[c(1,3,4)]
names(my_col)=c("CK","CF","BP")

kegg_anno=data.frame(KeggID=paste0("kegg",1:nrow(metabolite)),metabolite[,1:30])
metabolite_pct=sweep(metabolite[,31:48],colSums(metabolite[,31:48]),"/",MARGIN=2)



# OPLSDA
plsda_input=metabolite[,31:48] %>% t()
colnames(plsda_input)=metabolite$Compound_ID

plsda<-opls(plsda_input, metafile$Treatment,
            predI = 2, orthoI = 0,log10L = T,crossvalI =3,scaleC="pareto",permI=999)#,#orthoI = 1:opls-da; orthoI = 0:Pls-da
#log10L = F,
#crossvalI =10,
#scaleC="pareto",#pareto scaling
#fig.pdfC = c("none", "interactive", "outputfile/plsda.pdf")[3],
#permI=200)

# data and score data
## score data
plsdaScore <- data.frame(
  t1 =plsda@scoreMN[,1],
  to1 =plsda@scoreMN[,2] ) %>% 
  scale(center = T,scale = T) %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'SampleID') %>% 
  merge(metafile,by="SampleID")

t1Weight=sprintf("%.1f%%", plsda@modelDF[1,1] * 100);t1Weight
to1Weight=sprintf("%.1f%%", plsda@modelDF[2,1] * 100);to1Weight

R2X=plsda@modelDF[1,1]+plsda@modelDF[2,1]
R2Y=plsda@modelDF[1,3]+plsda@modelDF[2,3]
Q2Y=plsda@modelDF[1,6]+plsda@modelDF[2,6]

subTitle <- paste0("R2X=",R2X,"  R2Y=",R2Y,"  Q2Y=",Q2Y)
cent = aggregate(plsdaScore[, 2:3], by = list(Treatment=plsdaScore$Treatment),FUN = mean)
colnames(cent)[2:3] = c("cent1", "cent2")
plsdaScore=merge(plsdaScore,cent,by=c("Treatment"))
## opls-da 
F4A=ggplot(plsdaScore,aes(x=t1,y=to1,fill=Treatment))+
  geom_segment(aes(xend = cent1, yend = cent2,color=Treatment),show.legend = FALSE, size = 0.3, alpha = 0.8)+
  geom_point(size=1.5,alpha=.8,pch=21)+
  geom_hline(yintercept = 0, linetype = 4,color = "grey") +
  geom_vline(xintercept = 0, linetype = 4,color = "grey")+
  theme_zg()+
  scale_fill_manual(values=my_col)+
  scale_color_manual(values=my_col)+
  labs(#title = "OPLS-DA", 
    #subtitle = subTitle, 
    x = paste0("T score1 (",t1Weight,")"),
    y = paste0("T score2 (",to1Weight,")"))+
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1.2)), 
    plot.subtitle = element_text(hjust = 0.5, size = rel(0.6)) 
  )+guides(fill="none")
#ggsave("./results/structure/plsda.pdf",width = 5,height = 4)


##stat
set.seed(999)
adonis2(plsda_input~metafile$Treatment, by = "terms",method="euclidean")
set.seed(999)
vegdist(t(metabolite[,31:48]), method = "euclidean")  %>% 
  betadisper(
    group = metafile$Treatment,
    type = "median"
  ) %>% permutest(permutations = 999)


##composition####
metabolite_summary=aggregate(metabolite_pct,by=list(metabolite$ClassII),FUN=sum) 
metabolite_summary$Group.1[rowMeans(metabolite_summary[,-1])<=(rowMeans(metabolite_summary[,-1]) %>% sort(decreasing = T) %>% .[10])]="Others"
metabolite_summary=aggregate(metabolite_summary[,-1],by=list(metabolite_summary$Group.1),FUN=sum)
metabolite_summary_combine=data.frame(metafile,t(metabolite_summary[,-1]))
colnames(metabolite_summary_combine)[-c(1:2)]=metabolite_summary$Group.1
metabolite_summary_combine_mean=aggregate(metabolite_summary_combine[,-c(1:2)],by=list(metabolite_summary_combine$Treatment),FUN=mean)  #,metabolite_summary_combine$sex
metabolite_summary_long=gather(metabolite_summary_combine_mean,"classify","Rel",-c(Group.1))

metabolite_order=metabolite_summary$Group.1[rowMeans(metabolite_summary[,-1]) %>% order(.,decreasing=T)] %>%.[c(2:10,1)]
metabolite_summary_long$Group.1=factor(metabolite_summary_long$Group.1,levels =c("CK","CF","BP"))
metabolite_summary_long$classify=factor(metabolite_summary_long$classify,levels=metabolite_order)
library(ggalluvial)
classfy_col=color_scheme("Plan9",names = metabolite_order)


F4B=ggplot(metabolite_summary_long, aes(x = Group.1, y = Rel, fill = classify)) + 
  #facet_grid(~Group.1,scale="free")+
  geom_col(position = "stack", width = 0.8)+
  #geom_flow(alpha = 0.5, width = 0.1)+
  scale_y_continuous(expand = c(0,0),labels = scales::percent)+
  scale_x_discrete(expand = c(0.1,0))+
  theme_zg()+
  labs(x="",y="Relative abundance")+
  scale_fill_manual(values=classfy_col)+
  theme(legend.text = element_text(size = 8, face = "bold"), 
        strip.text = element_text(size = 8, face = "bold"), 
        axis.ticks.x = element_blank(), panel.background = element_rect(color = NA), 
        axis.line.y = element_line(size = 0.4, color = "black"), 
        axis.ticks.length.y = unit(0.4, "lines"), axis.ticks.y = element_line(color = "black", 
                                                                              size = 0.4))

#compunds stat####
phylum_stat=data.frame()
for(i in colnames(metabolite_summary_combine)[-c(1:2)]){
  y=metabolite_summary_combine[,i]
  stat=aov(y~Treatment,data=metabolite_summary_combine) %>% summary()%>% .[[1]] %>% as.data.frame()
  stat$Item=rownames(stat)
  stat$signif=case_when(
    stat$`Pr(>F)`>0.05 ~"ns",
    stat$`Pr(>F)`>0.01&stat$`Pr(>F)`<0.05 ~"*",
    stat$`Pr(>F)`>0.001&stat$`Pr(>F)`<0.01 ~"**",
    stat$`Pr(>F)`<0.001 ~"***",
    is.na( stat$`Pr(>F)`)==T~NA
  )
  stat$phylum=i
  phylum_stat=rbind(phylum_stat,stat)
}

phylum_stat$sign=case_when(phylum_stat$`Pr(>F)`>0.05~'ns',
                           phylum_stat$`Pr(>F)`<0.05&phylum_stat$`Pr(>F)`>0.01~'*',
                           phylum_stat$`Pr(>F)`<0.01&phylum_stat$`Pr(>F)`>0.001~'**',
                           phylum_stat$`Pr(>F)`<0.001~'***')

get_sign_letter <- function(v) {
  fit <- aov( as.numeric(metabolite_summary[(metabolite_summary$Group.1)==v,-1])~Treatment, data = metafile)
  lsd <- LSD.test(fit, "Treatment")
  df_let <- as.data.frame(lsd$groups)
  df_let$Treatment <- rownames(df_let)
  df_let$Index <- v
  df_let[, c("Treatment", "Index", "groups")] %>% rename(letter = groups)
}
letter_df_kegg <- map_dfr(metabolite_summary$Group.1, get_sign_letter)
letter_df_kegg

###
meta_anno=data.frame(ID=metabolite$Compound_ID,anno=paste(metabolite$ClassI,metabolite$ClassII,metabolite$ClassIII,metabolite$ClassIIImetabolite$ClassIIImetabolite$ClassIIImetabolite$ClassIII,sep=";"))

metabolome_obj=tax_summary(groupfile = metafile,inputtable = apply(metabolite[,31:48],MARGIN=2,FUN=round)/1e5,reads = T,taxonomytable = meta_anno,into = "standard",outputtax = "standard",sep = ";")
metabolome_obj_plan1=object_config(taxobj = metabolome_obj,treat_location = 2,treat_col = my_col,treat_order = c("CK","CF","BP")) 




###
CK_BP_deseq=Deseq_analysis(taxobj = metabolome_obj_plan1,taxlevel = "Base",comparison = c("CK","BP"),cutoff = 1,control_name = "CK")
CF_BP_deseq=Deseq_analysis(taxobj = metabolome_obj_plan1,taxlevel = "Base",comparison = c("CF","BP"),cutoff = 1,control_name = "CF")
CF_CK_deseq=Deseq_analysis(taxobj = metabolome_obj_plan1,taxlevel = "Base",comparison = c("CF","CK"),cutoff = 1,control_name = "CK")

table(CK_BP_deseq$tag)
table(CF_BP_deseq$tag)
table(CF_CK_deseq$tag)

CK_BP_vol=volcano_plot(inputframe = CK_BP_deseq,cutoff = 1,aes_col = my_col)
CF_BP_vol=volcano_plot(inputframe = CF_BP_deseq,cutoff = 1,aes_col = my_col)
CF_CK_vol=volcano_plot(inputframe = CF_CK_deseq,cutoff = 1,aes_col = my_col)

F4C=
(CK_BP_vol$FC_FDR|
CF_BP_vol$FC_FDR|
CF_CK_vol$FC_FDR)*guides(fill="none")

#core_up=intersect(CK_BP_deseq$ID[CK_BP_deseq$tag=="BP"],CF_BP_deseq$ID[CF_BP_deseq$tag=="BP"])
#core_up_frame=kegg_anno[kegg_anno$Compound_ID %in% core_up,]

##indicator####
indic=indicator_analysis(taxobj = metabolome_obj_plan1,taxlevel = "Class")
indic1=indic
indic1$tag[!indic1$tag %in% c("CK","CF","BP")]="None"

BP_enrich=indic1[indic1$tag %in% c("CK","CF","BP"),]

BP_enrich_frame=metabolome_obj_plan1$Class_percent[which(metabolome_obj_plan1$Class_taxonomy$ClassID %in% BP_enrich$ClassID),]
BP_enrich_frame=data.frame(BP_enrich_frame[,-1],row.names = BP_enrich$ClassID)
BP_enrich_anno=data.frame(ClassII=BP_enrich$Phylum,row.names =BP_enrich$ClassID )
BP_enrich_anno$ClassII[!BP_enrich_anno$ClassII %in% metabolite_order]="Others"
Treat_anno=data.frame(Treatment=metafile$Treatment,row.names = metafile$SampleID)
ann_col=list(
  Treatment=my_col,
  ClassII=classfy_col[1:10]
)

BP_enrich_frame_new=data.frame()
for(i in c("CK","CF","BP")){
  ID=BP_enrich$ClassID[BP_enrich$tag==i]
  for(j in names(classfy_col)[1:10]){
    ID2=rownames(BP_enrich_anno)[BP_enrich_anno$ClassII==j]
    core_ID=intersect(ID,ID2)
    temp_frame=BP_enrich_frame[core_ID,]
    BP_enrich_frame_new=rbind(BP_enrich_frame_new,temp_frame)
  }
}

pheatmap::pheatmap(BP_enrich_frame_new,scale = "row",annotation_row = BP_enrich_anno,annotation_col =Treat_anno,
                   show_rownames = F,show_colnames = F,annotation_colors = ann_col,cluster_rows = F,cellwidth = 4,cellheight = 1.4)



(F4A|F4B|F4C)*guides(fill="none")+plot_layout(widths = c(1,1,3))
ggsave("./Figure4/Figure4.pdf",width = 10,height = 2.2)

indic_kegg=indic1
##overall LSD####
library(agricolae)
n=224
L3_pct=metabolome_obj_plan1$Phylum_percent
L3_anno=metabolome_obj_plan1$Phylum_taxonomy

L3_all=data.frame(L3_pct,L3_anno)
L3_metab=L3_all[]
L3_metab_top=L3_metab[order(rowMeans(L3_metab[,2:19]),decreasing = T) %>% .[1:n],]
L3_metab_top_trans=t(L3_metab_top[,2:19]) %>% data.frame(metafile[,1:2])
colnames(L3_metab_top_trans)[1:n]=L3_metab_top$Phylum


var_names <- L3_metab_top$Phylum 

get_sign_letter <- function(v) {
  fit <- aov( L3_metab_top_trans[,v]~Treatment, data = L3_metab_top_trans)
  lsd <- LSD.test(fit, "Treatment")
  df_let <- as.data.frame(lsd$groups)
  df_let$Treatment <- rownames(df_let)
  df_let$Index <- v
  df_let[, c("Treatment", "Index", "groups")] %>% rename(letter = groups)
}
letter_df_kegg <- map_dfr(var_names, get_sign_letter)
letter_df_kegg=letter_df_kegg[letter_df_kegg$Index!="_",]

kept_list_metab=as.character()
for(i in unique(letter_df_kegg$Index) ){
  subdata=letter_df_kegg[letter_df_kegg$Index==i,]
  if(subdata[subdata$Treatment=="BP",3]=="a"&subdata[subdata$Treatment=="CF",3]!="a"&subdata[subdata$Treatment=="CK",3]!="a"){
    kept_list_metab=c(kept_list_metab,i)
  }
}

save(metabolome_obj_plan1,classfy_col,indic_kegg,kept_list_metab,file="./outdata/metabolome.rda")
