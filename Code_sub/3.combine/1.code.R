library(LorMe)
library(magrittr)
library(dplyr)


load("./input/kegg.rda")
load("./input/metabolome.rda")
load("./input/micro.rda")
load("./input/pheodata.rda")


#micro-metabolome network####
genus_pct=obj_plan1$Genus_percent
metabolite_pct=metabolome_obj_plan1$Class_percent
rownames(metabolite_pct)=metabolome_obj_plan1$Class_taxonomy$ClassID
metabolite_pct=metabolite_pct[,c(colnames(genus_pct)[-1])]

metabolite_pct_clean=metabolite_pct[which(rowMeans(metabolite_pct[,-1])>5e-5),]
genus_pct_clean=genus_pct[which(rowMeans(genus_pct[,-1])>5e-5),]


micro_metabo_newtwork=network_analysis2(input =genus_pct_clean,inputtype = 2,n = 9,
                                        threshold = 0.98,method = "pearson",display = T,
                                        input2 =metabolite_pct_clean,input2type = 3  )


network_adj=micro_metabo_newtwork$Adjacency_column_table
network_nodes=micro_metabo_newtwork$Nodes_info
network_nodes$group="Microbe"
network_nodes$group[grep("Class",network_nodes$nodes_id)]="Metabolites" 
network_nodes=left_join(data.frame(network_nodes,Genus=network_nodes$nodes_id),indic1[,8:14])

network_nodes_temp=left_join(data.frame(network_nodes[,1:14],ClassID=network_nodes$nodes_id),indic_kegg[,c(7,11)])
network_nodes$tag[network_nodes$group=="Metabolites"]=network_nodes_temp$tag[network_nodes_temp$group=="Metabolites"]

nodes_phylum=data.frame(Genus=network_nodes$nodes_id[network_nodes$group=="Microbe"]) %>% left_join(obj_plan1$Genus_taxonomy) %$% Phylum
main_phylum_col=community_plot(taxobj = obj_plan1,taxlevel = "Phylum",n = 10,nrow = 1) %$%
  filled_color
main_phylum_name= main_phylum_col%>% names()
nodes_phylum[!nodes_phylum %in% main_phylum_name]="Others"
network_nodes$Category[network_nodes$group=="Microbe"]=nodes_phylum
nodes_class=data.frame(ClassID=network_nodes$nodes_id[network_nodes$group=="Metabolites"]) %>% left_join(metabolome_obj_plan1$Class_taxonomy) %$% Phylum
nodes_class[!nodes_class %in% names(classfy_col)]="Others_metabo"
network_nodes$Category[network_nodes$group=="Metabolites"]=nodes_class
#col
network_nodes$col[network_nodes$group=="Microbe"]=main_phylum_col[nodes_phylum] %>% as.character()
network_nodes$col[network_nodes$group=="Metabolites"]=classfy_col[nodes_class] %>% as.character()
network_nodes$col[network_nodes$Category=="Others_metabo"]="#B09C85FF"
network_nodes$col2=case_when(
  network_nodes$tag=="BP"~"#6376A0",
  network_nodes$tag=="CF"~"#33B39F",
  network_nodes$tag=="CK"~"#FE5D5D",
  network_nodes$tag=="None"~"gray"
)
#size
network_nodes$size=case_when(
  network_nodes$node_degree==1 ~ 1,
  network_nodes$node_degree>1&network_nodes$node_degree<=5 ~ 1,
  network_nodes$node_degree>5&network_nodes$node_degree<=10 ~ 3,
  network_nodes$node_degree>10&network_nodes$node_degree<=20 ~ 5,
  network_nodes$node_degree>20 ~7
)
#shape
network_nodes$shape=case_when(
  network_nodes$group=="Microbe" ~"circle",
  network_nodes$group=="Metabolites" ~"square",
)
network_nodes$name=NA
network_nodes$name[network_nodes$tag %in% c("BP","CF","CK")]=network_nodes$nodes_id[network_nodes$tag%in% c("BP","CF","CK")]
rm_strange=unique(c(grep("0",network_nodes$name),
       grep("1",network_nodes$name),
       grep("2",network_nodes$name),
       grep("3",network_nodes$name),
       grep("4",network_nodes$name),
       grep("5",network_nodes$name),
       grep("6",network_nodes$name),
       grep("7",network_nodes$name),
       grep("8",network_nodes$name),
       grep("9",network_nodes$name))) %>% setdiff(grep("Class",network_nodes$name))
network_nodes$name[rm_strange]=NA

Meta_class_temp=metabolome_obj_plan1$Class_taxonomy
rownames(Meta_class_temp)=Meta_class_temp$ClassID

metabolome_obj_plan1$Class_taxonomy$ClassID %in% network_nodes$name[grep("Class",network_nodes$name)]
network_nodes$name[grep("Class",network_nodes$name)]=Meta_class_temp[network_nodes$name[grep("Class",network_nodes$name)],"Class"]

##condition####
network_nodes$name[which(network_nodes$name=="_")]=NA
network_nodes$name[which(network_nodes$node_degree<8&network_nodes$group=="Microbe")]=NA
network_nodes$name[which(network_nodes$node_degree<6&network_nodes$group=="Metabolites")]=NA

temp_stat=network_adj[network_adj$source %in% c("g__Lacisediminihabitans","g__Rhodanobacter","g__Homoserinimonas","g__Terrimesophilobacter"),] %$% 
  target %>% table() %>% sort()
key_names=temp_stat[temp_stat==4] %>% names()
temp=network_nodes[network_nodes$nodes_id %in% key_names,]


temp_stat=network_adj[network_adj$source %in% network_nodes$nodes_id[network_nodes$Order %in% c("o__Pedosphaerales","o__Gemmatimonadales","o__Burkholderiales")],] %$% 
  target %>% table() %>% sort()
key_names=temp_stat[temp_stat>0] %>% names()
temp=network_nodes[network_nodes$nodes_id %in% key_names,]

#network_nodes$name[network_nodes$node_degree>=5&network_nodes$group=="Microbe"]=network_nodes$nodes_id[network_nodes$node_degree>=5&network_nodes$group=="Microbe"]
#key_name=data.frame(Name=metabolome_obj_plan1$Class_taxonomy$Class,row.names = metabolome_obj_plan1$Class_taxonomy$ClassID) %>% .[network_nodes$nodes_id[network_nodes$node_degree>=5&network_nodes$group=="Metabolites"],]
#network_nodes$name[network_nodes$node_degree>=5&network_nodes$group=="Metabolites"]=key_name


library(igraph)
#visualization
adjfile=network_adj
verticefile=network_nodes
t_net_its <- graph_from_data_frame(adjfile,direct=F,vertices = verticefile)
#color
V(t_net_its)$color <- V(t_net_its)$col

#size
V(t_net_its)$size <- V(t_net_its)$size

#shape
V(t_net_its)$shape=V(t_net_its)$shape

#linetype
E(t_net_its)$lty <- ifelse(E(t_net_its)$value == 1, 1, 2) 


set.seed(999)
coords_t_its <- layout_(t_net_its,with_fr(niter=9999, grid="auto"))
pdf("./figure5/micro_metabo_network.pdf",width = 12,height = 8)
plot(t_net_its,vertex.label=V(t_net_its)$name,vertex.label.color=V(t_net_its)$col2,
     vertex.size=V(t_net_its)$size,
     layout=coords_t_its,
     vertex.shape=V(t_net_its)$shape,
     edge.lty = E(t_net_its)$lty)



legend("topright", legend=c("Microbe","Metabolite"),
       pch=c(21,22), pt.bg="grey70", pt.cex=2, bty="n", title="Group")

# 2) 节点大小
legend("topleft",
       legend=c("1","5","10","20"),
       pt.cex=c(1,1.5,2,2.5),
       pch=21, pt.bg="grey70", bty="n", title="Degree")

# 3) 节点颜色 (Phylum/Category)
# 假设你有一个颜色向量 main_phylum_col
legend("bottomleft",
       legend = names(main_phylum_col),
       inset = c(0, -0.1),
       pch = 21,                # 圆形
       pt.bg = main_phylum_col, # 填充色
       col = "black",           # 边框色
       pt.cex = 1.8,
       bty = "n",
       title = "Phylum")
legend("bottomright",
       legend = names(classfy_col)[1:10],
       inset = c(-0.3, 0.2),
       pch = 22,                # 圆形
       pt.bg = classfy_col[1:10], # 填充色
       col = "black",           # 边框色
       pt.cex = 1.8,
       bty = "n",
       title = "Phylum")

# 4) 边的线型 (相关性)
legend("bottomright",
       legend=c("Positive","Negative"),
       lty=c(1,2), col="black", bty="n", title="Correlation")
dev.off()

##alluvial####
###data
Module_frame=module_abundance$rowframe[,c(17:23)]
soil_frame=soil[,-c(1:2)]
colnames(soil_frame)=paste0("soil_",colnames(soil_frame))
tea_frame=tea[,-c(1:2)]
colnames(tea_frame)=paste0("tea_",colnames(tea_frame))
tea_frame_trans=t(tea_frame) %>% as.data.frame()

library(Hmisc)
ncyc_select$Phylum=paste0("Ncyc_",ncyc_select$Phylum)
func_prof=rbind(ncyc_select[,-1],
                kegg_obj_plan1$Class_percent[which(kegg_obj_plan1$Class_taxonomy$Class %in% kept_list_kegg),-1],
                metabolome_obj_plan1$Phylum_percent[which(metabolome_obj_plan1$Phylum_taxonomy$Phylum %in% kept_list_metab),c(colnames(genus_pct)[-1])])
rownames(func_prof)=c(ncyc_select$Phylum,
                      paste0("kegg_",kegg_obj_plan1$Class_taxonomy$ClassID[which(kegg_obj_plan1$Class_taxonomy$Class %in% kept_list_kegg)]),
                      paste0("meta_",metabolome_obj_plan1$Phylum_taxonomy$PhylumID[which(metabolome_obj_plan1$Phylum_taxonomy$Phylum %in% kept_list_metab)])
)

module_frame_trans=t(Module_frame) %>% as.data.frame()
rownames(module_frame_trans)=paste0("Micro_",rownames(module_frame_trans))




###calculate cor
make_corframe=function(x,y){
  micro_func_r=rcorr(t(x),t(y))
  micro_func_r$r[micro_func_r$P>0.05]=0
  
  rframe=micro_func_r$r[1:nrow(x),-c(1:nrow(x))] %>% 
    as.data.frame()
  rframe$from=rownames(rframe)
  rframe_trans=gather(rframe,"to","value",-from)
  rframe_trans=rframe_trans[rframe_trans$value!=0,]
  rframe_trans$Freq=1
  return(rframe_trans)
}

mic_func=make_corframe(module_frame_trans,func_prof)
mic_func$cor=ifelse(mic_func$value>0,"pos","neg")
mic_func$Category=mic_func$from
mic_func=mic_func[abs(mic_func$value)>0.8,]
func_tea=make_corframe(func_prof,tea_frame_trans)
func_tea$Category=func_tea$from
func_tea$cor=ifelse(func_tea$value>0,"pos","neg")
func_tea=func_tea[abs(func_tea$value)>0.8,]
rm=table(func_tea$from) %>% sort() %>% names %>% .[1:25]
rm=rm[grep("meta",rm)]
func_tea=func_tea[!func_tea$from %in% rm,]
core_name=intersect(mic_func$to,func_tea$from)
mic_func=mic_func[mic_func$to %in% core_name,]
func_tea=func_tea[func_tea$from %in% core_name,]

###combine####
library(ggsankeyfier)


sankey_data=rbind(mic_func,func_tea)

sankey_long <- pivot_stages_longer(
  data = sankey_data[,-3],
  stages_from = c( "from","to"),
  values_from = c("Freq"),  # 把 Phylum 当作附加美学变量带到长表，方便给 edge 着色
  additional_aes_from = c("Category","cor") 
)

sankey_long$stage=as.character(sankey_long$stage)
sankey_long$node=as.character(sankey_long$node)
sankey_long$stage=case_when(
  sankey_long$node %in% rownames(tea_frame_trans) ~ "Tea",
  sankey_long$node %in% rownames(func_prof) ~ "Function",
  sankey_long$node %in% rownames(module_frame_trans) ~ "Micro"
)


sankey_long$nodes_col=sankey_long$node
sankey_long$nodes_col[grep("Module",sankey_long$nodes_col)]="micro"
sankey_long$nodes_col[grep("kegg_",sankey_long$nodes_col)]="kegg"
#sankey_long$nodes_col[grep("kegg_",sankey_long$nodes_col)]%>% 
#  data.frame() %>% 
#  separate(col = ".",into=c("cate","ClassID"),sep = "_") %$% cate
sankey_long$nodes_col[grep("meta_",sankey_long$nodes_col)]="meta"
sankey_long$nodes_col[grep("Ncyc_",sankey_long$nodes_col)]="ncyc"
sankey_long$nodes_col[grep("tea_",sankey_long$nodes_col)]="tea"

sankey_long$Category1=sankey_long$Category
sankey_long$Category1[grep("Module",sankey_long$Category1)]="micro"
sankey_long$Category1[grep("kegg_",sankey_long$Category1)]="kegg"
sankey_long$Category1[grep("meta_",sankey_long$Category1)]="meta"
sankey_long$Category1[grep("Ncyc_",sankey_long$Category1)]="ncyc"


###config####
sankey_long$stage=factor(sankey_long$stage,levels = c("Micro","Function","Tea"))
traits_order=c(unique(sankey_long$node[grep("Ncyc",sankey_long$node)]),
               unique(sankey_long$node[grep("kegg",sankey_long$node)]),
               unique(sankey_long$node[grep("meta",sankey_long$node)]))
nodes_order=paste0("Module",c(3,6,10,1,2,4,8))
tea_order=sankey_long$node[grep("tea",sankey_long$node)] %>% unique()
sankey_long$node=factor(sankey_long$node,levels = c(traits_order,nodes_order,tea_order))
sankey_long$Category1=factor(sankey_long$Category1,levels = c("micro","ncyc","kegg","meta"))

nodes_col=c(color_scheme("Plan4")%>% .[1:5],"black","gray")
names(nodes_col)=c("micro","nycy","kegg","meta","tea","pos","neg")
pos <- position_sankey(order = "as_is", v_space =2,align = "center")

pos_text <- position_sankey(order = "as_is", v_space = 2, nudge_x = 0.08,align = "center")

sankey_long$nodeshow=as.character(sankey_long$node)
sankey_long$nodeshow[sankey_long$connector=="to"&sankey_long$Category1=="micro"]=NA
#---
sankey_long$nodeshow[grep("kegg_",sankey_long$nodeshow)]=
  sankey_long$nodeshow[grep("kegg_",sankey_long$nodeshow)]%>% 
  data.frame() %>% 
  separate(col = ".",into=c("cate","ClassID"),sep = "_")%>%
  left_join(kegg_obj_plan1$Class_taxonomy) %$% Class

sankey_long$nodeshow[grep("meta_",sankey_long$nodeshow)]=
  sankey_long$nodeshow[grep("meta_",sankey_long$nodeshow)]%>% 
  data.frame() %>% 
  separate(col = ".",into=c("cate","PhylumID"),sep = "_")%>%
  left_join(metabolome_obj_plan1$Phylum_taxonomy) %$% Phylum
###visual####
ggplot(sankey_long,
       aes(x = stage,
           y = Freq,
           group = node,
           connector = connector,
           edge_id = edge_id)) +
  # 边（用 Phylum 上色）-- 这里 legend 会只显示 Phylum，因为 geom_sankeyedge 的 aes(fill=Phylum)
  geom_sankeyedge(aes(fill = cor),
                  position = pos,
                  alpha = 0.8,
                  show.legend = TRUE) +
  # 节点（条形）-- 不显示 node 的图例，节点颜色按照 node_cols 映射（Phylum 节点会与 edge 颜色一致）
  geom_sankeynode(aes(fill = nodes_col),
                  position = pos,
                  show.legend = FALSE,alpha=.8) +
  # 标签（节点文字）
  geom_text(aes(label = nodeshow),
            stat = "sankeynode",
            position = pos_text,
            hjust = 0, size = 3) +
  # 填色映射：同时包含 phylum 名字 和其他 node 名字（需把 node_cols 与 phylum_cols 合并）
  scale_fill_manual(values = nodes_col) +
  # 美化 x 轴与主题
  scale_x_discrete(expand = expansion(add = c(0.2, 0.5)), position = "top") +
  theme_void() +
  theme(legend.position = "right",
        plot.margin = margin(8, 8, 8, 8))
ggsave("./figure5/figure5.pdf",width = 12,height = 5)






