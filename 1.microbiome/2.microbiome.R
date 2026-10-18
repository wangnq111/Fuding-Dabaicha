library(LorMe)
library(magrittr)
library(tidyverse)
library(agricolae)
library(vegan)
library(patchwork)
my_col=color_scheme("Plan2") %>% .[c(1,3,4)]
names(my_col)=c("CK","CF","BP")
metafile=read.table("./input/metafile.txt",header = T,sep="\t")
taxon=read.csv("./input/taxon.csv",header = T)

##pre-processing####
bac=taxon[grep("d__Bacteria",taxon$Tax),]
bac_clean=Filter_function(input = bac,threshold = 10,format = 2)

##packaging####
obj=tax_summary(groupfile = metafile,inputtable = bac_clean[,2:10],reads = T,taxonomytable = bac_clean[,c(1,11)],into = "standard",outputtax = "standard")
obj_plan1=object_config(taxobj = obj,treat_location = 2,treat_col = my_col,treat_order = c("CK","CF","BP")) 

##community feature####
alpha=Alpha_diversity_calculator(taxobj = obj_plan1,taxlevel = "Genus")
alpha$plotlist$Plotobj_Shannon$Barplot
alpha$plotlist$`Plotobj_Species number`$Barplot

F2A=(alpha$plotlist$Plotobj_Shannon$Barplot+alpha$plotlist$`Plotobj_Species number`$Barplot)*guides(fill="none")
beta=structure_plot(taxobj = obj_plan1,taxlevel = "Genus",diagram = "stick")
beta$PERMANOVA_statistics
set.seed(999)
vegdist(t(obj_plan1$Genus[,-1]), method = "bray")  %>% 
  betadisper(
    group = obj_plan1$groupfile$Treatment,
    type = "median"
  ) %>% permutest(permutations = 999)


F2B=beta$PCoA_Plot+guides(fill="none")
community=community_plot(taxobj = obj_plan1,taxlevel = "Phylum",n = 10)
F2C=community$barplot+guides(fill="none")

get_sign_letter <- function(v) {
  fit <- aov( as.numeric(community$Top10Phylum[community$Top10Phylum$Phylum==v,-1])~Treatment, data = metafile)
  lsd <- LSD.test(fit, "Treatment")
  df_let <- as.data.frame(lsd$groups)
  df_let$Treatment <- rownames(df_let)
  df_let$Index <- v
  df_let[, c("Treatment", "Index", "groups")] %>% rename(letter = groups)
}
letter_df_kegg <- map_dfr(community$Top10Phylum$Phylum, get_sign_letter)
letter_df_kegg

F2up=(F2A|F2B|F2C)+plot_layout(widths = c(3,2,4))

##indicator####

indic=indicator_analysis(taxobj = obj_plan1,taxlevel = "Genus")
indic1=indic
indic1$tag[!indic1$tag %in% c("CK","CF","BP")]="None"

indic_stat=table(indic1$Phylum[indic1$tag %in% "BP"]) %>% sort()
for(i in names(indic_stat)){
  orig_phylum_num=nrow(obj_plan1$Genus_taxonomy[obj_plan1$Genus_taxonomy$Phylum==i,])
  indic_stat[i]=indic_stat[i]/orig_phylum_num
}
sort(indic_stat)
##network####
network_results=network_analysis(taxobj = obj_plan1,taxlevel = "Genus",reads = T,n = 9,rel_threshold = 0.0001,threshold = 0.98,method = "pearson")
nodes_table=network_results$Nodes_info

set.seed(999)
network_v=network_withdiff(network_obj = network_results,diff_frame = indic1,aes_col =my_col )
pdf("./Figure2/network.pdf",width = 8,height = 6)
network_visual_re(network_visual_obj =network_v,module_num = c(3,6,10,1,2,4,8),module_paint = T,vertex.size = 4)
dev.off()

module_abundance=Module_abundance(network_obj = network_results,No.module = c(3,6,10,1,2,4,8))
F2E=(module_abundance$plotlist$Plotobj_Module3$Barplot/
  module_abundance$plotlist$Plotobj_Module6$Barplot/
  module_abundance$plotlist$Plotobj_Module10$Barplot)*guides(fill="none")

(module_abundance$plotlist$Plotobj_Module1$Barplot+guides(fill="none")|
  module_abundance$plotlist$Plotobj_Module2$Barplot+guides(fill="none"))/
  (module_abundance$plotlist$Plotobj_Module4$Barplot+guides(fill="none")|
  module_abundance$plotlist$Plotobj_Module8$Barplot+guides(fill="none"))

ggsave("./Supp/Sup_fig3.pdf",width = 6  ,height = 4)


module_pie=Module_composition(network_obj = network_results,No.module = c(3,6,10,1,2,4,8),taxlevel = "Order",mode = "most",top_n = 5,palette = "Plan6")

module_pie$Module3$Pie_plot_Module3/
  module_pie$Module6$Pie_plot_Module6/
  module_pie$Module10$Pie_plot_Module10
ggsave("./Figure2/Fi2F.pdf",width = 3,height = 7)

(module_pie$Module1$Pie_plot_Module1|
  module_pie$Module2$Pie_plot_Module2)/
 ( module_pie$Module4$Pie_plot_Module4|
  module_pie$Module8$Pie_plot_Module8)
ggsave("./Supp/Sup_fig3.pdf",width = 6  ,height = 4)

F2down=(plot_spacer()|F2E|plot_spacer())+plot_layout(widths = c(4.2,1.5,1))
(F2up/F2down)+plot_layout(heights = c(1,2.5))
ggsave("./Figure2/Fi2main.pdf",width = 9,height = 7)

save(obj_plan1,module_abundance,indic1,file="./outdata/micro.rda")

