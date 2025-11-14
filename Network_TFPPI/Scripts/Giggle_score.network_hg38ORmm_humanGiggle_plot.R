rm(list=ls())
args=commandArgs(trailingOnly = TRUE)
path=args[1]
species=args[2]
outpdf = args[3]
setwd(path)

library(openxlsx)
library(data.table)
library(dplyr)
library(showtext)

# Enable showtext for proper font embedding
showtext_auto()
font_add_google("Roboto", "roboto")  # or use Arial if you prefer

if(species=="mouse"){
  metafile="/workspace/rsrch1/tmp/DataBase_datatable0618/Network_TFGiggle1010/202501_CirTrom/mm10_tranfac_QC.txt"
}else{
    metafile="/workspace/rsrch1/tmp/DataBase_datatable0618/Network_TFGiggle1010/202501_CirTrom/hg38_tranfac_QC.txt"}

metaInfo = fread(metafile,header = T,fill = T) %>% select(2,3,5,6)

output=path
data <- "peakfile_GIGGLE_res_Top1K.xls"
col_name=c("file","file_size","overlaps","odds_ratio","fishers_two_tail","fishers_left_tail","fishers_right_tail","GIGGLE_score")

  df = as.data.frame(read.csv(data,sep = "\t",header = F) %>% select(1:8))
  df_raw=df[-1,] 
  names(df_raw)=col_name
  df_res = df_raw %>% mutate(External_id=gsub(".*\\/(.*)_sort.*.gz","\\1",file),
                             GIGGLE_score = as.numeric(GIGGLE_score)) %>% left_join(.,y=metaInfo)

  outfile= paste0("GIGGLE_res_",gsub("_GIGGLE_res_Top1K.xls","",basename(data)),"txt")
  if(species=="mouse"){
  res_tmp = df_res %>% 
  filter(grepl("[Ll]iver",Ontology)) %>% 
  group_by(Factor) %>% 
  arrange(desc(as.numeric(GIGGLE_score))) %>% dplyr::slice(1) %>% ungroup() %>% arrange(desc(GIGGLE_score))
  }else{
    res_tmp = df_res %>% 
      #filter(grepl("[Ll]iver",Ontology)) %>% 
      group_by(Factor) %>% 
      arrange(desc(as.numeric(GIGGLE_score))) %>% dplyr::slice(1) %>% ungroup() %>% 
      arrange(desc(GIGGLE_score)) %>% filter(GIGGLE_score >300)
    }
  #write.table(res_tmp,outfile,sep = "\t",quote = F,row.names = F)

  
  
geneList = fread("geneList.txt",header = F)  
genes_all <- c(res_tmp$Factor,geneList$V1[1:150]) %>% toupper()%>%unique() 


###### get the substract network######
network = fread("PPI_network.txt") %>% mutate(V1=toupper(v1),V2=toupper(v2)) %>% select(V1,V2)
#hub_TF=c("MYC","SP1","MAX","USF1","E2F1", "EGR1","JUN","TFAP2C","FOS","NFKB1","CREB1",  "MYCN", "MYB","RELA",   "STAT5A", "GATA1" ,"ESR1") 
sub_net = network %>% filter(
  V1 %in% genes_all,
  V2 %in% genes_all,
  (
    (V1 %in% toupper(res_tmp$Factor) & V2 %in% toupper(geneList$V1)) |
      (V2 %in% toupper(res_tmp$Factor) & V1 %in% toupper(geneList$V1))
  )
) #%>% filter(V1 %in% hub_TF | V2 %in% hub_TF )

#node = graph%>%activate(nodes)%>% as_tibble()
#edge = graph%>%activate(edges)%>% as_tibble()

### plot network by gggraph###

library(ggraph)
library(igraph)
library(tidygraph)
library(ggplot2)
graph = sub_net %>% as_tbl_graph() %>%
  mutate(degree = centrality_degree(),
         node_type = case_when(
           degree == 1 ~ "leaf",
           degree >= quantile(degree, 0.9) ~ "hub",  # top 10% = hubs
           TRUE ~ "other"
         )) 

#graph <- graph %>%
#  mutate(class = ifelse(name == "ESR1", "yellow",
#                        ifelse(name %in% toupper(res_tmp$Factor), "TF", "target"))) %>% na.omit()
graph <- graph %>%
  mutate(class = ifelse(toupper(name) %in% toupper(res_tmp$Factor), "TF", "target")) %>% 
 filter((class == "TF" & degree >15) | class=="target" |name=="ESR1") %>% filter(class=="target" | node_type == "hub" | name=="ESR1")


p=ggraph(graph, layout = "kk") +
  geom_edge_link(color = "grey80", edge_width = 0.1) +
  geom_node_point(aes(color = class,
                      shape = class,
                      size = degree)) +
  geom_node_text(
    aes(label = ifelse(class == "TF", name, "")),
    #aes(label = ifelse(node_type %in% c("hub", "leaf"), name, "")),
    repel = TRUE,
    size = 3,family = "roboto"
  )+
  scale_color_manual(values = c("TF" = "tomato", "target" = "steelblue"))+
                       scale_shape_manual(values = c("TF" = 19, "target" = 17))+
  #scale_color_manual(values = c("TF" = "tomato", "target" = "steelblue","yellow"="yellow")) +
  #scale_shape_manual(values = c("TF" = 19, "target" = 17,"yellow"=19)) +  #  circle & triangle
  scale_size_continuous(range = c(3, 6)) +
  theme_graph(base_family = "roboto")

ggsave(outpdf, device=cairo_pdf,plot = p, width = 8, height = 6,units = "in",family = "sans")
