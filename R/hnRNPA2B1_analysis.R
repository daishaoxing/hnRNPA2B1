setwd("/data/dsx/hnRNPA2B1")
getwd ()
rm(list = ls())
# edata <- read.delim(file = './file/fpkm.anno.xls',sep='\t',header=TRUE, stringsAsFactors=FALSE) #31553
# edatatmp <- edata[,c(1:10)]
# write.table(edatatmp, file = "./file/genes_fpkm.tab", quote = FALSE,sep="\t",row.names = TRUE)
edata <- read.delim(file = './file/genes_fpkm.tab',sep='\t',header=TRUE, stringsAsFactors=FALSE) #31553
protein_coding<-edata$name[edata$Biotype=='protein_coding']
edata <- edata[,c(1:10)]
rownames(edata) <- edata[,1]
edata <- edata[,-1]

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
edata <- apply(edata,2,fpkmToTpm)

#### protein coding
edata<-edata[protein_coding,]
edata<-edata[rowSums(edata)>0,]
write.table(edata, file = "./file/edata_TPM.tab", quote = FALSE,sep="\t",row.names = TRUE)

# 1.PCA check ----
library(ggord)
library(ggplot2)
edata <- read.table(file = './file/edata_TPM.tab',sep='\t',header=TRUE,stringsAsFactors=FALSE)
edata<-edata[rowSums(edata)>0,]
pca_group=factor(c(rep('Control',3),rep('PAC2',3),rep('PAC5',3)))
edata.pca <- prcomp(t(edata), scale. = TRUE)
p <- ggord(edata.pca, grp_in = pca_group, arrow=0, vec_ext =0,txt=NULL,cols=c('blue','#2ca25f','red'),
           poly = T,polylntyp='dashed', 
           ellipse_pro = 0.9,alpha_el = 0.2,alpha = 1) + theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.y = element_text(color = 'black',size = 14),
        axis.title.x = element_text(color = 'black',size = 14),
        axis.text.y = element_text(color = 'black',size = 10),
        axis.text.x = element_text(color = 'black',size = 12),
        axis.ticks = element_line(color = 'black'),
        axis.line = element_line(color = 'black', size = 0.5),
        legend.position = 'top',
        legend.key.height = unit(0.6,'cm'),
        legend.text = element_text(face = 'bold',color = 'black',size = 10)
  )
p
ggsave('./picture/group_PCA.pdf', plot=p, dpi = 600,width = 6, height = 6)

#### annotation  gene
rm(list=ls())
library( "clusterProfiler")
library("org.Hs.eg.db")
edata = read.table(file = './file/edata_TPM.tab',sep='\t',header=TRUE, stringsAsFactors=FALSE)
edata$ENSEMBL<-rownames(edata)
columns(org.Hs.eg.db)
geneinfo = select(org.Hs.eg.db, keys=keys(org.Hs.eg.db), columns = c('ENSEMBL','ENTREZID',"SYMBOL"))
all_data<-merge(x=edata,y=geneinfo,by="ENSEMBL")
all_data<-all_data[,c(-1,-11)]
edata_rmdup <- aggregate(all_data[,-10],by = list(all_data$SYMBOL), FUN = mean)
rownames(edata_rmdup) <- edata_rmdup[,1]
write.table(edata_rmdup[,-1], file = "./file/edata_rmdup.tab", quote = FALSE,sep="\t",row.names = TRUE)

# 2.所有基因及样本热图 -------------------------------------------------------------
rm(list=ls())
library(pheatmap)
dif_data = read.table(file = './file/edata_rmdup.tab',header = TRUE, sep = "\t", quote = "",stringsAsFactors=FALSE)
dif_data = dif_data[,c(1:3,7:9)]
dif_data<-dif_data[rowSums(dif_data)>0,]
Z_score <- (dif_data - apply(dif_data, 1, mean)) / apply(dif_data, 1, sd)
annotation_col <- data.frame(Group = c(rep('Control',3),rep('PAC5',3)))
rownames(annotation_col) <- colnames(dif_data)

ancols = c('blue','red')
names(ancols) = c("Control","PAC5")
ann_colors <- list(Group=ancols)
pheatmap(Z_score,border_color="#C5C5C5",color = colorRampPalette(c("blue", "white", "red"))(50),
         clustering_distance_cols = "correlation",clustering_distance_rows = "correlation",show_rownames=F,annotation_col=annotation_col,
         annotation_colors=ann_colors,fontsize=9,filename = "./picture/gene_heatmap_pac5.pdf")

# 3.差异分析 limma C-PAC5   ----
rm(list=ls())
library(limma)
dif_data = read.table(file = './file/edata_rmdup.tab',header = TRUE, sep = "\t", quote = "",stringsAsFactors=FALSE)
dif_data = dif_data[,c(1:3,7:9)]
dif_data<-dif_data[rowSums(dif_data)>0,]
group_list=c(rep('Control',3),rep('PAC5',3))
group_list <- factor(group_list,levels = c("Control","PAC5"))
#表达矩阵数据校正
boxplot(dif_data,outline=FALSE, notch=T,col=group_list, las=2)
exprSet=normalizeBetweenArrays(dif_data)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
#判断数据是否需要转换
exprSet <- log2(exprSet+1)
#差异分析：
dat <- exprSet
design=model.matrix(~factor( group_list ))
fit=lmFit(dat,design)
fit=eBayes(fit)
options(digits = 4)

bp=function(g){
  library(ggpubr)
  df=data.frame(gene=g,stage=group_list)
  p <- ggboxplot(df, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means()
}

deg=topTable(fit,coef=2,adjust='BH',number = Inf)
deg$gene<-rownames(deg)
head(deg) 

deg$change = as.factor(ifelse(deg$P.Value < 0.05,
                                ifelse( deg$logFC > 0 ,'UP','DOWN' ),
                                'NOT'))

deg$significance = as.factor(ifelse( deg$P.Value  <= 0.05,
                                                 ifelse( deg$P.Value> 0.01, 'TRUE', 'M_TRUE' ),
                                                 'FALSE' ))
table(deg$change)
# DOWN   NOT    UP 
# 568 14383   962
sum(deg$P.Value < 0.05)
# 1530
table(deg$significance)
# FALSE M_TRUE   TRUE 
# 14383    367   1163
dif_data$gene<-rownames(dif_data)
deg_matrix<-merge(dif_data,deg,by='gene')
write.table(deg_matrix, file = "./file/C_PAC5_deg_limma.tab", quote = FALSE,sep="\t",row.names = FALSE)  #用于后续分析的tab

# 4.volcano map_up down   ----
rm(list=ls())
library("ggplot2")
library("ggrepel")
nrDEG = read.table('./file/C_PAC5_deg_limma.tab',sep='\t',header=TRUE,stringsAsFactors=FALSE)
nrDEG$change = factor(nrDEG$change,levels = c('UP', 'DOWN','NOT'))

this_tile <- paste0( 'The number of up genes is ', nrow(nrDEG[ nrDEG$change == 'UP', ] ),"\n",
                     'The number of down genes is ', nrow(nrDEG[ nrDEG$change == 'DOWN', ] ))

upgene = nrDEG[nrDEG$change == 'UP',];upgene <- upgene[order(upgene$P.Value),]
top10_up = upgene$gene[1:10]

downgene = nrDEG[nrDEG$change == 'DOWN',];downgene <- downgene[order(downgene$P.Value),]
top10_down =downgene$gene[1:10]

volcano <-ggplot(data = nrDEG, aes( x = logFC, y = -log10(P.Value), color = change)) +
  scale_x_continuous(limits = c(-1.5,1.5),breaks = seq(-3,3,1))+
  scale_y_continuous(limits = c(0,4.5),breaks = seq(0,6,1))+
  geom_point( alpha = 0.4, size = 1.75) +  
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+  #y轴分界线
  theme_classic() + theme(legend.position = "none") +
  xlab( "log2 fold change" ) + ylab( "-log10 p-value" ) +
  labs(title = this_tile) + 
  theme_bw()+theme(plot.title = element_text(size = 18, hjust = 0.5),
                   panel.grid=element_blank(),
                   panel.border=element_blank(),
                   axis.line=element_line(size=1,colour='black'),
                   axis.ticks=element_line(size=1,colour='black'),
                   axis.text=element_text(size=16,colour='black'),
                   axis.text.x=element_text(angle=0,hjust=1,vjust=1),
                   axis.title=element_text(size=18))+
  scale_colour_manual( values = c('red','blue','black')) +
  geom_text_repel(data=subset(nrDEG, gene %in% top10_up | gene %in% top10_down ), aes(label=gene),col="black",alpha = 1)
print(volcano)
plotfile=paste('./picture/volcano','_C','_PAC5_limma.pdf',sep='')
ggsave(plotfile, plot=volcano, dpi = 600,width = 10, height = 8)

# 5.GSEA----
rm(list=ls())
library( "clusterProfiler")
library("org.Hs.eg.db")
library(ggpubr)
load(file = "/home/devdata/tongxin/enrichment_annotation/human_enrichment.RData")
Con_PAC5 = read.table('./file/C_PAC5_deg_limma.tab',sep='\t',header=TRUE,stringsAsFactors=FALSE)
Con_PAC5$FC<-2^(Con_PAC5$logFC)

##### go BP
term2gene <- gobp[,c(1,5)]
term2name <- gobp[,c(1,2)]
Con_PAC5 <- arrange(Con_PAC5, desc(FC))
glist<-Con_PAC5$FC
names(glist) <- Con_PAC5$gene
gsea <- GSEA(glist, TERM2GENE=term2gene, TERM2NAME=term2name,verbose=FALSE, minGSSize = 5,
             pvalueCutoff = 1)
gsea_result<-gsea@result
write.table(gsea_result, file = "./enrichment/control_pac5_gsea_BP.tab", quote = FALSE,sep="\t",row.names = FALSE)

#######top20 BP
# rm(list = ls())
gsea_result_0 <- read.delim(file = "/data/dsx/xuming_20.09.23/enrichment/control_pac5_gsea_BP.tab",
                          header = T,sep = "\t",stringsAsFactors = F)
gsea_result0 <- read.table(file = "/data/dsx/hnRNPA2B1_2021/enrichment/control_pac5_gsea_BP.tab",
                          header = T,sep = "\t",stringsAsFactors = F)
# gsea_result1 <- read.table(file = "/data/dsx/hnRNPA2B1_2021/enrichment/control_pac5_gsea_BP0.tab",
#                            header = T,sep = "\t",stringsAsFactors = F)
# 

df <- gsea_result_0[order(gsea_result_0$pvalue),][1:20,]
df$Description <- factor(df$Description,levels = rev(unique(df$Description)))
df$tag <- as.factor(ifelse( df$NES > 0, 'Activated', 'Suppressed')) 
breaks<-pretty(range(0,-log10(df$pvalue)), 6)
maximum<- breaks[length(breaks)]
p1 <- ggplot(df,aes(x=-log10(pvalue),y=Description,size=-log10(pvalue),color=NES,ylab=''))+
  geom_point()+
  scale_color_gradient(low = "blue", high = "red") +
  scale_x_continuous(limits = c(0,maximum+1),breaks = breaks)+
  theme_bw()+
  scale_size_continuous(range=c(3,10))+
  labs(y='',x='-log10(pvalue)',title='')+
  facet_grid(cols = vars(tag))+
  theme() +
  geom_vline(xintercept = -log10(0.05),linetype =2,colour = 'black')+
  theme(plot.title=element_text(size=22,color="black",hjust = 0.5),
        axis.text=element_text(size=15, color="black"),
        axis.text.x=element_text(angle = 0,hjust = 0.5),
        axis.title=element_text(size=20),
        axis.line=element_line(size=1,colour='black'),
        axis.ticks=element_line(size=1,colour='black'),
        strip.text = element_text(size = 18),
        legend.text=element_text(size=16),legend.title = element_text(size=20))
p1
ggsave(plot = p1,'./picture/con_pac5_limma_GSEA_goBP.pdf', width = 16, height = 8,dpi = 600)

#######interferon related BP
df <- gsea_result[grep('interferon',gsea_result$Description),]

df$Description <- factor(df$Description,levels = rev(unique(df$Description)))
df$tag <- as.factor(ifelse( df$NES > 0, 'Activated', 'Suppressed')) 

breaks<-seq(from=0, to=3, length.out=7)
maximum<- breaks[length(breaks)]
p1 <- ggplot(df,aes(x=-log10(pvalue),y=Description,size=-log10(pvalue),color=NES,ylab=''))+
  geom_point()+
  scale_color_gradient(low = "blue", high = "red") +
  scale_x_continuous(limits = c(0,maximum),breaks = breaks)+
  theme_bw()+
  scale_size_continuous(range=c(1,8))+
  labs(y='',x='-log10(pvalue)',title='')+
  facet_grid(cols = vars(tag))+
  theme() +
  geom_vline(xintercept = -log10(0.05),linetype =2,colour = 'black')+
  theme(plot.title=element_text(size=22,color="black",hjust = 0.5),
        axis.text=element_text(size=12, color="black"),
        axis.text.x=element_text(angle = 0,hjust = 0.5),
        axis.title=element_text(size=20),
        axis.line=element_line(size=1,colour='black'),
        axis.ticks=element_line(size=1,colour='black'),
        strip.text = element_text(size = 18),
        legend.text=element_text(size=16),legend.title = element_text(size=20))
p1
ggsave(plot = p1,'./picture/con_pac5_limma_interferon_GSEA_goBP.pdf', width = 14, height = 8,dpi = 600)

########################GSEA end

# 6.protein----
rm(list = ls())
library(ggrepel)
library(reshape2)
library(ggplot2)
pro_data <- read.delim(file = './protein/all_sample.xls',header = T, quote = "",sep = "\t",stringsAsFactors=FALSE)
pro_name <- read.delim(file = './protein/uniprot_name_map_ok.tab',header = T, quote = "",sep = "\t",stringsAsFactors=FALSE)
colnames(pro_name)<-c('Protein',"gene_name")
colnames(pro_data)
pro_data<-merge(pro_data,pro_name,by="Protein",all=T)
pro_data1<-pro_data[! is.na(pro_data$gene_name),]
pro_data2<-pro_data[is.na(pro_data$gene_name),]
pro_data2$gene_name<-pro_data2$Protein
pro_data<-rbind(pro_data1,pro_data2)

pro_data$negtive[pro_data$negtive=='-']<-0
pro_data$positive[pro_data$positive=='-']<-0

pro_data$negtive<-as.numeric(pro_data$negtive)
pro_data$positive<-as.numeric(pro_data$positive)
write.table(pro_data, file = "./file/pro_pos_neg_data_all.tab", quote = FALSE,sep="\t",row.names = FALSE)

# 6.1 geom_point ------------------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(ggplot2)
library("ggrepel")
pro_pos_neg <- read.table("/data/dsx/hnRNPA2B1_2021/file/pro_pos_neg_data_all.tab",sep = "\t",header = T,stringsAsFactors = F)

pro_pos_neg$log.neg <- log10(pro_pos_neg$negtive + 1)
pro_pos_neg$log.pos <- log10(pro_pos_neg$positive + 1)
pro_pos_neg$pos_neg <- pro_pos_neg$positive - pro_pos_neg$negtive
pro_pos_neg$log.pos_neg <- log((pro_pos_neg$positive - pro_pos_neg$negtive) + 1)

pro_pos_neg$order <- "low"

# only_pos
only_pos <- pro_pos_neg[pro_pos_neg$log.neg == 0,]
only_pos_top10 <- only_pos[order(only_pos$positive,decreasing = T),"gene_name"][1:10]
pro_pos_neg[pro_pos_neg$Protein %in% c("A0A024R6I7","P22626","Q9NXW2"),"gene_name"] <- c("A0A024R6I7","P22626","Q9NXW2")

pro_pos_neg[pro_pos_neg$gene_name %in% only_pos_top10 ,"order"] <- "only_pos"

# pos != 0,pos-neg
pos_neg <- pro_pos_neg[pro_pos_neg$negtive != 0,]
pos_neg_top10 <- pos_neg[order(pos_neg$pos_neg,decreasing = T),"gene_name"][1:10]
pro_pos_neg[pro_pos_neg$gene_name %in% pos_neg_top10 ,"order"] <- "pos_neg"

genebothzero<-pro_pos_neg$gene_name[pro_pos_neg$log.neg==0 & pro_pos_neg$log.pos==0]
pro_pos_neg[pro_pos_neg$gene_name %in% genebothzero ,"order"] <- "bothzero"
pro_pos_neg<-pro_pos_neg[pro_pos_neg$order %in% c("only_pos","pos_neg","low"),]

pro_pos_neg$order <- factor(pro_pos_neg$order,levels = c("only_pos","pos_neg","low"))
table(pro_pos_neg$order)
p <- ggplot(data = pro_pos_neg, aes( x = log.neg, y = log.pos, color = order)) +
  geom_point( alpha = 0.6, size = 2) +  
  scale_y_continuous(limits = c(-0.5,9),breaks = seq(0,9,9))+
  scale_x_continuous(limits = c(-0.5,9),breaks = seq(0,9,1))+
  theme_bw() + theme(legend.position = "none") +
  xlab("log10(Abundance) in NC" ) + ylab( "log10(Abundance) in Probe" ) +
  theme( axis.text = element_text(size = 18,colour = "black"),
         axis.title = element_text(size = 20,colour = "black"),
         panel.grid=element_blank(),
         panel.border=element_blank(),
         axis.line=element_line(size=0.8,colour='black'),
         axis.ticks = element_line(size = 0.8)) +
  scale_colour_manual( values = c("only_pos" = 'red',"pos_neg" = 'blue',"low" = 'grey')) +
  geom_text_repel(data=subset(pro_pos_neg, gene_name %in% only_pos_top10 |  gene_name %in% pos_neg_top10),
                  aes(label=gene_name),col="black",alpha = 1,size = 5 )
print(p)
ggsave(plot = p,'/data/dsx/hnRNPA2B1_2021/picture/bubble_protein_top10.pdf', width = 10, height = 10,dpi = 600)

# 6.2 bar_plot and enrich-----
rm(list = ls())
library(ggplot2)
setwd("/data/dsx/hnRNPA2B1_2021")
pro_pos_neg <- read.table("/data/dsx/hnRNPA2B1_2021/file/pro_pos_neg_data_all.tab",sep = "\t",header = T,stringsAsFactors = F)
pro_pos_neg$negtive

pro_pos_neg$count <- "Both"
pro_pos_neg[pro_pos_neg$Protein %in% 
              pro_pos_neg[pro_pos_neg$negtive == 0 & pro_pos_neg$positive == 0, "Protein"],"count"] <- "0"
pro_pos_neg[pro_pos_neg$Protein %in% 
              pro_pos_neg[pro_pos_neg$negtive == 0 & pro_pos_neg$positive != 0, "Protein"],"count"] <- "PAC5"
pro_pos_neg[pro_pos_neg$Protein %in% 
              pro_pos_neg[pro_pos_neg$positive == 0 & pro_pos_neg$negtive != 0, "Protein"],"count"] <- "NC"
pro_pos_neg_1 <- pro_pos_neg[pro_pos_neg$count != 0,]

bar_file <- as.data.frame(table(pro_pos_neg_1$count))
colnames(bar_file) <- c("cate","counts")
breaks<-pretty(range(0,bar_file$counts), 6)
maximum<- breaks[length(breaks)]

p <- ggplot(data = bar_file, aes( x = cate, y = counts, fill = cate)) +
  geom_bar(stat="identity",width=0.5) +  
  geom_text(aes(label = counts, vjust = -0.8, hjust = 0.5, color = cate), show.legend = TRUE,size = 6) + 
  theme_bw() + 
  xlab(" " ) + ylab( "Protein Number" ) +
  scale_y_continuous(limits = c(0,maximum),breaks = breaks)+
  # scale_x_discrete(labels=c("Both","NC","PAC5"))+
  theme( axis.text = element_text(size = 16,colour = "black"),
         axis.title.y = element_text(size=18),
         legend.text = element_text(size = 14,colour = "black"),
         legend.title = element_text(size = 17,colour = "black"),
         panel.grid=element_blank(),
         panel.border=element_blank(),
         axis.line=element_line(size=0.8,colour='black'),
         axis.ticks = element_line(size = 0.8),legend.position="none")
p
ggsave(plot = p,'/data/dsx/hnRNPA2B1_2021/picture/Protein_Number_barplot.pdf', width = 10, height = 10,dpi = 600)

# 6.3 enrich
library( "clusterProfiler")
library("org.Hs.eg.db")
### load(.RData)
load(file = "/home/devdata/tongxin/enrichment_annotation/human_enrichment.RData")
##### go ###
cate_file = pro_pos_neg_1
### only_pos ###
gene_list = pro_pos_neg_1$gene_name[pro_pos_neg_1$count=="PAC5"]
## BP
term2gene <- gobp[,c(1,5)]
term2name <- gobp[,c(1,2)]
ego_BP <- enricher(gene = gene_list,
                   pvalueCutoff = 0.05,pAdjustMethod = "BH",
                   minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                   TERM2GENE = term2gene,TERM2NAME = term2name)
ego_BP_result<-as.data.frame(ego_BP@result)

## KEGG
term2gene <- kegg[,c(1,5)]
term2name <- kegg[,c(1,2)]
ekegg <- enricher(gene = gene_list,
                  pvalueCutoff = 0.05,pAdjustMethod = "BH",
                  minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                  TERM2GENE = term2gene,TERM2NAME = term2name)
ekegg_dif<-as.data.frame(ekegg@result) 

write.table(ekegg_dif, file = "/data/dsx/hnRNPA2B1_2021/file/PAC5_only_pos_kegg_enrich.tab", quote = FALSE,sep="\t",row.names = FALSE)
write.table(ego_BP_result, file = "/data/dsx/hnRNPA2B1_2021/file/PAC5_only_pos_BP_enrich.tab", quote = FALSE,sep="\t",row.names = FALSE)

### only_neg ###
gene_list = pro_pos_neg_1$gene_name[pro_pos_neg_1$count=="NC"]
## BP
term2gene <- gobp[,c(1,5)]
term2name <- gobp[,c(1,2)]
ego_BP <- enricher(gene = gene_list,
                   pvalueCutoff = 0.05,pAdjustMethod = "BH",
                   minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                   TERM2GENE = term2gene,TERM2NAME = term2name)
ego_BP_result<-as.data.frame(ego_BP@result)


## KEGG
term2gene <- kegg[,c(1,5)]
term2name <- kegg[,c(1,2)]
ekegg <- enricher(gene = gene_list,
                  pvalueCutoff = 0.05,pAdjustMethod = "BH",
                  minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                  TERM2GENE = term2gene,TERM2NAME = term2name)
ekegg_dif<-as.data.frame(ekegg@result) 

write.table(ekegg_dif, file = "/data/dsx/hnRNPA2B1_2021/file/NC_only_pos_kegg_enrich.tab", quote = FALSE,sep="\t",row.names = FALSE)
write.table(ego_BP_result, file = "/data/dsx/hnRNPA2B1_2021/file/NC_only_pos_BP_enrich.tab", quote = FALSE,sep="\t",row.names = FALSE)

#### plot_bubble
rm(list=ls())
library(ggplot2)
ekegg_only_pos <- read.delim(file = '/data/dsx/hnRNPA2B1_2021/file/PAC5_only_pos_kegg_enrich.tab',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
ego_BP_only_pos <- read.delim(file = '/data/dsx/hnRNPA2B1_2021/file/PAC5_only_pos_BP_enrich.tab',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
ekegg_only_neg <- read.delim(file = '/data/dsx/hnRNPA2B1_2021/file/NC_only_pos_kegg_enrich.tab',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
ego_BP_only_neg <- read.delim(file = '/data/dsx/hnRNPA2B1_2021/file/NC_only_pos_BP_enrich.tab',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)

ekegg_only_pos$Group <- "PAC5"
ekegg_only_pos_top10 <- ekegg_only_pos[order(ekegg_only_pos$pvalue),][1 :10,]
ego_BP_only_pos$Group <- "PAC5"
ego_BP_only_pos_top10 <- ego_BP_only_pos[order(ego_BP_only_pos$pvalue),][1:10,]
ekegg_only_neg$Group <- "NC"
ekegg_only_neg_top10 <- ekegg_only_neg[order(ekegg_only_neg$pvalue),][1 :10,]
ego_BP_only_neg$Group <- "NC"
ego_BP_only_neg_top10 <- ego_BP_only_neg[order(ego_BP_only_neg$pvalue),][1:10,]

#####  kegg 
KEGG <- rbind(ekegg_only_pos_top10,ekegg_only_neg_top10)[,2]
KEGG_pos <- ekegg_only_pos[ekegg_only_pos$Description %in% KEGG,]
KEGG_neg <- ekegg_only_neg[ekegg_only_neg$Description %in% KEGG,]
KEGG_pos_neg <- rbind(KEGG_pos,KEGG_neg)
KEGG_pos_neg <-rbind(ekegg_only_pos_top10,ekegg_only_neg_top10)
table(KEGG_pos_neg$Group)
# NC PAC5 
# 10   10
KEGG_pos_neg$log.p <- -log10(KEGG_pos_neg$pvalue)

KEGG_pos_neg$Group <- factor(KEGG_pos_neg$Group,levels = c('PAC5','NC'))
KEGG_pos_neg$Description <- factor(KEGG_pos_neg$Description,levels = rev(unique(KEGG_pos_neg$Description)))
colnames(KEGG_pos_neg)
breaks<-pretty(range(0,KEGG_pos_neg$Count), 5)
maximum<- breaks[length(breaks)]
p1 <- ggplot(KEGG_pos_neg,aes(x=log.p,y=Description,size=Count,color=Group,ylab=''))+
  geom_point()+
  theme_bw()+
  scale_size_continuous(range=c(2,8),breaks=breaks)+
  scale_x_continuous(limits = c(1,6),breaks = c(1:6))+
  labs(y='',x='-log(pvalue)',title='KEGG',color = "Group")+
  theme() +
  geom_vline(xintercept = -log10(0.05),linetype =2)+
  theme(plot.title=element_text(size=23,hjust = 0.5),
        strip.text = element_text(size = rel(1.5)),
        axis.text.x =element_text(size=18,colour = "black"),
        axis.text.y =element_text(size=14,colour = "black"),
        axis.title.x = element_text(size=18,colour = "black"),
        legend.text=element_text(size=14),
        legend.title = element_text(size=18)) 
p1
ggsave(plot = p1,'/data/dsx/hnRNPA2B1_2021/picture/Protein_KEGG_pos_neg.pdf', width = 12, height = 10,dpi = 600)

#####  bp ego_BP_only_pos
BP <- rbind(ego_BP_only_pos_top10,ego_BP_only_neg_top10)[,2]
BP_pos <- ego_BP_only_pos[ego_BP_only_pos$Description %in% BP,]
BP_neg <- ego_BP_only_neg[ego_BP_only_neg$Description %in% BP,]
BP_pos_neg <- rbind(BP_pos,BP_neg)
BP_pos_neg <-  rbind(ego_BP_only_pos_top10,ego_BP_only_neg_top10)

table(BP_pos_neg$Group)
# NC PAC5 
# 10   10
BP_pos_neg$log.p <- -log10(BP_pos_neg$pvalue)

BP_pos_neg$Group <- factor(BP_pos_neg$Group,levels = c('PAC5','NC'))
BP_pos_neg$Description <- factor(BP_pos_neg$Description,levels = rev(unique(BP_pos_neg$Description)))
colnames(BP_pos_neg)

breaks<-pretty(range(0,BP_pos_neg$Count), 6)
maximum<- breaks[length(breaks)]

p1 <- ggplot(BP_pos_neg,aes(x=log.p,y=Description,size=Count,color=Group,ylab=''))+
  geom_point()+
  scale_size_area(name="genecounts")+
  theme_bw()+
  scale_size_continuous(range=c(4,10),breaks =breaks)+
  scale_x_continuous(limits = c(1,6),breaks = c(1:6))+
  labs(y='',x='-log(pvalue)',title='GO biological process',color = "Group")+
  theme() +
  geom_vline(xintercept = -log10(0.05),linetype =2)+
  theme(plot.title=element_text(size=23,hjust = 0.5),
        strip.text = element_text(size = rel(1.5)),
        axis.text.x =element_text(size=18,colour = "black"),
        axis.text.y =element_text(size=14,colour = "black"),
        axis.title.x = element_text(size=18,colour = "black"),
        legend.text=element_text(size=14),
        legend.title = element_text(size=18))
p1
ggsave(plot = p1,'/data/dsx/hnRNPA2B1_2021/picture/Protein_BP_pos_neg.pdf', width = 12, height = 10,dpi = 600)

