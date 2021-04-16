setwd("/data/dsx/xuming_20.09.23")
getwd ()
rm(list = ls())
edata <- read.delim(file = './file/fpkm.anno.xls',sep='\t',header=TRUE, stringsAsFactors=FALSE) #31553
protein_coding<-edata$name[edata$Biotype=='protein_coding']
edata <- edata[,c(1:10)]
rownames(edata) <- edata[,1]
edata <- edata[,-1]

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
edata <- apply(edata,2,fpkmToTpm)
# edata[1:3,]
# colSums(edata)
####protein coding
edata<-edata[protein_coding,]
edata<-edata[rowSums(edata)>0,]
write.table(edata, file = "./file/edata_TPM.tab", quote = FALSE,sep="\t",row.names = TRUE)

##############PCAcheck##############PCAcheck
library(ggord)
library(ggplot2)
edata<-edata[rowSums(edata)>0,]
pca_group=factor(c(rep('Control',3),rep('PAC2',3),rep('PAC5',3)))
edata.pca <- prcomp(t(edata), scale. = TRUE)
p <- ggord(edata.pca, grp_in = pca_group, arrow=0, vec_ext =0,txt=NULL,cols=c('blue','#2ca25f','red'),
           poly = T,polylntyp='dashed', # twodash, solid, longdash, dotted, dotdash, dashed, blank
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
        # legend.position = 'right',
        legend.position = 'top',
        legend.key.height = unit(0.6,'cm'),#定义图例中色块的高度
        legend.text = element_text(face = 'bold',color = 'black',size = 10)
  )
# + 
#   annotate(geom = 'segment', y = Inf,  x = -Inf, yend = Inf, xend = Inf, color = 'black', size = 1)+ 
#   annotate(geom = 'segment', y = -Inf, x = Inf, yend = Inf, xend = Inf, color = 'black', size = 1)
p
ggsave('./picture/group123_PCA_v2.pdf', plot=p, dpi = 600,width = 6, height = 6)

####only PAC5
edata<-edata[,c(1:3,7:9)]
edata<-edata[rowSums(edata)>0,]
pca_group=factor(c(rep('Control',3),rep('PAC5',3)))
edata.pca <- prcomp(t(edata), scale. = TRUE)
p <- ggord(edata.pca, grp_in = pca_group, arrow=0, vec_ext =0,txt=NULL,cols=c('blue','red'),
           poly = T,polylntyp='dashed', # twodash, solid, longdash, dotted, dotdash, dashed, blank
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
        # legend.position = 'right',
        legend.position = 'top',
        legend.key.height = unit(0.6,'cm'),#定义图例中色块的高度
        legend.text = element_text(face = 'bold',color = 'black',size = 10)
  )
p
ggsave('./picture/group_PAC5_PCA_v2.pdf', plot=p, dpi = 600,width = 6, height = 6)

####annotation  gene
rm(list=ls())
library( "clusterProfiler")
library("org.Hs.eg.db")
edata = read.table(file = './file/edata_TPM.tab',sep='\t',header=TRUE, stringsAsFactors=FALSE)
edata$ENSEMBL<-rownames(edata)
columns(org.Hs.eg.db)
geneinfo = select(org.Hs.eg.db, keys=keys(org.Hs.eg.db), columns = c('ENSEMBL','ENTREZID',"SYMBOL"))
all_data<-merge(x=edata,y=geneinfo,by="ENSEMBL")
all_data<-all_data[,c(-1,-11)]
edata_rmdup <- aggregate(all_data[,-10],by = list(all_data$SYMBOL), FUN = mean) #16203
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
library(pheatmap)
pheatmap(Z_score,border_color="#C5C5C5",color = colorRampPalette(c("blue", "white", "red"))(50),
         clustering_distance_cols = "correlation",clustering_distance_rows = "correlation",show_rownames=F,annotation_col=annotation_col,
         annotation_colors=ann_colors,fontsize=9,filename = "./picture/gene_heatmap_pac5.pdf")

# 5.差异分析  limma   C--PAC5   ----
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
# 639 17201  1001 
sum(deg$P.Value < 0.05)
# 1640
table(deg$significance)
# FALSE M_TRUE   TRUE 
# 17201    391   1249
dif_data$gene<-rownames(dif_data)
deg_matrix<-merge(dif_data,deg,by='gene')
write.table(deg_matrix, file = "./file/C_PAC5_deg_limma.tab", quote = FALSE,sep="\t",row.names = FALSE)  #用于后续分析的tab


# 5.差异分析  t.test   C--PAC5   ----
rm(list=ls())
library(limma)
dif_data = read.table(file = './file/edata_rmdup.tab',header = TRUE, sep = "\t", quote = "",stringsAsFactors=FALSE)
dif_data = dif_data[,c(1:3,7:9)]
dif_data<-dif_data[rowSums(dif_data)>0,]
pvalue<-c();tstat<-c();meansControl<-c();meansCase<-c();FC<-c();log2FC<-c()
Control_list<-1:3
Case_list<-4:6

for (i in 1:nrow(dif_data)){
  #t.test
  Control_value<-dif_data[i,Control_list]
  Case_value<-dif_data[i,Case_list]
  if (sum(Control_value)==0){
    Control_value<-rnorm(3,mean=0.0001,sd=0.00001)
    dif_data[i,Control_list]<-Control_value
  }
  if (sum(Case_value)==0){
    Case_value<-rnorm(3,mean=0.0001,sd=0.00001)
    dif_data[i,Case_list]<-Case_value
  }
  # print(Control_value)
  # print(Case_value)
  result<-t.test(as.numeric(Case_value), as.numeric(Control_value), paired=FALSE);
  pvalue[i]<-result[["p.value"]]
  tstat[i]<-result[["statistic"]][["t"]]
  meansControl[i]<-mean(as.numeric(Control_value))
  meansCase[i]<-mean(as.numeric(Case_value))
  FC[i]<-mean(as.numeric(Case_value))/mean(as.numeric(Control_value))
  log2FC[i]<-log2(FC[i])
}

p_bf = p.adjust(pvalue, method = "bonferroni")
dif_data1 = data.frame(rownames(dif_data),pvalue,tstat,meansControl,meansCase,FC,log2FC,p_bf,stringsAsFactors=FALSE)
colnames(dif_data1)[1] = 'SYMBOL';rownames(dif_data1) = dif_data1$SYMBOL
dif_data$SYMBOL<-rownames(dif_data)

#将表达矩阵同t.test结果合并
dif_data_all = merge(dif_data[,c(1:7)],dif_data1,by="SYMBOL") 
rownames(dif_data_all) = dif_data_all$SYMBOL

#根据P值和FC添加change_V1???
dif_data_all$change_V1 = as.factor( ifelse( dif_data_all$pvalue < 0.05,
                                            ifelse( dif_data_all$log2FC > 0 ,'UP','DOWN' ),
                                            'NOT' ))
dif_data_all$significance = as.factor( ifelse( dif_data_all$pvalue  <= 0.05,
                                               ifelse( dif_data_all$pvalue > 0.01, 'TRUE', 'M_TRUE' ),
                                               'FALSE' )) 

table(dif_data_all$change_V1)
# DOWN   NOT    UP 
# 525 17704   612 
sum(dif_data_all$pvalue < 0.05)
# 1137
table(dif_data_all$significance)
# FALSE M_TRUE   TRUE 
# 17704    221    91
write.table(dif_data_all, file = "./file/C_PAC5_deg_ttest.tab", quote = FALSE,sep="\t",row.names = FALSE)  #用于后续分析的tab


# volcanomap_up down   ----
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


# heatmap  ----
rm(list = ls())
library(pheatmap)
nrDEG = read.table('./file/C_PAC5_deg_limma.tab',sep='\t',header=TRUE,stringsAsFactors=FALSE)

dif_data_sig = nrDEG[nrDEG$significance == 'M_TRUE',]

row.names(dif_data_sig) = dif_data_sig$gene

upgene = dif_data_sig[dif_data_sig$change == 'UP',];upgene <- upgene[order(upgene$P.Value),]
rownames(upgene) = upgene$gene;top10_up = rownames(upgene)[1:10]

downgene = dif_data_sig[dif_data_sig$change == 'DOWN',];downgene <- downgene[order(downgene$P.Value),]
rownames(downgene) = downgene$gene;top10_down = rownames(downgene)[1:10]

dif_gene_top20 <- c(top10_up,top10_down)

dif_data_top20 <- dif_data_sig[dif_gene_top20,]

dif_data_top20 <- dif_data_top20[,-1] 
dif_data_top20 <- dif_data_top20[,1:6]

Z_score <- (dif_data_top20 - apply(dif_data_top20, 1, mean)) / apply(dif_data_top20, 1, sd)
annotation_col <- data.frame(Group = c(rep('Control',3),rep('PAC5',3)))
rownames(annotation_col) <- colnames(dif_data_top20)
# pheatmap(Z_score)
ancols = c('blue','red')
names(ancols) = c("Control","PAC5")
ann_colors <- list(Group=ancols)
library(pheatmap)
pheatmap(Z_score,border_color="#C5C5C5",
         color = colorRampPalette(c("blue", "white", "red"))(50),
         angle_col = 0,cluster_cols = F,cluster_rows = F,
         clustering_distance_cols = "correlation",
         clustering_distance_rows = "correlation",
         show_rownames=T,annotation_col=annotation_col,
         annotation_colors=ann_colors,fontsize=10,filename = "./picture/PAC5_heatmap_limma_top20.pdf",
         width = 10,height = 8)


## 差异基因热图 -------------------------------------------------------------
rm(list=ls())
library(pheatmap)
nrDEG = read.table('./file/C_PAC5_deg_limma.tab',sep='\t',header=TRUE,stringsAsFactors=FALSE)
dif_data_sig = nrDEG[nrDEG$significance == 'M_TRUE',]
row.names(dif_data_sig) = dif_data_sig$gene
dif_data<-dif_data_sig[,2:7]
Z_score <- (dif_data - apply(dif_data, 1, mean)) / apply(dif_data, 1, sd)
annotation_col <- data.frame(Group = c(rep('Control',3),rep('PAC5',3)))
rownames(annotation_col) <- colnames(dif_data)
# pheatmap(Z_score)
ancols = c('blue','red')
names(ancols) = c("Control","PAC5")
ann_colors <- list(Group=ancols)
library(pheatmap)
pheatmap(Z_score,border_color="#C5C5C5",color = colorRampPalette(c("blue", "white", "red"))(50),angle_col = 0,
         clustering_distance_cols = "correlation",clustering_distance_rows = "correlation",show_rownames=F,annotation_col=annotation_col,
         annotation_colors=ann_colors,fontsize=10,filename = "./picture/PAC5_heatmap_difgene_limma.pdf",
         width = 10,height = 8)

# 6.1 GO_BP/KEGG control_pac5 ----
getwd()
rm(list=ls())
library( "clusterProfiler")
library("org.Hs.eg.db")
### load(.RData)
load(file = "/home/devdata/tongxin/enrichment_annotation/human_enrichment.RData")
### go
Con_PAC5 = read.table('./file/C_PAC5_deg_limma.tab',sep='\t',header=TRUE,stringsAsFactors=FALSE)
Con_PAC5$FC<-2^(Con_PAC5$logFC)
dif_data <- Con_PAC5[which(Con_PAC5$change != "NOT"),]
# dif_data <- Con_PAC5[which(Con_PAC5$P.Value<0.01),]
gene_list = dif_data$gene
## BP
term2gene <- gobp[,c(1,5)]
term2name <- gobp[,c(1,2)]
ego_BP <- enricher(gene = gene_list,
                   pvalueCutoff = 0.05,pAdjustMethod = "BH",
                   minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                   TERM2GENE = term2gene,TERM2NAME = term2name)
ego_BP_result<-as.data.frame(ego_BP@result)
Con_PAC5 <- arrange(Con_PAC5, desc(FC))
glist<-Con_PAC5$FC
names(glist) <- Con_PAC5$gene
cnetplot(ego_BP, foldChange=glist,showCategory = 10)
enrichMap(ego_BP, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
dotplot(ego_BP, showCategory=30)
## KEGG
term2gene <- kegg[,c(1,5)]
term2name <- kegg[,c(1,2)]
ekegg <- enricher(gene = gene_list,
                  pvalueCutoff = 0.05,pAdjustMethod = "BH",
                  minGSSize = 10,maxGSSize = 500,qvalueCutoff = 0.2,
                  TERM2GENE = term2gene,TERM2NAME = term2name)
ekegg_dif<-as.data.frame(ekegg@result) 

write.table(ekegg_dif, file = "./enrichment/control_pac5_kegg_enrich_dif.tab", quote = FALSE,sep="\t",row.names = FALSE)
write.table(ego_BP_result, file = "./enrichment/control_pac5_BP_enrich.tab", quote = FALSE,sep="\t",row.names = FALSE)

# 6.2 GO_BP/KEGG control_pac2 ----
getwd()
rm(list=ls())
library( "clusterProfiler")
library("org.Hs.eg.db")
### load(.RData)
load(file = "/home/devdata/tongxin/enrichment_annotation/human_enrichment.RData")
### go
Con_PAC2 = read.table(file = './file/tab/Con_PAC2.tab',header = TRUE, sep = "\t", quote = "",stringsAsFactors=FALSE)
dif_data <- Con_PAC2[which(Con_PAC2$change_V1 != "NOT"),]
gene_list = dif_data$SYMBOL
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

write.table(ekegg_dif, file = "./enrichment/control_pac2_kegg_enrich_dif.tab", quote = FALSE,sep="\t",row.names = FALSE)
write.table(ego_BP_result, file = "./enrichment/control_pac2_BP_enrich.tab", quote = FALSE,sep="\t",row.names = FALSE)


#####GSEA
# devtools::install_github("ToledoEM/msigdf")
# library(msigdf)
###https://guangchuangyu.github.io/cn/2018/11/msigdf_clusterprofiler/
rm(list=ls())
library( "clusterProfiler")
library("org.Hs.eg.db")
library(ggpubr)
load(file = "/home/devdata/tongxin/enrichment_annotation/human_enrichment.RData")

Con_PAC5 = read.table('./file/C_PAC5_deg_limma.tab',sep='\t',header=TRUE,stringsAsFactors=FALSE)

Con_PAC5 = read.table('./file/C_PAC5_deg_limma.tab',sep='\t',header=TRUE,stringsAsFactors=FALSE)
Con_PAC5$FC<-2^(Con_PAC5$logFC)

##### go BP
term2gene <- gobp[,c(1,5)]
term2name <- gobp[,c(1,2)]
Con_PAC5 <- arrange(Con_PAC5, desc(FC))
glist<-Con_PAC5$FC
names(glist) <- Con_PAC5$gene
gsea <- GSEA(glist, TERM2GENE=term2gene, TERM2NAME=term2name,verbose=FALSE, minGSSize = 5,pvalueCutoff = 1)
gsea_result<-gsea@result
write.table(gsea_result, file = "./enrichment/control_pac5_gsea_BP.tab", quote = FALSE,sep="\t",row.names = FALSE)

#######top20 BP
df <- gsea_result[order(gsea_result$pvalue),][1:20,]
# num1<-unlist(lapply(df$GeneRatio, FUN = function(x) {return(strsplit(x, split = "/")[[1]][1])}))
# df$GeneRatio<- as.numeric(num1)/1301
df$Description <- factor(df$Description,levels = rev(unique(df$Description)))
df$tag <- as.factor(ifelse( df$NES > 0, 'Activated', 'Suppressed')) 
breaks<-pretty(range(0,-log10(df$pvalue)), 6)
maximum<- breaks[length(breaks)]
p1 <- ggplot(df,aes(x=-log10(pvalue),y=Description,size=-log10(pvalue),color=NES,ylab=''))+
  geom_point()+
  # scale_color_gradient2(low = "#bdbdbd", mid = "#ffeda0", high = "#de2d26",midpoint = max(-log10(dt$pvalue))/2) +
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
        #legend.position = 'top',
        strip.text = element_text(size = 18),
        legend.text=element_text(size=16),legend.title = element_text(size=20))
p1
ggsave(plot = p1,'./picture/con_pac5_limma_GSEA_goBP.pdf', width = 16, height = 8,dpi = 600)

#######interferon related BP
df <- gsea_result[grep('interferon',gsea_result$Description),]
# load(file = "/home/devdata/zhuzn/enrichment_annotation/human_enrichment.RData")
# NObp<-gobp$GO[grep('NO',gobp$relationship)]
# df<-df[! df$ID %in% NObp, ]
# df<-df[df$pvalue< 0.05,]
df<-df[1:20,]
df$Description <- factor(df$Description,levels = rev(unique(df$Description)))
df$tag <- as.factor(ifelse( df$NES > 0, 'Activated', 'Suppressed')) 
# breaks<-pretty(range(0,-log10(df$pvalue)), 6)
# maximum<- breaks[length(breaks)]
breaks<-seq(from=0, to=3, length.out=7)
maximum<- breaks[length(breaks)]
p1 <- ggplot(df,aes(x=-log10(pvalue),y=Description,size=-log10(pvalue),color=NES,ylab=''))+
  geom_point()+
  # scale_color_gradient2(low = "#bdbdbd", mid = "#ffeda0", high = "#de2d26",midpoint = max(-log10(dt$pvalue))/2) +
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
        #legend.position = 'top',
        strip.text = element_text(size = 18),
        legend.text=element_text(size=16),legend.title = element_text(size=20))
p1
ggsave(plot = p1,'./picture/con_pac5_limma_interferon_GSEA_goBP.pdf', width = 14, height = 8,dpi = 600)

########################GSEA end

# 6.3 bubble-plot GO/KEGG----
rm(list=ls())
library(ggplot2)
ekegg_pac5 <- read.delim(file = './enrichment/control_pac5_kegg_enrich_dif.tab',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
ego_BP_pac5 <- read.delim(file = './enrichment/control_pac5_BP_enrich.tab',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
ekegg_pac2 <- read.delim(file = './enrichment/control_pac2_kegg_enrich_dif.tab',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
ego_BP_pac2 <- read.delim(file = './enrichment/control_pac2_BP_enrich.tab',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)

#######top20 BP
df <- ego_BP_pac5[order(ego_BP_pac5$pvalue),][1:20,]
# num1<-unlist(lapply(df$GeneRatio, FUN = function(x) {return(strsplit(x, split = "/")[[1]][1])}))
# df$GeneRatio<- as.numeric(num1)/1301
df$Description <- factor(df$Description,levels = rev(unique(df$Description)))

p1 <- ggplot(df,aes(x=-log10(pvalue),y=Description,size=Count,color=-log10(pvalue),ylab=''))+
  geom_point()+
  # scale_color_gradient2(low = "#bdbdbd", mid = "#ffeda0", high = "#de2d26",midpoint = max(-log10(dt$pvalue))/2) +
  scale_color_gradient(low = "blue", high = "red") +
  # scale_color_manual(values = c('#e41a1c','#377eb8'))+ 
  # scale_size_area(name="genecounts")+
  theme_bw()+
  scale_size_continuous(range=c(4,10))+
  labs(y='',x='-log10(pvalue)',title='')+
  theme() +
  geom_vline(xintercept = -log10(0.05),linetype =2,colour = 'black')+
  theme(plot.title=element_text(size=22,color="black",hjust = 0.5),
        axis.text=element_text(size=15, color="black"),
        axis.text.x=element_text(angle = 0,hjust = 0.5),
        axis.title=element_text(size=20),
        axis.line=element_line(size=1,colour='black'),
        axis.ticks=element_line(size=1,colour='black'),
        #legend.position = 'top',
        legend.text=element_text(size=16),legend.title = element_text(size=20))
p1
ggsave(plot = p1,'./picture/con_pac5_limma_goBP.pdf', width = 12, height = 8,dpi = 600)

#######interferon related BP
df <- ego_BP_pac5[grep('interferon',ego_BP_pac5$Description),]
# num1<-unlist(lapply(df$GeneRatio, FUN = function(x) {return(strsplit(x, split = "/")[[1]][1])}))
# df$GeneRatio<- as.numeric(num1)/1301
df$Description <- factor(df$Description,levels = rev(unique(df$Description)))

p1 <- ggplot(df,aes(x=-log10(pvalue),y=Description,size=Count,color=-log10(pvalue),ylab=''))+
  geom_point()+
  # scale_color_gradient2(low = "#bdbdbd", mid = "#ffeda0", high = "#de2d26",midpoint = max(-log10(dt$pvalue))/2) +
  scale_color_gradient(low = "blue", high = "red") +
  # scale_color_manual(values = c('#e41a1c','#377eb8'))+ 
  # scale_size_area(name="genecounts")+
  theme_bw()+
  scale_size_continuous(range=c(4,10))+
  labs(y='',x='-log10(pvalue)',title='')+
  theme() +
  geom_vline(xintercept = -log10(0.05),linetype =2,colour = 'black')+
  theme(plot.title=element_text(size=22,color="black",hjust = 0.5),
        axis.text=element_text(size=15, color="black"),
        axis.text.x=element_text(angle = 0,hjust = 0.5),
        axis.title=element_text(size=20),
        axis.line=element_line(size=1,colour='black'),
        axis.ticks=element_line(size=1,colour='black'),
        #legend.position = 'top',
        legend.text=element_text(size=16),legend.title = element_text(size=20))
p1
ggsave(plot = p1,'./picture/con_pac5_limma_interferon_goBP.pdf', width = 14, height = 8,dpi = 600)

#######top20 KEGG
df <- ekegg_pac5[order(ekegg_pac5$pvalue),][1:20,]
# num1<-unlist(lapply(df$GeneRatio, FUN = function(x) {return(strsplit(x, split = "/")[[1]][1])}))
# df$GeneRatio<- as.numeric(num1)/1301
df$Description <- factor(df$Description,levels = rev(unique(df$Description)))

p1 <- ggplot(df,aes(x=-log10(pvalue),y=Description,size=Count,color=-log10(pvalue),ylab=''))+
  geom_point()+
  # scale_color_gradient2(low = "#bdbdbd", mid = "#ffeda0", high = "#de2d26",midpoint = max(-log10(dt$pvalue))/2) +
  scale_color_gradient(low = "blue", high = "red") +
  # scale_color_manual(values = c('#e41a1c','#377eb8'))+ 
  # scale_size_area(name="genecounts")+
  theme_bw()+
  scale_size_continuous(range=c(4,10))+
  labs(y='',x='-log10(pvalue)',title='')+
  theme() +
  geom_vline(xintercept = -log10(0.05),linetype =2,colour = 'black')+
  theme(plot.title=element_text(size=22,color="black",hjust = 0.5),
        axis.text=element_text(size=15, color="black"),
        axis.text.x=element_text(angle = 0,hjust = 0.5),
        axis.title=element_text(size=20),
        axis.line=element_line(size=1,colour='black'),
        axis.ticks=element_line(size=1,colour='black'),
        #legend.position = 'top',
        legend.text=element_text(size=16),legend.title = element_text(size=20))
p1
ggsave(plot = p1,'./picture/con_pac5_limma_kegg.pdf', width = 10, height = 8,dpi = 600)

ekegg_pac5$Group <- "c_pac5"
ekegg_pac5_top10 <- ekegg_pac5[order(ekegg_pac5$pvalue),][1 :10,]
ego_BP_pac5$Group <- "c_pac5"
ego_BP_pac5_top10 <- ego_BP_pac5[order(ego_BP_pac5$pvalue),][1:10,]
ekegg_pac2$Group <- "c_pac2"
ekegg_pac2_top10 <- ekegg_pac2[order(ekegg_pac2$pvalue),][1:10,]
ego_BP_pac2$Group <- "c_pac2"
ego_BP_pac2_top10 <- ego_BP_pac2[order(ego_BP_pac2$pvalue),][1:10,]

##  kegg
KEGG <- rbind(ekegg_pac5_top10,ekegg_pac2_top10)[,2]
KEGG_5 <- ekegg_pac5[ekegg_pac5$Description %in% KEGG,]
KEGG_2 <- ekegg_pac2[ekegg_pac2$Description %in% KEGG,]
KEGG_5_2 <- rbind(KEGG_5,KEGG_2)

table(KEGG_5_2$Group)
# c_pac2 c_pac5 
# 20     19 
KEGG_5_2$log.p <- -log10(KEGG_5_2$pvalue)

KEGG_5_2$Group <- factor(KEGG_5_2$Group,levels = c('c_pac5','c_pac2'))
KEGG_5_2$Description <- factor(KEGG_5_2$Description,levels = rev(unique(KEGG_5_2$Description)))
colnames(KEGG_5_2)
p1 <- ggplot(KEGG_5_2,aes(x=Group,y=Description,size=Count,color=log.p,ylab=''))+
  geom_point()+
  scale_color_gradient(low = "blue",high = "red")+
  # scale_size_area(name="genecounts")+
  theme_bw()+
  scale_size_continuous(range=c(4,7))+
  labs(y='',x='',title='kegg(c_pac5 vs c_pac2)',color = "-log.p")+
  theme() +
  geom_vline(xintercept = 0.05,linetype =2)+
  theme(plot.title=element_text(size=23,hjust = 0.5),
        axis.text.x =element_text(size=18),
        axis.text.y =element_text(size=16),
        legend.text=element_text(size=14),
        legend.title = element_text(size=18))#+
# facet_grid(~ Type,scales = 'free')
p1
ggsave(plot = p1,'./picture/kegg(c_pac5 vs c_pac2).pdf', width = 12, height = 10,dpi = 600)


## BP 
# BP <- rbind(ego_BP_pac5,ego_BP_pac2)
BP <- rbind(ego_BP_pac5_top10,ego_BP_pac2_top10)[,2]
BP_5 <- ego_BP_pac5[ego_BP_pac5$Description %in% BP,]
BP_2 <- ego_BP_pac2[ego_BP_pac2$Description %in% BP,]
BP_5_2 <- rbind(BP_5,BP_2)

table(BP_5_2$Group)
# c_pac2 c_pac5 
# 18     20
BP_5_2$log.p <- -log10(BP_5_2$pvalue)

BP_5_2$Group <- factor(BP_5_2$Group,levels = c('c_pac5','c_pac2'))
BP_5_2$Description <- factor(BP_5_2$Description,levels = rev(unique(BP_5_2$Description)))
colnames(BP_5_2)

p1 <- ggplot(BP_5_2,aes(x=Group,y=Description,size=Count,color=log.p,ylab=''))+
  geom_point()+
  scale_color_gradient(low = "blue",high = "red")+
  # scale_size_area(name="genecounts")+
  theme_bw()+
  scale_size_continuous(range=c(4,7))+
  labs(y='',x='',title='BP(c_pac5 vs c_pac2)',color = "-log.p")+
  theme() +
  geom_vline(xintercept = 0.05,linetype =2)+
  theme(plot.title=element_text(size=23,hjust = 0.5),
        axis.text.x =element_text(size=18),
        axis.text.y =element_text(size=16),
        legend.text=element_text(size=14),legend.title = element_text(size=18))#+
# facet_grid(~ Type,scales = 'free')
p1
ggsave(plot = p1,'./picture/BP(c_pac5 vs c_pac2).pdf', width = 14, height = 10,dpi = 600)


#######protein
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




data_filter<-pro_data[which(pro_data$positive >0 & pro_data$negtive ==0),]

data_filter <- data_filter[order(data_filter$positive,decreasing = T),]

top10<-data_filter$gene_name[1:10]

pos_list <- data_filter[,c(4,5,6)]
# pos_list <- pro_data[,c(4,5,6)]
pos_melt <- melt(pos_list,id.vars = "gene_name",variable.name = "sample",value.name = "value")
pos_melt$value<-log(pos_melt$value+1)
colnames(pos_melt)<-c("Protein","sample","value")
breaks<-pretty(range(0,18), 8)
maximum<- breaks[length(breaks)]
sort(unique(pos_melt$Protein))
pos_melt$Protein <- factor(pos_melt$Protein,levels = sort(unique(pos_melt$Protein)))

p1 <- ggplot(data = pos_melt, aes( x = Protein, y = value, color = sample)) +
  geom_point( alpha = 1, size = 2) +
  scale_colour_manual(values = c("blue","red"))+
  scale_x_discrete(breaks=NULL)+
  scale_y_continuous(limits = c(-1,maximum),breaks = breaks)+
  theme_classic() +
  xlab( "proteins" ) + ylab( "log10(Abundances)" ) +
  geom_hline(yintercept = 14.5,linetype =2)+
  theme(plot.title=element_text(size=22,color="black",hjust = 0.5),
        axis.text.y=element_text(size=15, color="black"),
        # axis.text.x=element_text(angle = 0,hjust = 0.5),
        #axis.ticks.x = element_blank(),
        axis.title=element_text(size=20),
        axis.line=element_line(size=1,colour='black'),
        axis.ticks=element_line(size=1,colour='black'),
        legend.position = 'top',
        legend.text=element_text(size=16),legend.title = element_text(size=20))+
   #geom_text_repel(data=subset(pos_melt, Protein %in% top10), aes(label=Protein),col="black",alpha = .8)
   geom_text(data=subset(pos_melt, Protein %in% top10 ), aes(label = Protein),col="black",size = 4,vjust = 0, nudge_y = 0.5)

p1
ggsave(plot = p1,'./picture/positive_proteinv2.pdf', width = 14, height = 8,dpi = 600)


pos_melt_top10<-pos_melt[which(pos_melt$Protein %in% top10),]
pos_melt_top10<-pos_melt_top10[order(pos_melt_top10$value,decreasing = T),]
pos_melt_top10$Protein<-factor(pos_melt_top10$Protein,levels = rev(pos_melt_top10$Protein[1:10]))

p1<-ggplot(pos_melt_top10, aes(x = sample, y=Protein)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = round(value, 3))) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 8, name="log10(Abundances)")+
  theme_classic()+
  theme(plot.title=element_text(size=22,color="black",hjust = 0.5),
        axis.text.y=element_text(size=15, color="black"),
        axis.text.x=element_text(size=15, color="black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title=element_text(size=20),
        axis.line=element_line(size=1,colour='black'),
        axis.ticks=element_line(size=1,colour='black'),
        # legend.position = 'top',
        legend.text=element_text(size=12),legend.title = element_text(size=12))
p1
ggsave(plot = p1,'./picture/positive_top_10v2.pdf', width = 14, height = 8,dpi = 600)


# cytokines ---------------------------------------------------------------
rm(list = ls())
setwd("/data/dsx/xuming_20.09.23")
library(VennDiagram)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(pheatmap)
# 细胞因子list导入
cytokines = read.csv('./file/cytokines.txt',sep=',',header=F,stringsAsFactors=FALSE)
colnames(cytokines) <- "cytokines"
Con_PAC5 = read.table('./file/C_PAC5_deg_limma.tab',sep='\t',header=TRUE,stringsAsFactors=FALSE)
diff_gene<-Con_PAC5$gene[Con_PAC5$P.Value < 0.05]
overlapgene<-intersect(diff_gene,cytokines$cytokines)
overlap<-Con_PAC5[Con_PAC5$gene %in% overlapgene,]

# 7.protein -----------------------------------------------------------------
rm(list = ls())
library(reshape2)
pro_pos <- read.delim(file = './file/pro_pos_top50.tab',header = F, quote = "",sep = "\t",stringsAsFactors=FALSE)
pro_pos <- pro_pos[order(pro_pos$V18,decreasing = T),]

pos_top10 <- pro_pos[1:10,][,16]

pos_list <- pro_pos[,c(16:18)]; pos_list$V17 <- 0
colnames(pos_list) <- c("Gene","Abundances..Normalized._negtive","Abundances..Normalized._positive")
pos_melt <- melt(pos_list,id.vars = "Gene",variable.name = "fact",value.name = "value")
pos_melt$log10 <- log(pos_melt$value)
pos_melt[1:50,"log10"] <- 0

library(ggrepel)
volcano = 
  ggplot(data = pos_melt, aes( x = Gene, y = log10, color = fact)) +
  geom_point( alpha = 3, size = 5) +  
  scale_y_continuous(limits = c(0,18),breaks = seq(0,18,9))+
  theme_classic() + theme(legend.position = "none") +
  xlab( "proteins" ) + ylab( "log10(Abundances)" ) +
  theme(axis.text.x  = element_blank(),
        axis.text.y  = element_text(size = 18),
        axis.title = element_text(size = 20)) +
  theme(axis.ticks.x = element_blank()) +
  scale_colour_manual( values = c("Abundances..Normalized._negtive" = 'blue',"Abundances..Normalized._positive" = 'red')) +
  # geom_text_repel(data=subset(pos_melt, Gene %in% pos_top10 ), aes(label = Gene),col="black",alpha = 1) +
  geom_text(data=subset(pos_melt, Gene %in% pos_top10 ), aes(label = Gene),
            col="black",size = 5,vjust = 0, nudge_y = 0.5)
print(volcano)
plotfile=paste('./picture/positive_top_10.pdf',sep='')
ggsave(plotfile, plot=volcano, dpi = 600,width = 15, height = 5)


# 8.positive_enrich -------------------------------------------------------
getwd()
rm(list=ls())
library( "clusterProfiler")
library("org.Hs.eg.db")
### load(.RData)
load(file = "/home/devdata/tongxin/enrichment_annotation/human_enrichment.RData")
### go
pro_pos = read.table(file = './file/pro_pos_top50.tab',header = F, sep = "\t", quote = "",stringsAsFactors=FALSE)
gene_list = pro_pos$V16
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

write.table(ekegg_dif, file = "./enrichment/positive_kegg_enrich_dif.tab", quote = FALSE,sep="\t",row.names = FALSE)
write.table(ego_BP_result, file = "./enrichment/positive_BP_enrich.tab", quote = FALSE,sep="\t",row.names = FALSE)

# ####
rm(list = ls())
kegg <- read.delim(file = "./enrichment/positive_kegg_enrich_dif.tab",header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
bp <- read.delim(file = "./enrichment/positive_BP_enrich.tab",header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)

p <- ggplot(kegg[1:20,],aes(x=-log(pvalue),y=reorder(Description,-pvalue),size=Count,color=-log(pvalue),ylab=''))+
  geom_point()+scale_colour_gradient(low="blue",high="red")+ 
  scale_size_area(name="count")+theme_bw()+
  labs(y='',x='-log(pvalue)',title='kegg_positive_top20')+
  theme(plot.title=element_text(size=23,hjust = 0.5),
        axis.text.x =element_text(size=18),
        axis.text.y =element_text(size=16),
        axis.title.x = element_text(size = 20),
        legend.text=element_text(size=14),legend.title = element_text(size=18))
plotfile = "./picture/kegg_positive_top20.pdf"
ggsave(plotfile, plot=p, dpi = 600,width = 12, height = 8)

p <- ggplot(bp[1:20,],aes(x=-log(pvalue),y=reorder(Description,-pvalue),size=Count,color=-log(pvalue),ylab=''))+
  geom_point()+scale_colour_gradient(low="blue",high="red")+ 
  scale_size_area(name="count")+theme_bw()+
  labs(y='',x='-log(pvalue)',title='bp_positive_top20')+
  theme(plot.title=element_text(size=23,hjust = 0.5),
        axis.text.x =element_text(size=18),
        axis.text.y =element_text(size=16),
        axis.title.x = element_text(size = 20),
        legend.text=element_text(size=14),legend.title = element_text(size=18))
plotfile = "./picture/bp_positive_top20.pdf"
ggsave(plotfile, plot=p, dpi = 600,width = 18, height = 8)




# 韦恩??? ---------------------------------------------------------------------
rm(list = ls())
library(VennDiagram)
library(RColorBrewer)
nrDEGC2 = read.table('./file/C_PAC2.tab',sep='\t',header=TRUE,stringsAsFactors=FALSE)
nrDEGC5 = read.table('./file/C_PAC5.tab',sep='\t',header=TRUE,stringsAsFactors=FALSE)
pro_pos <- read.delim(file = './file/pro_pos_top50.tab',header = F, quote = "",sep = "\t",stringsAsFactors=FALSE)
HBV <- read.delim(file = './file/HBV.xls',header = F, quote = "",sep = "\t",stringsAsFactors=FALSE)
nrDEGC2 <- nrDEGC2[nrDEGC2$change_V1 != "NOT",1] # 958
nrDEGC5 <- nrDEGC5[nrDEGC5$change_V1 != "NOT",1] # 997
pro_pos <- pro_pos$V16
HBV_gene <- HBV$V3
venn.diagram(list(nrDEGC2 = nrDEGC2,nrDEGC5 = nrDEGC5,positive = pro_pos,HBV_gene = HBV_gene),
             fill=c(brewer.pal(7,"Set1")[1:4]),
             alpha = c(0.5,0.5,0.5,0.5),cex = 2,
             cat.cex = 3,cat.fontface = 4,lty=2,fontfamily = 3,cat.pos=12,
             resolution = 300,filename = "./picture/HBV_positive_5_2.tiff",
             height = 3000, width = 3500)
venn.diagram(list(nrDEGC5 = nrDEGC5,positive = pro_pos),
             fill=c(brewer.pal(7,"Set1")[1:2]),
             alpha = c(0.5,0.5),cex = 2,
             cat.cex = 3,cat.fontface = 4,lty=2,fontfamily = 3,cat.pos=12,
             resolution = 300,filename = "./picture/positive_PAC5.tiff",
             height = 3000, width = 3000)
intersect(nrDEGC5,pro_pos)
# "ENO1"    "NDUFAF5" "SUMF1" 


# cytokines ---------------------------------------------------------------

rm(list = ls())
setwd("/data/tongxin/home/xuming_20.09.23")
library(VennDiagram)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(pheatmap)
# 细胞因子list导入
cytokines = read.csv('./file/cytokines.txt',sep=',',header=F,stringsAsFactors=FALSE)
colnames(cytokines) <- "cytokines"

# C_PAC5 and cytokines
nrDEGC5 = read.table('./file/C_PAC5.tab',sep='\t',header=TRUE,stringsAsFactors=FALSE)
sect_gene5 <- nrDEGC5[nrDEGC5$SYMBOL %in% intersect(cytokines$cytokines,nrDEGC5$SYMBOL),]

tmp<-setdiff(cytokines$cytokines,nrDEGC5$SYMBOL)
## volcano
cols <- c("UP" = "red","DOWN" = "blue","NOT" = "black")
volcano = 
  ggplot(data = sect_gene5, aes( x = log2FC, y = -log10(pvalue), color = change_V1)) +
  scale_x_continuous(limits = c(-2.5,2.5),breaks = seq(-3,3,1))+
  scale_y_continuous(limits = c(0,4.5),breaks = seq(0,5,1))+
  geom_point( alpha = 0.4, size = 1.75) +  
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+  #y轴分界线
  theme_classic() + theme(legend.position = "none") +
  xlab( "log2 fold change" ) + ylab( "-log10 p-value" ) +
  ggtitle( "PAC5_Cytokines.volcano") +
  theme( plot.title = element_text( size = 23, hjust = 0.5),
         axis.text = element_text(size = 18),
         axis.title = element_text(size = 20)) +
  scale_colour_manual( values = cols) +
  geom_text_repel(data=subset(sect_gene, change_V1 == "UP" | change_V1 == "DOWN" ), 
                  aes(label=SYMBOL),col="black",alpha = 1,size = 8)
print(volcano)
# plotfile=paste('./picture/volcano','_C','_PAC5.png',sep='')
plotfile=paste('./picture/volcano','_Cytokines','_PAC5.pdf',sep='')
ggsave(plotfile, plot=volcano, dpi = 600,width = 10, height = 8)


## heatmap
rownames(sect_gene5) <- sect_gene5[,1];sect_gene5 <- sect_gene5[-1]
sect_gene5 <- sect_gene5[1:6]
Z_score <- (sect_gene5 - apply(sect_gene5, 1, mean)) / apply(sect_gene5, 1, sd)
annotation_col <- data.frame(Group = c(rep('Control',3),rep('PAC5',3)))
rownames(annotation_col) <- colnames(sect_gene5)

ancols = c("#66c2a5","#fc8d62")
names(ancols) = c("Control","PAC5")
ann_colors <- list(Group=ancols)
library(pheatmap)
pheatmap(Z_score,border_color="#C5C5C5",color = colorRampPalette(c("blue", "white", "red"))(50),
         angle_col = 0,cluster_cols = F,cluster_rows = T,
         clustering_distance_cols = "correlation",
         clustering_distance_rows = "correlation",
         show_rownames=T,annotation_col=annotation_col,
         annotation_colors=ann_colors,fontsize=8,filename = "./picture/PAC5_Cytokines_heatmap.pdf",
         width = 10, height = 8)
dev.off()
## heatmap####mean
rownames(sect_gene5) <- sect_gene5[,1];sect_gene5 <- sect_gene5[-1]
sect_gene5 <- sect_gene5[1:6]
Z_score <- (sect_gene5 - apply(sect_gene5, 1, mean)) / apply(sect_gene5, 1, sd)
annotation_col <- data.frame(Group = c(rep('Control',3),rep('PAC5',3)))
rownames(annotation_col) <- colnames(sect_gene5)

ancols = c("#66c2a5","#fc8d62")
names(ancols) = c("Control","PAC5")
ann_colors <- list(Group=ancols)
library(pheatmap)
pheatmap(Z_score,border_color="#C5C5C5",color = colorRampPalette(c("blue", "white", "red"))(50),
         angle_col = 0,cluster_cols = F,cluster_rows = T,
         clustering_distance_cols = "correlation",
         clustering_distance_rows = "correlation",
         show_rownames=T,annotation_col=annotation_col,
         annotation_colors=ann_colors,fontsize=8,filename = "./picture/PAC5_Cytokines_heatmap.pdf",
         width = 10, height = 8)
dev.off()



# C_PAC2 and cytokines
nrDEGC2 = read.table('./file/C_PAC2.tab',sep='\t',header=TRUE,stringsAsFactors=FALSE)
sect_gene2 <- nrDEGC2[nrDEGC2$SYMBOL %in% intersect(cytokines$cytokines,nrDEGC2$SYMBOL),]

## volcano
cols <- c("UP" = "red","DOWN" = "blue","NOT" = "black")
volcano = 
  ggplot(data = sect_gene2, aes( x = log2FC, y = -log10(pvalue), color = change_V1)) +
  scale_x_continuous(limits = c(-2.5,2.5),breaks = seq(-3,3,1))+
  scale_y_continuous(limits = c(0,4.5),breaks = seq(0,5,1))+
  geom_point( alpha = 0.4, size = 1.75) +  
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8)+  #y轴分界线
  theme_classic() + theme(legend.position = "none") +
  xlab( "log2 fold change" ) + ylab( "-log10 p-value" ) +
  # ggtitle( this_tile ) + 
  theme( plot.title = element_text( size = 23, hjust = 0.5),
         axis.text = element_text(size = 18),
         axis.title = element_text(size = 20)) +
  scale_colour_manual( values = cols) +
  geom_text_repel(data=subset(sect_gene2, change_V1 == "UP" | change_V1 == "DOWN" ), 
                  aes(label=SYMBOL),col="black",alpha = 1,size = 8)
print(volcano)
# plotfile=paste('./picture/volcano','_C','_PAC5.png',sep='')
plotfile=paste('./picture/volcano','_Cytokines','_PAC2.pdf',sep='')
ggsave(plotfile, plot=volcano, dpi = 600,width = 10, height = 8)


## heatmap
rownames(sect_gene2) <- sect_gene2[,1];sect_gene2 <- sect_gene2[-1]
sect_gene2 <- sect_gene2[1:6]
Z_score <- (sect_gene2 - apply(sect_gene2, 1, mean)) / apply(sect_gene2, 1, sd)
annotation_col <- data.frame(Group = c(rep('Control',3),rep('PAC2',3)))
rownames(annotation_col) <- colnames(sect_gene2)

ancols = c("#66c2a5","#fc8d62")
names(ancols) = c("Control","PAC2")
ann_colors <- list(Group=ancols)
library(pheatmap)
pheatmap(Z_score,border_color="#C5C5C5",color = colorRampPalette(c("blue", "white", "red"))(50),
         angle_col = 0,cluster_cols = F,cluster_rows = T,
         clustering_distance_cols = "correlation",
         clustering_distance_rows = "correlation",
         show_rownames=T,annotation_col=annotation_col,
         annotation_colors=ann_colors,fontsize=8,filename = "./picture/PAC2_Cytokines_heatmap.pdf",
         width = 10, height = 8)

## heatmap Cytokines_PAC5_PAC2
dif_data = read.table(file = './file/edata_rm0.tab',header = TRUE, sep = "\t", quote = "",stringsAsFactors=FALSE)
# C5_C2.cytokines <- union(rownames(sect_gene5),rownames(sect_gene2))
C5_C2.cytokines <- dif_data[rownames(sect_gene5),-1]

Z_score <- (C5_C2.cytokines - apply(C5_C2.cytokines, 1, mean)) / apply(C5_C2.cytokines, 1, sd)
annotation_col <- data.frame(Group = c(rep('Control',3),rep('PAC2',3),rep('PAC5',3)))
rownames(annotation_col) <- colnames(C5_C2.cytokines)
# pheatmap(Z_score)
ancols = c("#66c2a5","#fc8d62","#8da0cb")
names(ancols) = c("Control","PAC2","PAC5")
ann_colors <- list(Group=ancols)
library(pheatmap)
pheatmap(Z_score,border_color="#C5C5C5",color = colorRampPalette(c("blue", "white", "red"))(50),
         angle_col = 0,cluster_cols = F,cluster_rows = T,
         clustering_distance_cols = "correlation",
         clustering_distance_rows = "correlation",
         show_rownames=T,annotation_col=annotation_col,
         annotation_colors=ann_colors,fontsize=8,filename = "./picture/C5_C2.cytokines.pdf",
         width = 10,height = 8)

# bubble_plot ------------------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(ggplot2)
library("ggrepel")
pro_pos_neg <- read.table("/data/dsx/xuming_20.09.23/file/pro_pos_neg_data_all.tab",sep = "\t",header = T,stringsAsFactors = F)

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
write.table(pro_pos_neg,file = "/data/dsx/xuming_20.09.23/file/pro_pos_neg.tab",col.names = TRUE,row.names = F,quote = FALSE,sep = ",")
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
ggsave(plot = p,'/data/dsx/xuming_20.09.23/picture/bubble_protein_top10.pdf', width = 10, height = 10,dpi = 600)


######################
#######################
########## barplot-----
rm(list = ls())
library(ggplot2)
setwd("/data/dsx/xuming_20.09.23")
pro_pos_neg <- read.table("/data/dsx/xuming_20.09.23/file/pro_pos_neg_data_all.tab",sep = "\t",header = T,stringsAsFactors = F)
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
ggsave(plot = p,'/data/dsx/xuming_20.09.23/picture/Protein_Number_barplot.pdf', width = 10, height = 10,dpi = 600)

######## enrich-----
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

write.table(ekegg_dif, file = "/data/dsx/xuming_20.09.23/file/PAC5_only_pos_kegg_enrich.tab", quote = FALSE,sep="\t",row.names = FALSE)
write.table(ego_BP_result, file = "/data/dsx/xuming_20.09.23/file/PAC5_only_pos_BP_enrich.tab", quote = FALSE,sep="\t",row.names = FALSE)

### only_neg ###
### only_pos ###
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

write.table(ekegg_dif, file = "/data/dsx/xuming_20.09.23/file/NC_only_pos_kegg_enrich.tab", quote = FALSE,sep="\t",row.names = FALSE)
write.table(ego_BP_result, file = "/data/dsx/xuming_20.09.23/file/NC_only_pos_BP_enrich.tab", quote = FALSE,sep="\t",row.names = FALSE)

#### plot_bubble
rm(list=ls())
library(ggplot2)
ekegg_only_pos <- read.delim(file = '/data/dsx/xuming_20.09.23/file/PAC5_only_pos_kegg_enrich.tab',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
ego_BP_only_pos <- read.delim(file = '/data/dsx/xuming_20.09.23/file/PAC5_only_pos_BP_enrich.tab',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
ekegg_only_neg <- read.delim(file = '/data/dsx/xuming_20.09.23/file/NC_only_pos_kegg_enrich.tab',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)
ego_BP_only_neg <- read.delim(file = '/data/dsx/xuming_20.09.23/file/NC_only_pos_BP_enrich.tab',header = TRUE, quote = "",sep = "\t",stringsAsFactors=FALSE)

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
# only_neg only_pos 
#    16       18
KEGG_pos_neg$log.p <- -log10(KEGG_pos_neg$pvalue)

KEGG_pos_neg$Group <- factor(KEGG_pos_neg$Group,levels = c('PAC5','NC'))
KEGG_pos_neg$Description <- factor(KEGG_pos_neg$Description,levels = rev(unique(KEGG_pos_neg$Description)))
colnames(KEGG_pos_neg)
breaks<-pretty(range(0,KEGG_pos_neg$Count), 5)
maximum<- breaks[length(breaks)]
p1 <- ggplot(KEGG_pos_neg,aes(x=log.p,y=Description,size=Count,color=Group,ylab=''))+
  geom_point()+
  # scale_color_gradient(low = "blue",high = "red")+
  # scale_size_area(name="genecounts")+
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
  # facet_grid(~ Group)
p1
ggsave(plot = p1,'/data/dsx/xuming_20.09.23/picture/Protein_KEGG_pos_neg.pdf', width = 12, height = 10,dpi = 600)

#####  bp ego_BP_only_pos
BP <- rbind(ego_BP_only_pos_top10,ego_BP_only_neg_top10)[,2]
BP_pos <- ego_BP_only_pos[ego_BP_only_pos$Description %in% BP,]
BP_neg <- ego_BP_only_neg[ego_BP_only_neg$Description %in% BP,]
BP_pos_neg <- rbind(BP_pos,BP_neg)
BP_pos_neg <-  rbind(ego_BP_only_pos_top10,ego_BP_only_neg_top10)

table(BP_pos_neg$Group)
# only_neg only_pos 
#   18       14
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
ggsave(plot = p1,'/data/dsx/xuming_20.09.23/picture/Protein_BP_pos_neg.pdf', width = 12, height = 10,dpi = 600)

