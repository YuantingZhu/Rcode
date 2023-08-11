if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")


library(tidyverse)
library(GEOquery)
library(limma)
library(tibble)
library(dplyr)

gset <- getGEO("GSE76701",getGPL = FALSE)
gset <- gset$GSE76701_series_matrix.txt.gz
sample_info <- pData(gset)%>%
  dplyr::select(geo_accession,title)

write.table(sample_info,file = "sample_info.txt",
            sep ="\t",
            quote = FALSE)
sample_info <- read_delim("sample_info2.txt", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)
sample_info <- column_to_rownames(sample_info,var = "1")

group_by(sample_info, group)%>%
  summarise(count= n())

gene_exp <- as.data.frame(exprs(gset))
write.table(gene_exp,file = "gene_exp.txt",
            sep ="\t",
            quote = FALSE)
boxplot(gene_exp,
        outline = FALSE,
        las=2)


gene_info<-gene_info<-GPL570_55999 <- read_delim("D:/ZYT/R/HF6/HF6/GPL570-55999.txt", 
                                                 delim = "\t", escape_double = FALSE, 
                                                 trim_ws = TRUE)%>%
  dplyr::select(ID,GB_ACC,ENTREZ_GENE_ID,Gene_Symbol = "Gene Symbol",Gene_title = "Gene Title")

save(gset,gene_exp,sample_info,gene_info,file = "D:/ZYT/R/HF6/HF6-3/HF6-3.rdata")
load("D:/ZYT/R/HF6/HF6-3/HF6-3.rdata")


##DEG
library(limma)
design <- model.matrix(~0 + sample_info$group)
colnames(design) <- levels(factor(sample_info$group))
rownames(design) <- rownames(sample_info)

contrasts <- makeContrasts(
  HN=HF-NF,
  levels = design)

fit <- lmFit(gene_exp,design)
fit <- contrasts.fit(fit,contrasts)
fit <- eBayes(fit)
de_result <- topTable(fit,
                      coef = "HN",
                      number = Inf)

gene_exp <- as.data.frame(gene_exp)

de_result <- rownames_to_column(de_result,var = "probe")%>%
  mutate(direction = if_else(abs(logFC)<0.6 | P.Value >0.05,"ns",
                             if_else(logFC >= 0.6,"up","down")))%>%
  left_join(gene_info,by = c("probe"= "ID"))%>%
  left_join(rownames_to_column(gene_exp,var = "probe"),by= "probe")%>%
  select(-t,-B)%>%
  arrange(desc(abs(logFC)))

write.table(de_result,file = "de_result.txt",
            sep ="\t",
            quote = FALSE)

deg<- 
  filter(de_result,!is.na(ENTREZ_GENE_ID))%>%
  mutate(ENTREZ_GENE_ID=str_split(ENTREZ_GENE_ID," /// ",simplify = TRUE)[,1])%>%
  group_by(ENTREZ_GENE_ID)%>%
  filter(abs(logFC) == max(abs(logFC)))%>%
  distinct(ENTREZ_GENE_ID,.keep_all = TRUE)

deg_diff <- filter(deg,direction != "ns")

write.table(deg_diff,file = "deg_diff.txt",
            sep ="\t",
            quote = FALSE)


##PCA
PCAgene_exp <- read_delim( "PCAgene_exp.txt", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)

PCAgene_exp <- column_to_rownames(PCAgene_exp,var = "1")
PCAgene_exp <- t(PCAgene_exp)

install.packages("FactoMineR")
install.packages("factoextra")
library("FactoMineR")
library("factoextra")
library(FactoMineR)

PCAgene_exp.pca <- PCA(PCAgene_exp, ncp = 2, scale.unit = TRUE, graph = FALSE)
plot(PCAgene_exp.pca) 


pca_sample <- data.frame(PCAgene_exp.pca$ind$coord[ ,1:2])
head(pca_sample)

pca_eig1 <- round(PCAgene_exp.pca$eig[1,2], 2)
pca_eig2 <- round(PCAgene_exp.pca$eig[2,2],2 )

group <- read.delim("group.txt",row.names = 1, sep = '\t', check.names = FALSE)
group <- group[rownames(pca_sample), ]
pca_sample <- cbind(pca_sample, group)
pca_sample$samples <- rownames(pca_sample)
head(pca_sample) 

library(ggplot2)

p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = group), size = 3) 
scale_color_manual(values = c('orange', 'green',"blue","red")) 
theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
      legend.key = element_rect(fill = 'transparent'))
labs(x =  paste('PCA1:', pca_eig1, '%'), y = paste('PCA2:', pca_eig2, '%'), color = '')  

p

library(ggrepel)

p <- p + 
  geom_text_repel(aes(label = group), size = 3, show.legend = FALSE, 
                  box.padding = unit(0.5, 'lines'))

p

p + stat_ellipse(aes(color = group), level = 0.95, show.legend = FALSE)

p + stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +
  scale_fill_manual(values =  c('orange', 'green',"blue","red"))

p


library(plyr)
cluster_border <- ddply(pca_sample, 'group', function(df) df[chull(df[[1]], df[[2]]), ])

p + geom_polygon(data = cluster_border, aes(color = group), fill = NA, show.legend = FALSE)

HF6PCA <- p + geom_polygon(data = cluster_border, aes(fill = group), alpha = 0.2, show.legend = FALSE) +
  scale_fill_manual(values =  c('orange', 'green',"blue","red"))

pdf(file ="HF6PCA.pdf",width=9,height=7)
print(HF6PCA)
dev.off()

tiff(filename = "HF6PCA.tiff",width = 2000,height = 2000,res = 300)
print(HF6PCA)
dev.off()

##volcano heatmap
library("ggplot2")
library("ggrepel")
hs_data <- read.delim("clipboard")

##HFvolcano
ggplot(data = de_result, aes(x = logFC, y = -log10(P.Value)))+geom_point() 
de_result$threshold = as.factor(ifelse(de_result$P.Value < 0.05 & abs(de_result$logFC) >= 0.6, ifelse(de_result$logFC> 0.6 ,'Up','Down'),'NoSignifi'))
HFvolcano <- ggplot(data = de_result, aes(x = logFC, y = -log10(P.Value), colour=threshold,label = probe)) +
  geom_point(alpha=0.4, size=3) +
  scale_color_manual(values=c("blue", "grey","red"))+
  xlim(c(-5, 5)) +
  geom_vline(xintercept=c(-0.6,0.6),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8)+
  labs(x="log2(fold change)",y="-log10 (p-value)",title="Differential genes") +
  theme(plot.title = element_text(hjust = 0.5), legend.position="right", legend.title = element_blank())

tiff(filename = "HFvolcano.tiff",width = 2000,height = 2000,res = 300)
print(HFvolcano)
dev.off()

##HFheatmap
HF_NF_sample_info <- filter(sample_info,group%in% c("HF","NF"))

temp<- slice(de_result,1:120)%>%
  select(Gene_Symbol,one_of(rownames(HF_NF_sample_info)))
temp <- as_tibble(temp)
temp <-  na.omit(temp)%>%
  distinct(Gene_Symbol,.keep_all =TRUE)
de_exp_top <-  column_to_rownames(temp,var = "Gene_Symbol")
write.table(de_exp_top,file = "de_exp_top.txt",
            sep ="\t",
            quote = FALSE)


cols <- list(group = c(HF= "red",
                       NF = "green"))


library(pheatmap)
HFheatmap <- pheatmap(de_exp_top, 
                      annotation=select(HF_NF_sample_info,group), 
                      annotation_colors = cols,
                      color = colorRampPalette(c("green","white","red"))(50),
                      cluster_cols =F,
                      show_colnames = F,
                      show_rownames = F,
                      scale="row",
                      fontsize = 8,
                      fontsize_row=7,
                      fontsize_col=10)

tiff(filename = "HFheatmap.tiff",width = 2000,height = 2000,res = 300)
print(HFheatmap)
dev.off()


##venn diagram
install.packages("VennDiagram")
install.packages("openxlsx")

library (VennDiagram)  
library(openxlsx)

set1<-read.xlsx('Venn.xlsx',sheet= "Sheet1",sep=',')
set2<-read.xlsx('Venn.xlsx',sheet= "Sheet2",sep=',')

set1=t(set1)
set2=t(set2)

venn.diagram(x=list(set1,set2),
             scaled = F, 
             alpha= 0.5, 
             lwd=1,lty=1,col=c('#FFFFCC','#CCFFFF'),
             label.col ='black' ,
             cex = 2, 
             fontface = "bold",
             fill=c('#FFFFCC','#CCFFFF'),
             category.names = c("DEGs", "Necroptosis genes") , 
             cat.dist = 0.02, 
             cat.pos = -180,
             cat.cex = 2, 
             cat.fontface = "bold", 
             cat.col='black' ,
             cat.default.pos = "outer",
             output=TRUE,
             filename='D:/ZYT/R/HF6/HF6-3/Venn.tiff',
             imagetype="tiff",
             resolution = 400,
             compression = "lzw")
grid.draw(data)

##violin
library(reshape2)
library(ggpubr)

geneviolin=read.table("geneviolin.txt",header=T,sep="\t",check.names=F,row.names=1)

geneviolin=as.data.frame(t(geneviolin))
x=colnames(geneviolin)[1]

dataviolin <- melt(geneviolin,id.vars = c("group"),variable.name ='Gene',
                   value.name = 'Expression')
dataviolin$Expression=as.numeric(dataviolin$Expression)

colnames(dataviolin)=c("Type","Gene","Expression")

pviolin=ggviolin(dataviolin, x="Gene", y="Expression", color = "Type", 
                 ylab="Gene expression",
                 legend.title=x,
                 add.params = list(fill="white"),
                 palette = c("#009393","#AE0000"),
                 width=1, add = "boxplot")
pviolin=pviolin+rotate_x_text(60)
tiff(filename = "pviolin.tiff",width = 2000,height = 2000,res = 300)
print(pviolin)
dev.off()


##NRDEGs volcano heatmap
library("ggplot2")
library("ggrepel")
venn_all_de_result <- read_delim("venn_all_result.txt", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)

ggplot(data =venn_all_de_result, aes(x = logFC, y = -log10(P.Value)))+geom_point() 
venn_all_de_result$threshold = as.factor(ifelse(venn_all_de_result$P.Value < 0.05 & abs(venn_all_de_result$logFC) >= 0.6, ifelse(venn_all_de_result$logFC> 0.6 ,'Up','Down'),'NoSignifi'))

p <- ggplot(data = venn_all_de_result, 
            aes(x = logFC, 
                y = -log10(P.Value),
                colour=threshold,
                label = venn_all_de_result$Gene_Symbol)) +
  geom_point(alpha=0.4, size=2) +
  scale_color_manual(values=c("blue", "grey","red"))+
  xlim(c(-3, 3)) +
  geom_vline(xintercept=c(-0.6,0.6),lty=2,col="black",lwd=0.4) +
  geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.4) +
  labs(x="log2(fold change)",
       y="-log10 (p-value)",
       title="Differential genes")  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

p
venn_all_de_result$label=ifelse(venn_all_de_result$P.Value < 0.05 & abs(venn_all_de_result$logFC) >= 0.6,venn_all_de_result$Gene_Symbol,"")
venn_de_resultVolcano <- p+geom_text_repel(data = venn_all_de_result, aes(x = venn_all_de_result$logFC, 
                                                                          y = -log10(venn_all_de_result$P.Value), 
                                                                          label = label),
                                           size = 4,box.padding =unit(0.1, "lines"),
                                           point.padding = FALSE, 
                                           segment.color = "black",
                                           max.overlaps=80,
                                           show.legend = FALSE)


tiff(filename = "venn_de_resultVolcano.tiff",width = 2000,height = 2000,res = 300)
print(venn_de_resultVolcano)
dev.off()



##NRDEGs heatmap
HF_NF_sample_info <- filter(sample_info,group%in%c("HF","NF"))

venn_diff <- read_delim("venn_diff.txt", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)
venn_diff <- column_to_rownames(venn_diff,var = "Gene_Symbol")

cols <- list(group = c(HF= "red",
                       NF = "green"))

library(pheatmap)
venn_diffheatmap <- pheatmap(venn_diff, 
                             annotation=select(HF_NF_sample_info,group), 
                             annotation_colors = cols,
                             color = colorRampPalette(c("green","white","red"))(50),
                             cluster_cols =F,
                             show_colnames = F,
                             show_rownames = F,
                             scale="row",
                             fontsize = 8,
                             fontsize_row=7,
                             fontsize_col=10)


tiff(filename = "venn_diffheatmap.tiff",width = 2000,height = 2000,res = 300)
print(venn_diffheatmap)
dev.off()





##GSEA
library(data.table)
library(tidyverse)
library(magrittr)
library(knitr)
library(clusterProfiler)
library(enrichplot)

de_resultgsea <- read.table("de_resultgsea.txt",header = T,sep = "\t")

de_resultgsea=de_resultgsea[,c(3,1)]
colnames(de_resultgsea)=c("SYMBOL","foldChange")

a=na.omit(de_resultgsea)
a=a[!duplicated(a$SYMBOL),]

library(clusterProfiler)
library(DOSE)
library(stringr)
library("org.Hs.eg.db")
library(ggplot2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(data.table)
library(tidyverse)
library(magrittr)
library(knitr)
library(clusterProfiler)
library(enrichplot)



df.id<-bitr(a$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
head(df.id)
dim(df.id)
easy.df<-merge(a,df.id,by="SYMBOL",all=F)

sortdf<-easy.df[order(easy.df$foldChange, decreasing = T),]
head(sortdf)

gene.expr = sortdf$foldChange
names(gene.expr) <-sortdf$ENTREZID
head(gene.expr)

require(enrichplot)
require(clusterProfiler)

R.utils::setOption( "clusterProfiler.download.method",'auto' )
kk <- gseKEGG(gene.expr, organism = "hsa",pvalueCutoff  = 1)

kk <- DOSE::setReadable(kk, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
sortkk<-kk[order(kk$enrichmentScore, decreasing = T),]
write.table(sortkk,"gseaKEGG_output.txt",sep = "\t",quote = F)


kkGO  <- gseGO(geneList     = gene.expr,
               ont          = "BP",  # "BP"???"MF"???"CC"???"ALL"
               OrgDb        = org.Hs.eg.db,#浜虹被org.Hs.eg.db 榧爋rg.Mm.eg.db
               keyType      = "ENTREZID",
               pvalueCutoff = 1) 
kkGO <- DOSE::setReadable(kkGO, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kkGO <-kkGO [order(kkGO$enrichmentScore, decreasing = T),]
write.table(kkGO,"gseaGO_output.txt",sep = "\t",quote = F)

####GSEA plot
apo <- gseaplot2(kk, "hsa04210",
                 color = "red", 
                 base_size = 20, 
                 rel_heights = c(2, 0.5, 1),
                 subplots = 1:3,
                 ES_geom = "line",
                 pvalue_table = T)

tiff(filename = "apo.tiff",width = 3000,height = 2000,res = 300)
print(apo)
dev.off()

necr <- gseaplot2(kk, "hsa04217",
                  color = "red", 
                  base_size = 20, 
                  rel_heights = c(2, 0.5, 1),
                  subplots = 1:3,
                  ES_geom = "line",
                  pvalue_table = T)

tiff(filename = "necr.tiff",width = 3000,height = 2000,res = 300)
print(necr)
dev.off()


##predefined gesa
library(data.table)
library(tidyverse)
library(magrittr)
library(knitr)
library(clusterProfiler)
library(enrichplot)

gene_List<-as.numeric(a$foldChange)
names(gene_List)<-as.character(a$SYMBOL)
gene_List<-sort(gene_List,decreasing = T)
print(head(gene_List))

geneset=read.gmt("c5.go.v7.5.1.symbols.gmt")

library(GSEABase) 
egmt <- GSEA(geneList=gene_List,pvalueCutoff = 0.05, TERM2GENE=geneset, nPermSimple = 10000,verbose=FALSE)
head(egmt)
gsea_results_df <- egmt@result

write.table(gsea_results_df,"GO_gsea.xls",sep = "\t",quote=
              F, row.names =F)

geneset2=read.gmt("c2.cp.kegg.v7.5.1.symbols.gmt")
library(GSEABase) 

egmt2 <- GSEA(geneList=gene_List,pvalueCutoff = 0.05, TERM2GENE=geneset2, nPermSimple = 10000,verbose=FALSE)
head(egmt2)

gsea_results_df2 <- egmt2@result
write.table(gsea_results_df2,"KEGG_gsea.xls",sep = "\t",quote=
              F, row.names =F)



WNT <- gseaplot2(egmt2, "KEGG_WNT_SIGNALING_PATHWAY",
                 color = "red", 
                 base_size = 20, 
                 rel_heights = c(2, 0.5, 1),
                 subplots = 1:3,
                 ES_geom = "line",
                 pvalue_table = T)

tiff(filename = "WNT.tiff",width = 3000,height = 2000,res = 300)
print(WNT)
dev.off()

NOTCH <- gseaplot2(egmt2, "KEGG_NOTCH_SIGNALING_PATHWAY",
                   color = "red", 
                   base_size = 20, 
                   rel_heights = c(2, 0.5, 1),
                   subplots = 1:3,
                   ES_geom = "line",
                   pvalue_table = T)

tiff(filename = "NOTCH.tiff",width = 3000,height = 2000,res = 300)
print(NOTCH)
dev.off()


###
install.packages("hrbrthemes")
library(ggplot2)
library(hrbrthemes)
KEGG_gsea <- read.table("KEGG_gsea.xls", 
                        header = T,sep = "\t")
KEGG_gsea <- KEGG_gsea[order(KEGG_gsea$NES, decreasing = T),]

dat<-data.frame(group=KEGG_gsea$NES, 
                val= KEGG_gsea$ID)

ggplot(dat, aes(x = group, y = reorder(val,group)))+ 
  geom_bar(stat = "identity",   
           show.legend = FALSE,   
           width =0.8) + 
  xlab("NES") + 
  ylab("Kegg pathways")

library(data.table)
setDT(dat)
dat[,col:=ifelse(group>0,1,2)]
dat$col<-as.factor(dat$col)

ggplot(dat, aes(x = group, y = reorder(val,group),fill=col))+ 
  geom_bar(stat = "identity",   
           show.legend = FALSE,   
           width =0.8) + 
  xlab("NES") + 
  ylab("Kegg pathways")


Keggbar <- ggplot(dat, aes(x = group, y = reorder(val,group),fill=col))+ 
  geom_bar(stat = "identity",   
           show.legend = FALSE,   
           width =0.8,
           aes(fill=group)) + 
  xlab("NES") + 
  ylab("Kegg pathways")+
  scale_fill_gradient2(low = "skyblue",         
                       mid = "aliceblue",          
                       high = "#56B1F7")


tiff(filename = "Keggbar.tiff",width = 4000,height = 2000,res = 300)
print(Keggbar)
dev.off()


##NRDEGs GO KEGG
install.packages("colorspace")
install.packages("stringi")
install.packages("ggplot2")
install.packages("digest")
install.packages("GOplot")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DOSE")
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")

library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library(GOplot)

pvalueFilter=0.05    
qvalueFilter=0.05  

colorSel="p.adjust"
if( qvalueFilter>0.05){
  colorSel="pvalue"
}


rt=read.table("gene.txt", header=T, sep="\t", check.names=F)
genes=as.vector(rt[,1])
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]   
gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

for(i in c("BP")){
  
  kk=enrichGO(gene=gene,OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont=i, readable=T)
  GO=as.data.frame(kk)
  GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
  
  write.table(GO, file=paste0(i,".GO.txt"), sep="\t", quote=F, row.names=F)
  
  showNum=10
  if(nrow(GO)<10){
    showNum=nrow(GO)
  }
  
  
  pdf(file=paste0(i,".barplot.pdf"), width=9, height=7)
  bar=barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, color=colorSel)
  print(bar)
  dev.off()
  
  
  pdf(file=paste0(i,".bubble.pdf"), width=12, height=7)
  bub=dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, color=colorSel)
  print(bub)
  dev.off()
  
  tiff(file=paste0(i,".bubble.tiff"),width = 3000,height = 2000,res = 300)
  print(bub)
  dev.off()
}



colorSel="pvalue"
if( qvalueFilter>0.05){
  colorSel="pvalue"
}


for(i in c("CC")){
  kk=enrichGO(gene=gene,OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont=i, readable=T)
  GO=as.data.frame(kk)
  GO=GO[(GO$pvalue<pvalueFilter),]
  write.table(GO, file=paste0(i,".GO.txt"), sep="\t", quote=F, row.names=F)
  showNum=10
  if(nrow(GO)<10){
    showNum=nrow(GO)
  }
  
  pdf(file=paste0(i,".barplot.pdf"), width=9, height=7)
  bar=barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, color=colorSel)
  print(bar)
  dev.off()
  
  pdf(file=paste0(i,".bubble.pdf"), width=9, height=7)
  bub=dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, color=colorSel)
  print(bub)
  dev.off()
  
  tiff(file=paste0(i,".bubble.tiff"),width = 3000,height = 2000,res = 300)
  print(bub)
  dev.off()
}


for(i in c("MF")){
  kk=enrichGO(gene=gene,OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont=i, readable=T)
  GO=as.data.frame(kk)
  GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
  write.table(GO, file=paste0(i,".GO.txt"), sep="\t", quote=F, row.names=F)
  showNum=10
  if(nrow(GO)<10){
    showNum=nrow(GO)
  }
  pdf(file=paste0(i,".barplot.pdf"), width=9, height=7)
  bar=barplot(kk, drop=TRUE, showCategory=showNum, label_format=130, color=colorSel)
  print(bar)
  dev.off()
  
  pdf(file=paste0(i,".bubble.pdf"), width=9, height=7)
  bub=dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=130, color=colorSel)
  print(bub)
  dev.off()
  
  tiff(file=paste0(i,".bubble.tiff"),width = 3000,height = 2000,res = 300)
  print(bub)
  dev.off()
}


##NRDEGs KEGG
library(ggplot2)
library(tidyverse)
Necr_diff_de_result <- read_delim("venn_diff_de_result.txt", 
                                  delim = "\t", escape_double = FALSE, 
                                  trim_ws = TRUE)
deg_Necr_diff_deresult <- dplyr::select(Necr_diff_de_result,
                                        probe,
                                        ENTREZ_GENE_ID,
                                        logFC,
                                        P.Value,
                                        adj.P.Val,
                                        direction
)%>%
  filter(!is.na(ENTREZ_GENE_ID))%>%
  mutate(ENTREZ_GENE_ID=str_split(ENTREZ_GENE_ID," /// ",simplify = TRUE)[,1])%>%
  group_by(ENTREZ_GENE_ID)%>%
  filter(abs(logFC) == max(abs(logFC)))%>%
  distinct(ENTREZ_GENE_ID,.keep_all = TRUE)

diffgeneNecr_diff <- filter(deg_Necr_diff_deresult,direction != "ns")%>%
  pull(ENTREZ_GENE_ID)

geneList <- deg_Necr_diff_deresult$logFC
names(geneList) <- deg_Necr_diff_deresult$ENTREZ_GENE_ID
geneList <- sort(geneList,decreasing = TRUE)


library(clusterProfiler)
library(org.Hs.eg.db)

install.packages("R.utiils")
R.utils::setOption("clusterProfiler.download.method","auto")

de_ekp_Necr_diff <- enrichKEGG(gene = diffgeneNecr_diff,
                               organism = 'hsa',
                               keyType = "kegg",
                               pvalueCutoff = 0.05)
de_ekp_Necr_diff  <- setReadable(de_ekp_Necr_diff,OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
de_ekp_Necr_diff_df <- de_ekp_Necr_diff@result

write.table(de_ekp_Necr_diff_df,file = "de_ekp_Necr_diff_df.txt",
            sep ="\t",
            quote = FALSE)


install.packages("ggnewscale")
library(ggnewscale)
barplot(de_ekp_Necr_diff,showCategory = 10)
kegg <- dotplot(de_ekp_Necr_diff,showCategory = 10)


tiff(filename = "kegg.tiff",width = 2000,height = 2000,res = 300)
print(kegg)
dev.off()


cnetplot <- cnetplot(de_ekp_Necr_diff,foldChange = geneList,
                     showCategory = 20,
                     cex_label_gene =1, 
                     cex_label_category =1,
                     circular = TRUE,
                     colorEdge = TRUE,
                     node_label = "all")

tiff(filename = "cnetplot.tiff",width = 3500,height = 3000,res = 300)
print(cnetplot)
dev.off()


##corheatmap
library(tidyverse)
library(corrplot)
install.packages("ggcorrplot")
library(ggcorrplot)

venn_diff <- read.table("venn_diff.txt",row.names = 1, header=T, sep="\t", check.names=F)
cormtcars <- round(cor(t(venn_diff),method='spearman'), 3)
p.mtcars <- cor_pmat(cormtcars)


cormtcars <- as.data.frame(cormtcars)
ggcorrplot(cormtcars,method = "circle")

corheatmap <- ggcorrplot(cormtcars,method = "circle",
                         hc.order = TRUE,hc.method = "ward.D",
                         outline.color = "white",ggtheme = theme_bw(),
                         type = "upper",
                         colors=c("#6D9EC1","white","#E46726"),
                         lab = TRUE,lab_size = 2,
                         p.mat = p.mtcars,pch = 4)



tiff(filename = "corheatmap.tiff",width = 2000,height = 2000,res = 300)
print(corheatmap)
dev.off()




##GSE21610

#引用包
library(ggplot2)
library(reshape2)
library(ggpubr)
inputFile="de_result26210DDX58.txt"       #输入
outFile="barplot.pdf"       #输出

#读取文件
rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)

x=colnames(rt)[1]

data <- melt(rt,id.vars = c("group"),variable.name ='Gene',
             value.name = 'Expression')
data$Expression=as.numeric(data$Expression)

colnames(data)=c("group","Gene","Expression")

DDX58=ggviolin(data, x="Gene", y="Expression", color = "group", 
               ylab="Gene expression",
               add.params = list(fill="white"),
               palette = c("#009393","#AE0000"),
               width=0.5, add = "mean_ci")+xlab(NULL)


tiff(filename = "DDX58.tiff",width = 800,height = 1000,res = 300)
print(DDX58)
dev.off()



#CASP1
inputFile="de_result26210CASP1.txt"   

rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)

x=colnames(rt)[1]

data <- melt(rt,id.vars = c("group"),variable.name ='Gene',
             value.name = 'Expression')
data$Expression=as.numeric(data$Expression)

colnames(data)=c("group","Gene","Expression")

CASP1=ggviolin(data, x="Gene", y="Expression", color = "group", 
               ylab="Gene expression",
               add.params = list(fill="white"),
               palette = c("#009393","#AE0000"),
               width=0.5, add = "mean_ci")+xlab(NULL)


tiff(filename = "CASP1.tiff",width = 800,height = 1000,res = 300)
print(CASP1)
dev.off()


#JAK2
inputFile="de_result26210JAK2.txt"   

rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)

x=colnames(rt)[1]

data <- melt(rt,id.vars = c("group"),variable.name ='Gene',
             value.name = 'Expression')
data$Expression=as.numeric(data$Expression)

colnames(data)=c("group","Gene","Expression")

JAK2=ggviolin(data, x="Gene", y="Expression", color = "group", 
              ylab="Gene expression",
              add.params = list(fill="white"),
              palette = c("#009393","#AE0000"),
              width=0.5, add = "mean_ci")+xlab(NULL)


tiff(filename = "JAK2.tiff",width = 800,height = 1000,res = 300)
print(JAK2)
dev.off()

#STAT4
inputFile="de_result26210STAT4.txt"   
rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)

x=colnames(rt)[1]

data <- melt(rt,id.vars = c("group"),variable.name ='Gene',
             value.name = 'Expression')
data$Expression=as.numeric(data$Expression)

colnames(data)=c("group","Gene","Expression")

STAT4=ggviolin(data, x="Gene", y="Expression", color = "group", 
               ylab="Gene expression",
               add.params = list(fill="white"),
               palette = c("#009393","#AE0000"),
               width=0.5, add = "mean_ci")+xlab(NULL)


tiff(filename = "STAT4.tiff",width = 800,height = 1000,res = 300)
print(STAT4)
dev.off()



#HSP90AB1
inputFile="de_result26210HSP90AB1.txt"   
rt=read.table(inputFile,header=T,sep="\t",check.names=F,row.names=1)

x=colnames(rt)[1]

data <- melt(rt,id.vars = c("group"),variable.name ='Gene',
             value.name = 'Expression')
data$Expression=as.numeric(data$Expression)

colnames(data)=c("group","Gene","Expression")

HSP90AB1=ggviolin(data, x="Gene", y="Expression", color = "group", 
                  ylab="Gene expression",
                  add.params = list(fill="white"),
                  palette = c("#009393","#AE0000"),
                  width=0.5, add = "mean_ci")+xlab(NULL)


tiff(filename = "HSP90AB1.tiff",width = 800,height = 1000,res = 300)
print(HSP90AB1)
dev.off()









