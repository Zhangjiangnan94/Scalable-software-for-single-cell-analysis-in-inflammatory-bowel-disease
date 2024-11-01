setwd("./UC-master")
source('analysis.r')
cell_subsets = read.table('cell_subsets.txt', sep='\t', header=F, stringsAsFactors=F)
epi.seur = readRDS('train.Epi.seur.rds')
fib.seur = readRDS('train.Fib.seur.rds')
imm.seur = readRDS('train.Imm.seur.rds')

#
epi.counts = epi.seur@assays[['RNA']]@counts
fib.counts = fib.seur@assays[['RNA']]@counts
imm.counts = imm.seur@assays[['RNA']]@counts
data = sparse_cbind(list(epi.counts, fib.counts, imm.counts))
numi = colSums(data)
data = scaleMargins(data, cols=1e4/numi)
idents = c(as.character(epi.seur@active.ident), as.character(fib.seur@active.ident), as.character(imm.seur@active.ident))
groups = data.frame(cell_subsets, row.names=1)[idents,1]
samples = c(epi.seur@meta.data$Sample, fib.seur@meta.data$Sample, imm.seur@meta.data$Sample)
epi.freq = as.matrix(as.data.frame.matrix(table(epi.seur@meta.data$Sample, epi.seur@meta.data$Cluster)))
fib.freq = as.matrix(as.data.frame.matrix(table(fib.seur@meta.data$Sample, fib.seur@meta.data$Cluster)))
imm.freq = as.matrix(as.data.frame.matrix(table(imm.seur@meta.data$Sample, imm.seur@meta.data$Cluster)))
all.freq = sparse_cbind(list(epi.freq, fib.freq, imm.freq))
reps = gsub('[12]*[ab]*$', '', rownames(all.freq))
temp = as.matrix(data.frame(aggregate(as.matrix(all.freq), list(reps), sum), row.names=1))
colnames(temp) = colnames(all.freq)
all.freq = temp[,colSums(temp) > 0]
ep.ident = levels(epi.seur@active.ident)
lp.ident = c(levels(fib.seur@active.ident), levels(imm.seur@active.ident))
ep.freq = all.freq[grep('Epi', rownames(all.freq)), ep.ident]
lp.freq = all.freq[grep('LP', rownames(all.freq)), lp.ident]
meta = read.table('G:/ulcerative_colitis/all.meta2.txt', sep='\t', header=T, stringsAsFactors=F)
sample2health = data.frame(unique(data.frame(sample=gsub('[12]*[ab]*$', '', meta[,'Sample']), health=meta[,'Health'])), row.names=1)
ep.cov = data.frame(condition=factor(sample2health[rownames(ep.freq),1], levels=c('Healthy', 'Inflamed', 'Non-inflamed')), row.names=rownames(ep.freq))
lp.cov = data.frame(condition=factor(sample2health[rownames(lp.freq),1], levels=c('Healthy', 'Inflamed', 'Non-inflamed')), row.names=rownames(lp.freq))
ep.pvals = dirichlet_regression(counts=ep.freq, covariates=ep.cov, formula=counts ~ condition)$pvals
colnames(ep.pvals) = colnames(ep.freq)
lp.pvals = dirichlet_regression(counts=lp.freq, covariates=lp.cov, formula=counts ~ condition)$pvals
colnames(lp.pvals) = colnames(lp.freq)
ep.pct = 100*ep.freq/rowSums(ep.freq)
#matrix_barplot(ep.pct, group_by=ep.cov$condition, pvals=ep.pvals, colors=set.colors)
matrix_barplot(ep.pct, group_by=ep.cov$condition, pvals=ep.pvals, colors=c("#7BABC9","#FF8F70","#63C1A2"))

lp.pct = 100*lp.freq/rowSums(lp.freq)
#matrix_barplot(lp.pct, group_by=lp.cov$condition, pvals=lp.pvals, colors=set.colors)
matrix_barplot(lp.pct, group_by=lp.cov$condition, pvals=lp.pvals, colors=c("#7BABC9","#FF8F70","#63C1A2"))

fib.ident = levels(fib.seur@active.ident)
imm.ident = levels(imm.seur@active.ident)
fib.freq = all.freq[grep('LP', rownames(all.freq)), fib.ident]
imm.freq = all.freq[grep('LP', rownames(all.freq)), imm.ident]

fib.cov = data.frame(condition=factor(sample2health[rownames(fib.freq),1], levels=c('Healthy', 'Inflamed', 'Non-inflamed')), row.names=rownames(fib.freq))
imm.cov = data.frame(condition=factor(sample2health[rownames(imm.freq),1], levels=c('Healthy', 'Inflamed', 'Non-inflamed')), row.names=rownames(imm.freq))

fib.pvals = dirichlet_regression(counts=fib.freq, covariates=fib.cov, formula=counts ~ condition)$pvals
colnames(fib.pvals) = colnames(fib.freq)
fib.pct = 100*fib.freq/rowSums(fib.freq)

imm.pvals = dirichlet_regression(counts=imm.freq, covariates=imm.cov, formula=counts ~ condition)$pvals
colnames(imm.pvals) = colnames(imm.freq)
imm.pct = 100*imm.freq/rowSums(imm.freq)

matrix_barplot(fib.pct, group_by=fib.cov$condition, pvals=fib.pvals, colors=c("#7BABC9","#FF8F70","#63C1A2"))
matrix_barplot(imm.pct, group_by=imm.cov$condition, pvals=imm.pvals, colors=c("#7BABC9","#FF8F70","#63C1A2"))

#
#ep.ident = levels(epi.seur@active.ident)
#fib.ident = levels(fib.seur@active.ident)
#imm.ident = levels(imm.seur@active.ident)
#ep.freq = all.freq[grep('Epi', rownames(all.freq)), ep.ident]
#fib.freq = all.freq[grep('LP', rownames(all.freq)), fib.ident]
#imm.freq = all.freq[grep('LP', rownames(all.freq)), imm.ident]


epi.expr <- epi.seur@assays$RNA@data
gene_name <- c("CLPP")
epi.gene_expression <- epi.expr %>% .[gene_name,] %>% as.data.frame()
colnames(epi.gene_expression) <- paste0(gene_name)
epi.seur$CLPP <- epi.gene_expression[,paste0(gene_name)]
epi.meta<-epi.seur@meta.data

fib.expr <- fib.seur@assays$RNA@data
fib.gene_expression <- fib.expr %>% .[gene_name,] %>% as.data.frame()
colnames(fib.gene_expression) <- paste0(gene_name)
fib.seur$CLPP <- fib.gene_expression[,paste0(gene_name)]
fib.meta<-fib.seur@meta.data

imm.expr <- imm.seur@assays$RNA@data
imm.gene_expression <- imm.expr %>% .[gene_name,] %>% as.data.frame()
colnames(imm.gene_expression) <- paste0(gene_name)
imm.seur$CLPP <- imm.gene_expression[,paste0(gene_name)]
imm.meta<-imm.seur@meta.data

epi.CLPP = as.matrix(as.data.frame.matrix(aggregate(epi.meta$CLPP, by=list(epi.meta$Cluster, epi.meta$Health),mean)))
fib.CLPP = as.matrix(as.data.frame.matrix(aggregate(fib.meta$CLPP, by=list(fib.meta$Cluster, fib.meta$Health),mean)))
imm.CLPP = as.matrix(as.data.frame.matrix(aggregate(imm.meta$CLPP, by=list(imm.meta$Cluster, imm.meta$Health),mean)))
epi.CLPP2 = as.matrix(as.data.frame.matrix(aggregate(epi.meta$CLPP, by=list(epi.meta$Sample, epi.meta$Cluster),mean)))
fib.CLPP2 = as.matrix(as.data.frame.matrix(aggregate(fib.meta$CLPP, by=list(fib.meta$Sample, fib.meta$Cluster),mean)))
imm.CLPP2 = as.matrix(as.data.frame.matrix(aggregate(imm.meta$CLPP, by=list(imm.meta$Sample, imm.meta$Cluster),mean)))

write.table(epi.CLPP2, file='epi.CLPP.csv', row.names = TRUE,quote = FALSE,sep = ',')
write.table(fib.CLPP2, file='fib.CLPP.csv', row.names = TRUE,quote = FALSE,sep = ',')
write.table(imm.CLPP2, file='imm.CLPP.csv', row.names = TRUE,quote = FALSE,sep = ',')

epi.CLPP3 = read.table('clipboard', sep='\t', header=T, stringsAsFactors=F)
rownames(epi.CLPP3) = epi.CLPP3$X
epi.CLPP3 = epi.CLPP3[,-1]
epiCLPP.cov = data.frame(condition=factor(sample2health[rownames(epi.CLPP3),1], levels=c('Healthy', 'Inflamed', 'Non-inflamed')), row.names=rownames(epi.CLPP3))
epiCLPP.pvals = dirichlet_regression(counts=epi.CLPP3, covariates=epiCLPP.cov, formula=counts ~ condition)$pvals
colnames(epiCLPP.pvals) = colnames(epi.CLPP3)
matrix_barplot(epi.CLPP3, group_by=epiCLPP.cov$condition, pvals=epiCLPP.pvals, colors=c("#7BABC9","#FF8F70","#63C1A2"))

fib.CLPP3 = read.table('clipboard', sep='\t', header=T, stringsAsFactors=F)
rownames(fib.CLPP3) = fib.CLPP3$X
fib.CLPP3 = fib.CLPP3[,-1]
fibCLPP.cov = data.frame(condition=factor(sample2health[rownames(fib.CLPP3),1], levels=c('Healthy', 'Inflamed', 'Non-inflamed')), row.names=rownames(fib.CLPP3))
fibCLPP.pvals = dirichlet_regression(counts=fib.CLPP3, covariates=fibCLPP.cov, formula=counts ~ condition)$pvals
colnames(fibCLPP.pvals) = colnames(fib.CLPP3)
matrix_barplot(fib.CLPP3, group_by=fibCLPP.cov$condition, pvals=fibCLPP.pvals, colors=c("#7BABC9","#FF8F70","#63C1A2"))

imm.CLPP3 = read.table('clipboard', sep='\t', header=T, stringsAsFactors=F)
rownames(imm.CLPP3) = imm.CLPP3$X
imm.CLPP3 = imm.CLPP3[,-1]
immCLPP.cov = data.frame(condition=factor(sample2health[rownames(imm.CLPP3),1], levels=c('Healthy', 'Inflamed', 'Non-inflamed')), row.names=rownames(imm.CLPP3))
immCLPP.pvals = dirichlet_regression(counts=imm.CLPP3, covariates=immCLPP.cov, formula=counts ~ condition)$pvals
colnames(immCLPP.pvals) = colnames(imm.CLPP3)
matrix_barplot(imm.CLPP3, group_by=immCLPP.cov$condition, pvals=immCLPP.pvals, colors=c("#7BABC9","#FF8F70","#63C1A2"))

#visualization
library(ggplot2)
library(resharp2)
epi.CLPP4 = as.matrix(as.data.frame.matrix(aggregate(epi.meta$CLPP, by=list(epi.meta$Sample, epi.meta$Cluster, epi.meta$Health),mean)))
epi.CLPP4=as.data.frame(epi.CLPP4)
str(epi.CLPP4)
epi.CLPP4$x <- as.numeric(epi.CLPP4$x)
ggplot(epi.CLPP4, aes(x = Group.3, y = x , fill = Group.3)) +
  geom_violin(alpha = 0.5, trim = FALSE, scale = "width", aes(linetype=NA)) +
  geom_jitter(shape=21,aes(fill=Group.3),position = position_jitter(width = 0.2))+
  xlab("Subset")+ylab("CLPP")+
  theme_bw()+theme(legend.position = "none")+facet_wrap(~Group.2)+
  theme(axis.text.x = element_text(angle = 315, hjust = 0, vjust = 0))+
  scale_fill_manual(values = c("#7BABC9","#FF8F70","#63C1A2"))

fib.CLPP4 = as.matrix(as.data.frame.matrix(aggregate(fib.meta$CLPP, by=list(fib.meta$Sample, fib.meta$Cluster, fib.meta$Health),mean)))
fib.CLPP4=as.data.frame(fib.CLPP4)
str(fib.CLPP4)
fib.CLPP4$x <- as.numeric(fib.CLPP4$x)
ggplot(fib.CLPP4, aes(x = Group.3, y = x , fill = Group.3)) +
  geom_violin(alpha = 0.5, trim = FALSE, scale = "width", aes(linetype=NA)) +
  geom_jitter(shape=21,aes(fill=Group.3),position = position_jitter(width = 0.2))+
  xlab("Subset")+ylab("CLPP")+
  theme_bw()+theme(legend.position = "none")+facet_wrap(~Group.2)+
  theme(axis.text.x = element_text(angle = 315, hjust = 0, vjust = 0))+
  scale_fill_manual(values = c("#7BABC9","#FF8F70","#63C1A2"))

imm.CLPP4 = as.matrix(as.data.frame.matrix(aggregate(imm.meta$CLPP, by=list(imm.meta$Sample, imm.meta$Cluster, imm.meta$Health),mean)))
imm.CLPP4=as.data.frame(imm.CLPP4)
str(imm.CLPP4)
imm.CLPP4$x <- as.numeric(imm.CLPP4$x)
ggplot(imm.CLPP4, aes(x = Group.3, y = x , fill = Group.3)) +
  geom_violin(alpha = 0.5, trim = FALSE, scale = "width", aes(linetype=NA)) +
  geom_jitter(shape=21,aes(fill=Group.3),position = position_jitter(width = 0.2))+
  xlab("Subset")+ylab("CLPP")+
  theme_bw()+theme(legend.position = "none")+facet_wrap(~Group.2)+
  theme(axis.text.x = element_text(angle = 315, hjust = 0, vjust = 0))+
  scale_fill_manual(values = c("#7BABC9","#FF8F70","#63C1A2"))

#
epi.CLPP5<-transform(epi.CLPP4,dist_cat_n=as.numeric(as.factor(Group.2)),
                     scat_adj=ifelse(Group.3 == "Uninflamed",0.2,ifelse(Group.3 == "Inflamed",0,-0.2)))
ggplot(epi.CLPP5, aes(x = Group.2, y = x , fill = Group.3)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar")+
  geom_jitter(shape=21,aes(scat_adj+dist_cat_n, x, fill=Group.3),position = position_jitter(width = 0.1))+
  xlab("Subset")+ylab("CLPP")+
  theme_bw()+theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 315, hjust = 0, vjust = 0))+
  scale_fill_manual(values = c("#7BABC9","#FF8F70","#63C1A2"))

fib.CLPP5 = fib.CLPP4[grep('LP', fib.CLPP4$Group.1), ]
fib.CLPP5<-transform(fib.CLPP5,dist_cat_n=as.numeric(as.factor(Group.2)),
                     scat_adj=ifelse(Group.3 == "Uninflamed",0.2,ifelse(Group.3 == "Inflamed",0,-0.2)))
ggplot(fib.CLPP5, aes(x = Group.2, y = x , fill = Group.3)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar")+
  geom_jitter(shape=21,aes(scat_adj+dist_cat_n, x, fill=Group.3),position = position_jitter(width = 0.1))+
  xlab("Subset")+ylab("CLPP")+
  theme_bw()+theme(legend.position="right")+
  theme(axis.text.x = element_text(angle = 315, hjust = 0, vjust = 0))+
  scale_fill_manual(values = c("#7BABC9","#FF8F70","#63C1A2"))

imm.CLPP5 = imm.CLPP4[grep('LP', imm.CLPP4$Group.1), ]
imm.CLPP5<-transform(imm.CLPP5,dist_cat_n=as.numeric(as.factor(Group.2)),
                     scat_adj=ifelse(Group.3 == "Uninflamed",0.2,ifelse(Group.3 == "Inflamed",0,-0.2)))
ggplot(imm.CLPP5, aes(x = Group.2, y = x , fill = Group.3)) +
   geom_boxplot() +
   stat_boxplot(geom = "errorbar")+
   geom_jitter(shape=21,aes(scat_adj+dist_cat_n, x, fill=Group.3),position = position_jitter(width = 0.1))+
   xlab("Subset")+ylab("CLPP")+
   theme_bw()+theme(legend.position="right")+
   theme(axis.text.x = element_text(angle = 315, hjust = 0, vjust = 0))+
   scale_fill_manual(values = c("#7BABC9","#FF8F70","#63C1A2"))

#
epi.CLPP5<-transform(epi.CLPP4,dist_cat_n=as.numeric(as.factor(Group.3)))
ggplot(epi.CLPP5, aes(x = Group.3, y = x , fill = Group.3)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar")+
  geom_jitter(shape=21,aes(dist_cat_n, x, fill=Group.3),position = position_jitter(width = 0.2))+
  xlab("Subset")+ylab("CLPP")+
  theme_bw()+theme(legend.position="none")+facet_wrap(~Group.2)+
  theme(axis.text.x = element_text(angle = 315, hjust = 0, vjust = 0))+
  scale_fill_manual(values = c("#7BABC9","#FF8F70","#63C1A2"))

fib.CLPP5<-transform(fib.CLPP4,dist_cat_n=as.numeric(as.factor(Group.3)))
ggplot(fib.CLPP5, aes(x = Group.3, y = x , fill = Group.3)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar")+
  geom_jitter(shape=21,aes(dist_cat_n, x, fill=Group.3),position = position_jitter(width = 0.2))+
  xlab("Subset")+ylab("CLPP")+
  theme_bw()+theme(legend.position="none")+facet_wrap(~Group.2)+
  theme(axis.text.x = element_text(angle = 315, hjust = 0, vjust = 0))+
  scale_fill_manual(values = c("#7BABC9","#FF8F70","#63C1A2"))

imm.CLPP5<-transform(imm.CLPP4,dist_cat_n=as.numeric(as.factor(Group.3)))
ggplot(imm.CLPP5, aes(x = Group.3, y = x , fill = Group.3)) +
  geom_boxplot() +
  stat_boxplot(geom = "errorbar")+
  geom_jitter(shape=21,aes(dist_cat_n, x, fill=Group.3),position = position_jitter(width = 0.2))+
  xlab("Subset")+ylab("CLPP")+
  theme_bw()+theme(legend.position="none")+facet_wrap(~Group.2)+
  theme(axis.text.x = element_text(angle = 315, hjust = 0, vjust = 0))+
  scale_fill_manual(values = c("#7BABC9","#FF8F70","#63C1A2"))

#data_melt <- melt(t(ep.pct))
#ggplot(data_melt, aes(x = Var1 , y = value , fill = Var1)) +
#  geom_violin(alpha = 0.5, scale = "width", aes(linetype=NA)) + 
#  geom_jitter(shape=21,aes(fill=Var1),position = position_jitter(width = 0.2))+
#  xlab("Subset")+ylab("Ratio")+
#  theme_bw()+theme(legend.position = "none")+facet_wrap(~Var1)



#rubbish
ep.cov = data.frame(condition=factor(sample2health[rownames(ep.freq),1], levels=c('Healthy', 'Inflamed', 'Non-inflamed')), row.names=rownames(ep.freq))
lp.cov = data.frame(condition=factor(sample2health[rownames(lp.freq),1], levels=c('Healthy', 'Inflamed', 'Non-inflamed')), row.names=rownames(lp.freq))
ep.pvals = dirichlet_regression(counts=ep.freq, covariates=ep.cov, formula=counts ~ condition)$pvals
colnames(ep.pvals) = colnames(ep.freq)
lp.pvals = dirichlet_regression(counts=lp.freq, covariates=lp.cov, formula=counts ~ condition)$pvals
colnames(lp.pvals) = colnames(lp.freq)
ep.pct = 100*ep.freq/rowSums(ep.freq)

ggplot(epi.CLPP2)+geom_violin()
ggplot(epi.CLPP)+geom_violin()
save.image("E:/Intra and inter/ulcerative_colitis-master/ulcerative_colitis-master/variant20230214.RData")
