library(Seurat)
library(dplyr)
library(scCustomize)
library(patchwork)
library(glue)
library(limma)
library(ggplot2)
library(ggrepel)

##################
## Load samples ##
##################
base_path <- '/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/09_mosiacism_analysis/'
load(glue('{base_path}/all.female.cortex.parsed.RData'))
mosaic.cortex <- subset(x = all.female.cortex, subset = Condition == 'MUTANT')
cluster <- subset(mosaic.cortex, idents = c("L2_3_IT", "L4", "L5", "L6","Pvalb", "Vip", "Sst","Sncg","Lamp5", "Oligo","Astro","Non-neuronal"))

###################################################################
## Limma DEG analysis of mutant vs wt cells from the het females ##
###################################################################
# Get expression info 
expr <- as.matrix(GetAssayData(cluster))
# Filter out genes that are 0 for every cell in this cluster
bad <- which(rowSums(expr) == 0)
expr <- expr[-bad,]
logcpm <- cpm(expr, prior.count=2, log=TRUE)
mm <- model.matrix(~0 + Mecp2_allele, data = cluster@meta.data)
y <- voom(expr, mm, plot = TRUE)
fit <- lmFit(y, mm)  
head(coef(fit)) # means in each sample for each gene
contr<- makeContrasts(c(Mecp2_alleleMecp2_MUT) - c(Mecp2_alleleMecp2_WT), levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contrasts = contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "M", n = Inf) # top 20 DE genes
length(which(top.table$adj.P.Val < 0.05))
DEG<-as.matrix(topTable(tmp, sort.by = "adj.P.Val", n = 100))
summary(decideTests(tmp))

###########################
## Prepare and Visualize ##
###########################
top.table$Gene <- rownames(top.table)
top.table$diffexpressed <- 'NO'
top.table$diffexpressed[top.table$logFC > 0 & top.table$adj.P.Val < 0.05] <- 'UP'
top.table$diffexpressed[top.table$logFC < 0 & top.table$adj.P.Val < 0.05] <- 'DOWN'
top.table$diffexpressed[top.table$adj.P.Val > 0.05] <- 'Not Sig'
top.table$delabel <- NA
thresh = head(arrange(top.table, adj.P.Val), 10)$adj.P.Val[10]
top.table$delabel[top.table$adj.P.Val <=thresh] <-(top.table$Gene[top.table$adj.P.Val<=thresh])

# Volcano Plot
ggplot(data = top.table, aes(x = logFC, y = -log(adj.P.Val), col = diffexpressed, label = delabel))+
  geom_point(size=3)+
  theme_minimal()+
  geom_text_repel()+
  scale_color_manual(values = c('blue', 'black', 'red'))+
  theme(text = element_text(size=16))




                                                   