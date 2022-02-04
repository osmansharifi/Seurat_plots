---
#title: "alleler visualization"
#author: "Osman Sharifi"
---
library(ggplot2)
library(reshape2)

counts_table <- read.table("Desktop/miscellaneous/sequence.alleler", sep="\t", header=FALSE)
names(counts_table) <- c("Barcode", "UMI", "WT", "MUT", "Body")
counts_table$count <- seq.int(nrow(counts_table))
Condition <- c("WT", " MUT", "Body")
counts_table$Conditions <- Condition
counts_table.m<- melt(counts_table,id.vars='count', measure.vars=c('WT','MUT','Body'))
ggplot(counts_table.m, aes(x=variable, y=value, fill = variable))+
    geom_bar(stat="identity")+
    #geom_errorbar(aes(ymin=value-2, ymax=value+2))+
    scale_fill_brewer(palette="Dark2")+
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(face = "bold", size = 14),
          plot.title = element_text(face = "bold", size = 18, hjust = 0.5)) +
    ggtitle("WT vs MUT cells") +
    theme_minimal()

  