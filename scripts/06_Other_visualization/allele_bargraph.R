---
#title: "alleler visualization"
#author: "Osman Sharifi"
---
library(ggplot2)
library(reshape2)

counts_table <- read.table("/Users/osman/Desktop/miscellaneous/sequence.alleler", sep="\t", header=FALSE)
names(counts_table) <- c("Barcode", "UMI", "WT", "MUT", "Body")
counts_table$count <- seq.int(nrow(counts_table))
Condition <- c("WT", " MUT", "Body")
counts_table$sums <- Condition
counts_table.m<- melt(counts_table,id.vars='count', measure.vars=c('WT','MUT','Body'))
final_table <- subset(counts_table.m, value > 0)
# dotplot of cell parsing
p1 <- ggplot(final_table, aes(x=variable, y=value, fill = variable))+
  geom_dotplot(
  binaxis = "y", stackdir = "center", binpositions = "all", stackgroups = TRUE, binwidth = 1/7) +
  scale_fill_brewer(palette="Dark2")+
  xlab("Condition of nuclei") + ylab("Mecp2 transcript counts per nucleus") +
  labs(fill = "Condition") +
  ggtitle("Parsing WT and MUT nuclei") +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(face = "bold", size = 14),
        plot.title = element_text(face = "bold", size = 18, hjust = 0.5))+
  theme_linedraw()
p1
