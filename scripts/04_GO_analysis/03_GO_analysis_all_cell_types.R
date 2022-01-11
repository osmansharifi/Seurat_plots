library(ggplot2)
library(dplyr)
# By Osman Sharifi & Viktoria Haghani

################################################################################
## Variables

# Paths
#csv_path <- "~/GitHub/snRNA-seq-pipeline/GO_data/GO_term_tables/master_go_data_males.csv"
#figure_path <- "~/GitHub/snRNA-seq-pipeline/figures/go_analysis/"
csv_path <- "~/Documents/GitHub/snRNA-seq-pipeline/GO_data/GO_term_tables/master_go_data_males_top3.csv"
figure_path <- "~/Documents/GitHub/snRNA-seq-pipeline/figures/go_analysis/"
# Names
plot_title <- "GO Analysis"
plot_subtitle <- "All Cell Types, Time Points, and Ontologies for Male Mice"

################################################################################

# Read in GO Data from master CSV file
go_data <- read.csv(file = csv_path)

#subset and sort GO terms
E18_GO <- go_data[which(go_data$Time.Point=='E18'),]
E18_GO <- E18_GO[order(E18_GO$Fisher),]
#for (celltype in E18_GO$Cell.Type){
 # order(E18_GO$Fisher)
 # return(E18_GO)
#}
p.adjust <- round(p.adjust(E18_GO$Fisher,method="BH"),digits = 4)
E18_GO=cbind(E18_GO,p.adjust)
E18_GO <- E18_GO[order(E18_GO$p.adjust),]
E18_GO <- filter(E18_GO, p.adjust < 0.05)


P30_GO <- go_data[which(go_data$Time.Point=='P30'),]
P30_GO <- P30_GO[order(P30_GO$Fisher),]
P60_GO <- go_data[which(go_data$Time.Point=='P60'),]
P60_GO <- P60_GO[order(P60_GO$Fisher),]
P120_GO <- go_data[which(go_data$Time.Point=='P120'),]
P120_GO <- P120_GO[order(P120_GO$Fisher),]

# Visualize GO Data
p1<-ggplot(data = E18_GO, aes(x = Cell.Type, y = Term, 
                        color = p.adjust, size = Significant)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1))+
  ggtitle("TopGO E18 male GO enrichment analysis")
ggsave(p1,
       filename = "/Users/osman/Documents/GitHub/snRNA-seq-pipeline/figures/go_analysis/E18_male_GO.pdf",
       height = 15, 
       width = 10)

p2<-ggplot(data = P30_GO, aes(x = Cell.Type, y = Term, 
                              color = -log10(`Fisher`), size = Significant)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("TopGO P30 male GO enrichment analysis")
ggsave(p2,
       filename = "/Users/osman/Documents/GitHub/snRNA-seq-pipeline/figures/go_analysis/P30_male_GO.pdf",
       height = 20, 
       width = 10)

p3<-ggplot(data = P60_GO, aes(x = Cell.Type, y = Term, 
                              color = -log10(`Fisher`), size = Significant)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("TopGO P60 male GO enrichment analysis")
ggsave(p3,
       filename = "/Users/osman/Documents/GitHub/snRNA-seq-pipeline/figures/go_analysis/P60_male_GO.pdf",
       height = 20, 
       width = 10)

p4<-ggplot(data = P120_GO, aes(x = Cell.Type, y = Term, 
                              color = -log10(`Fisher`), size = Significant)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("TopGO P120 male GO enrichment analysis")
ggsave(p4,
       filename = "/Users/osman/Documents/GitHub/snRNA-seq-pipeline/figures/go_analysis/P120_male_GO.pdf",
       height = 20, 
       width = 10)


#########################################################################################################
ggplot(go_data, aes(x = Metadata, y = Term, size = -log10(go_data[9]))) +
  expand_limits(y = 1) +
  geom_point(shape = 21) +
  scale_size(range = c(2.5,12.5)) +
  scale_fill_continuous(low = 'royalblue', high = 'red4') +
  xlab('') +
  ylab('Enrichment score') +
  labs(title = plot_title, subtitle = plot_subtitle) +
  geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
             linetype = c("dotted", "longdash", "solid"),
             colour = c("black", "black", "black"),
             size = c(0.5, 1.5, 3)) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = 'right',
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
    plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
    plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
    axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
    axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
    axis.title = element_text(size = 12, face = 'bold'),
    axis.title.x = element_text(size = 12, face = 'bold'),
    axis.title.y = element_text(size = 12, face = 'bold'),
    axis.line = element_line(colour = 'black'),
    # Legend
    legend.key = element_blank(), # removes the border
    legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
    legend.text = element_text(size = 14, face = "bold"), # Text size
    title = element_text(size = 14, face = "bold")) +
  coord_flip()

#ggplot2::ggsave(glue('{figure_path}master_go_data_Fisher.pdf'),
#                device = NULL,
#                height = 8.5,
#                width = 12)