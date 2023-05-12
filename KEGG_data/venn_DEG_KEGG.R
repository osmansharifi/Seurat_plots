####################
## Load libraries ##
####################
suppressWarnings(suppressMessages(library(VennDiagram)))
suppressWarnings(suppressMessages(library(gridExtra)))

###############
## Load Data ##
###############
mouse <- read.csv("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/KEGG_data/mouse_total_kegg.csv", header = TRUE)
female_mouse <- filter(mouse, Sex == "males")
female_mouse <- filter(female_mouse, P.value <= 0.05)

female_human <- read.csv("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/KEGG_data/human_total_kegg.csv", header = TRUE)
female_human <- filter(female_human, P.value <= 0.05)
# Create a list to store the plots
plot_list <- list()

# Prepare data for each Cell_Type in male_mouse and female_human data frames
for (Cell_Type in unique(c(female_mouse$Cell_Type, female_human$Cell_Type))) {
  venn_data <- list(Human=female_human$Term[female_human$Cell_Type == Cell_Type],
                    Mouse=female_mouse$Term[female_mouse$Cell_Type == Cell_Type])
  
  # Calculate the number of overlapping genes
  overlap <- intersect(venn_data$Human, venn_data$Mouse)
  n <- length(overlap)
  
  # Calculate the total number of genes
  N1 <- length(venn_data$Human)
  N2 <- length(venn_data$Mouse)
  N <- N2 + N1
  
  options(scipen=999)
  # Perform the hypergeometric test
  pval <- phyper(n-1, N1, N2, N, lower.tail=FALSE)
  # Set the number of significant digits for p-values

  
  # Print the p-value
  print(pval)
  # Set the significance level
  alpha <- 0.05
  
  # Determine if the overlap is significant
  if (pval < alpha) {
    sig <- "Yes"
  } else {
    sig <- "No"
  }
  
  # Make plot for each Cell_Type and add it to the plot list
  venn.plot <- venn.diagram( venn_data, filename = NULL,
                             main=paste0(Cell_Type), 
                             main.fontface="bold",
                             sub=paste0("p-val=", round(pval, 3), ", Significant=", sig),
                             col=c("red","blue"),
                             fill=c("red","blue"),
                             alpha=0.3,
                             hyper.test = TRUE,
                             width = 2)
  plot_list[[Cell_Type]] <- venn.plot
}

# Arrange the plots in a grid
f1 <- grid.arrange(grobs = plot_list, ncol = 3, top=textGrob(expression(bold("KEGG terms of Human and Male Mouse corticies p-value <=0.05"))))
dev.off()
ggsave(file="/Users/osman/Documents/GitHub/snRNA-seq-pipeline/KEGG_data/human_male_mouse_kegg_venn.tiff",
       f1,
       width = 15, 
       height = 12)
