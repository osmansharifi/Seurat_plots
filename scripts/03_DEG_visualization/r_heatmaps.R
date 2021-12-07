#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("ComplexHeatmap"))

parser <- ArgumentParser()

parser$add_argument("-f", "--frame", type="character", metavar="<string>",
	required=FALSE, help="feather frame to plot")

args <- parser$parse_args()
print(args$frame)
print(file.exists(args$frame))

df <- read.csv(file=args$frame)
print(df)
# m <- as.matrix(df[,4:dim(df)[2]])
# rownames(m) <- df$gene
# 
# Heatmap(
# 	m,
# 	row_names_gp=gpar(fontsize=9),
# 	column_names_gp=gpar(fontsize=9),
# 	height=unit(15, "cm")
# )