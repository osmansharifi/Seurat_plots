#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("argparse"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("circlize"))
suppressPackageStartupMessages(library("ComplexHeatmap"))

parser <- ArgumentParser()

parser$add_argument("-l", "--logfc", type="character", metavar="<string>",
	required=TRUE, help="csv of log fold changes")
parser$add_argument("-p", "--pv", type="character", metavar="<string>",
	required=TRUE, help="csv of p-values")	
parser$add_argument("-s", "--save", type="character", metavar="<string>",
	required=TRUE, help="save location to save pdf of plot")
####
#parser$add_argument("-x", "--sex", type="character", metavar="<string>",
                    #required=TRUE, help="sex metadata on graph")
####
parser$add_argument("-t", "--title", type="character", metavar="<string>",
	required=TRUE, nargs='+', help="title for plot")
parser$add_argument("-r", "--rotate", action="store_true", default=FALSE,
	help="whether to rotate plot")

args <- parser$parse_args()
print(args$title)
print(paste(args$title, sep=""))
logfc_df <- read.csv(file=args$logfc)
pv_df <- read.csv(file=args$pv)
####
#x_df <- read.csv(file=args$sex)
####
m <- as.matrix(logfc_df[,4:dim(logfc_df)[2]])
rownames(m) <- logfc_df$SYMBOL
pvm <- as.matrix(pv_df[,4:dim(pv_df)[2]])
rownames(pvm) <- pv_df$SYMBOL
####
#x <- as.matrix(x_df[,4:dim(x_df)[2]])
#rownames(x) <- x_df$SYMBOL
####
col_fun = colorRamp2(c(min(m), 0.0, max(m)), c("#2166AC", "#EEEEEE", "#B2182B")) 


pdf(file=args$save)
map = grid.grabExpr(
	draw(
		Heatmap(
			m,
 			col = col_fun,
 			row_names_gp=gpar(fontsize=5),
 			column_names_gp=gpar(fontsize=6),
			row_names_max_width = max_text_width(rownames(m)),
 			heatmap_legend_param = list(title="logFC"),
			#column_split = m$tp,
 			cell_fun = function(j, i, x, y, width, height, fill) {
 				if( pvm[i, j] <= 0.05 ) {
 					grid.text(print("*"), x, y-height/2, gp = gpar(fontsize=8)) #grid.text(print("*"), x, y, gp = gpar(fontsize=9)) 
 				}
			},
			column_title = paste(args$title, sep="", collapse=" ")
			) 
	)
)
grid.newpage()
if( args$rotate ) {
	pushViewport(viewport(angle = -90))
} else {
	pushViewport(viewport(angle = 0))
}
grid.draw(map)
popViewport()
dev.off()