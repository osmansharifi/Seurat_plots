###################
##Load Libraries ##
###################
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)

##########################
##Load and prepare data ##
##########################
male_matrix <- read.csv("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/DEG_data/LimmaVoomCC/male_mouse_logfc.csv")
rownames(male_matrix) = male_matrix$SYMBOL
male_matrix$SYMBOL <- NULL
male_matrix$X <- NULL
male_matrix = as.matrix(male_matrix)
pv_male <- read.csv("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/DEG_data/LimmaVoomCC/male_mouse_pv.csv")
rownames(pv_male) = pv_male$SYMBOL
pv_male$SYMBOL <- NULL
pv_male$X <- NULL
pv_male = as.matrix(pv_male)
male_metadata <- read.csv("/Users/osman/Documents/GitHub/snRNA-seq-pipeline/DEG_data/LimmaVoomCC/male_mouse_meta.csv")

####################
## Create Heatmap ##
####################

#Set color palette
polychrome_palette <- c("#5A5156FF","#E4E1E3FF","#F6222EFF","#FE00FAFF","#16FF32FF","#3283FEFF","#FEAF16FF","#B00068FF","#1CFFCEFF","#90AD1CFF","#2ED9FFFF","#DEA0FDFF","#AA0DFEFF","#F8A19FFF","#325A9BFF","#C4451CFF","#1C8356FF","#85660DFF","#B10DA1FF","#FBE426FF","#1CBE4FFF","#FA0087FF","#FC1CBFFF","#F7E1A0FF","#C075A6FF","#782AB6FF","#AAF400FF","#BDCDFFFF","#822E1CFF","#B5EFB5FF","#7ED7D1FF","#1C7F93FF","#D85FF7FF","#683B79FF","#66B0FFFF", "#3B00FBFF")
col_fun = colorRamp2(c(min(male_matrix), 0.0, max(male_matrix)), c("#2166AC", "#EEEEEE", "#B2182B")) 

# Set column and row annotations

row_ha = rowAnnotation(Genes = rownames(male_matrix))
column_ha = HeatmapAnnotation(`Cell Type` = male_metadata$Cell.Type, 
                              `Time Point` = male_metadata$Time.Point,
                              Sex = male_metadata$Sex,
                              col = list(`Time Point` = c("P30" = "#E6B8BFFF","P60"= "#CC7A88FF","P120"= "#B33E52FF", "P150" = "#990F26FF"), 
                                         Sex = c("Male" = "#0F8299FF", "Female" = "#3D0F99FF"),
                                         `Cell Type` = c("L2_3_IT" = polychrome_palette[1], "L4" = polychrome_palette[2], "L5"= polychrome_palette[3], "Pvalb"= polychrome_palette[4],"Vip"= polychrome_palette[5], "Sst"= polychrome_palette[6],"Sncg"= polychrome_palette[7], "Lamp5"= polychrome_palette[8], "Peri"= polychrome_palette[9], "Astro"= polychrome_palette[13])),
                              annotation_name_gp = gpar(fontsize = 20, fontface = 'bold', fontfamily = 'Times'))
column_ha@anno_list$`Cell Type`@color_mapping@levels <- c("L2_3_IT", "L4", "L5", "Pvalb","Vip", "Sst","Sncg", "Lamp5", "Peri", "Astro")
column_ha@anno_list$`Time Point`@color_mapping@levels<- c("P30", "P60", "P120")
# Create heatmap

pdf(file="Top Postnatal Male DEGs.pdf", height = 15, width = 16)
map = grid.grabExpr(
  draw(
    Heatmap(male_matrix, 
            name = "logFC", 
            top_annotation = column_ha, 
            col = col_fun, 
            row_names_gp=gpar(fontsize=18, fontface = 'bold', fontfamily = 'Times'),
            column_names_gp=gpar(fontsize=20, fontface = 'bold'),
            cluster_columns = FALSE, 
            column_title = "Top Postnatal Male DEGs", 
            heatmap_legend_param = list(title="logFC", 
                                        title_gp = gpar(fontsize = 18, fontface = 'bold', fontfamily = 'Times'), 
                                        labels_gp = gpar(fontsize = 12, fontface = 'bold', fontfamily = 'Times')), 
            cell_fun = function(j, i, x, y, width, height, fill) {
  if( pv_male[i, j] <= 0.05 ) {
    grid.text(print("*"), x, y-height/3, gp = gpar(fontsize=18, fontface = 'bold')) #grid.text(print("*"), x, y, gp = gpar(fontsize=9)) 
  }
})))
grid.newpage()
grid.draw(map)
dev.off()

