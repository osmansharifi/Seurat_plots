library(dplyr)
library(glue)

#######################
## load parsed table ##
#######################
counts_table <- read.csv('/Users/osman/Desktop/LaSalle_lab/Rett_Data/E18/cell_parsing_outputs/cell_parsing_df_bodyadded.csv')

################################
## remove duplicated barcodes ##
################################
p120_M_counts <- filter(counts_table, Sex == 'M' & Timepoint == 'P120')
p120_M_counts <- p120_M_counts[!duplicated(p120_M_counts$Barcode) & !duplicated(p120_M_counts$Barcode, fromLast = TRUE), ]
p150_F_counts <- filter(counts_table, Sex == 'F' & Timepoint == 'P150')
p150_F_counts <- p150_F_counts[!duplicated(p150_F_counts$Barcode) & !duplicated(p150_F_counts$Barcode, fromLast = TRUE), ]
p60_M_counts <- filter(counts_table, Sex == 'M' & Timepoint == 'P60')
# Remove non-unique rows from the original p60_M_counts
p60_M_counts <- p60_M_counts[!duplicated(p60_M_counts$Barcode) & !duplicated(p60_M_counts$Barcode, fromLast = TRUE), ]
p60_F_counts <- filter(counts_table, Sex == 'F' & Timepoint == 'P60')
# Remove non-unique rows from the original p60_F_counts
p60_F_counts <- p60_F_counts[!duplicated(p60_F_counts$Barcode) & !duplicated(p60_F_counts$Barcode, fromLast = TRUE), ]
p30_M_counts <- filter(counts_table, Sex == 'M' & Timepoint == 'P30')
# Remove non-unique rows from the original p30_M_counts
p30_M_counts <- p30_M_counts[!duplicated(p30_M_counts$Barcode) & !duplicated(p30_M_counts$Barcode, fromLast = TRUE), ]
p30_F_counts <- filter(counts_table, Sex == 'F' & Timepoint == 'P30')
# Remove non-unique rows from the original p30_F_counts
p30_F_counts <- p30_F_counts[!duplicated(p30_F_counts$Barcode) & !duplicated(p30_F_counts$Barcode, fromLast = TRUE), ]
e18_M_counts <- filter(counts_table, Sex == 'M' & Timepoint == 'E18')
# Remove non-unique rows from the original e18_M_counts
e18_M_counts <- e18_M_counts[!duplicated(e18_M_counts$Barcode) & !duplicated(e18_M_counts$Barcode, fromLast = TRUE), ]
e18_F_counts <- filter(counts_table, Sex == 'F' & Timepoint == 'E18')
# Remove non-unique rows from the original e18_F_counts
e18_F_counts <- e18_F_counts[!duplicated(e18_F_counts$Barcode) & !duplicated(e18_F_counts$Barcode, fromLast = TRUE), ]

###################
## Export counts ##
###################
base_path <- '/Users/osman/Documents/GitHub/snRNA-seq-pipeline/scripts/09_mosiacism_analysis'
write.csv(p120_M_counts, file = glue('{base_path}/p120_M_counts.csv'))
write.csv(p150_F_counts, file = glue('{base_path}/p150_F_counts.csv'))
write.csv(p60_F_counts, file = glue('{base_path}/p60_F_counts.csv'))
write.csv(p60_M_counts, file = glue('{base_path}/p60_M_counts.csv'))
write.csv(p30_F_counts, file = glue('{base_path}/p30_F_counts.csv'))
write.csv(p30_M_counts, file = glue('{base_path}/p30_M_counts.csv'))
write.csv(e18_F_counts, file = glue('{base_path}/e18_F_counts.csv'))
write.csv(e18_M_counts, file = glue('{base_path}/e18_M_counts.csv'))

# Males
p120_M_counts <- read.csv(glue('{base_path}/p120_M_counts.csv'))
p60_M_counts <- read.csv(glue('{base_path}/p60_M_counts.csv'))
p30_M_counts <- read.csv(glue('{base_path}/p30_M_counts.csv'))
combined_male <- rbind(p120_M_counts, p60_M_counts, p30_M_counts)
# Identify non-unique rows based on Barcode
non_unique_rows <- combined_male[duplicated(combined_male$Barcode) | duplicated(combined_male$Barcode, fromLast = TRUE), ]
# Remove non-unique rows from the original e18_F_counts
combined_male <- combined_male[!duplicated(combined_male$Barcode) & !duplicated(combined_male$Barcode, fromLast = TRUE), ]
write.csv(combined_male, file = glue('{base_path}/combined_male_unique.csv'))

# Females
p150_F_counts <- read.csv(glue('{base_path}/p150_F_counts.csv'))
p60_F_counts <- read.csv(glue('{base_path}/p60_F_counts.csv'))
p30_F_counts <- read.csv(glue('{base_path}/p30_F_counts.csv'))
combined_female <- rbind(p150_F_counts, p60_F_counts, p30_F_counts)
# Identify non-unique rows based on Barcode
non_unique_rows <- combined_female[duplicated(combined_female$Barcode) | duplicated(combined_female$Barcode, fromLast = TRUE), ]
# Remove non-unique rows from the original e18_F_counts
combined_female <- combined_female[!duplicated(combined_female$Barcode) & !duplicated(combined_female$Barcode, fromLast = TRUE), ]
write.csv(combined_male, file = glue('{base_path}/combined_female_unique.csv'))
