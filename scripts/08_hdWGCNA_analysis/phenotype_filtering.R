#Read phenotype file
phenotype <- readxl::read_xlsx("/Users/osman/Downloads/Mecp2e1_sn-RNA seq samples_disease score_weights.xlsx", col_names = TRUE)

#Keep only cortex data
phenotype <- phenotype[phenotype$region != "HYP",]

# create an empty vector to store the values of Time_point
sample_name <- vector(mode = "character", length = nrow(phenotype))

# loop through each row in the phenotype dataframe
for(i in 1:nrow(phenotype)) {
  if(phenotype$mouse[i] == "19-0106") {
    sample_name[i] <- "MUT_F_P150_CORT1"
  } else if(phenotype$mouse[i] == "19-0107") {
    sample_name[i] <- "WT_F_P150_CORT1"
  } else if(phenotype$mouse[i] == "19-0108") {
    sample_name[i] <- "MUT_F_P150_CORT2"
  } else if(phenotype$mouse[i] == "19-0111") {
    sample_name[i] <- "WT_F_P150_CORT2"
  } else if(phenotype$mouse[i] == "19-0112") {
    sample_name[i] <- "WT_M_P120_CORT1"
  } else if(phenotype$mouse[i] == "19-0113") {
    sample_name[i] <- "WT_M_P120_CORT2"
  } else if(phenotype$mouse[i] == "19-0114") {
    sample_name[i] <- "MUT_M_P120_CORT1"
  } else if(phenotype$mouse[i] == "19-0117") {
    sample_name[i] <- "MUT_M_P120_CORT2"
  } else if(phenotype$mouse[i] == "19-0132") {
    sample_name[i] <- "MUT_F_P150_CORT3"
  } else if(phenotype$mouse[i] == "19-0133") {
    sample_name[i] <- "MUT_F_P150_CORT4"
  } else if(phenotype$mouse[i] == "19-0134") {
    sample_name[i] <- "WT_F_P150_CORT3"
  } else if(phenotype$mouse[i] == "19-0136") {
    sample_name[i] <- "WT_F_P150_CORT4"
  } else if(phenotype$mouse[i] == "19-1982") {
    sample_name[i] <- "MUT_F_P60_CORT1"
  } else if(phenotype$mouse[i] == "19-1992") {
    sample_name[i] <- "WT_F_P60_CORT1"
  } else if(phenotype$mouse[i] == "19-1993") {
    sample_name[i] <- "WT_F_P60_CORT2"
  } else if(phenotype$mouse[i] == "19-1994") {
    sample_name[i] <- "MUT_F_P60_CORT2"
  } else if(phenotype$mouse[i] == "19-2652") {
    sample_name[i] <- "MUT_F_P30_CORT1"
  } else if(phenotype$mouse[i] == "19-2653") {
    sample_name[i] <- "WT_F_P30_CORT1"
  } else if(phenotype$mouse[i] == "19-2654") {
    sample_name[i] <- "WT_M_P30_CORT1"
  } else if(phenotype$mouse[i] == "19-2655") {
    sample_name[i] <- "WT_F_P30_CORT2"
  } else if(phenotype$mouse[i] == "19-2656") {
    sample_name[i] <- "MUT_F_P30_CORT2"
  } else if(phenotype$mouse[i] == "19-2658") {
    sample_name[i] <- "MUT_M_P30_CORT1"
  } else if(phenotype$mouse[i] == "19-2659") {
    sample_name[i] <- "MUT_M_P30_CORT2"
  } else if(phenotype$mouse[i] == "19-2662") {
    sample_name[i] <- "WT_M_P30_CORT2"
  } else if(phenotype$mouse[i] == "19-4133") {
    sample_name[i] <- "MUT_M_P60_CORT1"
  } else if(phenotype$mouse[i] == "19-4134") {
    sample_name[i] <- "WT_M_P60_CORT1"
  } else if(phenotype$mouse[i] == "19-4135") {
    sample_name[i] <- "MUT_M_P60_CORT2"
  } else if(phenotype$mouse[i] == "19-4136") {
    sample_name[i] <- "WT_M_P60_CORT2"
  }
}

# add the time_points vector as a new column to the phenotype dataframe
phenotype$Sample_name <- sample_name

phenotype[phenotype == "?"] <- NA

write.csv(phenotype, "phenotype.csv", row.names = FALSE)
