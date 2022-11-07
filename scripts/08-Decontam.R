#........................................................................
# ANALYSIS OF AMPLICON DATA: Decontamination of LCA results
#........................................................................

# Currently this script cannot be automated and is voyage specific
# Usage: Rscript scripts/Decontam.R -v RSV5 -a 16S -o custom

suppressPackageStartupMessages(library(tidyverse)) 
suppressPackageStartupMessages(library(RColorBrewer)) 
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))

# Define options for command line
option_list = list(
  make_option(c("-v", "--voyage"), action="store", default=NA, type='character',
              help="voyage identifier code"),
  make_option(c("-a", "--assay"), action="store", default=NA, type='character',
              help="assay, e.g. '16S' or 'MiFish"),
  make_option(c("-o", "--option"), action="store", default=NA, type='character',
              help="nt or custom blast database"))  

opt = parse_args(OptionParser(option_list=option_list))

voyage <- opt$voyage
assay  <- opt$assay
option <- opt$option



#......................................................................................
# WE CALL THE SCRIPT WE NEED BASED ON THE OPTIONS INPUT
#......................................................................................

# Run the analysis by executing the function above
if(option == "nt"){
  
  # Read in filtered LCA results
  lca_tab <- read_csv(paste0("05-taxa/LCA_out/LCA_filtered_", voyage, "_", assay, "_nt.csv"))
  
  # Mark all potential contaminant ASV sequences in new column
  lca_tab$Contam <- "False"
  
  # All ASV sequences identified in any control samples -excluding bleach controls
  # lca_tab$Contam[lca_tab$MT1_BC_1 >0] <- "True"
  # lca_tab$Contam[lca_tab$MT1_BC_2 >0] <- "True"
  # lca_tab$Contam[lca_tab$MT1_BC_3 >0] <- "True"
  # lca_tab$Contam[lca_tab$MT1_BC_4 >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_EB_1 >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_EB_2 >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_EB_3 >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_EB_4 >0] <- "True"
  lca_tab$Contam[lca_tab$NTC1 >0] <- "True"
  lca_tab$Contam[lca_tab$NTC2 >0] <- "True"
  lca_tab$Contam[lca_tab$NTC3 >0] <- "True"
  # lca_tab$Contam[lca_tab$RS1_CL_BC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_CL_EB >0] <- "True"
  # lca_tab$Contam[lca_tab$RS1_IM_BC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_IM_EB >0] <- "True"
  # lca_tab$Contam[lca_tab$RS1_ME_BC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_ME_EB >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_10_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_11_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_12_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_13_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_14_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_15_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_16_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_17_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_18_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_19_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_20_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_21_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_22_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_23_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_24_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_1_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_2_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_3_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_4_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_5_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_6_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_7_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_8_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_9_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_CL_L1_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_CL_L2_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_CL_L3_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_CL_S1_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_CL_S2_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_CL_S3_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_CL_S4_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_IM_L1_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_IM_L2_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_IM_L3_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_IM_S1_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_IM_S2_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_IM_S3_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_IM_S4_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_ME_L1_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_ME_L2_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_ME_L3_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_ME_S1_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_ME_S2_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_ME_S3_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_ME_S4_WC >0] <- "True"
  
  lca_tab <- lca_tab %>%
    relocate(OTU, .before = domain) %>%
    rename(ASV = OTU)
  
  # Write final output with contam labels
  write_csv(lca_tab, file = paste0("05-taxa/", voyage, "_", assay, "_contam_table_nt.csv"))
  
  # Create final file with no contaminants
  nocontam <- read_csv(file = paste0("05-taxa/", voyage, "_", assay, "_contam_table_nt.csv"))
  nocontam$Contam
  nocontam <- subset(nocontam, Contam=="FALSE")
  
  nocontam <- nocontam %>% 
    select(where(~ any(. != 0)))
  
  write_csv(nocontam, file = paste0("05-taxa/", voyage, "_", assay, "_nocontam_nt.csv"))
  
  
  
}

if(option == "custom"){
  
  # Read in filtered LCA results
  lca_tab <- read_delim(paste0("05-taxa/LCA_out/", voyage, "_", assay, "_LCA.tsv"))
  
  # Mark all potential contaminant ASV sequences in new column
  lca_tab$Contam <- "False"
  
  # All ASV sequences identified in any control samples -excluding bleach controls
  # lca_tab$Contam[lca_tab$MT1_BC_1 >0] <- "True"
  # lca_tab$Contam[lca_tab$MT1_BC_2 >0] <- "True"
  # lca_tab$Contam[lca_tab$MT1_BC_3 >0] <- "True"
  # lca_tab$Contam[lca_tab$MT1_BC_4 >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_EB_1 >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_EB_2 >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_EB_3 >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_EB_4 >0] <- "True"
  lca_tab$Contam[lca_tab$NTC1 >0] <- "True"
  lca_tab$Contam[lca_tab$NTC2 >0] <- "True"
  lca_tab$Contam[lca_tab$NTC3 >0] <- "True"
  # lca_tab$Contam[lca_tab$RS1_CL_BC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_CL_EB >0] <- "True"
  # lca_tab$Contam[lca_tab$RS1_IM_BC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_IM_EB >0] <- "True"
  # lca_tab$Contam[lca_tab$RS1_ME_BC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_ME_EB >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_10_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_11_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_12_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_13_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_14_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_15_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_16_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_17_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_18_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_19_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_20_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_21_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_22_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_23_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_24_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_1_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_2_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_3_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_4_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_5_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_6_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_7_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_8_WC >0] <- "True"
  lca_tab$Contam[lca_tab$MT1_9_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_CL_L1_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_CL_L2_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_CL_L3_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_CL_S1_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_CL_S2_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_CL_S3_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_CL_S4_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_IM_L1_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_IM_L2_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_IM_L3_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_IM_S1_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_IM_S2_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_IM_S3_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_IM_S4_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_ME_L1_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_ME_L2_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_ME_L3_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_ME_S1_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_ME_S2_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_ME_S3_WC >0] <- "True"
  lca_tab$Contam[lca_tab$RS1_ME_S4_WC >0] <- "True"
  
  lca_tab <- lca_tab %>%
    relocate(OTU, .before = domain) %>%
    rename(ASV = OTU)
  
  # Write final output with contam labels
  write_csv(lca_tab, file = paste0("05-taxa/", voyage, "_", assay, "_contam_table.csv"))
  
  # Create final file with no contaminants
  nocontam <- read_csv(file = paste0("05-taxa/", voyage, "_", assay, "_contam_table.csv"))
  nocontam$Contam
  nocontam <- subset(nocontam, Contam=="FALSE")
  
  nocontam <- nocontam %>% 
    select(where(~ any(. != 0)))
  
  write_csv(nocontam, file = paste0("05-taxa/", voyage, "_", assay, "_nocontam.csv"))
  
}

### Done!
print(paste0(voyage, " ", assay, " decontamination done!"))
