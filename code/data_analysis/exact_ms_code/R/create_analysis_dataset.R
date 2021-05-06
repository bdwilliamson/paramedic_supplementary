## create analysis dataset
library("tidyverse")
library("cowplot")
theme_set(theme_cowplot())
## -----------------------------------------------------------------------------------
## load in the data, clean it up
## -----------------------------------------------------------------------------------
data_dir <- "data/"
full_data <- read_csv(file = paste0(data_dir, "RSMp01_counts.csv"))

## check each relative abundance to make sure it's a valid taxon (not genus, for example)
full_data %>% names
get_nm_indices <- function(x, dat) dat %>% names %>% startsWith(x) %>% which 
taxon_nms <- c("Lacto", "Gardnerella", "Prev", "Lepto", "Veillon", "Sneath", "Dialist",
               "BVAB", "Mega", "Clostridiales", "Anaero", "Aero", "Parv", "Anaero",
               "Pepton", "Fineg", "Gemella", "Bifi", "Strep", "Egger", "Mycopl", 
               "Atop", "Eubac", "Mobilun", "Fusobac", "Porphyrom", "Actin", "Coryne", 
               "Mobil", "Varibac", "Bacteroid", "Campy", "Staph", "Moryella", "Bacill", "Jonque", 
               "Rumino", "Haemo", "Blauti", "Faecal", "Kocur", "Beta", "Gamma", "Brevi",
               "Lachno", "Neisser", "Escheri", "Clostru", "Klebi", "Urea", "Syner") %>% unique

all_taxon_nm_indices <- as.matrix(taxon_nms) %>% apply(1, function(x) get_nm_indices(x, full_data))
lapply(all_taxon_nm_indices, function(x) (full_data %>% names)[x])
(full_data %>% names)[all_taxon_nm_indices[[1]]]

# Hash table: qPCR name = BR16S name(s)
# aero = Aerococcuschristensenii
# atop = Atopobiumvaginae
# bvab2 = BVAB2                     
# dial1 = Dialistermicraerophilus                    
# dial2 = Dialister prop/sp type 2 (Dialisterpropionicifaciens + Dialistersp_type2)                     
# egg = Eggerthellasp_type1                  
# gemella = Gemella asaccharolytica              
# gvag  = Gardnerellavaginalis
# lcrisp = Lactobacilluscrispatus                       
# liners = Lactobacillusiners                        
# ljens = Lactobacillusjensenii                       
# lepto = Lepto/Sneathia (Leptotrichiaamnionii + Sneathiasanguinegens + Sneathia)                      
# mega = Megasphaera (Megasphaerasp_type1 + Megasphaerasp_type2)
# myco = Mycoplasumhominis
# parv1 + parv2 = Parvamonas sp type 1 + Parvamonas sp type 2(Parvimonasmicra)
# porphasacch = P. asacch/uenonis (Porphyromonasasaccharolyticaueno)             
# porphtype1 = Porphyromonassp_type1                    
# porphbenn =  Porphyromonasbennonis                   
# prevo = Prevotella

# remove any qPCRs that match to multiple br16Ss or zero br16s; 
# sum up any qPCRs that match to the same 16S
# this means removing dial2, lepto, mega, porhpasacch, prevo (multiple) and gemella (zero)
# summing up parv1 and parv2 to match to parvimonasmicra
# and remove identifiers
slim_data <- full_data %>% 
  select(-id, -label, -specimen, -case, -cohort, -group, -lactotest, -viral_load_quantifiable, -viral_load, -br16s) %>% 
  mutate(parv = parv1 + parv2) %>% 
  select(-dial2, -lepto, -mega, -parv1, -parv2, -gemella, -porphasacch, -prevo)
slim_data %>% names
# reorder so that qPCR goes first
input_data <- slim_data %>% 
  select(aero:porphbenn, parv, Lactobacillusiners:Synergistaceae)
qpcr_inds <- 1:13
br_inds <- 14:length(input_data %>% names)
pcr_to_16s_nm_match <- tibble(pcr_nm = names(input_data)[qpcr_inds],
                              br16s_nm = c("Aerococcuschristensenii", "Atopobiumvaginae", "BVAB2", 
                                           "Dialistermicraerophilus",
                                           "Eggerthellasp_type1", "Gardnerellavaginalis",
                                           "Lactobacilluscrispatus", "Lactobacillusiners", "Lactobacillusjensenii",
                                           "Mycoplasmahominis",
                                           "Porphyromonassp_type1", "Porphyromonasbennonis", "Parvimonasmicra"))
## may have to do some pre-processing to get the correct column indices
pcr_plus_br_inds <- unlist(lapply(as.list(pcr_to_16s_nm_match$br16s_nm), function(x) grep(paste0("^", x, "$"), names(input_data))))
## process the data, yielding qpcr and br16s matrices
qpcr <- input_data %>% 
  select(tidyselect::all_of(qpcr_inds)) 
br <- input_data %>% 
  select(tidyselect::all_of(qpcr_inds), tidyselect::all_of(pcr_plus_br_inds), !tidyselect::all_of(pcr_plus_br_inds)) %>% 
  select(!tidyselect::all_of(qpcr_inds)) 

## ----------------------------------------------------------------------------------------
## filter based on prevalence
## ----------------------------------------------------------------------------------------
## check how many of each taxon there are
colSums(br)
## compute the number of samples in which a taxon appears at least once
prevalence <- apply(X = br, MARGIN = 2, FUN = function(x) sum(x > 0))
## make a data frame with taxonomy and total read counts
prev_df <- tibble(Prevalence = prevalence, TotalAbundance = colSums(br), Taxon = names(br))
## plot
prev_df %>% 
  ggplot(aes(x = TotalAbundance, y = Prevalence / nrow(br), color = Taxon)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +
  geom_hline(yintercept = 0.1, alpha = 0.5, linetype = 2) +
  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +
  xlab("Total Abundance") +
  ylab("Prevalence [Frac. Samples]") +
  guides(color = FALSE)
## define prevalence filter as 5% of total samples
prev_number <- 0.05
prev_threshold <- prev_number * nrow(br)
prev_threshold
## execute prevalence filter
keep_taxa <- ifelse(prev_df$Prevalence >= prev_threshold, TRUE, FALSE)
sum(keep_taxa)
pruned_br <- br %>% 
  select(all_of(prev_df$Taxon[keep_taxa]))
## define prevalence filter as 10% of total samples
prev_threshold_1 <- 0.1 * nrow(br)
prev_threshold_1
## execute prevalence filter
keep_taxa_1 <- ifelse(prev_df$Prevalence >= prev_threshold_1, TRUE, FALSE)
sum(keep_taxa_1)
## ----------------------------------------------------------------------------------------
## save off the datasets
## ----------------------------------------------------------------------------------------
analysis_data <- list(qpcr = qpcr, br16S = br, filtered_br16S = pruned_br, case_control = full_data %>% select(case))
saveRDS(analysis_data, paste0(data_dir, "analysis_data.rds"))
