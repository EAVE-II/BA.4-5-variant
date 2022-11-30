##########################################################
# Name of file: 02c_comparison_sequenced.R
# Data release (if applicable):
# Original author(s): Chris Robertson chris.robertson@phs.scot
# Original date: 30 June 2022
# Latest update author (if not using version control) - Chris Robertson chris.robertson@phs.scot
# Latest update date (if not using version control) - 
# Latest update description (if not using version control)
# Type of script: Descriptive stats
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 3.6.1
# Description of content: Runs expected hospitalisations analysis
#                         run 01b_get_all_cases.R
#                         
# Approximate run time: Unknown
##########################################################

# Libraries
library("finalfit")

source('00_functions.R')

remove(list=ls(pa="^z"))
output_list <- list()
output_list$a_begin_analysis <- as.Date("2022-04-01")
output_list$a_begin <- a_begin
output_list$a_end_wgs <- a_end_wgs
output_list$a_end_smr01 <- a_end_smr01
output_list$a_end <- as.Date("2022-09-30")  #specified end date for the letter

df_pos <- readRDS("output/temp/All_positive.rds") #positive tests
df_pos <- df_pos %>% filter(date_ecoss_specimen <= output_list$a_end)  #used for paper
df_cohort <- EAVE_cohort %>% left_join(dplyr::select(rg, EAVE_LINKNO, n_risk_gps), by="EAVE_LINKNO")
df_cohort <- df_cohort %>% mutate(n_risk_gps = fct_explicit_na(n_risk_gps, na_level = "Unknown"))
#get vaccination status at a_begin_analysis
a_begin_analysis <- output_list$a_begin_analysis
df_cohort <- df_cohort %>%  
  left_join(Vaccinations, by="EAVE_LINKNO" )
df_cohort <- df_cohort %>%  mutate(vs1 = case_when(is.na(date_vacc_1) | date_vacc_1 > a_begin_analysis ~ "uv",
                                         date_vacc_1 <= a_begin_analysis &  date_vacc_1 > a_begin_analysis - 28 ~ "v1_0:3",
                                         TRUE ~ "v1_4+"),
                         vs2 = case_when(is.na(date_vacc_2) | date_vacc_2 > a_begin_analysis ~ "uv",
                                         date_vacc_2 <= a_begin_analysis &  date_vacc_2 > a_begin_analysis - 174 ~ "v2_0:24",
                                         TRUE ~ "v2_25+"),
                         vs3 = case_when(is.na(date_vacc_3) | date_vacc_3 > a_begin_analysis ~ "uv",
                                         date_vacc_3 <= a_begin_analysis &  date_vacc_3 > a_begin_analysis - 83 ~ "v3_0:11",
                                         date_vacc_3 <= a_begin_analysis-84 &  date_vacc_3 > a_begin_analysis - 111 ~ "v3_12:15",
                                         date_vacc_3 <= a_begin_analysis-112 &  date_vacc_3 > a_begin_analysis - 139 ~ "v3_16:19",
                                         date_vacc_3 <= a_begin_analysis-140 &  date_vacc_3 > a_begin_analysis - 181 ~ "v3_20:25",
                                         TRUE ~ "v3_26+"),
                         vs4 = case_when(is.na(date_vacc_4) | date_vacc_4 > a_begin_analysis ~ "uv",
                                         date_vacc_4 <= a_begin_analysis &  date_vacc_4 > a_begin_analysis - 13 ~ "v4_0:1",
                                         date_vacc_4 <= a_begin_analysis-14 &  date_vacc_4 > a_begin_analysis - 55~ "v4_2:7",
                                         TRUE ~ "v4_8+") ) %>% 
  mutate(vs = case_when(vs4 !="uv" ~ vs4,
                        vs3 !="uv" ~ vs3,
                        vs2 !="uv" ~ vs2,
                        TRUE ~ vs1))
z_levs <- c("uv","v1_0:3","v1_4+","v2_0:24","v2_25+",
            "v3_0:11","v3_12:15","v3_16:19", "v3_20:25", "v3_26+", "v4_0:1","v4_2:7","v4_8+")
df_cohort <- df_cohort %>% mutate(vs=factor(vs, levels=z_levs))
#merge_levels
df_cohort <- df_cohort %>% mutate(vs=fct_recode(vs, "v1" = "v1_0:3", "v1" = "v1_4+"))

df_cohort <- df_cohort %>% 
  mutate(vacc_type = if_else(is.na(vacc_type) , "uv", as.character(vacc_type)) ) %>% 
  filter(vacc_type != "UNK") %>% 
  mutate(vacc_type = factor(vacc_type, levels=c("uv","AZ","Mo","PB") ) ) %>% 
  mutate(vt = fct_cross(vs,vacc_type, sep="_")) %>% 
  mutate(vt = fct_recode(vt, "uv" ="uv_uv", "uv" = "uv_AZ", "uv" = "uv_Mo","uv" = "uv_PB"))
df_cohort <- df_cohort %>% mutate(age_gp = cut(ageYear, breaks = c(-1, 11, 19, 39, 59, 74, 130),
                                     labels=c("0-11", "12-19","20-39","40-59","60-74","75+")))
df_cohort <- df_cohort %>% dplyr::rename(sex=Sex)


#comparison of the cohort with tested positive and sequenced since a_begin
explanatory <- c("sex","age_gp","n_risk_gps", "simd2020_sc_quintile", "vs", "eave_weight")
z_df_pos <- df_pos %>% dplyr::select_at(all_of(explanatory)) %>% 
  filter(!is.na(simd2020_sc_quintile)) %>% 
  mutate(n_risk_gps = fct_drop(n_risk_gps)) %>% 
  mutate(group="positive")
z_df_cohort <- df_cohort %>% dplyr::select_at(all_of(explanatory))%>% 
  mutate(group="population") %>% 
  mutate(sex=case_when(sex=="F" ~ "Female",
                       sex=="M" ~ "Male"))
z_df_seq <- df_pos %>%   filter(variant != "not_sequenced") %>% 
  dplyr::select_at(all_of(explanatory)) %>% 
  filter(!is.na(simd2020_sc_quintile)) %>% 
  mutate(n_risk_gps = fct_drop(n_risk_gps)) %>% 
  mutate(group="sequenced")

z_df <- z_df_cohort %>% bind_rows(z_df_pos) %>% bind_rows(z_df_seq)
z_df$Total <- "Total"

explanatory <- c("Total","sex","age_gp", "simd2020_sc_quintile", "n_risk_gps", "vs")
summary_tbl_wt_chrt <- summary_factorlist_wt(z_df, "group", explanatory = explanatory, wght="eave_weight") 
names(summary_tbl_wt_chrt)[1:2] <- c('Characteristic', 'Levels')
summary_tbl_wt_chrt$Characteristic[duplicated(summary_tbl_wt_chrt$Characteristic)] <- ''
summary_tbl_wt_chrt[1, 'Levels'] <- ''
write.csv(summary_tbl_wt_chrt , paste0("./output/comparison_table_population.csv"), row.names = F)


#comparison of the tested positive and sequenced in different periods
explanatory <- c("sex","age_gp","n_risk_gps", "simd2020_sc_quintile", "vs", "test_type", "lab", "variant", "eave_weight")
#period 1
z_df_pos <- df_pos %>% filter(date_ecoss_specimen < a_begin_analysis) %>% 
  dplyr::select_at(all_of(explanatory)) %>% 
  filter(!is.na(simd2020_sc_quintile)) %>% 
  mutate(n_risk_gps = fct_drop(n_risk_gps)) %>% 
  mutate(group="Period 1 positive")
z_df_seq <- df_pos %>%   filter(variant != "not_sequenced") %>% filter(date_ecoss_specimen < a_begin_analysis) %>% 
  dplyr::select_at(all_of(explanatory)) %>% 
  filter(!is.na(simd2020_sc_quintile)) %>% 
  mutate(n_risk_gps = fct_drop(n_risk_gps)) %>% 
  mutate(group="Period 1 sequenced")
z_df <- z_df_pos %>% bind_rows(z_df_seq)
#period 2
z_df_pos <- df_pos %>% filter(date_ecoss_specimen >= a_begin_analysis) %>% 
  dplyr::select_at(all_of(explanatory)) %>% 
  filter(!is.na(simd2020_sc_quintile)) %>% 
  mutate(n_risk_gps = fct_drop(n_risk_gps)) %>% 
  mutate(group="Period 2 positive")
z_df_seq <- df_pos %>%   filter(variant != "not_sequenced") %>% filter(date_ecoss_specimen >= a_begin_analysis) %>% 
  dplyr::select_at(all_of(explanatory)) %>% 
  filter(!is.na(simd2020_sc_quintile)) %>% 
  mutate(n_risk_gps = fct_drop(n_risk_gps)) %>% 
  mutate(group="Period 2 sequenced")
z_df <- z_df %>% bind_rows(z_df_pos) %>% bind_rows(z_df_seq)
z_df$Total <- "Total"

explanatory <- c("Total","sex","age_gp", "simd2020_sc_quintile", "n_risk_gps", "vs", "test_type", "lab", "variant")
summary_tbl_wt_chrt <- summary_factorlist_wt(z_df, "group", explanatory = explanatory, wght="eave_weight") 
names(summary_tbl_wt_chrt)[1:2] <- c('Characteristic', 'Levels')
summary_tbl_wt_chrt$Characteristic[duplicated(summary_tbl_wt_chrt$Characteristic)] <- ''
summary_tbl_wt_chrt[1, 'Levels'] <- ''
write.csv(summary_tbl_wt_chrt , paste0("./output/comparison_table_pos_seq.csv"), row.names = F)
