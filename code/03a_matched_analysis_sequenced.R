##########################################################
# Name of file: 02a_analysis_sequenced.R
# Data release (if applicable):
# Original author(s): Chris Robertson chris.robertson@phs.scot
# Original date: 30 June 2022
# Latest update author (if not using version control) - Chris Robertson chris.robertson@phs.scot
# Latest update date (if not using version control) - 
# Latest update description (if not using version control)
# Type of script: Descriptive stats
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 3.6.1
# Description of content: Runs matched analysis on ba2/ba5 only
#                         run 01b_get_all_cases.R
#                         
# Approximate run time: Unknown
##########################################################

# Libraries
library("finalfit")
library(gtsummary)

source('00_functions.R')

remove(list=ls(pa="^z"))
output_list <- list()
output_list$a_begin <- a_begin #earliest start of data analysis
#output_list$a_end <- a_end
output_list$a_end <- as.Date("2022-09-30")  #specified end date for the letter
output_list$a_end_hosp <- a_end_hosp #RAPID admissions
output_list$a_end_death <- a_end_death #NRS registrations
output_list$a_end_wgs <- a_end_wgs
output_list$a_end_smr01 <- a_end_smr01
output_list$endpoint <- "covid_admission" # "covid_admission", "hosp_pos_test_emerg" , "covid_death_cert"
output_list$nsims <- 25

df_pos <- readRDS("output/temp/All_positive.rds") 
df_pos_orig <- df_pos

z_df <- df_pos %>% filter(date_ecoss_specimen <= output_list$a_end) %>%  #use this for the paper
  filter(variant %in% c("omicron_ba2","omicron_ba5") ) %>% 
  filter(in_eave==1) %>% 
  filter(lab %in% c("lh","nhs"))

if (output_list$endpoint == "covid_admission")  {
  z_df <- z_df %>%   filter(covid_admission_status != "covid_admission_admitted_before_test")
  z_df <- z_df %>%  mutate(event=covid_admission)
  output_list$a_end_analysis <- output_list$a_end_smr01
  }
if (output_list$endpoint == "hosp_pos_test_emerg")  {
  z_df <- z_df %>%   filter(In_Hosp_At_Test == "no")
  z_df <- z_df %>%  mutate(event=hosp_pos_test_emerg)
  output_list$a_end_analysis <- output_list$a_end_hosp-14 #subtract 14 to give everyone 14 days for an admission
  }
if (output_list$endpoint == "covid_death_cert")  {
  z_df <- z_df %>%  mutate(event=covid_death_cert)
  output_list$a_end_analysis <- output_list$a_end_death - 28  #subtract 28 to give everyone 28 days for a death
}

#start of data analysis when ba5 detected
output_list$a_begin_analysis <- z_df %>% filter(variant=="omicron_ba5") %>% dplyr::summarize(min_date=min(date_ecoss_specimen)) %>%  pull(min_date) 
z_df %>% group_by(variant) %>% dplyr::summarize(max_date=max(date_ecoss_specimen))  #latest dates


z_df <- z_df %>% filter(date_ecoss_specimen >= output_list$a_begin_analysis) %>% 
  filter(date_ecoss_specimen <=  output_list$a_end_analysis)

df_pos <- z_df
df_pos$Total = 'Total'

#get the date of any positive test before a_begin_analysis
z_pcr <- Positive_Tests %>% filter(date_ecoss_specimen < output_list$a_begin_analysis) %>% dplyr::select(-test_id)
z_lft <- lft %>% filter(date_ecoss_specimen < output_list$a_begin_analysis) %>% dplyr::select(-test_result, -test_type)
z <- bind_rows(z_pcr,z_lft) %>% arrange(EAVE_LINKNO, desc(date_ecoss_specimen))
z <- z %>% filter(!duplicated(EAVE_LINKNO))

z_df <- df_pos %>%  left_join(z, by="EAVE_LINKNO", suffix=c("","_prev")) %>% 
  mutate(days_from_last_pos = as.numeric(date_ecoss_specimen - date_ecoss_specimen_prev)) %>% 
  mutate(prev_pos_gp = case_when(is.na(days_from_last_pos) ~ "not_prev_pos",
                                 !is.na(days_from_last_pos) & days_from_last_pos <= 28 ~ "< 29 days",
                                 !is.na(days_from_last_pos) & days_from_last_pos <= 90 ~ "29-90 days",
                                 !is.na(days_from_last_pos) & days_from_last_pos <= 180 ~ "91-180 days",
                                 !is.na(days_from_last_pos) & days_from_last_pos <= 360 ~ "181-360 days",
                                 TRUE ~ "361+ days")) %>% 
  mutate(prev_pos_gp = factor(prev_pos_gp, levels=c("not_prev_pos","< 29 days","29-90 days","91-180 days","181-360 days","361+ days" )))
z_df <- z_df %>% mutate(prev_pos_variant = case_when(is.na(date_ecoss_specimen_prev) ~ "not_prev_pos",
                          !is.na(date_ecoss_specimen_prev) & date_ecoss_specimen_prev < as.Date("2020-12-15") ~ "wild_type",
                          !is.na(date_ecoss_specimen_prev) & date_ecoss_specimen_prev < as.Date("2021-05-15") ~ "alpha", 
                          !is.na(date_ecoss_specimen_prev) & date_ecoss_specimen_prev < as.Date("2021-12-20") ~ "delta",
                          TRUE ~ "omicron") ) %>% 
  mutate(prev_pos_variant = factor(prev_pos_variant, levels = c("not_prev_pos", "wild_type", "alpha","delta","omicron")))
df_pos <- z_df

saveRDS(df_pos, "output/temp/BA2_BA5_matched_analysis.rds") 
saveRDS(output_list, "output/temp/BA2_BA5_matched_output_list.rds") 

df_pos <- df_pos %>% mutate(age_match = case_when(age_year <= 4 ~ 4,
                                                  age_year >= 90 ~ 90,
                                                  TRUE ~ age_year))
df_pos <- df_pos %>% mutate(age_match = cut(age_match, breaks=seq(0,90, by=5)) )
df_pos <- df_pos  %>% mutate(prev_pos_gp = fct_recode(prev_pos_gp, "<=90 days" = "< 29 days", "<=90 days" = "29-90 days"))                           
#start of the matching process

results_list_tab <- list()
results_list_coefs <- list()
results_list_var <- list()

for (isims in 1:output_list$nsims) {
#isims <- 1

z_pos <- df_pos %>% filter(variant=="omicron_ba5") %>% dplyr::select(EAVE_LINKNO, date_ecoss_specimen, age_match,  lab)
z_neg <- df_pos %>% filter(variant=="omicron_ba2") %>% dplyr::select(EAVE_LINKNO, date_ecoss_specimen, age_match,  lab)

z_merge <- z_pos %>% inner_join(z_neg, by=c("date_ecoss_specimen", "age_match", "lab"), suffix=c("_case","_cont")) %>% 
  filter(!is.na(EAVE_LINKNO_cont))
z_merge <- z_merge %>%   mutate(random_id = runif(nrow(z_merge))) %>% 
  arrange(EAVE_LINKNO_case, random_id) %>% 
  filter(!duplicated(EAVE_LINKNO_case)) %>% dplyr::select(-random_id)
z_merge <- z_merge %>%   mutate(random_id = runif(nrow(z_merge))) %>% 
  arrange(EAVE_LINKNO_cont, random_id) %>% 
  filter(!duplicated(EAVE_LINKNO_cont)) %>% dplyr::select(-random_id)

z_pos <- z_merge %>% mutate(EAVE_LINKNO = EAVE_LINKNO_case) %>%  dplyr::select(EAVE_LINKNO, EAVE_LINKNO_case)
z_neg <- z_merge %>% mutate(EAVE_LINKNO = EAVE_LINKNO_cont) %>%  dplyr::select(EAVE_LINKNO, EAVE_LINKNO_case)
z_df <- bind_rows(z_pos,z_neg) %>% arrange(EAVE_LINKNO_case) %>% left_join(df_pos, by="EAVE_LINKNO")
z_df <- z_df %>% mutate(n_risk_gps = fct_drop(n_risk_gps),
                        prev_pos_gp = fct_drop(prev_pos_gp))

z_tab <- z_df %>% group_by(across(c(variant, output_list$endpoint))) %>% dplyr::summarise(N=n())
results_list_tab[[isims]] <- z_tab

#z <- clogit(event ~  variant + strata(EAVE_LINKNO_case) + simd2020_sc_quintile + n_risk_gps + immuno + shielding  + vs, data=z_df, method="efron")
z <- clogit(event ~  variant + strata(EAVE_LINKNO_case) + sex + simd2020_sc_quintile + n_risk_gps + immuno + shielding + prev_pos_gp + vs, data=z_df, method="efron")
#summary(z)
print(summary(z)$conf.int[1,])
z_names <- names(z$coefficients)
results_list_coefs[[isims]] <- data.frame(names = z_names, coefs =z$coefficients)
results_list_var[[isims]] <- data.frame(names = 1:(length(z_names)^2), var=as.numeric(z$var))


}

z <- ldply(results_list_tab, data.frame)
z_tab <- z %>% group_by(across(c(variant, output_list$endpoint))) %>% dplyr::summarise(N=round(mean(N)))

z <- ldply(results_list_coefs, data.frame)
z_coefs <- z %>% group_by(names) %>% dplyr::summarise(coef = mean(coefs, na.rm=T), var_coef = var(coefs,na.rm=T))
#z_names <- z$names[1:nrow(z_coefs)]

z <- ldply(results_list_var, data.frame)
z_var <- z %>% group_by(names) %>% dplyr::summarise(var = mean(var, na.rm=T)) 
z_var <- matrix(z_var$var, nrow=nrow(z_coefs), ncol=nrow(z_coefs)) 
dimnames(z_var) <- list(z_names,z_names)
z_var <- diag(z_var)
z_var <- data.frame(names=names(z_var),var=z_var)

z_coefs <- z_coefs %>% left_join(z_var, by="names") %>% 
  mutate(se = sqrt(var_coef+var)) %>% 
  mutate(lcl = coef - 1.96*se, ucl = coef+1.96*se) %>%
  mutate(OR = exp(coef), LCL=exp(lcl), UCL=exp(ucl))

saveRDS(list(table=z_tab, coefficients=z_coefs), paste0("./output/matched_analysis_",output_list$endpoint,".rds") )

        