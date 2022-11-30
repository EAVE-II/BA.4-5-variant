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
# Description of content: Runs expected hospitalisations analysis
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
output_list$a_begin_analysis <- as.Date("2022-04-01")
output_list$a_begin <- a_begin
#output_list$a_end <- a_end
output_list$a_end <- as.Date("2022-09-30")  #specified end date for the letter
output_list$a_end_hosp_expected <- a_end_smr01 + 7 #add 7 days to give a bit more follow up time
output_list$a_end_hosp_14days_expected <- a_end_hosp-7 #a_end_hosp is max date from all_hospitalisation (RAPID) - go back 7 days for completness
output_list$a_end_deaths <- a_end_death-7 #max date  - go back 7 days for completness
output_list$a_end_wgs <- a_end_wgs
output_list$a_end_smr01 <- a_end_smr01

df_pos <- readRDS("output/temp/All_positive.rds") 
df_pos_orig <- df_pos
df_pos <- df_pos_orig %>%  filter(date_ecoss_specimen >= output_list$a_begin_analysis) #select the data in the study period
df_pos <- df_pos %>%  filter(date_ecoss_specimen <= output_list$a_end) #select the data in the study period
df_pos$Total = 'Total'
df_pos$days <- df_pos$days - min(df_pos$days) + 1
z_gps <- 1:(trunc(max(df_pos$days)/14)+1)
df_pos <- df_pos %>% 
  mutate(days_gp = cut(days, breaks = c(-1, (z_gps*14 - 1)) , labels=paste0("week_",z_gps)))
#line below only done for the letter
df_pos <- df_pos %>% mutate(days_gp = fct_recode(days_gp, "week_13" = "week_14"))

table(df_pos$days_gp, df_pos$covid_admission, exclude=NULL)
table(df_pos$days_gp, df_pos$variant, exclude=NULL)
df_pos %>% group_by(days_gp) %>% dplyr::summarise(mindate=min(date_ecoss_specimen), maxdate=max(date_ecoss_specimen)) %>% as.data.frame()

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

saveRDS(df_pos, "output/temp/All_positive_analysis.rds") 
saveRDS(output_list, "output/temp/output_list_analysis.rds") 


#plot out the trends in the endpoints by date_specimen
z <- df_pos %>% group_by(date_ecoss_specimen) %>% 
  dplyr::summarise(covid_admission=sum(covid_admission), hosp_pos_test_emerg = sum(hosp_pos_test_emerg), covid_death=sum(covid_death), covid_death_cert=sum(covid_death_cert)) %>% 
  pivot_longer(cols= -date_ecoss_specimen)
z %>% filter(name %in% c("covid_admission","hosp_pos_test_emerg")) %>% ggplot(aes(x=date_ecoss_specimen, y=value, colour=name)) +geom_point() +
  geom_vline(xintercept=output_list$a_end_smr01, size=0.5, linetype=3) + geom_vline(xintercept=output_list$a_end_hosp_expected, size=0.5, linetype=3) 
ggsave(paste0("./output/covid_admissions_day.png"), width=14, height=10, unit="cm")

z %>% filter(name %in% c("covid_death","covid_death_cert")) %>% ggplot(aes(x=date_ecoss_specimen, y=value, colour=name)) +geom_point() 
ggsave(paste0("./output/covid_deaths_day.png"), width=14, height=10, unit="cm")

#comparing ba4/5 to ba2 in the period post a_begin_analysis
z_df <- df_pos %>% filter(date_ecoss_specimen >=  output_list$a_begin_analysis ) %>% 
  filter(variant %in% c("omicron_ba2","omicron_ba4","omicron_ba5") ) %>% 
  filter(in_eave==1) %>% 
  filter(lab != "lft") %>% 
  mutate(days_gp = fct_recode(days_gp, week_2="week_1"))

explanatory <- c("Total","sex","age_gp","n_risk_gps", "simd2020_sc_quintile", "lab", "vs", "prev_pos_gp","prev_pos_variant","days_gp")
summary_tbl_wt_chrt <- summary_factorlist_wt(z_df, "variant", explanatory = explanatory, wght="wght_all") 
names(summary_tbl_wt_chrt)[1:2] <- c('Characteristic', 'Levels')
summary_tbl_wt_chrt$Characteristic[duplicated(summary_tbl_wt_chrt$Characteristic)] <- ''
summary_tbl_wt_chrt[1, 'Levels'] <- ''
write.csv(summary_tbl_wt_chrt , paste0("./output/summary_table_variant.csv"), row.names = F)

z_df_sub <- z_df %>% mutate(event=if_else(variant=="omicron_ba5",1L,0L)) %>% 
  filter(variant %in% c("omicron_ba2","omicron_ba5")) %>% 
  mutate(vs = fct_relevel(vs, "v3_26+")) #%>% filter(age_year >= 18)
z <- glm(event ~ age_gp + sex + simd2020_sc_quintile + n_risk_gps + lab+ vs + prev_pos_gp + days_gp, 
         data=z_df_sub, family=binomial)
#z <- mgcv::gam(variant=="omicron_ba5" ~ s(days) + s(age_year) + sex + simd2020_sc_quintile + n_risk_gps + lab+ vs + prev_pos_gp, data=z_df, family=binomial,
#subset = variant %in% c("omicron_ba2","omicron_ba5"))

saveRDS(z , paste0("./output/ba5_glm.RDS"))

summary(z)
#drop1(z)
z1 <- update(z, ~ . - vs)
#saveRDS(z1 , paste0("./output/ba5_glm_vs.RDS"))
anova(z1,z,test="Chisq") %>% saveRDS(paste0("./output/ba5_glm_vs_anova.RDS"))
z_tab <- gtsummary::tbl_regression(z, exponentiate=TRUE) %>% 
  add_n(location="level") %>% add_nevent(location="level") 
saveRDS(as_kable_extra(z_tab) , paste0("./output/ba5_glm_tab.RDS"))

z_df_sub <- z_df %>% mutate(event=if_else(variant=="omicron_ba4",1L,0L)) %>% 
  filter(variant %in% c("omicron_ba2","omicron_ba4")) %>% 
  mutate(vs = fct_relevel(vs, "v3_26+"))
z <- glm(variant=="omicron_ba4" ~ age_gp + sex + simd2020_sc_quintile + n_risk_gps + lab+ vs+ prev_pos_gp + days_gp,
         data=z_df_sub, family=binomial )
saveRDS(z , paste0("./output/ba4_glm.RDS"))
z1 <- update(z, ~ . - vs)
#saveRDS(z1 , paste0("./output/ba4_glm_vs.RDS"))
anova(z1,z,test="Chisq") %>% saveRDS(paste0("./output/ba4_glm_vs_anova.RDS"))
summary(z)

z_tab <- gtsummary::tbl_regression(z, exponentiate=TRUE) %>% 
  add_n(location="level") %>% add_nevent(location="level")
saveRDS(as_kable_extra(z_tab) , paste0("./output/ba4_glm_tab.RDS"))

##########################################################################################
#Confirmed Covid Admission to hospital Person years to event by variant
#Individuals not in hospital and then get covid.
z.rv <- "covid_admission" 
z.rv.time <- "time_to_hosp" 

fmla.plot <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  variant"))
z_df <- df_pos %>% #filter(date_ecoss_specimen >=  output_list$a_begin_analysis ) %>% 
  filter(variant %in% c("omicron_ba2","omicron_ba4","omicron_ba5") ) %>% 
  filter(in_eave==1) %>% 
  filter(lab != "lft") %>%  
  filter(covid_admission_status != "covid_admission_admitted_before_test") %>% 
  filter(date_ecoss_specimen <= output_list$a_end_hosp_expected)%>% 
  filter(lab %in% c("lh","nhs"))
z.tab <- pyears(fmla.plot, data=z_df, data.frame=TRUE)$data %>% 
  mutate(rate_100=event/pyears*100)
write.csv(z.tab, paste0("./output/pyears_by_variant.csv"), row.names = F)

# cumulative incidence curves   - uses z_df
z.survfit <- survfit(fmla.plot, data=z_df)

png("./output/cuminc_variant.png", width=8, height=6, unit="in", res=72)
plot(z.survfit, fun="event", col=1:3, xlab="Days from test to confirmed covid hospital admission",
     ylab="Risk")
z_names <- gsub("variant=","",names(z.survfit$strata))
legend("topright",col=1:3, lty=1, legend=z_names )
dev.off()

#calculated the expected hospitalisations for ba4/5 based upon ba2   - uses z_df
fmla.coxph <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  pspline(days) + pspline(age_year) + sex + lab +
                               simd2020_sc_quintile + n_risk_gps +  immuno + shielding + prev_pos_gp + vs"))
z_df$n_risk_gps <- fct_drop(z_df$n_risk_gps)
z_df$days_gp <- fct_drop(z_df$days_gp)
z <- coxph(fmla.coxph, data=z_df, subset=variant=="omicron_ba2")
summary(z)
saveRDS(z,"./output/coxph_model_ba2.rds")
z_plot1 <- plot_HR(z, "pspline(days)")
z_plot2 <- plot_HR(z, "pspline(age_year)")
z_g <- gridExtra::arrangeGrob(z_plot1,z_plot2,ncol=2, nrow=1)
grid::grid.draw(z_g)
ggsave(paste0("./output/coxph_model_ba2.png"), plot=z_g, width=14, height=10, unit="cm")

z_pred <- predict(z, newdata=z_df, type="expected")
summary(z_pred) #z_pred can be greater than 1 so calculate the probability of an admission and scale
z_df$pred <- 1-exp(-z_pred)  #expected prob on an admissions within time period 
z_correction <-  sum(z_df$pred[z_df$variant=="omicron_ba2"])/sum(z_df[z_df$variant=="omicron_ba2",z.rv])
z_df$pred <- z_df$pred/z_correction # correct so that the sum of expected in comparator group is equal to observed
z_tab_pred <- z_df %>% group_by(variant) %>% 
  dplyr::summarise(N=n(), Person_Years = round(sum(get(z.rv.time))/365.25,1), Covid_Hosp= sum(get(z.rv)),
                   Expected = sum(pred)) %>% 
  as.data.frame() 
z_ci <- epitools::pois.byar(z_tab_pred$Covid_Hosp, z_tab_pred$Expected)
z_tab_pred <- bind_cols(z_tab_pred, z_ci[,c("rate","lower","upper")])
write.csv(z_tab_pred, paste0("./output/pyears_by_variant_exp_hosp_eave.csv"), row.names = F)

#Hospitalisation rate by age group  - uses z_df
fmla.plot <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  variant + age_gp"))
z.tab <- pyears(fmla.plot, data=z_df, data.frame=TRUE)$data
z.tab <- z.tab %>% mutate(rate_100py=event/(pyears/100))
z_ci <- epitools::pois.byar(z.tab$event, z.tab$pyears/100)
z.tab <- bind_cols(z.tab,z_ci) %>% mutate(lower=if_else(lower<0,0,lower))
write.csv(z.tab, paste0("./output/pyears_by_variant_age.csv"), row.names = F)

z_position <- position_dodge(width=0.2)
z.tab %>%  filter(variant != "omicron_ba4") %>% 
  ggplot(aes(x=age_gp, colour=variant)) + geom_point(aes(y=rate_100py), position=z_position) + 
  geom_errorbar(aes(ymin=lower, ymax=upper, width=0.2), position = z_position) +
  labs(x="Age Group", y="Rate Per 100 person years", colour="variant") + ylim(0,500)
ggsave(paste0("./output/rate_variant_age.png"), width=14, height=10, unit="cm")


#######################################################################
# comparing ba5/ba2
#  Changes z_df

z_df <- df_pos %>% filter(date_ecoss_specimen >=  as.Date("2022-04-14") ) %>% 
  filter(variant %in% c("omicron_ba2","omicron_ba5") ) %>% 
  filter(in_eave==1) %>% 
  filter(lab != "lft") %>% 
  filter(covid_admission_status != "covid_admission_admitted_before_test") %>% 
  filter(date_ecoss_specimen <= output_list$a_end_hosp_expected) %>% 
  filter(lab %in% c("lh","nhs"))
z_df$n_risk_gps <- fct_drop(z_df$n_risk_gps)
z_df$days_gp <- fct_drop(z_df$days_gp)

fmla.plot <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  variant"))
z.tab <- pyears(fmla.plot, data=z_df, data.frame=TRUE)$data %>% 
  mutate(rate_100=event/pyears*100)
write.csv(z.tab, paste0("./output/pyears_by_variant_ba2_ba5.csv"), row.names = F)

fmla.coxph <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(days) + pspline(age_year) + sex + lab +
                               simd2020_sc_quintile + n_risk_gps +  immuno + shielding + prev_pos_gp + variant + vs"))
z1 <- coxph(fmla.coxph, data=z_df)
summary(z1)
saveRDS(z1,"./output/coxph_model_ba2_ba5.rds")
z <- z1
z_plot1 <- plot_HR(z1, "pspline(days)")
z_plot2 <- plot_HR(z1, "pspline(age_year)")
z_g <- gridExtra::arrangeGrob(z_plot1,z_plot2,ncol=2, nrow=1)
grid::grid.draw(z_g)
ggsave(paste0("./output/coxph_model_ba2_ba5.png"), plot=z_g, width=14, height=10, unit="cm")
fmla.coxph <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(days) + pspline(age_year) + sex + lab +
                               simd2020_sc_quintile + n_risk_gps +  immuno + shielding + prev_pos_gp + variant + variant:vs"))
z1 <- coxph(fmla.coxph, data=z_df)
summary(z1)
saveRDS(z1,"./output/coxph_model_ba2_ba5_vs.rds")
z_int_test <- anova(z,z1,test="Chisq")
saveRDS(z_int_test,"./output/coxph_model_ba2_ba5_vs_interaction.rds")
z_plot1 <- plot_HR(z1, "pspline(days)")
z_plot2 <- plot_HR(z1, "pspline(age_year)")
z_g <- gridExtra::arrangeGrob(z_plot1,z_plot2,ncol=2, nrow=1)
grid::grid.draw(z_g)
ggsave(paste0("./output/coxph_model_ba2_ba5_vs.png"), plot=z_g, width=14, height=10, unit="cm")

#######################################################################
# comparing ba5/ba2
#  Changes z_df
#sensitivity using time selections

z_df <- df_pos %>% filter(days_gp %in% c("week_4","week_5") ) %>% 
  filter(variant %in% c("omicron_ba2","omicron_ba5") ) %>% 
  filter(in_eave==1) %>% 
  filter(lab != "lft") %>% 
  filter(covid_admission_status != "covid_admission_admitted_before_test") %>% 
  filter(date_ecoss_specimen <= output_list$a_end_hosp_expected) %>% 
  filter(lab %in% c("lh","nhs"))
z_df$n_risk_gps <- fct_drop(z_df$n_risk_gps)
z_df$days_gp <- fct_drop(z_df$days_gp)

fmla.plot <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  variant"))
z.tab <- pyears(fmla.plot, data=z_df, data.frame=TRUE)$data %>% 
  mutate(rate_100=event/pyears*100)
write.csv(z.tab, paste0("./output/pyears_by_variant_ba2_ba5_wk45.csv"), row.names = F)

fmla.coxph <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(days) + pspline(age_year) + sex + lab +
                               simd2020_sc_quintile + n_risk_gps +  immuno + shielding + prev_pos_gp + variant + vs"))
z1 <- coxph(fmla.coxph, data=z_df)
summary(z1)
saveRDS(z1,"./output/coxph_model_ba2_ba5_wk45.rds")
z_plot1 <- plot_HR(z1, "pspline(days)")
z_plot2 <- plot_HR(z1, "pspline(age_year)")
z_g <- gridExtra::arrangeGrob(z_plot1,z_plot2,ncol=2, nrow=1)
grid::grid.draw(z_g)
ggsave(paste0("./output/coxph_model_ba2_ba5_wk45.png"), plot=z_g, width=14, height=10, unit="cm")

z_df <- df_pos %>% filter(days_gp %in% c("week_6","week_7", "week_8") ) %>% 
  filter(variant %in% c("omicron_ba2","omicron_ba5") ) %>% 
  filter(in_eave==1) %>% 
  filter(lab != "lft") %>% 
  filter(covid_admission_status != "covid_admission_admitted_before_test") %>% 
  filter(lab %in% c("lh","nhs"))
z_df$n_risk_gps <- fct_drop(z_df$n_risk_gps)
z_df$days_gp <- fct_drop(z_df$days_gp)

fmla.plot <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  variant"))
z.tab <- pyears(fmla.plot, data=z_df, data.frame=TRUE)$data %>% 
  mutate(rate_100=event/pyears*100)
write.csv(z.tab, paste0("./output/pyears_by_variant_ba2_ba5_wk67.csv"), row.names = F)

fmla.coxph <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(days) + pspline(age_year) + sex + lab +
                               simd2020_sc_quintile + n_risk_gps +  immuno + shielding + prev_pos_gp + variant + vs"))
z1 <- coxph(fmla.coxph, data=z_df)
summary(z1)
saveRDS(z1,"./output/coxph_model_ba2_ba5_wk67.rds")
z_plot1 <- plot_HR(z1, "pspline(days)")
z_plot2 <- plot_HR(z1, "pspline(age_year)")
z_g <- gridExtra::arrangeGrob(z_plot1,z_plot2,ncol=2, nrow=1)
grid::grid.draw(z_g)
ggsave(paste0("./output/coxph_model_ba2_ba5_wk67.png"), plot=z_g, width=14, height=10, unit="cm")


###########################################################################################
#
#Repeat of analysis using admission within 14 days of a positive test  - copy of code with some changes to selections graphs and output files
#
###########################################################################################

#14 days of a positive test  Covid Admission to hospital Person years to event by variant
#Individuals not in hospital and then get covid.
z.rv <- "hosp_pos_test_emerg" 
z.rv.time <- "Time.To.Hosp" 

fmla.plot <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  variant"))
z_df <- df_pos %>% #filter(date_ecoss_specimen >=  output_list$a_begin_analysis ) %>% 
  filter(variant %in% c("omicron_ba2","omicron_ba4","omicron_ba5") ) %>% 
  filter(in_eave==1) %>% 
  filter(In_Hosp_At_Test == "no") %>% 
  filter(date_ecoss_specimen <= output_list$a_end_hosp_14days_expected) %>% 
  filter(lab %in% c("lh","nhs"))
z.tab <- pyears(fmla.plot, data=z_df, data.frame=TRUE)$data %>% 
  mutate(rate_100=event/pyears*100)
write.csv(z.tab, paste0("./output/pyears_by_variant_14days.csv"), row.names = F)

# cumulative incidence curves   - uses z_df
z.survfit <- survfit(fmla.plot, data=z_df)

png("./output/cuminc_variant_14days.png", width=8, height=6, unit="in", res=72)
plot(z.survfit, fun="event", col=1:3, xlab="Days from test to confirmed covid hospital admission",
     ylab="Risk")
z_names <- gsub("variant=","",names(z.survfit$strata))
legend("topright",col=1:3, lty=1, legend=z_names )
dev.off()

#calculated the expected hospitalisations for ba4/5 based upon ba2   - uses z_df
fmla.coxph <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  pspline(days) + pspline(age_year) + sex + lab +
                               simd2020_sc_quintile + n_risk_gps +  immuno + shielding + prev_pos_gp + vs"))
z_df$n_risk_gps <- fct_drop(z_df$n_risk_gps)
z_df$days_gp <- fct_drop(z_df$days_gp)
z <- coxph(fmla.coxph, data=z_df, subset=variant=="omicron_ba2")
summary(z)
saveRDS(z,"./output/coxph_model_ba2_14days.rds")
z_plot1 <- plot_HR(z, "pspline(days)")
z_plot2 <- plot_HR(z, "pspline(age_year)")
z_g <- gridExtra::arrangeGrob(z_plot1,z_plot2,ncol=2, nrow=1)
grid::grid.draw(z_g)
ggsave(paste0("./output/coxph_model_ba2_14days.png"), plot=z_g, width=14, height=10, unit="cm")

z_pred <- predict(z, newdata=z_df, type="expected")
summary(z_pred) #z_pred can be greater than 1 so calculate the probability of an admission and scale
z_df$pred <- 1-exp(-z_pred)  #expected prob on an admissions within time period 
z_correction <-  sum(z_df$pred[z_df$variant=="omicron_ba2"])/sum(z_df[z_df$variant=="omicron_ba2",z.rv])
z_df$pred <- z_df$pred/z_correction # correct so that the sum of expected in comparator group is equal to observed
z_tab_pred <- z_df %>% group_by(variant) %>% 
  dplyr::summarise(N=n(), Person_Years = round(sum(get(z.rv.time))/365.25,1), Covid_Hosp= sum(get(z.rv)),
                   Expected = sum(pred)) %>% 
  as.data.frame() 
z_ci <- epitools::pois.byar(z_tab_pred$Covid_Hosp, z_tab_pred$Expected)
z_tab_pred <- bind_cols(z_tab_pred, z_ci[,c("rate","lower","upper")])
write.csv(z_tab_pred, paste0("./output/pyears_by_variant_exp_hosp_eave_14days.csv"), row.names = F)

#Hospitalisation rate by age group  - uses z_df
fmla.plot <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  variant + age_gp"))
z.tab <- pyears(fmla.plot, data=z_df, data.frame=TRUE)$data
z.tab <- z.tab %>% mutate(rate_py=event/(pyears))
z_ci <- epitools::pois.byar(z.tab$event, z.tab$pyears)
z.tab <- bind_cols(z.tab,z_ci) %>% mutate(lower=if_else(lower<0,0,lower))
write.csv(z.tab, paste0("./output/pyears_by_variant_age_14days.csv"), row.names = F)

z_position <- position_dodge(width=0.2)
z.tab %>%  filter(variant != "omicron_ba4") %>% 
  ggplot(aes(x=age_gp, colour=variant)) + geom_point(aes(y=rate_py), position=z_position) + 
  geom_errorbar(aes(ymin=lower, ymax=upper, width=0.2), position = z_position) +
  labs(x="Age Group", y="Rate Per 1 person years", colour="variant") + ylim(0,10)
ggsave(paste0("./output/rate_variant_age_14days.png"), width=14, height=10, unit="cm")

# comparing ba5/ba2
#  Changes z_df

z_df <- df_pos %>% filter(date_ecoss_specimen >=  as.Date("2022-04-14") ) %>% filter(days_gp != "week_10") %>% 
  filter(variant %in% c("omicron_ba2","omicron_ba5") ) %>% 
  filter(in_eave==1) %>% 
  filter(In_Hosp_At_Test == "no") %>% 
  filter(date_ecoss_specimen <= output_list$a_end_hosp_14days_expected) %>% 
  filter(lab %in% c("lh","nhs"))
z_df$n_risk_gps <- fct_drop(z_df$n_risk_gps)
z_df$days_gp <- fct_drop(z_df$days_gp)

fmla.plot <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  variant"))
z.tab <- pyears(fmla.plot, data=z_df, data.frame=TRUE)$data %>% 
  mutate(rate_100=event/pyears*100)
write.csv(z.tab, paste0("./output/pyears_by_variant_ba2_ba5_14days.csv"), row.names = F)

fmla.coxph <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(days) + pspline(age_year) + sex + lab +
                               simd2020_sc_quintile + n_risk_gps +  immuno + shielding + prev_pos_gp + variant + vs"))
z1 <- coxph(fmla.coxph, data=z_df)
summary(z1)
saveRDS(z1,"./output/coxph_model_ba2_ba5_14days.rds")
z <- z1
z_plot1 <- plot_HR(z1, "pspline(days)")
z_plot2 <- plot_HR(z1, "pspline(age_year)")
z_g <- gridExtra::arrangeGrob(z_plot1,z_plot2,ncol=2, nrow=1)
grid::grid.draw(z_g)
ggsave(paste0("./output/coxph_model_ba2_ba5_14days.png"), plot=z_g, width=14, height=10, unit="cm")

fmla.coxph <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(days) + pspline(age_year) + sex + lab +
                               simd2020_sc_quintile + n_risk_gps +  immuno + shielding + prev_pos_gp + variant + variant:vs"))
z1 <- coxph(fmla.coxph, data=z_df)
summary(z1)
saveRDS(z1,"./output/coxph_model_ba2_ba5_vs_14days.rds")
z_int_test <- anova(z,z1,test="Chisq")
saveRDS(z_int_test,"./output/coxph_model_ba2_ba5_vs_interaction_14days.rds")
z_plot1 <- plot_HR(z1, "pspline(days)")
z_plot2 <- plot_HR(z1, "pspline(age_year)")
z_g <- gridExtra::arrangeGrob(z_plot1,z_plot2,ncol=2, nrow=1)
grid::grid.draw(z_g)
ggsave(paste0("./output/coxph_model_ba2_ba5_vs_14days.png"), plot=z_g, width=14, height=10, unit="cm")

###########################################################################################
#
#Repeat of analysis using confirmed covid death certificate death  - copy of code with some changes to selections graphs and output files
#
###########################################################################################

#14 days of a positive test  Covid Admission to hospital Person years to event by variant
#Individuals not in hospital and then get covid.
z.rv <- "covid_death_cert" #any position in death certificate
z.rv.time <- "Time.To.Death" 

fmla.plot <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  variant"))
z_df <- df_pos %>% #filter(date_ecoss_specimen >=  output_list$a_begin_analysis ) %>% 
  filter(variant %in% c("omicron_ba2","omicron_ba4","omicron_ba5") ) %>% 
  filter(in_eave==1) %>% 
  filter(In_Hosp_At_Test == "no") %>% 
  filter(date_ecoss_specimen <= output_list$a_end_deaths) %>% 
  filter(lab %in% c("lh","nhs")) %>% 
  mutate(age_gp = fct_recode(age_gp, "0-39"="0-11", "0-39"="12-19", "0-39"="20-39"))
z.tab <- pyears(fmla.plot, data=z_df, data.frame=TRUE)$data %>% 
  mutate(rate_100=event/pyears*100)
write.csv(z.tab, paste0("./output/pyears_by_variant_death.csv"), row.names = F)

# cumulative incidence curves   - uses z_df
z.survfit <- survfit(fmla.plot, data=z_df)

png("./output/cuminc_variant_death.png", width=8, height=6, unit="in", res=72)
plot(z.survfit, fun="event", col=1:3, xlab="Days from test to confirmed covid hospital admission",
     ylab="Risk")
z_names <- gsub("variant=","",names(z.survfit$strata))
legend("topright",col=1:3, lty=1, legend=z_names )
dev.off()

#calculated the expected hospitalisations for ba4/5 based upon ba2   - uses z_df
fmla.coxph <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  pspline(days) + pspline(age_year) + sex + lab +
                               simd2020_sc_quintile + n_risk_gps +  immuno + shielding + prev_pos_gp + vs"))
z_df$n_risk_gps <- fct_drop(z_df$n_risk_gps)
z_df$days_gp <- fct_drop(z_df$days_gp)
z <- coxph(fmla.coxph, data=z_df, subset=variant=="omicron_ba2")
summary(z)
saveRDS(z,"./output/coxph_model_ba2_death.rds")
z_plot1 <- plot_HR(z, "pspline(days)")
z_plot2 <- plot_HR(z, "pspline(age_year)")
z_g <- gridExtra::arrangeGrob(z_plot1,z_plot2,ncol=2, nrow=1)
grid::grid.draw(z_g)
ggsave(paste0("./output/coxph_model_ba2_death.png"), plot=z_g, width=14, height=10, unit="cm")

z_pred <- predict(z, newdata=z_df, type="expected")
summary(z_pred) #z_pred can be greater than 1 so calculate the probability of an admission and scale
z_df$pred <- 1-exp(-z_pred)  #expected prob on an admissions within time period 
z_correction <-  sum(z_df$pred[z_df$variant=="omicron_ba2"])/sum(z_df[z_df$variant=="omicron_ba2",z.rv])
z_df$pred <- z_df$pred/z_correction # correct so that the sum of expected in comparator group is equal to observed
z_df$pred <- z_pred
z_tab_pred <- z_df %>% group_by(variant) %>% 
  dplyr::summarise(N=n(), Person_Years = round(sum(get(z.rv.time))/365.25,1), Covid_Death= sum(get(z.rv)),
                   Expected = sum(pred)) %>% 
  as.data.frame() 
z_ci <- epitools::pois.byar(z_tab_pred$Covid_Death, z_tab_pred$Expected)
z_tab_pred <- bind_cols(z_tab_pred, z_ci[,c("rate","lower","upper")])
write.csv(z_tab_pred, paste0("./output/pyears_by_variant_exp_death_eave.csv"), row.names = F)

#Hospitalisation rate by age group  - uses z_df
fmla.plot <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  variant + age_gp"))
z.tab <- pyears(fmla.plot, data=z_df, data.frame=TRUE)$data
z.tab <- z.tab %>% mutate(rate_100py=event/(pyears/100))
z_ci <- epitools::pois.byar(z.tab$event, z.tab$pyears/100)
z.tab <- bind_cols(z.tab,z_ci) %>% mutate(lower=if_else(lower<0,0,lower))
write.csv(z.tab, paste0("./output/pyears_by_variant_age_death.csv"), row.names = F)

z_position <- position_dodge(width=0.2)
z.tab %>%  filter (variant != "omicron_ba4") %>% 
  ggplot(aes(x=age_gp, colour=variant)) + geom_point(aes(y=rate_100py), position=z_position) + 
  geom_errorbar(aes(ymin=lower, ymax=upper, width=0.2), position = z_position) +
  labs(x="Age Group", y="Rate Per 100 person years", colour="variant") + ylim(0,20)
ggsave(paste0("./output/rate_variant_age_death.png"), width=14, height=10, unit="cm")

# comparing ba5/ba2
#  Changes z_df

z_df <- df_pos %>% filter(date_ecoss_specimen >=  as.Date("2022-04-14") )  %>% filter(days_gp != "week_10") %>% 
  filter(variant %in% c("omicron_ba2","omicron_ba5") ) %>% 
  filter(in_eave==1) %>% 
  filter(In_Hosp_At_Test == "no") %>% 
  filter(date_ecoss_specimen <= output_list$a_end_deaths) %>% 
  filter(lab %in% c("lh","nhs")) %>% 
  mutate(age_gp = fct_recode(age_gp, "0-39"="0-11", "0-39"="12-19", "0-39"="20-39"))
z_df$n_risk_gps <- fct_drop(z_df$n_risk_gps)
z_df$days_gp <- fct_drop(z_df$days_gp)

fmla.plot <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~  variant"))
z.tab <- pyears(fmla.plot, data=z_df, data.frame=TRUE)$data %>% 
  mutate(rate_100=event/pyears*100)
write.csv(z.tab, paste0("./output/pyears_by_variant_ba2_ba5_death.csv"), row.names = F)

fmla.coxph <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(days) + pspline(age_year) + sex + lab +
                               simd2020_sc_quintile + n_risk_gps +  immuno + shielding + prev_pos_gp + variant + vs"))
z1 <- coxph(fmla.coxph, data=z_df)
summary(z1)
saveRDS(z1,"./output/coxph_model_ba2_ba5_death.rds")
z_plot1 <- plot_HR(z1, "pspline(days)")
z_plot2 <- plot_HR(z1, "pspline(age_year)")
z_g <- gridExtra::arrangeGrob(z_plot1,z_plot2,ncol=2, nrow=1)
grid::grid.draw(z_g)
ggsave(paste0("./output/coxph_model_ba2_ba5_death.png"), plot=z_g, width=14, height=10, unit="cm")
#z_tab <- gtsummary::tbl_regression(z1, exponentiate=TRUE) %>% 
#  add_n(location="level") %>% add_nevent(location="level") 
#termplot(z1,terms="pspline(days)",se=T)
#termplot(z1,terms="pspline(age_year)",se=T)
#fmla.coxph <- as.formula(paste("Surv(",z.rv.time,",",z.rv,") ~   pspline(days) + pspline(age_year) + sex +
#simd2020_sc_quintile + n_risk_gps +  immuno + shielding + prev_pos_gp + variant + variant:vs"))
#z1 <- coxph(fmla.coxph, data=z_df)
#summary(z1)
#saveRDS(z1,"./output/coxph_model_ba2_ba5_vs_death.rds")
