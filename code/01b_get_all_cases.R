##########################################################
# Name of file: 01b_get_all_cases.R
# Data release (if applicable):
# Original author(s): Chris Robertson chrisobertson@chris.robertson@phs.scot
# Original date: 24 June 2022
# Latest update author (if not using version control) - Chris Robertson chris.robertson@phs.scot
# Latest update date (if not using version control) - 
# Latest update description (if not using version control)
# Type of script: Descriptive stats
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 3.6.1
# Description of content: sets up analysis for all sequenced positive tests
#                         run 01a_Input_Data_VarSeq.R
#                         
# Approximate run time: Unknown
##########################################################

remove(list=ls(pa="^z"))
#Need to get Positive Tests going right back
a_begin <- as_date("2021-11-01")  #start of omicron

#link all positive tests to all sequenced by eave_linkno and then select one observation per person 
#picking the test closest to the sequence when there is a sequence and the first test when there is no sequence

z_df <- Positive_Tests %>% filter (date_ecoss_specimen >= a_begin) %>% 
  #dplyr::select(EAVE_LINKNO, date_ecoss_specimen) %>% 
  mutate(test_type="pcr")
z <- lft %>% filter (date_ecoss_specimen >= a_begin) %>% 
  dplyr::select(EAVE_LINKNO, date_ecoss_specimen) %>% 
  mutate(test_type="lft")
z_df <- bind_rows(z_df,z) %>% 
  arrange(EAVE_LINKNO, date_ecoss_specimen, desc(test_type)) %>% #pick pcr test if both done on the same day
  filter(!duplicated(paste(EAVE_LINKNO, date_ecoss_specimen)))
df_pos <- z_df
z <-  df_pos %>% group_by(date_ecoss_specimen, test_type) %>% 
  dplyr::summarise(N=n())
z %>% ggplot(aes(x=date_ecoss_specimen, y=N, colour=test_type)) +geom_point()


z_df <- wgs %>% dplyr::select(EAVE_LINKNO, Collection_Date, variant)
table(table(z_df$EAVE_LINKNO)) #check for duplicates

z <- df_pos %>% full_join(z_df, by="EAVE_LINKNO") %>% 
  mutate(date_ecoss_specimen = if_else(is.na(date_ecoss_specimen), Collection_Date, date_ecoss_specimen)) %>% #a few do not link
  mutate(days_spec_collection = as.numeric(Collection_Date - date_ecoss_specimen))
#z_id <- z %>% filter(duplicated(EAVE_LINKNO)) %>% pull(EAVE_LINKNO)
#z1 <- z %>% filter(EAVE_LINKNO %in% z_id)
#
z <- z %>% arrange(EAVE_LINKNO, abs(days_spec_collection), date_ecoss_specimen ) %>% filter(!duplicated(EAVE_LINKNO))
#z1 <- z %>% filter(!is.na(days_spec_collection) & days_spec_collection > 14)
#there are not many tests a long way from the date_ecoss_specimen and it is likely the variant is the same
#table(z$days_spec_collection, exclude=NULL)

df_pos <- z %>% dplyr::select(-days_spec_collection)

z <- smr01_covid %>% filter(ADMISSION_DATE >= a_begin)
table(table(z$EAVE_LINKNO)) #check for duplicates
z_df <- df_pos %>% left_join(z, by="EAVE_LINKNO")
table(table(z_df$EAVE_LINKNO)) #check for duplicates
z_df <- z_df %>% mutate(days = as.numeric(ADMISSION_DATE - date_ecoss_specimen))

z_id <- z_df %>% filter(duplicated(EAVE_LINKNO)) %>% pull(EAVE_LINKNO) %>% unique()
z0 <- z_df %>% filter(!(EAVE_LINKNO %in% z_id)) # unique records
z1 <- z_df %>% filter(EAVE_LINKNO %in% z_id) #duplicated records
#arrange so that the Admission closest to the specimen date is first and pick this record
z1 <- z1 %>% arrange(EAVE_LINKNO, abs(days) ) %>%  
  filter(!duplicated(EAVE_LINKNO)) 
z_df <- bind_rows(z0,z1)  

z_df <- z_df %>% mutate(covid_admission = if_else(!is.na(days) & days >= -13 & days <= 28, 1L,0L)) %>% 
  mutate(covid_admission_status = case_when(is.na(days) ~ "no_covid_admission",
                                             !is.na(days) & days < -13 ~"covid_admission_admitted_before_test",
                                             !is.na(days) & days >= -13 & days <= 28 ~ "covid_admission",
                                             !is.na(days) & days > 28 ~ "covid_admission_28_plus_post_test",
                                             TRUE ~ "other")) %>% 
  mutate(variant = if_else(is.na(variant), "not_sequenced", variant)) 
df_pos <- z_df %>% dplyr::select(-days)

a_end <- max(df_pos$date_ecoss_specimen)

#add in deaths for censoring
a_end_death <- max(all_deaths$NRS.Date.Death)
z <- df_pos %>% left_join(dplyr::select(all_deaths, EAVE_LINKNO, NRS.Date.Death, covid_death_cert ), by="EAVE_LINKNO") %>% 
  mutate(death = if_else(is.na(NRS.Date.Death), 0L, 1L)) %>% 
  mutate(days = as.numeric(NRS.Date.Death - date_ecoss_specimen))
z <- z %>%  filter(is.na(NRS.Date.Death) | NRS.Date.Death >= a_begin) #remove those dead before study started (despite spec date OK, likely chi mismatch)
z <- z %>%  filter(is.na(days) | days >= -7) %>% #remove those dead long time before test date , likely chi mismatch
  mutate(date_ecoss_specimen = if_else(!is.na(days) & days < 0, NRS.Date.Death, date_ecoss_specimen)) %>%  # change test to death date for those within a week
  filter(date_ecoss_specimen >= a_begin)
df_pos <- z %>% dplyr::select(-days)

df_pos <- df_pos %>% mutate(time_to_hosp = as.numeric(if_else(covid_admission==1, ADMISSION_DATE, a_end) - date_ecoss_specimen)) %>% 
  mutate(time_to_hosp = if_else(death==1 & !is.na(NRS.Date.Death), as.numeric(NRS.Date.Death - date_ecoss_specimen), time_to_hosp )) %>%
  mutate(time_to_hosp = case_when(time_to_hosp < 0 ~ 0,
                                  time_to_hosp > 28 ~ 28,
                                  TRUE ~ time_to_hosp)) %>% 
  mutate(test_type = if_else(is.na(test_type), "unknown", test_type))

z <- df_pos %>% group_by(date_ecoss_specimen) %>% dplyr::summarise(N=n(), R=sum(covid_admission))
z %>% ggplot(aes(x=date_ecoss_specimen))  + geom_point(aes(y=R), colour="red") + geom_vline(xintercept=a_end_smr01)

z <- df_pos %>% group_by(ADMISSION_DATE) %>% dplyr::summarise(N=n(), R=sum(covid_admission))
z %>% ggplot(aes(x=ADMISSION_DATE))  + geom_point(aes(y=R), colour="red") + geom_vline(xintercept=a_end_smr01)


##################
z <- cdw_full %>% filter(test_result == "POSITIVE") %>% filter(date_ecoss_specimen >= a_begin)
#Ecoss is the nhs labs
z <- z %>% mutate(lab = if_else(test_result_record_source == "ECOSS", "nhs","lh") )
#  dplyr::select(EAVE_LINKNO, ecossid, sex, age_year, specimen_date, lab)
z <- z%>%  
  dplyr::select(EAVE_LINKNO, test_id, subject_sex, age, date_ecoss_specimen, lab, test_result_record_source, flag_covid_symptomatic) %>% 
  dplyr::rename(age_year=age) 

z_df_tid <- df_pos %>% filter(!is.na(test_id))
z_df_no_tid <- df_pos %>% filter(is.na(test_id))

#for those with a test_id - lighthouse - match on this
z1 <- z %>% filter(!is.na(test_id)) %>% filter(!duplicated(test_id))
z_df_tid <- z_df_tid %>% left_join(dplyr::select(z1, -date_ecoss_specimen, -EAVE_LINKNO), by="test_id")
#for those with no test_id - match on EAVE_LINKNO and date_ecoss_specimen - only pcr tests will link
z_df_no_tid <- z_df_no_tid %>% left_join(dplyr::select(z, -test_id), by=c("EAVE_LINKNO","date_ecoss_specimen") )

print(nrow(z_df_tid) + nrow(z_df_no_tid) == nrow(df_pos))
table(z_df_tid$EAVE_LINKNO %in% z_df_no_tid$EAVE_LINKNO)
z_df <- bind_rows(z_df_tid,z_df_no_tid)

z_df <- z_df %>%  mutate(lab = if_else(is.na(lab) & test_type=="lft", "lft", lab))
df_pos <- z_df

#plot 
z1 <- df_pos %>% group_by(date_ecoss_specimen) %>% dplyr::summarise(N=n())
z2 <- df_pos %>% group_by(date_ecoss_specimen, lab) %>% dplyr::summarise(R=n())
z <- z2 %>% left_join(z1, by="date_ecoss_specimen") %>% mutate(P=R/N*100) 
z %>% filter(lab != "unknown") %>% ggplot(aes(x=date_ecoss_specimen, y=P, colour=lab)) +geom_point() +
  labs(x="Specimen Date", y="Percentage", colour="Test Type/Lab", title="Proportion of positive tests by day")
ggsave(paste0("./output/prop_lab_tests_day.png"), width=14, height=10, unit="cm")


#link in hospitalisations - use all as some will be in hospital at time of test

#keep all admissions post a_begin and those in hospital at a_begin
z_h <- all_hospitalisations %>% dplyr::select(-validchi) %>% 
  filter(!(!is.na(discharge_date) & discharge_date <= a_begin))
#summary(filter(z_h, admission_date < a_begin))  # all admitted before are discharged after a_begin or no discharge date

# Steven: I have taken the end date of hospitalisation records before we remove any records
# for people who have had multiple hospitalisations.
a_end_hosp <- max(c(max(z_h$admission_date, na.rm=T), max(z_h$discharge_date, na.rm=T)))


z <- df_pos %>%
  left_join(z_h, by="EAVE_LINKNO") %>%
  # Create variable for if they were in hospital at time of test
  # I think Chris's original code can miss if they were in hospital at time of test for people
  # with multiple admissions, since it ony keeps the first admission post test.
  mutate(In_Hosp_At_Test = case_when(date_ecoss_specimen > a_end_hosp ~ 'unknown',
                                     date_ecoss_specimen >= admission_date & date_ecoss_specimen <= discharge_date ~ 'yes',
                                     TRUE ~ 'no')) %>%
  group_by(EAVE_LINKNO) %>%
  mutate(In_Hosp_At_Test = case_when( any(In_Hosp_At_Test == 'yes') ~ 'yes',
                                      TRUE ~ In_Hosp_At_Test)) %>%
  ungroup()


z <- z %>% mutate(ad=admission_date, dd=discharge_date, em = emergency) %>% 
  mutate(ad = as.Date(ifelse(!is.na(dd) & dd < date_ecoss_specimen, NA, ad), origin="1970-01-01"),
         em = ifelse(!is.na(dd) & dd < date_ecoss_specimen, NA, em), 
         dd = as.Date(ifelse(!is.na(dd) & dd < date_ecoss_specimen, NA, dd), origin="1970-01-01") )
#multiple admissions - take the first admission in the study period for those with multiple
#this will be the first after specimen date or if before will have no discharge
z_id <- z %>% filter(duplicated(EAVE_LINKNO)) %>% pull(EAVE_LINKNO) %>% unique()
z_01 <- z %>% filter(!(EAVE_LINKNO %in% z_id)) # 0 or 1 admission
z_m <- z %>% filter((EAVE_LINKNO %in% z_id))  #multiple admissions
z_m <- z_m %>%  arrange(EAVE_LINKNO, ad) %>% 
  filter(!duplicated(EAVE_LINKNO))  # NA on ad go to the end, pick the first admission post specimen

z <- bind_rows(z_01, z_m) %>% 
  dplyr::select(-admission_date, -discharge_date, -emergency) %>% 
  dplyr::rename(admission_date=ad, discharge_date=dd, emergency=em)

z_df <- z

#now modify the admission dates for those admitted a long time ago but with no discharge and
#evidence of a recent test in the community 
#impute discharge dates for people admitted a long time ago with no discharge
#so different rules for those who are tested in teh community as opposed to thos tests in the nhs labs (ECOSS)
#z <- z_df %>% filter(!is.na(admission_date) & admission_date < date_ecoss_specimen - 7 & is.na(discharge_date) & test_result_record_source != "ECOSS")
#if there is an admission date more than 7 days before the specimen date and the test was a lighhouse
#or more than 30 days before and the test was in NHS
#- asssume the person was discharged and set admission date to NA so it won't count as a covid admission (in hosp at time of test)
z_df <- z_df %>% mutate(discharge_date_orig = discharge_date, admission_date_orig = admission_date)
z_df <- z_df %>%
  mutate(admission_date = case_when(
    !is.na(admission_date_orig) & admission_date_orig < date_ecoss_specimen - 7 & is.na(discharge_date_orig) & test_result_record_source != "ECOSS" ~ NA_Date_,
    !is.na(admission_date_orig) & admission_date_orig < date_ecoss_specimen - 30 & is.na(discharge_date_orig) & test_result_record_source == "ECOSS" ~ NA_Date_,
    TRUE ~ admission_date_orig))
#z_df <- z_df %>% 
#  mutate(discharge_date = case_when(
#    !is.na(admission_date_orig) & admission_date_orig < a_begin - 7 & is.na(discharge_date_orig) & test_result_record_source != "ECOSS" ~ admission_date_orig+6,
#    !is.na(admission_date_orig) & admission_date_orig < a_begin - 30 & is.na(discharge_date_orig) & test_result_record_source == "ECOSS" ~ admission_date_orig+14,
#    TRUE ~ discharge_date_orig))
#change the admission and discharge dates to NA if the imputed discharge date is before the speciment date
#use discharge_date as that has the imputed discharges - use admission_data_orig as that is unchanged
#z_df <- z_df %>%  mutate(admission_date = if_else(
#    !is.na(admission_date_orig) & !is.na(discharge_date) & discharge_date < date_ecoss_specimen,  NA_Date_, admission_date_orig)) %>% 
#  mutate(discharge_date = if_else(
#     !is.na(discharge_date) & discharge_date < date_ecoss_specimen,  NA_Date_, discharge_date))

#a_end_hosp <- max(c(max(z_df$admission_date, na.rm=T), max(z_df$discharge_date, na.rm=T)))

z_df <- z_df %>% mutate(Time.To.Hosp = case_when(
  !is.na(admission_date) & !is.na(discharge_date) & discharge_date <= date_ecoss_specimen ~ as.numeric(a_end_hosp - date_ecoss_specimen),
  !is.na(admission_date) & is.na(discharge_date)  ~  as.numeric(admission_date - date_ecoss_specimen),
  !is.na(admission_date) & !is.na(discharge_date) & discharge_date > date_ecoss_specimen ~ as.numeric(admission_date - date_ecoss_specimen),
  is.na(admission_date) ~ as.numeric(a_end_hosp - date_ecoss_specimen),
  TRUE ~ NA_real_) ) %>% 
# change time fro those who die during follow up  
  mutate(Time.To.Hosp = if_else(!is.na(NRS.Date.Death) & is.na(admission_date), as.numeric(NRS.Date.Death - date_ecoss_specimen),Time.To.Hosp)) %>% 
 # mutate(In_Hosp_At_Test = ifelse(Time.To.Hosp <= -3  ~ "yes", In_Hosp_At_Test) %>% 
  mutate(hosp_pos_test = if_else(!is.na(admission_date) & Time.To.Hosp <= 14,  1L, 0L ),
         Time.To.Hosp = case_when(Time.To.Hosp < 0 ~ 0,
                                  Time.To.Hosp >= 15 ~ 15,
                                  TRUE ~ Time.To.Hosp))
z_df <- z_df %>% mutate(hosp_pos_test_emerg = if_else(!is.na(emergency) & emergency ,hosp_pos_test, 0L ) )
z_df <- z_df %>% mutate(days = as.numeric(date_ecoss_specimen- min(date_ecoss_specimen)) )

#add in the vaccinations and risk groups
z_df <- z_df %>%  
  left_join(dplyr::select(Vaccinations, -age, -patient_sex), by="EAVE_LINKNO" )   #use age, sex from testing data
z_df <- z_df %>%  mutate(vs1 = case_when(is.na(date_vacc_1) | date_vacc_1 > date_ecoss_specimen ~ "uv",
                                         date_vacc_1 <= date_ecoss_specimen &  date_vacc_1 > date_ecoss_specimen - 28 ~ "v1_0:3",
                                         TRUE ~ "v1_4+"),
                         vs2 = case_when(is.na(date_vacc_2) | date_vacc_2 > date_ecoss_specimen ~ "uv",
                                         date_vacc_2 <= date_ecoss_specimen &  date_vacc_2 > date_ecoss_specimen - 174 ~ "v2_0:24",
                                         TRUE ~ "v2_25+"),
                         vs3 = case_when(is.na(date_vacc_3) | date_vacc_3 > date_ecoss_specimen ~ "uv",
                                         date_vacc_3 <= date_ecoss_specimen &  date_vacc_3 > date_ecoss_specimen - 83 ~ "v3_0:11",
                                         date_vacc_3 <= date_ecoss_specimen-84 &  date_vacc_3 > date_ecoss_specimen - 111 ~ "v3_12:15",
                                         date_vacc_3 <= date_ecoss_specimen-112 &  date_vacc_3 > date_ecoss_specimen - 139 ~ "v3_16:19",
                                         date_vacc_3 <= date_ecoss_specimen-140 &  date_vacc_3 > date_ecoss_specimen - 181 ~ "v3_20:25",
                                         TRUE ~ "v3_26+"),
                         vs4 = case_when(is.na(date_vacc_4) | date_vacc_4 > date_ecoss_specimen ~ "uv",
                                         date_vacc_4 <= date_ecoss_specimen &  date_vacc_4 > date_ecoss_specimen - 13 ~ "v4_0:1",
                                         date_vacc_4 <= date_ecoss_specimen-14 &  date_vacc_4 > date_ecoss_specimen - 55~ "v4_2:7",
                                         TRUE ~ "v4_8+") ) %>% 
  mutate(vs = case_when(vs4 !="uv" ~ vs4,
                        vs3 !="uv" ~ vs3,
                        vs2 !="uv" ~ vs2,
                        TRUE ~ vs1))
z_levs <- c("uv","v1_0:3","v1_4+","v2_0:24","v2_25+",
            "v3_0:11","v3_12:15","v3_16:19", "v3_20:25", "v3_26+", "v4_0:1","v4_2:7","v4_8+")
z_df <- z_df %>% mutate(vs=factor(vs, levels=z_levs))
#merge_levels
z_df <- z_df %>% mutate(vs=fct_recode(vs, "v1" = "v1_0:3", "v1" = "v1_4+"))

z_df <- z_df %>% 
  mutate(vacc_type = if_else(is.na(vacc_type) , "uv", as.character(vacc_type)) ) %>% 
  filter(vacc_type != "UNK") %>% 
  mutate(vacc_type = factor(vacc_type, levels=c("uv","AZ","Mo","PB") ) ) %>% 
  mutate(vt = fct_cross(vs,vacc_type, sep="_")) %>% 
  mutate(vt = fct_recode(vt, "uv" ="uv_uv", "uv" = "uv_AZ", "uv" = "uv_Mo","uv" = "uv_PB"))

#add in demographics
z_df <- z_df %>%  
  left_join(dplyr::select(EAVE_cohort, -NRS.Date.Death), by="EAVE_LINKNO" )
z_df <- z_df %>% mutate(in_eave = if_else(is.na(eave_weight), 0L,1L))
#update age for those missing age/sex (lft test)
z_df <- z_df %>% mutate(age_year = if_else(is.na(age_year), ageYear,age_year),
                        subject_sex = ifelse(is.na(subject_sex), Sex, subject_sex))

z_df <- z_df %>% mutate(age_gp = cut(age_year, breaks = c(-1, 11, 19, 39, 59, 74, 130),
                                     labels=c("0-11", "12-19","20-39","40-59","60-74","75+")))
z_df <- z_df %>% dplyr::rename(sex=subject_sex)
z_df <- z_df %>% filter(!is.na(age_year))
z_df <- z_df %>% dplyr::select(-Sex, - ageYear)
df_pos <- z_df
df_pos <- df_pos %>%  mutate(sex = case_when(sex=="F" ~ "Female",
                                             sex=="M" ~ "Male",
                                             TRUE ~ sex))
#add in risk groups as at Dec 2020
z_df <- df_pos %>% left_join(rg, by="EAVE_LINKNO")
z_df <- z_df %>% mutate(n_risk_gps = fct_explicit_na(n_risk_gps, na_level = "Unknown"))
df_pos <- z_df

#add in icu/death  - icu dates depend upon the endpoints linkage

#covid death derived from all deaths - calculate Time.To.Death
z <- df_pos %>% 
  mutate(Time.To.Death = if_else(is.na(NRS.Date.Death), as.numeric(a_end_death - date_ecoss_specimen),
                                       as.numeric(NRS.Date.Death - date_ecoss_specimen))) %>% 
  mutate(Time.To.Death = if_else(Time.To.Death < 0, 0, Time.To.Death)) %>% 
  mutate(covid_death_cert = if_else(is.na(covid_death_cert), 0, covid_death_cert)) %>% 
  mutate(covid_death = if_else(is.na(covid_death_cert), 0, covid_death_cert)) %>% 
  mutate(covid_death = if_else(!is.na(NRS.Date.Death) & (NRS.Date.Death - date_ecoss_specimen <= 28) , 1, covid_death)) %>% 
  mutate(exclude_death_analysis = if_else(date_ecoss_specimen > a_end_death, 1,0)) # exclude those testing after latest death date
z_df <- z

#covid icu derived from EAVE_Endpoints
a_end_icu <- max(covid_icu$covid_icu_date)
z_icu <- covid_icu %>% filter(covid_icu_date > a_begin) %>% 
  arrange(EAVE_LINKNO, covid_icu_date) %>% 
  filter(!duplicated(EAVE_LINKNO))
z <- z_df %>% left_join(dplyr::select(z_icu, -SpecimenDate), by="EAVE_LINKNO") %>% 
  mutate(Time.To.ICU = if_else(is.na(covid_icu_date), as.numeric(a_end_icu - date_ecoss_specimen),
                                 as.numeric(covid_icu_date - date_ecoss_specimen))) %>% 
  mutate(Time.To.ICU = if_else(!is.na(NRS.Date.Death) & is.na(covid_icu_date) , as.numeric(NRS.Date.Death - date_ecoss_specimen), Time.To.ICU)) %>% 
  mutate(Time.To.ICU = if_else(Time.To.ICU < 0, 0, Time.To.ICU)) %>% 
  mutate(covid_icu = if_else(is.na(covid_icu_date), 0, 1)) %>% 
  mutate(exclude_icu_analysis = if_else(date_ecoss_specimen > a_end_icu, 1,0)) # exclude those testing after latest death date
z_df <- z

df_pos <- z_df

#z_df <- z_df %>% mutate(prev_pos = case_when(is.na(days_from_previous_pos) ~ "not_prev_pos",
#                                             days_from_previous_pos >= 1 & days_from_previous_pos <= 28 ~ "pos_1:28",
#                                             days_from_previous_pos >= 29 & days_from_previous_pos <= 90 ~ "pos_29:90",
#                                             TRUE ~ "pos_91+"))

z.df <- df_pos
z.df <- mutate(z.df, Time.To.Hosp = if_else(Time.To.Hosp==0,0.1,Time.To.Hosp))  #add 0.1 to 0 days so survival models work
z.df <- mutate(z.df, time_to_hosp = if_else(time_to_hosp==0,0.1,time_to_hosp))  #add 0.1 to 0 days so survival models work
z.df <- mutate(z.df, Time.To.Death = if_else(Time.To.Death==0,0.1,Time.To.Death))
z.df <- mutate(z.df, Time.To.ICU = if_else(Time.To.ICU==0,0.1,Time.To.ICU))



z.df <- z.df %>% mutate(bmi.gp = cut(bmi_impute, breaks=c(-1, 20,25,30,35,40,51), labels=FALSE))
z.df <- z.df %>% mutate(bmi_gp = case_when(age_year <= 17 ~ "NA_too_young",
                                           age_year >= 18 & !is.na(bmi.gp) ~ as.character(bmi.gp),
                                           TRUE ~ "Unknown")) %>% 
  mutate(bmi_gp = factor(bmi_gp, levels=c("1","2","3","4","5","6","NA_too_young","Unknown"),
                         labels=c("<20","20-24","25-29","30-34","35-39","40+","NA_too_young","Unknown")))

df_pos <- z.df
df_pos$eave_weight[is.na(df_pos$eave_weight)] <- 0
df_pos$wght_all <- 1

df_pos <- df_pos %>% mutate(shielding = if_else(EAVE_LINKNO %in% shielding$EAVE_LINKNO, 1, 0))
df_pos <- df_pos %>% mutate(immuno = if_else(EAVE_LINKNO %in% immuno$EAVE_LINKNO, 1, 0))


rgs <- colnames(df_pos)[startsWith(colnames(df_pos), "Q")]
df_pos <- df_pos %>% mutate(across(all_of(rgs), ~ as.factor(.)))

saveRDS(df_pos, "output/temp/All_positive.rds") 


