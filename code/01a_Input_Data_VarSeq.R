##########################################################
# Name of file: 01_Input_Data_Omicron.R
# Data release (if applicable):
# Original author(s): Chris Robertson chris.robertson@phs.scot
# Original date: 24 June 2022
# Latest update author (if not using version control) - Chris Robertson chris.robertson@phs.scot
# Latest update date (if not using version control) - 
# Latest update description (if not using version control)
# Type of script: Descriptive stats
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 3.6.1
# Description of content: reads in the cohort and merges in the risk groups 
#                         
# Approximate run time: Unknown
##########################################################

# 01 Setup ####
#Libraries
library(plyr)
library(tidyverse)
library(survival)
library(lubridate)
#Load data

Location <- "/conf/"  # Server
#Location <- "//isdsf00d03/"  # Desktop
project_path <- paste0(Location,"EAVE/GPanalysis/progs/CR/Variants_Sequencing")
project_path_vaccine <- paste0(Location,"EAVE/GPanalysis/progs/CR/Vaccine")

#just use the demographics from here
EAVE_cohort <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/Cohort_Demog_Endpoints_Times2021-07-28.rds"))
EAVE_cohort <- filter(EAVE_cohort, !duplicated(EAVE_LINKNO)) %>% dplyr::select(EAVE_LINKNO:ur6_2016_name) %>% 
  mutate(ageYear=ageYear+2)  # make age as at March 2022

a_begin <- as.Date("2021-11-01")  #start date for the analysis
#read in all deaths and then omit any who have died prior to a_begin
all_deaths  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/all_deaths.rds"))
summary(all_deaths)
EAVE_cohort <- EAVE_cohort %>% left_join(dplyr::select(all_deaths, EAVE_LINKNO, NRS.Date.Death), by="EAVE_LINKNO")
EAVE_cohort <- EAVE_cohort %>% filter(is.na(NRS.Date.Death) | !is.na(NRS.Date.Death) & NRS.Date.Death > a_begin)
EAVE_cohort <- filter(EAVE_cohort, !duplicated(EAVE_LINKNO)) #duplicate deaths

EAVE_Weights <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/CR_Cohort_Weights.rds"))
EAVE_cohort  <- EAVE_cohort %>% left_join(EAVE_Weights, by="EAVE_LINKNO")
EAVE_cohort$eave_weight[is.na(EAVE_cohort$eave_weight)] <- mean(EAVE_cohort$eave_weight, na.rm=T)

#read in the Previous Tests data
cdw_full  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/CDW_full.rds"))
cdw_full <- cdw_full %>% mutate(date_ecoss_specimen = as_date(date_ecoss_specimen))
cdw_full <- cdw_full %>%  
  arrange(EAVE_LINKNO, date_ecoss_specimen, desc(test_result)) %>%
  filter(!duplicated(paste(EAVE_LINKNO, date_ecoss_specimen)))  #get one test per person per day - preferentially positive test 
cdw_full <- filter(cdw_full, date_ecoss_specimen <= Sys.Date())  
summary(cdw_full)

Positive_Tests <- cdw_full %>% filter(test_result=="POSITIVE") %>% 
  dplyr::select(EAVE_LINKNO, test_id, date_ecoss_specimen) 
summary(Positive_Tests)
length(unique(paste(Positive_Tests$EAVE_LINKNO, Positive_Tests$date_ecoss_specimen))) # no duplicates on the same date

#read in the Previous Tests data
lft  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/lft_positives.rds"))
lft <- lft %>% 
  arrange(EAVE_LINKNO, date_ecoss_specimen) %>% 
  filter(!duplicated(paste(EAVE_LINKNO, date_ecoss_specimen))) #keep one positive test per day

#read in the GP vaccination data
#source("../00_Read_GP_Vaccinations.R")
source("../00_Read_DV_Vaccinations_Dose4.R")

#get covid death certificate deaths
z <- all_deaths %>%  
  dplyr::mutate(across(UNDERLYING_CAUSE_OF_DEATH:CAUSE_OF_DEATH_CODE_9, ~if_else(. %in% c("U071","U072"), 1,0)))
z_names <- names(z)[grepl("CAUSE", names(z))]
z <- z %>% mutate(rowsum = apply(z[,z_names],1,sum)) %>% 
  mutate(covid_death_cert = if_else(rowsum>=1,1,0)) 
all_deaths <- z %>% dplyr::select(-rowsum) 
summary(dplyr::select(all_deaths, NRS.Date.Death, NRS.Reg.Date, covid_death_cert, UNDERLYING_CAUSE_OF_DEATH))
#plot to check
z <- all_deaths %>% filter(NRS.Date.Death >= a_begin ) %>% group_by(NRS.Date.Death) %>% 
  dplyr::summarise(N=n())
z %>% ggplot((aes(x=NRS.Date.Death, y=N))) + geom_point() + labs(title="All Deaths")


all_hospitalisations  <- readRDS(paste0(Location,"EAVE/GPanalysis/data/automated_any_hospitalisation_post_01022020.rds")) #RAPID
summary(all_hospitalisations)
#plot to check
z <- all_hospitalisations %>% filter(admission_date >= a_begin ) %>% group_by(admission_date) %>% 
  dplyr::summarise(N=n())
z %>% ggplot((aes(x=admission_date, y=N))) + geom_point() + labs(title="All Admissions to Hospital - RAPID")

#cohort + risk groups
rg <- readRDS(paste0(project_path_vaccine,"/output/temp/Qcovid_all.rds"))
rg <- rg %>% dplyr::select(-(Sex:ur6_2016_name), -Q_BMI)
rg <- filter(rg, !duplicated(EAVE_LINKNO))

z <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/CR_Cohort_RG_EAVE_BP_Smoke.rds"))
z <- filter(z, !duplicated(EAVE_LINKNO))
z <- z %>% dplyr::select(EAVE_LINKNO, EAVE_Smoking_Status_Worst, EAVE_BP) %>% 
  dplyr::rename(EAVE_Smoke = EAVE_Smoking_Status_Worst)

rg <- rg %>% left_join(z, by="EAVE_LINKNO")
rg <- rg %>% mutate(EAVE_Smoke = if_else(!is.na(EAVE_Smoke), as.character(EAVE_Smoke), "Unknown"),
                    EAVE_BP = if_else(!is.na(EAVE_BP), as.character(EAVE_BP), "No Investigation"))

#update weights
#those with a pis records over the last 12 months before March 2020
bnf <- readRDS(paste0(Location,"EAVE/GPanalysis/data/BNF_paragraphs.rds"))
sicsag <- readRDS(paste0(Location,"EAVE/GPanalysis/data/SICSAG_episode_level_.rds"))
smr01 <- readRDS(paste0(Location,"EAVE/GPanalysis/data/SMR01_allstays.rds"))
pis <- readRDS(paste0(Location,"EAVE/GPanalysis/data/PIS_2019_2021_EAVELink.rds"))
#z <- readRDS(paste0(Location,"EAVE/GPanalysis/data/patient_file.rds"))
z_ids <- c(Vaccinations$EAVE_LINKNO, all_deaths$EAVE_LINKNO,  bnf$EAVE_LINKNO,
           cdw_full$EAVE_LINKNO, all_hospitalisations$EAVE_LINKNO, sicsag$EAVE_LINKNO, 
           smr01$EAVE_LINKNO, lft$EAVE_LINKNO, pis$EAVE_LINKNO) %>% unique()
#summary(filter(EAVE_cohort, !(EAVE_LINKNO %in% z_ids))$eave_weight)
z_N <- round(sum(EAVE_cohort$eave_weight) )
z_k <- sum(EAVE_cohort$EAVE_LINKNO %in% z_ids)
z_m <- round(sum(filter(EAVE_cohort, (EAVE_LINKNO %in% z_ids))$eave_weight))
z <- EAVE_cohort %>% mutate(ew = if_else(EAVE_LINKNO %in% z_ids, 1, eave_weight*(z_N - z_k)/(z_N - z_m)) )
EAVE_cohort <- z %>% dplyr::select(-eave_weight) %>% dplyr::rename(eave_weight=ew)

#remove data sets only used for identification of patients who exist
rm(bnf, sicsag, pis)

z <- read_csv(paste0(Location,"/EAVE/GPanalysis/data/restored/map_files/Datazone2011Lookup.csv")) %>% 
  dplyr::select(DataZone, InterZone, Council, HB)
EAVE_cohort <- EAVE_cohort %>% left_join(z, by="DataZone") %>% 
  mutate(HB = if_else(is.na(HB),"Unknown", HB),
         InterZone = if_else(is.na(InterZone),"Unknown", InterZone),
         Council = if_else(is.na(Council),"Unknown", Council))


wgs <- readRDS(paste0(Location,"EAVE/GPanalysis/data/WGS_latest.rds")) %>% 
  mutate_at(c("Collection_Date","Sequencing_Date","Alignment_Date"), ~ as.Date(. , format="%d/%m/%Y")) %>% 
  filter(Collection_Date >= a_begin)

#z_id <- filter(wgs, duplicated(EAVE_LINKNO)) %>% pull(EAVE_LINKNO)
#z <- filter(wgs, EAVE_LINKNO %in% z_id) %>% arrange(EAVE_LINKNO)
#duplicates all have the same lineage
wgs <- wgs %>% arrange(EAVE_LINKNO, Collection_Date) # %>% filter(!duplicated(EAVE_LINKNO))
a_end_wgs <- max(wgs$Collection_Date) - 2  #table(wgs$Collection_Date)  #check each time
z <- wgs %>% mutate(variant = case_when(is.na(VariantofInterest) ~ "not_sequenced",
                       VariantofInterest=="VOC-22APR-04" ~ "omicron_ba5",
                       VariantofInterest=="VOC-22APR-03" ~ "omicron_ba4",
                       VariantofInterest=="V-22SEP-01" ~ "omicron_ba4",  #4.6 grouped with 4
                       VariantofInterest=="V-21APR-02" ~ "delta",
                       VariantofInterest=="V-21OCT-01" ~ "delta_ay42",
                       VariantofInterest=="V-22APR-02" ~ "omicron_xe",
                       VariantofInterest=="VOC-21NOV-01" ~ "omicron_ba1",
                       VariantofInterest=="VOC-22JAN-01" ~ "omicron_ba2",
                       VariantofInterest=="V-22JUL-01" ~ "omicron_ba2", #2.75 grouped with 2
                       VariantofInterest=="Omicron_Unassigned" ~ "omicron",
                       VariantofInterest=="SIM-BA3" ~ "omicron",
                       VariantofInterest=="V-22OCT-01" ~ "omicron_bq1",
                       TRUE ~ "other") ) 
wgs <- z

shielding <- readRDS(paste0(Location,"EAVE/GPanalysis/data/Shielding_list.rds"))
immuno <- readRDS(paste0(Location,"EAVE/GPanalysis/data/cleaned_data/Imm_supp_cohort_Nov2021.rds"))

#mostly these people are dead and so come out of the analysis at the time of death
#not_resident <- readRDS(paste0(Location,"EAVE/GPanalysis/data/not_resident_list.rds"))


#get covid deaths in Omicron Period
#all with covid on death certificate
covid_death <- all_deaths %>% filter(covid_death_cert == 1) %>% filter(NRS.Date.Death >= a_begin) %>%
  dplyr::rename(covid_uc_death = UNDERLYING_CAUSE_OF_DEATH) %>% 
  dplyr::select(EAVE_LINKNO, NRS.Date.Death, covid_death_cert, covid_uc_death)

z <- covid_death %>% group_by(NRS.Date.Death) %>% 
  dplyr::summarise(N=n(), N_dc = sum(covid_death_cert), N_uc = sum(covid_uc_death) )
z %>% ggplot(aes(x=NRS.Date.Death, y=N)) + geom_point() + geom_point(aes(y=N_dc), colour="blue") +
  geom_point(aes(y=N_uc), colour="red") + labs(title="covid deaths")
z %>% ggplot(aes(x=NRS.Date.Death, y=N)) + geom_smooth(colour="black",se=F, span=0.3) + 
  geom_smooth(aes(y=N_dc), colour="blue",se=F, span=0.3) +
  geom_smooth(aes(y=N_uc), colour="red",se=F, span=0.3) + labs(title="covid deaths")
summary(covid_death)

#endpoints
EAVE_endpoints <- readRDS(paste0(Location,"EAVE/GPanalysis/outputs/temp/severe_endpoints2022-06-23.rds"))
EAVE_endpoints <- EAVE_endpoints %>% filter(!duplicated(paste(EAVE_LINKNO, KEYDATE))) %>% 
  arrange(EAVE_LINKNO, KEYDATE)
#z_id <- EAVE_endpoints %>% filter(duplicated(paste(EAVE_LINKNO, KEYDATE))) %>% pull(EAVE_LINKNO)
#z <- EAVE_endpoints %>% filter(EAVE_LINKNO %in% z_id) %>% arrange(EAVE_LINKNO, KEYDATE)

table(EAVE_endpoints$dead28, is.na(EAVE_endpoints$NRS.Date.Death), exclude=NULL) 
table(EAVE_endpoints$covid_cod, is.na(EAVE_endpoints$NRS.Date.Death), exclude=NULL)
table(EAVE_endpoints$covid_ucod, is.na(EAVE_endpoints$NRS.Date.Death), exclude=NULL)
table(EAVE_endpoints$covid_cod, EAVE_endpoints$covid_ucod, is.na(EAVE_endpoints$NRS.Date.Death), exclude=NULL)
table(EAVE_endpoints$covid_cod, EAVE_endpoints$dead28, is.na(EAVE_endpoints$NRS.Date.Death), exclude=NULL)
#for a covid death use death28==1 | covid_cod==1
#dead28 death within 28 days of a pos test
#9 for covid_cod means no NRS info on cause of death
#for event date use NRS.Date.Death

table(EAVE_endpoints$icu28, is.na(EAVE_endpoints$ICU_admit_date), exclude=NULL)
table(EAVE_endpoints$covid_mcoa_icu, is.na(EAVE_endpoints$ICU_admit_date), exclude=NULL)
table(EAVE_endpoints$icu28,EAVE_endpoints$covid_mcoa_icu, is.na(EAVE_endpoints$ICU_admit_date), exclude=NULL)
#for covid_icu use icu28==1 | covid_mcao_icu==1
#for event date use ICU_admit_date - present for all ICU admissions
#for the combined covid_death_icu use death_covid==1 | covid_icu==1 
#for event date use ICU_admit_date for covid_icu==1 or NRS.Date.Death for the others

table(EAVE_endpoints$hosp28, EAVE_endpoints$hosp14, is.na(EAVE_endpoints$hosp_admit_date), exclude=NULL)
table(EAVE_endpoints$hosp14, is.na(EAVE_endpoints$hosp_admit_date), exclude=NULL)
table(EAVE_endpoints$covid_mcoa_hosp, is.na(EAVE_endpoints$hosp_admit_date), exclude=NULL)
table(EAVE_endpoints$hosp14,EAVE_endpoints$covid_mcoa_hosp, is.na(EAVE_endpoints$hosp_admit_date), exclude=NULL)
#covid hosp is hosp14==1 | covid_mcoa_hosp==1
#covid_mcoa_hosp = 9 for RAPID admissions
#for event date use hosp_admit_date - present for all admissions
#emergency for emergency admissions

#adjust inconsistencies in the endpoints and times - all hosp have an admission date
#single endpoints
EAVE_endpoints <- EAVE_endpoints %>% 
  mutate(covid_death = case_when(!is.na(dead28) & dead28==1 ~ 1L,
                                 !is.na(covid_cod) & covid_cod==1 ~ 1L,
                                 TRUE ~ 0L) ) %>%
  mutate(covid_icu = case_when(!is.na(icu28) & icu28==1 ~ 1L,
                               !is.na(covid_mcoa_icu) & covid_mcoa_icu==1 ~ 1L,
                               TRUE ~ 0L) ) %>%
  mutate(covid_hosp = case_when(!is.na(hosp14) & hosp14==1 ~ 1L,
                                !is.na(covid_mcoa_hosp) & covid_mcoa_hosp==1 ~ 1L,
                                TRUE ~ 0L) ) %>% 
  mutate(covid_hosp = if_else(covid_hosp==0 & covid_icu==1, 1L, covid_hosp))
#dates for single endpoints
EAVE_endpoints <- EAVE_endpoints %>% 
  mutate(covid_death_date = if_else(covid_death==1, NRS.Date.Death, NA_Date_),
         covid_icu_date = if_else(covid_icu==1, ICU_admit_date, NA_Date_),
         covid_hosp_date = if_else(covid_hosp==1 & !is.na(hosp_admit_date), hosp_admit_date, NA_Date_ ))


table(EAVE_endpoints$covid_hosp, EAVE_endpoints$covid_icu, exclude=NULL)
table(EAVE_endpoints$covid_hosp, is.na(EAVE_endpoints$covid_hosp_date), exclude=NULL)
table(EAVE_endpoints$covid_icu, is.na(EAVE_endpoints$covid_icu_date), exclude=NULL)
table(EAVE_endpoints$covid_death, is.na(EAVE_endpoints$covid_death_date), exclude=NULL)

#get individual endpoints separately
#covid_death <- EAVE_endpoints %>% filter(covid_death==1) %>% 
#  dplyr::select(EAVE_LINKNO, SpecimenDate, covid_death_date, covid_cod)
#table(table(covid_death$EAVE_LINKNO))  #all unique
#table(covid_death$covid_cod, exclude=NULL)

covid_icu <- EAVE_endpoints %>% filter(covid_icu==1) %>% 
  dplyr::select(EAVE_LINKNO, SpecimenDate, covid_icu_date, covid_mcoa_icu) %>% 
  arrange(paste(EAVE_LINKNO, covid_icu_date)) #%>% 
#filter(!duplicated(paste(EAVE_LINKNO, covid_icu_date)))
table(table(covid_icu$EAVE_LINKNO))  #a few have 2+ ICU admissions

smr01_covid <- smr01 %>% filter(covid_main_diag_admit==1 |covid_main_other_ep==1 ) %>% 
  dplyr::select(EAVE_LINKNO, covid_main_diag_admit, covid_main_other_ep, ADMISSION_DATE) %>% 
  #filter(ADMISSION_DATE >= a_begin) %>% 
  mutate(ADMISSION_DATE = as_date(ADMISSION_DATE))
z <- smr01_covid %>% group_by(ADMISSION_DATE) %>% 
  dplyr::summarise(N=n(), N_other = sum(covid_main_other_ep), N_uc = sum(covid_main_diag_admit) )
z %>% ggplot(aes(x=ADMISSION_DATE, y=N)) + geom_point() + geom_point(aes(y=N_other), colour="blue") +
  geom_point(aes(y=N_uc), colour="red") + labs(title="smr01_covid_hosp")


#update lft to omit positive lft's followed by a negative PCR within 2 days

z <- cdw_full %>% filter(test_result=="NEGATIVE") %>% 
  filter(date_ecoss_specimen >= min(lft$date_ecoss_specimen)) %>% 
  dplyr::select(EAVE_LINKNO, date_ecoss_specimen) %>% 
  arrange(EAVE_LINKNO, date_ecoss_specimen) %>% 
  filter(!duplicated(paste(EAVE_LINKNO, date_ecoss_specimen))) #keep one negative test per day

z_lft <- lft %>% dplyr::mutate(id=row_number()) 
z1 <- z_lft %>% 
  left_join(z, by="EAVE_LINKNO", suffix=c("_lft","_neg")) %>% 
  filter((date_ecoss_specimen_neg >= date_ecoss_specimen_lft) &  (date_ecoss_specimen_neg <= date_ecoss_specimen_lft))
z1 <- z1 %>% pull(id) %>% unique()
z_lft <- z_lft %>% mutate(pcr_neg = ifelse(id %in% z1, 1L,0L)) %>%
  filter(pcr_neg==0) %>%  #drop those followed by a negative PCR
  dplyr::select(-id, -pcr_neg)
lft <- z_lft  # these are true lft positives

#get a date for completeness of SMR01

z <- smr01 %>% filter(ADMISSION_DATE >= a_begin) %>% group_by(ADMISSION_DATE) %>% 
  dplyr::summarise(N=n())
z %>% ggplot(aes(x=ADMISSION_DATE, y=N)) + geom_point()
#check but this is a rough rule
a_end_smr01 <- z %>% filter(N >= 1000) %>% filter(ADMISSION_DATE==max(ADMISSION_DATE)) %>% pull(ADMISSION_DATE) %>% as.Date()

remove(list=ls(pa="^z"))
