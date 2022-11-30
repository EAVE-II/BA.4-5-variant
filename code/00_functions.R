##########################################################
## Title: 1st dose COVID-19 vaccine waning
## Code author(s): Rachel Mulholland <rachel.mulholland@ed.ac.uk> 
##                 Chris Robertson <chrisrobertson@nhs.net>
## Description: 00_functions - Unique functions for analysis
##########################################################

#### Libraries ####
library("spatstat")

# Table defaults ---------------------------------------------------------------------------
# This makes table resize or continue over multiple pages in all output types
# PDF powered by kableExtra, Word by flextable
mytable = function(x, caption = "", row.names = FALSE, longtable = TRUE,
                   latex_options = c("hold_position"), font_size = 7.0, ...){
  
  # if not latex or html then else is Word
  if (is_latex_output()) {
    knitr::kable(x, row.names = row.names, align = c("l", "l", "r", "r", "r", "r", "r", "r", "r", "r", "r"),
                 booktabs = TRUE, caption = caption, #longtable = longtable,
                 linesep = "", ...) %>%
      kableExtra::kable_styling(font_size = font_size,
                                latex_options = latex_options)
  } else if(is_html_output()) {
    knitr::kable(x, row.names = row.names, align = c("l", "l", "r", "r", "r", "r", "r", "r", "r", "r", "r"),
                 booktabs = TRUE, caption = caption, longtable = longtable,
                 linesep = "", ...) %>%
      kableExtra::kable_styling(latex_options = c("scale_down", "hold_position"))
  } else {
    flextable::flextable(x) %>%
      flextable::autofit() %>%
      flextable::width(j = 1, width = 1.5) %>%
      flextable::height(i = 1, height = 0.5, part = "header")
  }
}


#### Summary table using weights ####
# Creates a table of cohort summaries using weights
# For categorical variables, the sum of the weights and the % is calculated
# For numerical variables, the weighed mean and weighted sd are calculated, as well as
# the weighted median and weighted IQR

# Input:
# - data = the dataset (must have weights by default named eave_weight, but can be changed by assigning wght-)
# - dependent = a character of the dependent variables name
# - explanatory = a string of characters of the explanatory variables

# Output:
# A table with each explanatory variable as a row (multiple rows for each category if categorical)
# with two columns of the weighted summaries for the levels in the dependent variable

summary_factorlist_wt <- function(data, dependent, explanatory, wght="eave_weight"){
  # Create list to put in summaries into each element
  summary_tbl_list <- list()
  data$wght <- data[,wght]  # create a variable called wght
  for(i in 1:length(explanatory)){
    
    # Extract variable
    n <- data %>%
      pull(!!sym(explanatory[i]))
    
    # If numeric then make weighted mean
    if(is.numeric(n)) {
      z_mean <- data %>%
        group_by(!!sym(dependent)) %>%
        summarise(mean = round(weighted.mean(!!sym(explanatory[i]), w = wght, na.rm = TRUE),1),
                  sd = round(sqrt(spatstat.geom::weighted.var(!!sym(explanatory[i]), w = wght)),1)) %>%
        mutate(mean.sd = paste0(mean, " (",sd,")")) %>%
        select(-mean, -sd) %>%
        mutate("characteristic" = explanatory[i]) %>%
        pivot_wider(names_from = !!sym(dependent), values_from = mean.sd) %>%
        relocate(characteristic) %>%
        mutate(levels = "mean.sd")
      
      
      z_median <- data %>%
        group_by(!!sym(dependent)) %>%
        summarise(median = spatstat.geom::weighted.median(!!sym(explanatory[i]), w = wght),
                  q1 = spatstat.geom::weighted.quantile(!!sym(explanatory[i]), w = wght, probs = 0.25),
                  q3 = spatstat.geom::weighted.quantile(!!sym(explanatory[i]), w = wght, probs = 0.75)) %>%
        mutate("characteristic" = explanatory[i]) %>%
        mutate(iqr = q3 -q1) %>%
        mutate(median.iqr = paste0(median, " (",iqr,")")) %>%
        select(-q1, -q3, -median, -iqr) %>%
        pivot_wider(names_from = !!sym(dependent), values_from = median.iqr) %>%
        relocate(characteristic) %>%
        mutate(levels = "median.iqr")
      
      # Combine!!
      summary_tbl_list[[i]] <- full_join(z_mean, z_median)
      
      
      # Else get sum of weights of each level
    } else if (length(unique(data %>% pull(!!sym(dependent) ) ) ) ==1) {
      
      # This is for when there is only one level in the dependent variable
      summary_tbl_list[[i]] <- data %>%
        group_by(!!sym(explanatory[i])) %>%
        summarise(n = sum(wght)) %>%
        ungroup() %>%
        mutate(perc = sprintf("%.1f",round(n/sum(n)*100,1))) %>%
        mutate_if(is.numeric, ~formatC(round(.,0), format = "f", big.mark = ",", drop0trailing = TRUE)) %>%
        mutate(n_perc := paste0(n, " (", perc,"%)")) %>%
        select(-n, -perc) %>%
        rename("levels"=explanatory[i], !!dependent := n_perc) %>%
        mutate("characteristic" = explanatory[i]) %>%
        relocate(characteristic)
      
    } else {

      summary_tbl_list[[i]] <- data %>%
        group_by(!!sym(explanatory[i]), !!sym(dependent)) %>%
        summarise(n = sum(wght)) %>%
        ungroup() %>%
        group_by(!!sym(dependent)) %>%
        mutate(perc = sprintf("%.1f",round(n/sum(n)*100,1))) %>%
        mutate_if(is.numeric, ~formatC(round(.,0), format = "f", big.mark = ",", drop0trailing = TRUE)) %>%
        mutate(n_perc := paste0(n, " (", perc,"%)")) %>%
        select(-n, -perc) %>%
        pivot_wider(names_from = !!sym(dependent), values_from = n_perc) %>%
        rename("levels"=explanatory[i]) %>%
        mutate("characteristic" = explanatory[i]) %>%
        relocate(characteristic)
      
    }
  }
  
  # Combine list together to make dataset
  summary_tbl_wt <- summary_tbl_list %>%
    reduce(full_join)
  
  summary_tbl_wt
}

##### Table of events and event rates by vaccine category
# Calculates the standardised mean differences (smd) between the uv and vacc for each of 
# the categorical explanatory variables for a vaccine type.

# Input:
# - data = the cohort descriptive dataset - z_chrt_desc

# Output:
# A table with weighted millions of person years spent with in each vaccination category in
# the cohort, and weighted count of events and event rates per million person years for
# hospitaliation, death, and hospitalisation or death post vaccination, and more than 14
# days post vaccination

# Table output to be used to plot comparisons between the matched and overall population 
# (before matching - crude)

event_summary_wt <- function(data, wght=eave_weight){
  
  summary_tbl_list <- list()
  
  # First row is person years spent with each vaccination status in cohort
  first_row <-  t(select(data, starts_with('days'))) %*% pull(data, wght)/(365.21 * 1000) 
  
  dependent <- grep('vacc_at', names(data), value = TRUE)
  
  for (i in 1:length(dependent)){
    summary_tbl_list[[i]] <- data %>%
      group_by(!!sym(dependent[i]) ) %>%
      summarise(n = sum(wght)) %>%
      na.omit() %>%
      mutate(rate = sprintf('%.2f',n/first_row))  %>% 
      # format numbers with commas every 3 digits,  
      mutate_if(is.numeric, ~formatC(., format = "f", big.mark = ",", drop0trailing = TRUE)) %>%
      mutate(n = paste0( n, ' (', rate, ')') )  %>%
      select(-rate) %>%
      pivot_wider(names_from = !!sym(dependent[i]), values_from = n) %>%
      mutate(Event = dependent[i])  
  }
  
  # Combine list together to make dataset
  summary_tbl_wt <- summary_tbl_list %>% 
    reduce(full_join) %>%
    mutate(Event = c('Hospitalisation',
                     '14 days prior to hospitalisation',
                     'Death',
                     '14 days prior to death',
                     'Hospitliation or death',
                     '14 days prior to hospitalisation or death')) 
  
  first_row <- formatC(sprintf('%.2f',first_row), format = "f", big.mark = ",", drop0trailing = TRUE)
  
  names(first_row) <- names(summary_tbl_wt)[1:5]
  
  first_row <- data.frame(as.list(first_row), stringsAsFactors = FALSE) %>% 
    mutate(Event = 'Person years (thousands)') %>%
    relocate(Event)
  
  summary_tbl_wt <- bind_rows(first_row, summary_tbl_wt) %>% relocate(uv, .after = Event)
  
  names(summary_tbl_wt) <- c('Event', 'Unvaccinated', 'First dose ChAdOx1', 
                             'Second dose ChAdOx1', 'First dose BNT162b2', 
                             'Second dose BNT162b2')
  
  summary_tbl_wt
}


fun_ve_glm <- function(z_raw){
  z_coef <- cbind(z_raw$coefficients, confint.default(z_raw)) %>% as.data.frame()
  names(z_coef) <- c("est","lcl","ucl")
  z_coef$names <- names(z_raw$coefficients)
  z_coef
}

fun_extract_ests <- function(z_res, factor_name) {
  z <- filter(z_res, grepl(factor_name, names)) %>% 
    dplyr::select(names, est, lcl,ucl) %>% 
    mutate(x_val = gsub(factor_name, "", names)) %>% 
    dplyr::select(-names)
  z <- z %>% mutate(across(est:ucl, ~exp(.)))
  z <- bind_rows(data.frame(est=1,lcl=1,ucl=1,x_val="Ref"), z)
  z
}


fun_plot_rr <- function(z, x_title="", ul=10, ll = 0){
  
  z <- z %>% mutate(plot_order=row_number() ) %>% 
    mutate( x_val=fct_reorder(x_val, plot_order)) 
  
  z <- z %>% mutate(across(est:ucl, ~ if_else(. > ul, ul, .)))
  z <- z %>% mutate(across(est:ucl, ~ if_else(. < ll, ll, .)))
  
  g1 <- z %>% ggplot(aes(x=x_val)) + geom_point(aes(y=est)) + geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.4) + 
    labs(x=x_title, y="Rate Ratio") + theme(axis.text.x=element_text(angle=45))
  g1
}

fun_plot_rr_pub <- function(z, x_title="", ul=10, ll = 0){
  
  z <- z %>% mutate(plot_order=row_number() ) %>% 
    mutate( x_val=fct_reorder(x_val, plot_order)) 
  
  z <- z %>% mutate(across(est:ucl, ~ if_else(. > ul, ul, .)))
  z <- z %>% mutate(across(est:ucl, ~ if_else(. < ll, ll, .)))
  
  g1 <- z %>% ggplot(aes(x=x_val)) + theme_Publication()+ geom_point(aes(y=est)) + geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.4) + 
    labs(x=x_title, y="Rate Ratio") + theme(axis.text.x=element_text(angle=45))
  g1
}

fun_summ_agg <- function(df_res, var){
  z <- df_res %>% group_by_at(var) %>% 
    dplyr::summarise(pyears=sum(pyears), n=sum(n), event=sum(event)) %>% 
    mutate(rate=event/pyears) %>% 
    mutate(RR=rate/first(rate) ) %>% 
    ungroup() %>% as.data.frame()
  z
}

fun_extract_clogit <- function(z) {
  #z is clogit fit
  z_se <-  sqrt(diag(z$var))
  z_out <- cbind.data.frame(RR = z$coefficients, LCL= z$coefficients - 1.96*z_se, UCL= z$coefficients + 1.96*z_se )
  z_out <- exp(z_out)
  z_out
}

fun_plot_tnd_cc <- function(z_est, z_vt, ul=1.5, ll = 0){
  #z_est is output of fun_extract_clogit
  #z_vt <- "Mo"
  z_title <- case_when(z_vt=="PB" ~ "BNT162b2", z_vt=="AZ" ~ "ChAdOx1", z_vt=="Mo" ~ "mRNA-1273", z_vt=="" ~ "Any Vaccine")
  z1 <- z_est %>% filter(grepl(z_vt, rownames(z_est))) %>% 
    mutate(name = gsub(paste0("\\_",z_vt),"", rownames(.))) %>% 
    mutate(name = gsub("^vt","",name)) 
  z1 <- bind_rows(data.frame(RR=1, LCL=1,UCL=1, name="uv"), z1)
  z1 <- z1 %>% mutate(across(RR:UCL, ~ if_else(. > ul, ul, .)))
  z1 <- z1 %>% mutate(across(RR:UCL, ~ if_else(. < ll, ll, .)))
  
  z1 <- z1 %>% mutate(name=factor(name,levels=name))
  g1 <- z1 %>% ggplot(aes(x=name, y=RR)) + geom_point() + 
    geom_errorbar(aes(ymin=LCL, ymax=UCL), width=0.4) +
    theme(axis.text.x=element_text(angle=45)) + 
    labs(x="Vaccine Status", y="Risk Ratio", title=z_title)
  g1
  
}

fun_plot_vs_tnd_cc <- function(z_est, z_vt, z_ref, z_levs, ul=1.5, ll = 0){
  #z_est is output of fun_extract_clogit, z_ref the reference level
  #z_vt <- "vs_d3b"
  z_title <- case_when(z_vt=="vs" ~ "Any Vaccine", z_vt=="vs_d3b" ~ "Any Vaccine  - Booster Moderna or Pfizer" )
  z1 <- z_est %>% filter(grepl(z_vt, rownames(z_est))) %>% 
    mutate(name = gsub(paste0("^",z_vt),"", rownames(.))) 
  z1 <- bind_rows(data.frame(RR=1, LCL=1,UCL=1, name=z_ref), z1)
  z1 <- z1 %>% mutate(across(RR:UCL, ~ if_else(. > ul, ul, .)))
  z1 <- z1 %>% mutate(across(RR:UCL, ~ if_else(. < ll, ll, .)))
  
  z1 <- z1 %>% mutate(name=factor(name,levels=z_levs)) %>% 
    arrange(name)
  if (z_vt=="vs_d3b") z1 <- z1 %>% filter(!grepl("AZ$", name))
    g1 <- z1 %>% ggplot(aes(x=name, y=RR)) + geom_point() + 
    geom_errorbar(aes(ymin=LCL, ymax=UCL), width=0.4) +
    theme(axis.text.x=element_text(angle=45)) + 
    labs(x="Vaccine Status", y="Risk Ratio", title=z_title)
  g1
  
}

fun.extract.coxph <- function(z.fit) {
  #takes a coxph filt using penalised splines and drops off the ps terms
  #make sure no variable begins ps
  z <- summary(z.fit)
  z <- data.frame(z$conf.int)
  z <- z %>% mutate(names = row.names(z)) %>% 
    filter(!(grepl("^ps", names))) %>% 
    dplyr::relocate(names, .before=1) %>% 
    dplyr::select(-exp..coef.)
  names(z) <- c("names","HR","LCL","UCL")
  z
}

fun_print_hr <- function(z,z_var){
  # z_var <- "vsv2_14+"
  z1 <- filter(z,names==z_var)
  z1[c("HR","LCL","UCL")] <- round(z1[c("HR","LCL","UCL")], 2)
  str_c(z1["HR"], ", 95% CI (", z1["LCL"], ", ", z1["UCL"],")")
}

plot_HR <- function(model_fit, term){
  # plots hazard ratios for a single term in a fitted model
  
  hr <- termplot(model_fit, term = term, se = T, plot = F)
  
  var <- names(hr)
  
  hr <- hr[[var]]
  
  hr <- mutate(hr, ucl = y + 1.96*se,
               lcl = y - 1.96*se) %>%
    mutate_at(c('y', 'ucl', 'lcl'), exp)
  
  hr <- do.call(data.frame,lapply(hr, function(x) replace(x, is.infinite(x),NA)))
  
  output <- ggplot(data=hr, aes(x=x, y=y)) + geom_line() +
    geom_ribbon(aes(ymin=lcl, ymax=ucl), linetype=2, alpha=0.1, fill = 'steelblue')  + 
    ylab("Hazard Ratio")
  
  if (var == 'age_year'){
    output <- output + xlab("Age")
  } else if (var == 'days'){
    output <- output + xlab("Calendar Time - days since beginning of study")
  }
  
  output
}
