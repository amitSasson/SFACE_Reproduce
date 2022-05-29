library(tidyverse)
df_read <- read_csv("/udd/nhasn/amit_data/smoking/data_Amit_20211201_smoking.csv")
##
df_raw <- df_read



# Func to calc who is removed when each filter is applied
stats_func <- function(df, filter_used)
{
  n <- df %>% pull(id) %>% n_distinct()

  msin_missing <- df %>%
    filter(crc == 1) %>%
    mutate(msin_missing = ifelse(is.na(msin), 1, 0)) %>%
    group_by(id) %>%
    mutate(msin_missing_max = max(msin_missing)) %>%
    select(id, msin_missing_max) %>%
    distinct() %>%
    ungroup() %>%
    pull(msin_missing_max) %>%
    sum()

  msin <- df %>%
    filter(!is.na(msin)) %>%
    group_by(id) %>%
    mutate(msin_max = max(msin)) %>%
    select(id, msin_max) %>%
    distinct() %>%
    ungroup() %>%
    pull(msin_max) %>%
    table()

  ever_smoke <- df %>%
    drop_na(packyr) %>%
    filter(ageyr > (baseline-10) & ageyr <= baseline ) %>%
    mutate(ever_smoke = ifelse(packyr > 0, 1, 0)) %>%
    group_by(id) %>%
    mutate(ever_smoke_max = max(ever_smoke)) %>%
    select(id, ever_smoke_max) %>%
    distinct() %>%
    ungroup() %>%
    pull(ever_smoke_max) %>%
    mean()

  stats <- tibble("msin == 0" = msin[1],
                  "msin == 1" = msin[2],
                  "msin_missing" = msin_missing,
                  "n" = n,
                  "ever_smoke_max" = ever_smoke,
                  filter_used = filter_used)

  return(stats)
}

clean_data <- function(df_raw, baseline, outcome_age)
{
  stats <- tibble("msin == 0" = numeric(),
                  "msin == 1" = numeric(),
                  "msin_missing" = numeric(),
                  "n" = numeric(),
                  "ever_smoke_max" = numeric(),
                  filter_used = character())

  ##add irt in months since 1900 (colname - irt_amit)
  df_raw <- df_raw %>% mutate(irt_amit = case_when(interval == 1 ~ irt80,
                                                   interval == 2 ~ irt82,
                                                   interval == 3 ~ irt84,
                                                   interval == 4 ~ irt86,
                                                   interval == 5 ~ irt88,
                                                   interval == 6 ~ irt90,
                                                   interval == 7 ~ irt92,
                                                   interval == 8 ~ irt94,
                                                   interval == 9 ~ irt96,
                                                   interval == 10 ~ irt98,
                                                   interval == 11 ~ irt00,
                                                   interval == 12 ~ irt02,
                                                   interval == 13 ~ irt04,
                                                   interval == 14 ~ irt06,
                                                   interval == 15 ~ irt08,
                                                   interval == 16 ~ irt10,
                                                   interval == 17 ~ irt12))

  # some irt rows are in text format (e.g. "jun 1992"), translating them all to number of months since 1900
  df_raw <- df_raw %>%
    mutate(irt_amit_numeric = as.numeric(irt_amit)) %>%
    separate(irt_amit, c("irt_amit_month", "irt_amit_year")) %>%
    mutate(irt_amit_year = as.numeric(irt_amit_year)) %>%
    mutate(irt_amit_month = tolower(irt_amit_month))

  df_raw <- df_raw %>% mutate(irt_amit_month_numeric = case_when(irt_amit_month == "jan" ~ 1,
                                                                 irt_amit_month == "feb" ~ 2,
                                                                 irt_amit_month == "mar" ~ 3,
                                                                 irt_amit_month == "apr" ~ 4,
                                                                 irt_amit_month == "may" ~ 5,
                                                                 irt_amit_month == "jun" ~ 6,
                                                                 irt_amit_month == "jul" ~ 7,
                                                                 irt_amit_month == "aug" ~ 8,
                                                                 irt_amit_month == "sep" | irt_amit_month == "sept" ~ 9,
                                                                 irt_amit_month == "oct" ~ 10,
                                                                 irt_amit_month == "nov" ~ 11,
                                                                 irt_amit_month == "dec" ~ 12))
  df_raw <- df_raw %>%
    mutate(irt_amit_numeric_calc = (irt_amit_year - 1900)*12 + irt_amit_month_numeric) %>%
    mutate(irt_amit = case_when(!is.na(irt_amit_numeric) ~ irt_amit_numeric,
                                !is.na(irt_amit_numeric_calc) ~ irt_amit_numeric_calc))

  df_raw$irt_amit %>% is.na() %>% sum() #make sure it is 0 - every row has irt

  ## select columns
  df_raw <- df_raw %>%
  #  select(id, irt_amit, gender, msin, regasp, interval, agedx, ageyr, dt_crc,
  #         time_crc, dtdth, dtdx, bmicum, mvt, nsaid, actmcum, packyr,
  #         caloriecum, alcocum, endo, caloriecum, crcfh) %>%
    mutate(msin = na_if(msin, "999")) %>%
    mutate(packyr = na_if(packyr, "999=Missing")) %>%
    mutate(age_of_death = dtdth/12 - irt_amit/12  + ageyr) %>% #age of death
    group_by(id) %>%
    mutate(crc = max(dt_crc)) %>% # crate crc = this person had/have/will ahve crc at some point
    ungroup()

  df_raw <- df_raw %>% filter(!(crc == 0 & !is.na(msin)))

  #remove data after death (just to be sure - there isn't supposed to be any data after death)
  df_raw <-  df_raw %>% filter(dtdth - irt_amit  >= 0 | is.na(dtdth) | is.na(irt_amit))


  #remove data after diagnosis (just to be sure - there isn't supposed to be any data after diagnosis)
  df_raw <-  df_raw %>% filter(dtdx - irt_amit >= 0 | is.na(dtdx) | is.na(irt_amit))

  stats <- rbind(stats, stats_func(df_raw, "before any filter"))


  ##### age 50 #######

  #keep baseline-10 < age < baseline
  df_baseline <- df_raw %>% filter(ageyr > (baseline-10) & ageyr <= baseline )
  stats <- rbind(stats, stats_func(df_baseline, "no data for baseline"))



  #keep the oldest row up until baseline
  df_baseline <- df_baseline %>% arrange(id, desc(ageyr)) %>% distinct(id, .keep_all = T)

  #Filter those who were diagnosed before baseline
  df_baseline <- df_baseline %>% filter(agedx >= baseline | is.na(agedx))
  stats <- rbind(stats, stats_func(df_baseline, "diagnosed before baseline"))

  #Filter those who were diagnosed before baseline
  df_baseline <- df_baseline %>% filter(!is.na(packyr))
  stats <- rbind(stats, stats_func(df_baseline, "missing smoking at baseline data"))

  #Filter those who died before outcome_age-5 without being diagnosed
  df_baseline <- df_baseline %>% filter(!(age_of_death < outcome_age-5 & crc == 0) | is.na(age_of_death))
  stats <- rbind(stats, stats_func(df_baseline, "died before outcome_age - 5, without crc"))

  #change to undiagnosed those who were diagnosed after outcome_age
  df_baseline <- df_baseline %>%
                        mutate(crc = ifelse(agedx < outcome_age | is.na(agedx), crc, 0)) %>%
                        mutate(msin = ifelse(agedx < outcome_age | is.na(agedx), msin, NA))   %>%
                        mutate(agedx = ifelse(agedx < outcome_age | is.na(agedx), agedx, NA))

  stats <- rbind(stats, stats_func(df_baseline, "change diagnosed after outcome_age"))
  stats

  outcome_age_data <- df_raw %>%
                         filter(id %in% df_baseline$id) %>%
                         filter(ageyr <= outcome_age & ageyr >= outcome_age - 2) %>%
                         pull(id)


  # who don't have data at outcome_age+
  df_baseline <- df_baseline %>% mutate(outcome_age_data = (id %in% outcome_age_data) | crc == 1)


  stats <- rbind(stats, stats_func(df_baseline %>% filter(outcome_age_data == TRUE), "has outcome_age data"))
  print(stats)

  #a more relaxed version - filter who does not follow it - only one!
  ten_years_before_outcome_age_data <- df_raw %>%
    filter(id %in% df_baseline$id) %>%
    filter(ageyr <= outcome_age & ageyr >= outcome_age - 10) %>%
    pull(id)
  df_baseline <- df_baseline %>% mutate(ten_years_before_outcome_age_data = (id %in% ten_years_before_outcome_age_data) | crc == 1)
  df_baseline <- df_baseline %>% filter(ten_years_before_outcome_age_data)

  return(df_baseline)
}

baseline <- 60
outcome_age <- 70
df_baseline  <- clean_data(df_read, baseline, outcome_age)

write_csv(df_baseline, "data_after_pre_process_for_weights_60_70.csv")

