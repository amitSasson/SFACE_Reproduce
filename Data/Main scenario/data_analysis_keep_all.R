library(tidyverse)
library(patchwork)
library(nnet)
library(boot)
library(TheSFACE)
#### setup
set.seed(31415)
q <- qnorm(0.975)
M <- 100000
df <- read_csv("/udd/nhasn/data_after_weights_for_analysis_keep_all_60_70.csv")

df$gender <- as.factor(df$gender)
df$regasp <- as.factor(df$regasp)
df$crcfh <- as.factor(df$crcfh)
df$endo <- as.factor(df$endo)

#### Different estimands under S-Monotonicity ####
#func to calculate sface1 or sface2 using standardization for a given df
stand <- function(df, sface = "sface1")
{
  subtype <- ifelse(sface == "sface1", 1, 2)
  val <- sface(stand_formula = y ~ A + X1 + X2,
               iptw_formula = A ~ X1 + X2,
               exposure = "A",
               outcome = "y",
               df = df,
               subtype = subtype,
               scale = "diff",
               method = "stand",
               weight = "weight",
               MultPer = M)
  val$sface %>% unlist()

  }

#func to calculate sface1 or sface2 using IPTW for a given df
iptw <- function(df, sface = "sface1")
{
  subtype <- ifelse(sface == "sface1", 1, 2)
  val <- sface(stand_formula = y ~ A + X1 + X2,
               iptw_formula = A ~ X1 + X2,
               exposure = "A",
               outcome = "y",
               df = df,
               subtype = subtype,
               scale = "diff",
               method = "IPTW",
               weight = "weight",
               MultPer = M)
  val$sface %>% unlist()
}


#func to calculate sface1 or sface2 using the total naive estimand using standartization
stand_total_naive <- function(df, sface = "sface1")
{
  fit_y_by_A_X <- multinom(y ~  A + gender + bmicum_cat + regasp + alcocum + caloriecum + crcfh + endo + actmcum,
                           df,
                           trace = FALSE,
                           weights = weight)

  df_treat <- df_untr <- df
  df_treat$A <- 1
  df_untr$A <- 0

  pred_treat <- as.data.frame(predict(fit_y_by_A_X, newdata = df_treat, type = "probs"))
  pred_untr <- as.data.frame(predict(fit_y_by_A_X, newdata = df_untr, type = "probs"))

  self <- ifelse(sface == "sface1", "1", "2")
  other <- ifelse(sface == "sface1", "2", "1")

  n <- nrow(df)
  return(M*(sum(df$relative_weight*pred_treat[,self]) - sum(df$relative_weight*pred_untr[,self])))
}


#func to calculate sface1 or sface2 using a Doubly-robust estimator for a given df
DR <- function(df, sface = "sface1")
{
  subtype <- ifelse(sface == "sface1", 1, 2)
  val <- sface(stand_formula = y ~ A + X1 + X2,
               iptw_formula = A ~ X1 + X2,
               exposure = "A",
               outcome = "y",
               df = df,
               subtype = subtype,
               scale = "diff",
               method = "DR",
               weight = "weight",
               MultPer = M)
  val$sface %>% unlist()
}


#func to calculate sface1-sface2 for a given sample / population and estimand
diff <- function(df, FUN)
{
  return(FUN(df, "sface1") - FUN(df, "sface2"))
}



#### BOOTSTRAP function
#func to calculate the variance of an estimetor in a specific sample, using BS
boot_func <- function(df, estimand_func, sface = "sface1", BS = 200)
{
  boot_estimand_func <- function(df, index, estimand_func, sface = "sface1")
  {
    df <- df[index,]
    return(estimand_func(df, sface))
  }

  boot_estimators <- boot(data = df,
                          statistic = boot_estimand_func,
                          estimand_func = estimand_func,
                          sface = sface,
                          R = BS,
                          parallel = "multicore",
                          ncpus = 8)
  return(c(sd(boot_estimators$t), quantile(boot_estimators$t, 0.0275), quantile(boot_estimators$t, 0.975)))
}



#### Estimating sface1 and sface2 in multiple dfs

#func to caculate sface1 or sface2 for each sample and compute stats
compute_stats_sface <- function(df, sface = "sface1")
{
  #for each estimand
  estimand_compute_stats <- function(estimand_func, df, sface)
  {
    estimator <- estimand_func(df, sface = sface)
    boot_sample <- boot_func(df, estimand_func, sface = sface)
    se_boot <- boot_sample[1]
    quantile_low <-  boot_sample[2]
    quantile_high <-  boot_sample[3]
    return(data.frame(estimator = estimator,
                      se_boot = se_boot,
                      quantile_low = quantile_low,
                      quantile_high = quantile_high))
  }

  estimands <- list("stand" = stand,
                    "iptw" = iptw,
                    "stand_total_naive" = stand_total_naive,
                    "DR" = DR)

  stats <- sapply(estimands, estimand_compute_stats, df, sface, USE.NAMES = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("value") %>%
    unnest(cols = c(stand, iptw, stand_total_naive, DR)) %>%
    mutate_if(is.numeric, round, 6) %>%
    pivot_longer(cols = -c(value), names_to = "estimand", values_to = "res") %>%
    pivot_wider(id_cols = estimand, names_from = value, values_from = res) %>%
    mutate(CI_low = estimator-q*se_boot,
           CI_high = estimator+q*se_boot)

  return(stats)
}

#calculate the estimators and variances using BS
# #
# stand(df, "sface1")
# iptw(df, "sface1")
# stand_total_naive(df, "sface1")
# stand_cond_naive(df, "sface1")
#
# stand(df, "sface2")
# iptw(df, "sface2")
# stand_total_naive(df, "sface2")
# stand_cond_naive(df, "sface2")

ans <- compute_stats_sface(df,"sface1")
save(ans, file = str_c("data_diff_sface1_keep_all_60_70", ".RData"))

ans <- compute_stats_sface(df,"sface2")
save(ans, file = str_c("data_diff_sface2_keep_all_60_70", ".RData"))


###### RR

stand_RR <- function(df, sface = "sface1")
{
  subtype <- ifelse(sface == "sface1", 1, 2)
  val <- sface(stand_formula = y ~ A + X1 + X2,
               iptw_formula = A ~ X1 + X2,
               exposure = "A",
               outcome = "y",
               df = df,
               subtype = subtype,
               scale = "RR",
               method = "stand",
               weight = "weight",
               MultPer = M)
  val$sface %>% unlist()
}

#func to calculate sface1 or sface2 using IPTW for a given df
iptw_RR <- function(df, sface = "sface1")
{
  subtype <- ifelse(sface == "sface1", 1, 2)
  val <- sface(stand_formula = y ~ A + X1 + X2,
               iptw_formula = A ~ X1 + X2,
               exposure = "A",
               outcome = "y",
               df = df,
               subtype = subtype,
               scale = "RR",
               method = "IPTW",
               weight = "weight",
               MultPer = M)
  val$sface %>% unlist()
}


#func to calculate sface1 or sface2 using the total naive estimand using standartization
stand_total_naive_RR <- function(df, sface = "sface1")
{
  fit_y_by_A_X <- multinom(y ~  A + gender + bmicum_cat + regasp + alcocum + caloriecum + crcfh + endo + actmcum,
                           df,
                           trace = FALSE,
                           weights = weight)

  df_treat <- df_untr <- df
  df_treat$A <- 1
  df_untr$A <- 0

  pred_treat <- as.data.frame(predict(fit_y_by_A_X, newdata = df_treat, type = "probs"))
  pred_untr <- as.data.frame(predict(fit_y_by_A_X, newdata = df_untr, type = "probs"))

  self <- ifelse(sface == "sface1", "1", "2")
  other <- ifelse(sface == "sface1", "2", "1")

  return((sum(df$weight*pred_treat[,self]) / sum(df$weight*pred_untr[,self])))
}

DR_RR <- function(df, sface = "sface1")
{
  subtype <- ifelse(sface == "sface1", 1, 2)
  val <- sface(stand_formula = y ~ A + X1 + X2,
               iptw_formula = A ~ X1 + X2,
               exposure = "A",
               outcome = "y",
               df = df,
               subtype = subtype,
               scale = "RR",
               method = "DR",
               weight = "weight",
               MultPer = M)
  val$sface %>% unlist()
}


#### Estimating sface1 and sface2 in multiple dfs

#func to caculate sface1 or sface2 for each sample and compute stats
compute_stats_sface_RR <- function(df, sface = "sface1")
{
  #for each estimand
  estimand_compute_stats <- function(estimand_func, df, sface)
  {
    estimator <- estimand_func(df, sface = sface)
    boot_sample <- boot_func(df, estimand_func, sface = sface)
    se_boot <- boot_sample[1]
    quantile_low <-  boot_sample[2]
    quantile_high <-  boot_sample[3]
    return(data.frame(estimator = estimator,
                      se_boot = se_boot,
                      quantile_low = quantile_low,
                      quantile_high = quantile_high))
  }

  estimands <- list("stand" = stand_RR,
                    "iptw" = iptw_RR,
                    "stand_total_naive" = stand_total_naive_RR,
                    "DR" = DR_RR)

  stats <- sapply(estimands, estimand_compute_stats, df, sface, USE.NAMES = TRUE) %>%
    as.data.frame() %>%
    rownames_to_column("value") %>%
    unnest(cols = c(stand, iptw, stand_total_naive, DR)) %>%
    mutate_if(is.numeric, round, 6) %>%
    pivot_longer(cols = -c(value), names_to = "estimand", values_to = "res") %>%
    pivot_wider(id_cols = estimand, names_from = value, values_from = res) %>%
    mutate(CI_low = estimator-q*se_boot,
           CI_high = estimator+q*se_boot)

  return(stats)
}

#ans <- compute_stats_sface_RR(df,"sface1")
#save(ans, file = str_c("data_RR_sface1_keep_all_60_70",".RData"))

#ans <- compute_stats_sface_RR(df,"sface2")
#save(ans, file = str_c("data_RR_sface2_keep_all_60_70", ".RData"))


