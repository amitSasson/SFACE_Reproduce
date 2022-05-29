library(tidyverse)
library(patchwork)
library(nnet)
library(boot)

#### setup 
set.seed(31415)
q <- qnorm(0.975)
M <- 100000
df <- read_csv("/udd/nhasn/data_after_weights_for_analysis_has_outcome_data_60_70.csv") 


#### Different estimands under S-Monotonicity ####
#func to calculate sface1 or sface2 using standardization for a given df
stand <- function(df, sface = "sface1") 
{
  #model y ~ A + X
  fit_y_by_A_X <- multinom(y ~ A + gender + bmicum_cat + regasp + alcocum + caloriecum + crcfh + endo + actmcum,
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
  
  n_w <- sum(df$weight)
  p_Y11_A1 <- sum(df$weight*pred_treat[,self])/n_w
  p_Y11_A0 <- sum(df$weight*pred_untr[,self])/n_w
  p_Y12_A1 <- sum(df$weight*pred_treat[,other])/n_w
  
  #p_stand <- c(p_Y11_A1, p_Y11_A0, p_Y12_A1)
  
  return(M*(p_Y11_A1 - p_Y11_A0)/(1-p_Y12_A1)) 
}

#func to calculate sface1 or sface2 using IPTW for a given df
iptw <- function(df, sface = "sface1") 
{
  self <- ifelse(sface == "sface1", "y1", "y2")
  other <- ifelse(sface == "sface1", "y2", "y1")
  
  #model A ~ X
  fit_A_by_X <- glm(A ~ gender + bmicum_cat + regasp + alcocum + caloriecum + crcfh + endo + actmcum,
                    df, 
                    family = "binomial", 
                    weights = weight)
  
  pred_A <-  predict(fit_A_by_X, type = "response")
  pr_A_1 <- weighted.mean(x = df$A, w = df$weight)
  n_w <- sum(df$weight)
  
  df$w_A <- ifelse(df$A == 1, pr_A_1/pred_A, (1-pr_A_1)/(1-pred_A) ) #Stabilized weights
  
  q99 <- quantile(df$w_A, .99)
  df$w_A <- ifelse(df$w_A > q99, q99, df$w_A)
  
  p_Y11_A1 <- sum(df$weight*df$w_A*df$A*df[,self])/sum(df$weight*df$A)
  p_Y11_A0 <- sum(df$weight*df$w_A*(1-df$A)*df[,self])/sum(df$weight*(1-df$A))
  p_Y12_A1 <- sum(df$weight*df$w_A*df$A*df[,other])/sum(df$weight*df$A)
  
  #p_iptw <- c(p_Y11_A1, p_Y11_A0, p_Y12_A1)
  
  return(M*(p_Y11_A1 - p_Y11_A0)/(1-p_Y12_A1)) 
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
  self <- ifelse(sface == "sface1", "y1", "y2")
  other <- ifelse(sface == "sface1", "y2", "y1")  
  
  #model A ~ X
  fit_A_by_X <- glm(A ~ gender + bmicum_cat + regasp + alcocum + caloriecum + crcfh + endo + actmcum,
                    df, 
                    family = "binomial", 
                    weights = weight)
  
  pred_A <-  predict(fit_A_by_X, type = "response")
  pr_A_1 <- mean(df$A)
  
  #model y ~ A + X
  fit_y_by_A_X <- multinom(y ~ A + gender + bmicum_cat + regasp + alcocum + caloriecum + crcfh + endo + actmcum,
                           df, 
                           trace = FALSE, 
                           weights = weight)
  
  df_treat <- df_untr <- df 
  df_treat$A <- 1 
  df_untr$A <- 0 
  
  pred_treat <- as.data.frame(predict(fit_y_by_A_X, newdata = df_treat, type = "probs"))  
  colnames(pred_treat) <-c("0","y1", "y2")
  pred_untr <- as.data.frame(predict(fit_y_by_A_X, newdata = df_untr, type = "probs"))
  colnames(pred_untr) <-c("0","y1", "y2")
  
  n_w <- sum(df$weight)
  
  p_Y11_A1 <- sum(df$weight*(df$A*df[,self]/(pred_A) - ((df$A-pred_A)*pred_treat[,self])/pred_A))/n_w
  p_Y11_A0 <- sum(df$weight*((1-df$A)*df[,self]/(1-pred_A) + ((df$A-pred_A)*pred_untr[,self])/(1-pred_A)))/n_w 
  p_Y12_A1 <- sum(df$weight*(df$A*df[,other]/(pred_A) - ((df$A-pred_A)*pred_treat[,other])/pred_A))/n_w 
  
  #p_DR <- c(p_Y11_A1, p_Y11_A0, p_Y12_A1)
  
  return(M*(p_Y11_A1 - p_Y11_A0)/(1-p_Y12_A1))     
}


#func to calculate sface1-sface2 for a given sample / population and estimand 
diff <- function(df, FUN)
{
  return(FUN(df, "sface1") - FUN(df, "sface2"))
}



#### BOOTSTRAP function 
#func to calculate the variance of an estimetor in a specific sample, using BS
boot_func <- function(df, estimand_func, sface = "sface1", BS = 100)
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
save(ans, file = str_c("data_diff_sface1_has_outcome_data_60_70", ".RData"))

ans <- compute_stats_sface(df,"sface2")
save(ans, file = str_c("data_diff_sface2_has_outcome_data_60_70", ".RData"))


###### RR 

stand_RR <- function(df, sface = "sface1") 
{
  #model y ~ A + X
  fit_y_by_A_X <- multinom(y ~ A + gender + bmicum_cat + regasp + alcocum + caloriecum + crcfh + endo + actmcum,
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

#func to calculate sface1 or sface2 using IPTW for a given df
iptw_RR <- function(df, sface = "sface1") 
{
  #model A ~ X
  fit_A_by_X <- glm(A ~ gender + bmicum_cat + regasp + alcocum + caloriecum + crcfh + endo + actmcum,
                    df, 
                    family = "binomial")
  pred_A <-  predict(fit_A_by_X, type = "response")
  pr_A_1 <- mean(df$A)
  
  self <- ifelse(sface == "sface1", "y1", "y2")
  other <- ifelse(sface == "sface1", "y2", "y1")
  
  df$w_A <- ifelse(df$A == 1, pr_A_1/pred_A, (1-pr_A_1)/(1-pred_A) ) #Stabilized weights
  q99 <- quantile(df$w_A, .99)
  df$w_A <- ifelse(df$w_A > q99, q99, df$w_A)
  
  p_Y11_A1 <- sum(df$weight*df$w_A*df$A*df[,self])/sum(df$weight*df$A)
  p_Y11_A0 <- sum(df$weight*df$w_A*(1-df$A)*df[,self])/sum(df$weight*(1-df$A))
  
  return(p_Y11_A1/ p_Y11_A0) 
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
  self <- ifelse(sface == "sface1", "y1", "y2")
  other <- ifelse(sface == "sface1", "y2", "y1")  
  
  #model A ~ X
  fit_A_by_X <- glm(A ~ gender + bmicum_cat + regasp + alcocum + caloriecum + crcfh + endo + actmcum,
                    df, 
                    family = "binomial", 
                    weights = weight)
  
  pred_A <-  predict(fit_A_by_X, type = "response")
  pr_A_1 <- mean(df$A)
  
  #model y ~ A + X
  fit_y_by_A_X <- multinom(y ~ A + gender + bmicum_cat + regasp + alcocum + caloriecum + crcfh + endo + actmcum,
                           df, 
                           trace = FALSE, 
                           weights = weight)
  
  df_treat <- df_untr <- df 
  df_treat$A <- 1 
  df_untr$A <- 0 
  
  pred_treat <- as.data.frame(predict(fit_y_by_A_X, newdata = df_treat, type = "probs"))  
  colnames(pred_treat) <-c("0","y1", "y2")
  pred_untr <- as.data.frame(predict(fit_y_by_A_X, newdata = df_untr, type = "probs"))
  colnames(pred_untr) <-c("0","y1", "y2")
  
  p_Y11_A1 <- sum(df$weight*(df$A*df[,self]/(pred_A) - ((df$A-pred_A)*pred_treat[,self])/pred_A)) 
  p_Y11_A0 <- sum(df$weight*((1-df$A)*df[,self]/(1-pred_A) + ((df$A-pred_A)*pred_untr[,self])/(1-pred_A))) 
  
  return(p_Y11_A1 /p_Y11_A0)     
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


ans <- compute_stats_sface_RR(df,"sface1")
save(ans, file = str_c("data_RR_sface1_has_outcome_data_60_70",".RData"))

ans <- compute_stats_sface_RR(df,"sface2")
save(ans, file = str_c("data_RR_sface2_has_outcome_data_60_70", ".RData"))


