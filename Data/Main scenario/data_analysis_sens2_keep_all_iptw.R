library(tidyverse)
library(patchwork)
library(nnet)
library(boot)
library(latex2exp)

#### setup 
set.seed(31415)
q <- qnorm(0.975)
M <- 100000
effect_scale <- "RR"

df <- read_csv("/udd/nhasn/data_after_weights_for_analysis_keep_all_60_70.csv") 

df$gender <- as.factor(df$gender)
df$regasp <- as.factor(df$regasp)
df$crcfh <- as.factor(df$crcfh)
df$endo <- as.factor(df$endo)

#### sensitivity analysis #####
stand_sesnsitivity_ps_ps_calc_p <- function(df, sface)
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
  n_w <- sum(df$weight)
  
  df$w_A <- ifelse(df$A == 1, pr_A_1/pred_A, (1-pr_A_1)/(1-pred_A) ) #Stabilized weights
  q99 <- quantile(df$w_A, .99)
  df$w_A <- ifelse(df$w_A > q99, q99, df$w_A)
  n_w <- sum(df$weight)
  
  p_Y11_A1 <- sum(df$weight*df$w_A*df$A*df[,self])/sum(df$weight*df$A)
  p_Y11_A0 <- sum(df$weight*df$w_A*(1-df$A)*df[,self])/sum(df$weight*(1-df$A))
  p_Y21_A0 <- sum(df$weight*df$w_A*(1-df$A)*df[,other])/sum(df$weight*(1-df$A))
  p_Y20_A1 <- 1- sum(df$weight*df$w_A*df$A*df[,other])/sum(df$weight*df$A)
  
  return(c(p_Y11_A1, p_Y11_A0, p_Y21_A0, p_Y20_A1))
}

probs1 <- stand_sesnsitivity_ps_ps_calc_p(df, "sface1")
probs2 <- stand_sesnsitivity_ps_ps_calc_p(df, "sface2")

#fun for the effect of the subtype under S-Mono (usually sface1) with sensitivity analysis 
stand_sesnsitivity_ps_ps <- function(lambda2 = 0, lambda1 = 0, effect_scale, probs) 
{
  p_Y11_A1 <- probs[1]
  p_Y11_A0 <- probs[2]
  p_Y21_A0 <- probs[3]
  p_Y20_A1 <- probs[4]

  if (effect_scale == "diff") 
  {return(M*(p_Y11_A1-lambda2*p_Y21_A0+(lambda1-1)*p_Y11_A0)/(p_Y20_A1-lambda2*p_Y21_A0))}
  if (effect_scale == "RR") 
  {return((p_Y11_A1-lambda2*p_Y21_A0)/((1-lambda1)*p_Y11_A0))}
}


#choose lambda values for sensitivity analysis
max_lambda_func <- function(df, sface)
{
  self <- ifelse(sface == "sface1", "y1", "y2")
  other <- ifelse(sface == "sface1", "y2", "y1") 
  
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
  
  ratio <- sum(df$relative_weight*pred_treat[,self])/sum(df$relative_weight*pred_untr[,other])
  
  return(min(1, ratio))
}

boot_func_sface <- function(df, sface, BS = 200)
{
  boot_estimand_func <- function(df, index, sface) 
  {
    df <- df[index,]
    return(stand_sesnsitivity_ps_ps_calc_p(df,  sface))
  }
  
  boot_estimators <- boot(data = df, 
                          statistic = boot_estimand_func, 
                          sface = sface, 
                          R = BS,
                          parallel = "multicore", 
                          ncpus = 8)
  return(boot_estimators$t)
}

boot_p_sface1 <- boot_func_sface(df, "sface1")
boot_p_sface2 <- boot_func_sface(df, "sface2")
max_lambda2 <- max_lambda_func(df, "sface1")
max_lambda1 <- max_lambda_func(df, "sface2")

lambda2_vals <- seq(0, max_lambda2, by=0.02)
lambda1_vals <- seq(0, max_lambda1, by=0.02)

sface_sensitiviey <- tibble(lambda2 = numeric(), 
                            lambda1 = numeric(),  
                            estimator_sface1 = numeric(),  
                            estimator_sface2 = numeric(),  
                            estimator_theta = numeric(), 
                            sd_sface1 = numeric(), 
                            sd_sface2 = numeric(), 
                            sd_theta = numeric())


for (lambda2_val in lambda2_vals) {
  for (lambda1_val in lambda1_vals) {
    
    sface1 <- stand_sesnsitivity_ps_ps(lambda2_val, lambda1_val, effect_scale, probs = probs1) 
    sface2 <- stand_sesnsitivity_ps_ps(lambda1_val, lambda2_val, effect_scale, probs = probs2) 
    theta <- sface1-sface2
    
    sface1_boot <- apply(boot_p_sface1, 1, stand_sesnsitivity_ps_ps, 
                         lambda1 = lambda2_val, 
                         lambda2 = lambda1_val, 
                         effect_scale = effect_scale)
    sd_sface1 <- sd(sface1_boot)
    
    sface2_boot <- apply(boot_p_sface2, 1, stand_sesnsitivity_ps_ps, 
                         lambda1 = lambda1_val, 
                         lambda2 = lambda2_val, 
                         effect_scale = effect_scale)
    sd_sface2 <- sd(sface2_boot)
    
    theta_boot <- sface1_boot-sface2_boot
    sd_theta <- sd(theta_boot)
    
    
    sface_sensitiviey <- sface_sensitiviey %>% add_row(tibble_row(lambda2 = lambda2_val, 
                                                                  lambda1 = lambda1_val, 
                                                                  estimator_sface1 = sface1, 
                                                                  estimator_sface2 = sface2, 
                                                                  estimator_theta = theta, 
                                                                  sd_sface1 = sd_sface1, 
                                                                  sd_sface2 = sd_sface2, 
                                                                  sd_theta = sd_theta))
  }
}


save(sface_sensitiviey, file = str_c(as.character(effect_scale), "_sens2_iptw_",Sys.time(),".RData"))


