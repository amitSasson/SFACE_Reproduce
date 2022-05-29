args = commandArgs(trailingOnly=TRUE)
args <- as.numeric(args)

library(tidyverse)
library(patchwork)
library(nnet)
library(parallel)
library(boot)
library(TheSFACE)

#### setup
set.seed(31415)
M <- 100000

#func to generate the data (Y1,Y2 ~ U + X + A) under S-Monotonicity
generate_df_covariates <- function(n.sample = 1000,
                                   y1_par = y1_par,
                                   y2_par = y2_par,
                                   A_par = A_par)
{
  U <-
    rnorm(n = n.sample, 0, 1)
  X1 <- rbinom(n = n.sample, size = 1, prob = 0.5)
  X2 <- rnorm(n = n.sample, mean = 0, sd = 1)
  X <- matrix(c(X1,X2), nrow = n.sample, byrow = FALSE)

  #func to create y (equals 0/1/2) using A, U and parameters
  po_probs <- function(U = U, X = X, A = A, y1_par, y2_par)
  {
    e.y1 <- exp(y1_par$gamma0 + y1_par$gammaU*U + y1_par$gammaA*A + c(y1_par$gammaX%*%(t(X))))
    e.y2 <- exp(y2_par$gamma0 + y2_par$gammaU*U + y2_par$gammaA*A + c(y2_par$gammaX%*%(t(X))))

    p.y1 <- e.y1 / (1 + e.y1 + e.y2)
    p.y2 <- e.y2 / (1 + e.y1 + e.y2)

    probs <- data.frame(p.y1, p.y2, p.y0 = 1-p.y1-p.y2)

    return(probs)
  }

  #create y, y1, y2 under A = 0
  y_A0_probs <- po_probs(U = U, X = X, A = 0, y1_par, y2_par)
  y_A0 <- apply(y_A0_probs, 1, function(i) {sample(c(1,2,0), 1, prob = i)})
  y1_A0 <- ifelse(y_A0 == 1, 1, 0)
  y2_A0 <- ifelse(y_A0 == 2, 1, 0)

  #create y under A = 1
  #before adj
  y_A1_probs <- po_probs(U = U, X = X, A = 1, y1_par, y2_par)
  #asjust Using Law of total probability: P(y(A=1)=1) = P(y(A=0)=1) + P(y(A=0)=0)*p
  y_A1_probs <- (y_A1_probs - y_A0_probs)/(y_A0_probs$p.y0)
  #replace negative values with zero
  y_A1_probs[y_A1_probs < 0] <- 0
  y_A1_probs$p.y0 <- 1 - y_A1_probs$p.y1 - y_A1_probs$p.y2

  y_A1 <- apply(y_A1_probs, 1, function(i) {sample(c(1,2,0), 1, prob = i)})

  #Monotonicity -
  #if y(A=0) = 1 or y(A=0) = 2 we set y(A=1) = y(A=0), else we use the sampled Y(A=1)
  y_A1 <- ifelse(y_A0 == 0, y_A1, y_A0)

  y1_A1 <- ifelse(y_A1 == 1, 1, 0)
  y2_A1 <- ifelse(y_A1 == 2, 1, 0) #assert y1A1 >= y1A0. same for y2

  #create observed data, A as a function of X:
  p <- exp(A_par$gamma0 + A_par$gammaX%*%(t(X))) / (1 + exp(A_par$gamma0 + A_par$gammaX%*%(t(X))))
  p <- c(p)
  A <- rbinom(n.sample, 1, prob = p)

  y1 <- y1_A0*(1-A) + y1_A1*A
  y2 <- y2_A0*(1-A) + y2_A1*A
  y <- case_when(y1 == 1 ~ 1, y2 == 1 ~ 2, TRUE ~ 0)

  df <- data.frame(A, X, y, y1, y2, y1_A0, y1_A1, y2_A0, y2_A1, U)
  df$miss <- 1
  df$weight <- 1

  return(df)
}

#func to create n.sim samples of size n.sample
generate_samples <- function(n.sim = 1000, n.sample = 500,
                             y1_par = y1_par, y2_par = y2_par, A_par = A_par)
{
  dfs <- list()
  for (i in 1:n.sim)
  {
    dfs[[i]] <- generate_df_covariates(n.sample,
                                       y1_par = y1_par, y2_par = y2_par, A_par = A_par)
  }
  return(dfs)
}

#### Different estimands
#func to calculate sface1 or sface2 using standartization for a given df
stand <- function(df, sface)
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
iptw <- function(df, sface)
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
stand_total_naive <- function(df, sface)
{
  df <- df %>% filter(miss == 1)
  fit_y_by_A_X <- multinom(y ~ A + X1 + X2 ,df, trace = FALSE,
                           weights = weight)

  df_treat <- df_untr <- df
  df_treat$A <- 1
  df_untr$A <- 0

  pred_treat <- as.data.frame(predict(fit_y_by_A_X, newdata = df_treat, type = "probs"))
  colnames(pred_treat) <-c("0","y1", "y2")
  pred_untr <- as.data.frame(predict(fit_y_by_A_X, newdata = df_untr, type = "probs"))
  colnames(pred_untr) <-c("0","y1", "y2")

  self <- ifelse(sface == "sface1", "y1", "y2")
  other <- ifelse(sface == "sface1", "y2", "y1")

  n <- nrow(df)


  return(M*(mean(pred_treat[,self]) - mean(pred_untr[,self])))
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

#func to calculate sface1 or sface2 in terms of PO
po_sface_func <- function(df, sface)
{
  if(sface == "sface1")
  {
    df <- subset(df, y2_A0 == 0 & y2_A1 == 0)
    return(M*mean(df$y1_A1 - df$y1_A0))
  }

  if(sface == "sface2")
  {
    df <- subset(df, y1_A0 == 0 & y1_A1 == 0)
    return(M*mean(df$y2_A1 - df$y2_A0))
  }
}


#func to calculate sface1 or sface2 using the naive estimand that conditions on y2=0 using standartization
stand_cond_naive <- function(df, sface)
{
  if (sface == "sface1")
  {
    fit_y1_by_A_X_y2 <- glm(y1 ~ A + X1 + X2 + y2 ,df, family = "binomial")

    df_treat <- df_untr <- df
    df_treat$A <- 1
    df_treat$y2 <- 0
    df_untr$A <- 0
    df_untr$y2 <- 0

    pred_treat <- predict(fit_y1_by_A_X_y2, newdata = df_treat, type="response")
    pred_untr <- predict(fit_y1_by_A_X_y2, newdata = df_untr, type="response")

    #model y ~ A + X
    fit_y2_by_A_X <- glm(y2 ~ A + X1 + X2 ,df, family = "binomial")


    pred_treat_y2 <- predict(fit_y2_by_A_X, newdata = df_treat, type = "response")
    n_star <- sum(1-pred_treat_y2)
    return(M*sum(pred_treat - pred_untr)/n_star)
  }
  if (sface == "sface2")
  {
    fit_y2_by_A_X_y1 <- glm(y2 ~ A + X1 + X2 + y1 ,df, family = "binomial")

    df_treat <- df_untr <- df
    df_treat$A <- 1
    df_treat$y1 <- 0
    df_untr$A <- 0
    df_untr$y1 <- 0

    pred_treat <- predict(fit_y2_by_A_X_y1, newdata = df_treat, type="response")
    pred_untr <- predict(fit_y2_by_A_X_y1, newdata = df_untr, type="response")

    #model y ~ A + X
    fit_y1_by_A_X <- glm(y1 ~ A + X1 + X2 ,df, family = "binomial")

    pred_treat_y1 <- predict(fit_y1_by_A_X, newdata = df_treat, type = "response")
    n_star <- sum(1-pred_treat_y1)
    return(M*sum(pred_treat - pred_untr)/n_star)
  }
}


#func to calculate sface1-sface2 for a given sample / population and estimand
diff <- function(df, FUN)
{
  return(FUN(df, "sface1") - FUN(df, "sface2"))
}



#### BOOTSTRAP function
#func to calculate the variance of an estimetor in a specific sample, using BS
boot_func <- function(df, estimand_func_vec, sface = "sface1", BS = 200)
{
  boot_estimand_func <- function(df, index, estimand_func_vec, sface = "sface1")
  {
    df <- df[index,]
    return(estimand_func_vec(df, sface))
  }

  boot_estimators <- boot(data = df,
                          statistic = boot_estimand_func,
                          estimand_func_vec = estimand_func_vec,
                          sface = sface,
                          R = BS,
                          parallel = "multicore",
                          ncpus = 8)
  return(apply(boot_estimators$t, 2, sd))
}

#func to calculate the variance of an estimetor in a specific sample, using BS, for diff
boot_func_diff <- function(df, estimand_func_vec, sface = "sface1", BS = 200)
{
  boot_estimand_func <- function(df, index, diff, estimand_func_vec, sface = "sface1")
  {
    df <- df[index,]
    return(diff(df, estimand_func_vec))
  }

  boot_estimators <- boot(data = df,
                          statistic = boot_estimand_func,
                          diff = diff,
                          estimand_func = estimand_func_vec,
                          sface = sface,
                          R = BS,
                          parallel = "multicore",
                          ncpus = 8)
  return(apply(boot_estimators$t, 2, sd))
}


#### Estimating sface1 and sface2 in multiple dfs

#func to caculate sface1 or sface2 for each sample and compute stats
compute_stats_sface <- function(val_real, dfs, sface = "sface1", boot_func)
{
  estimand_func_vec <- function(df, sface)
  {
    return(c(stand(df, sface),
             iptw(df, sface),
             DR(df, sface),
             stand_total_naive(df, sface),
             stand_cond_naive(df, sface)))
  }
  estimands <- list("stand" = stand,
                    "iptw" = iptw,
                    "DR" = DR,
                    "stand_total_naive" = stand_total_naive,
                    "stand_cond_naive" = stand_cond_naive)

  boot_sample <- sapply(dfs, boot_func, estimand_func_vec, sface = sface)
  boot_sample <- t(boot_sample)
  colnames(boot_sample) <- names(estimands)

  #for each estimand
  estimand_compute_stats <- function(estimand_func, name_estimand, val_real, dfs, sface, boot_sample)
  {
    boot_sample <- boot_sample[,name_estimand]
    se_boot <- mean(boot_sample)
    q <- qnorm(0.975)
    sample <- sapply(dfs, estimand_func, sface = sface)
    mean_estimator <- mean(sample)
    bias <- mean_estimator - val_real
    sd_estimator <- sd(sample)
    cov_rate <- mean(sample - q*boot_sample <= val_real &
                       sample + q*boot_sample >= val_real)*100


    return(data.frame(mean_estimator = mean_estimator,
                      bias = bias,
                      sd_estimator = sd_estimator,
                      se_boot = se_boot,
                      cov_rate = cov_rate))
  }

  stats <- mapply(estimand_compute_stats,
                  estimands,
                  names(estimands),
                  MoreArgs = list(val_real, dfs, sface, boot_sample)) %>%
    as.data.frame() %>%
    rownames_to_column("value") %>%
    unnest(cols = c(stand, iptw, DR, stand_total_naive, stand_cond_naive)) %>%
    mutate_if(is.numeric, round, 6)
  return(stats)
}


compute_stats_diff <- function(val_real, dfs, boot_func_diff)
{
  estimand_func_vec <- function(df, sface)
  {
    return(c(stand(df, sface),
             iptw(df, sface),
             DR(df, sface),
             stand_total_naive(df, sface),
             stand_cond_naive(df, sface)))
  }
  estimands <- list("stand" = stand,
                    "iptw" = iptw,
                    "DR" = DR,
                    "stand_total_naive" = stand_total_naive,
                    "stand_cond_naive" = stand_cond_naive)

  boot_sample <- sapply(dfs, boot_func_diff, estimand_func_vec)
  boot_sample <- t(boot_sample)
  colnames(boot_sample) <- names(estimands)

  #for each estimand
  estimand_compute_stats <- function(estimand_func, name_estimand, val_real, dfs, boot_sample)
  {
    boot_sample <- boot_sample[,name_estimand]
    q <- qnorm(0.975)
    sample <- sapply(dfs, diff, estimand_func)
    mean_estimator <- mean(sample)
    bias <- mean_estimator - val_real
    sd_estimator <- sd(sample)
    se_boot <- mean(boot_sample)
    cov_rate <- mean(sample - q*boot_sample <= val_real &
                       sample + q*boot_sample >= val_real)*100

    return(data.frame(mean_estimator = mean_estimator,
                      bias = bias,
                      sd_estimator = sd_estimator,
                      se_boot = se_boot,
                      cov_rate = cov_rate))
  }


  estimands <- list("stand" = stand,
                    "iptw" = iptw,
                    "DR" = DR,
                    "stand_total_naive" = stand_total_naive,
                    "stand_cond_naive" = stand_cond_naive)


  stats <- mapply(estimand_compute_stats,
                  estimands,
                  names(estimands),
                  MoreArgs = list(val_real, dfs, boot_sample)) %>%
    as.data.frame() %>%
    rownames_to_column("value") %>%
    unnest(cols = c(stand, iptw, DR, stand_total_naive, stand_cond_naive)) %>%
    mutate_if(is.numeric, round, 6)

  return(stats)
}

#using all previos functions to create population, samples, and calculate all stats:
compute_stats_sface_sample_sizes <- function(sizes = c(10000),
                                             y1_par = y1_par, y2_par = y2_par, A_par = A_par)
{
  #create population and calculate "real" sface1 and sface2
  df_pop <- generate_df_covariates(n.sample = 1000000,
                                   y1_par = y1_par, y2_par = y2_par, A_par = A_par)

  sface1_real <- po_sface_func(df_pop, "sface1")
  sface2_real <- po_sface_func(df_pop, "sface2")
  diff_real <- diff(df_pop, po_sface_func)

  #loop: each time using a different n.sample, calculatng n.sim samples of this size
  #and calculate all stats (both for sface1 and for sface2)
  stats_sface1 <- list()
  stats_sface2 <- list()
  stats_diff <- list()

  for(i in 1:length(sizes))
  {
    samples <- generate_samples(n.sample = sizes[i],  y1_par = y1_par, y2_par = y2_par, A_par = A_par)
    stats_sface1[[i]] <- compute_stats_sface(sface1_real, samples, sface = "sface1", boot_func) %>%
      mutate(size = sizes[i])

    stats_sface2[[i]] <- compute_stats_sface(sface2_real, samples, sface = "sface2", boot_func) %>%
      mutate(size = sizes[i])
    stats_diff[[i]] <- compute_stats_diff(diff_real, samples, boot_func_diff) %>%
      mutate(size = sizes[i])

  }

  stats_sface1 <- bind_rows(stats_sface1) %>%
    select(value, size, everything()) %>%
    arrange(value, size)

  stats_sface2 <- bind_rows(stats_sface2) %>%
    select(value, size, everything()) %>%
    arrange(value, size)

  stats_diff <- bind_rows(stats_diff) %>%
    select(value, size, everything()) %>%
    arrange(value, size)


  sface <- list(stats_sface1, stats_sface2, stats_diff)
  sface[[4]] <- df_pop %>% summarise(prop_y1 = mean(y1),
                                     prop_y2 = mean(y2),
                                     prop_both =  mean(y1) + mean(y2),
                                     prop_A = mean(A)
  )

  names(sface) <- c(str_c("sface1 real value: ", as.character(round(sface1_real, 6))),
                    str_c("sface2 real value: ", as.character(round(sface2_real, 6))),
                    str_c("sface1-sface2 real value: ", as.character(round(diff_real, 6))),
                    "Y1 Y2 A pervelance"
  )
  return(sface)
}

ans <- compute_stats_sface_sample_sizes(sizes = c(10000),
                                               y1_par = list(gamma0 = log(0.05),
                                                             gammaA = log(2),
                                                             gammaU = log(2),
                                                             gammaX = c(log(0.25),
                                                                        log(2))),
                                               y2_par = list(gamma0 = log(0.005),
                                                             gammaA = log(2),
                                                             gammaU = log(args),
                                                             gammaX = c(log(2),
                                                                        log(2))),
                                               A_par = list(gamma0 = log(0.7),
                                                            gammaX = c(log(2),
                                                                       log(2))) )


save(ans, file = str_c("Upar_log", as.character(args),".RData"))

