library(tidyverse)
library(patchwork)
library(nnet)
library(boot)
library(latex2exp)
library(TheSFACE)

#### setup
set.seed(3141)
q <- qnorm(0.975)
M <- 100000
effect_scale <- "RR"

setwd("/Users/amitsasson/git/SFACE_Reproduce/Data")
df <- read_csv("pseudo_dataset_for_analysis_SFACE.csv")

df$gender <- as.factor(df$gender)
df$regasp <- as.factor(df$regasp)
df$crcfh <- as.factor(df$crcfh)
df$endo <- as.factor(df$endo)

stand_formula <- y ~ A + gender + bmicum_cat + regasp + alcocum + caloriecum + crcfh + endo + actmcum
iptw_formula <- A ~ gender + bmicum_cat + regasp + alcocum + caloriecum + crcfh + endo + actmcum

#### sensitivity analysis #####

max_lambda <- max_lambda(stand_formula = stand_formula,
                       exposure = "A",
                       df,
                       weight = "weight")

lambda2_vals <- seq(0, max_lambda[2], by=0.02)
lambda1_vals <- seq(0, max_lambda[1], by=0.02)


val <- sface(stand_formula = stand_formula,
             iptw_formula = iptw_formula,
             exposure = "A",
             outcome = "y",
             df = df,
             scale = effect_scale,
             weight = "weight",
             MultPer = M,
             lambda1 = lambda1_vals,
             lambda2 = lambda2_vals)

print(val)

plot(val)
