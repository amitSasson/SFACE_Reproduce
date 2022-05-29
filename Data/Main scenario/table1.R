library(furniture)
library(stargazer)
#### setup 
set.seed(31415)
q <- qnorm(0.975)
M <- 100000
df <- read_csv("/udd/nhasn/data_after_weights_for_analysis_keep_all_60_70.csv") 
df$gender <- as.factor(df$gender)
df$regasp <- as.factor(df$regasp)
df$crcfh <- as.factor(df$crcfh)
df$endo <- as.factor(df$endo)
df$A <- as.factor(df$A)

df <- df %>% mutate(y = case_when(y == 0 ~ "CRC free",
                                  y == 1 ~ "Non-MSI-high",
                                  y == 2 ~ "MSI-high"),
                    gender = case_when(gender == 0 ~ "Women (NHS)",
                                       gender == 1 ~ "Men (HPFS)"),
                    regasp = case_when(regasp == 0 ~ "No",
                                       regasp == 1 ~ "Yes"),                   
                    A = case_when(A == 0 ~ "No",
                                  A == 1 ~ "Yes"),    
                    crcfh = case_when(crcfh == 0 ~ "No",
                                      crcfh == 1 ~ "Yes"),   
                    endo = case_when(endo == 0 ~ "No",
                                     endo== 1 ~ "Yes"))   
df <- df %>% rename(Gender = gender,
                    Smoking = A)

table1(df, Smoking, 
              Gender, 
              bmicum_cat , 
              regasp, 
              crcfh,
              endo ,
              alcocum,
              caloriecum,
              actmcum, splitby = ~y,
            total = TRUE, 
            output = "latex")



fit_y_by_A_X <- multinom(y ~ A + gender + bmicum_cat + regasp + alcocum + caloriecum + crcfh + endo + actmcum,
                         df, 
                         trace = FALSE, 
                         weights = weight)
OR <- exp(coef(fit_y_by_A_X))
stargazer(fit_y_by_A_X, coef = list(OR), ci = TRUE, t.auto = FALSE, single.row = TRUE, digits = 2, header = F)
        
fit_A_by_X <- glm(A ~ gender + bmicum_cat + regasp + alcocum + caloriecum + crcfh + endo + actmcum,
                  df, 
                  family = "binomial", 
                  weights = weight)
stargazer(fit_A_by_X, ci = TRUE, t.auto = FALSE)
OR <- exp(coef(fit_A_by_X))
# pass in coef directly
stargazer(fit_A_by_X, coef = list(OR), ci = TRUE, t.auto = FALSE, single.row = TRUE, digits = 2, header = F)


ggplot(df) + aes(weight_msin_avaliable) + geom_histogram() + scale_x_continuous(breaks = seq(1, 10, by=1)) + xlab("Weight") + theme_bw()
