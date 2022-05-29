library(tidyverse)
df <- read_csv("/udd/nhasn/data_after_pre_process_for_weights_60_70.csv") 

df %>% 
  select(gender, packyr, bmicum, regasp, alcocum, caloriecum, actmcum, endo, crcfh) %>% 
  summarise_all(funs(mean(is.na(.))))

### Handle missing values
df <- df %>% replace_na(list(bmicum = mean(df$bmicum, na.rm = T), 
                             actmcum = mean(df$actmcum, na.rm = T))) 

## turn to cat
df  <- df %>% 
  mutate(grade = ifelse(grade == 5, NA, grade)) %>% 
  mutate(dtdx_cat = case_when(dtdx/12 + 1900 < 1996 ~ "before_1996", 
                              dtdx/12 + 1900 > 2000 ~ "after_2000", 
                              TRUE ~ "1996_2000"), 
         bmicum_cat = case_when(bmicum < 25 ~ "low", 
                                bmicum >=30 ~ "high", 
                                TRUE ~ "normal"),
         gender = factor(gender), 
         crcfh = factor(crcfh))

## create y 
df <- df %>% mutate(y = case_when(is.na(msin) ~ 0, 
                                  msin == 0 ~ 1, 
                                  msin == 1 ~ 2))

## change factor labels 
df <- df %>% mutate(subsite = case_when(subsite == 1 ~ "prox",
                                        subsite == 2 ~ "distal",
                                        subsite == 3 ~ "rectum", 
                                        is.na(subsite) & crc == 1 ~ "missing"), 
                    differentiation = case_when(grade == 1 ~ "well",
                                                grade == 2 ~ "mderate",
                                                grade == 3 ~ "poor",
                                                grade == 4 ~ "unspecified", 
                                                is.na(grade) & crc == 1 ~ "missing")) %>% 
  mutate(subsite = factor(subsite), 
         stageall = factor(stageall), 
         differentiation = factor(differentiation), 
         dtdx_cat = factor(dtdx_cat), 
         y = factor(y), 
         bmicum_cat = factor(bmicum_cat)) 


#df <- df %>% filter(!is.na(packyr))


### Weights for missing msin value 
df_diagnosed <- df %>% filter(crc != 0)
df_diagnosed %>% 
  select(gender, packyr, bmicum_cat, regasp, alcocum, caloriecum, actmcum, endo, crcfh, subsite , stageall , agedx, differentiation) %>% 
  summarise_all(funs(mean(is.na(.))))

df_diagnosed <- df_diagnosed %>% mutate(msin_avaliable = !is.na(msin))


fit_msin_avaliable <- glm(msin_avaliable ~ gender + subsite + stageall + differentiation + agedx + dtdx_cat + crcfh, 
                          df_diagnosed, 
                          family = "binomial")

summary(fit_msin_avaliable)

# pass in coef directly
stargazer(fit_msin_avaliable, ci = TRUE, t.auto = FALSE, single.row = TRUE, digits = 2, header = F)


pred_msin_avaliable <-  predict(fit_msin_avaliable, type = "response")
df_diagnosed$pred_msin_avaliable <- pred_msin_avaliable

#library(pROC)
#g <- roc(msin_avaliable ~ pred_msin_avaliable, data = df_diagnosed)
#plot(g)    
#auc(df_diagnosed$msin_avaliable, df_diagnosed$pred_msin_avaliable)

df <- df %>% 
  left_join(df_diagnosed %>% select(id, pred_msin_avaliable)) %>% 
  mutate(weight_msin_avaliable = ifelse(is.na(pred_msin_avaliable), 1, 1/pred_msin_avaliable)) %>% 
  filter(!(crc == 1 & is.na(msin))) #filter after weights 

df %>% filter(weight_msin_avaliable != 1) %>% pull(weight_msin_avaliable) %>% hist()
df %>% filter(weight_msin_avaliable != 1) %>% pull(weight_msin_avaliable) %>% summary()
df %>%  pull(weight_outcome_missing) %>% summary()

cap_weight_msin_avaliable <- quantile(df %>% filter(weight_msin_avaliable != 1) %>% pull(weight_msin_avaliable), .99)

df <- df %>% mutate(weight_msin_avaliable = case_when(weight_msin_avaliable > cap_weight_msin_avaliable ~ cap_weight_msin_avaliable,
                                                      TRUE ~ weight_msin_avaliable))

df <- df %>% 
  mutate(A = ifelse(packyr > 0, 1, 0)) %>% 
  mutate(y1 = ifelse(y == 1, 1, 0), 
         y2 = ifelse(y == 2, 1, 0)) 

df <- df %>% 
  mutate(weight = weight_msin_avaliable)

sum_weights <- sum(df$weight)
df <- df %>% mutate(relative_weight=weight/sum_weights)

#write_csv(df, "data_after_weights_for_analysis_keep_all_60_70.csv")

