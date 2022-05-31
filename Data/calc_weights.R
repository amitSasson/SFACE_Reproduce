library(tidyverse)
setwd("/Users/amitsasson/git/SFACE_Reproduce/Data")
df <- read_csv("pseudo_dataset_for_weights_SFACE.csv")

### Weights for missing msin value
df_diagnosed <- df %>% filter(y != 0)

fit_msin_avaliable <- glm(msin_avaliable ~ gender + subsite + stageall + differentiation + agedx + dtdx_cat + crcfh,
                          df_diagnosed,
                          family = "binomial")

summary(fit_msin_avaliable)

# create table of coef
#stargazer(fit_msin_avaliable, ci = TRUE, t.auto = FALSE, single.row = TRUE, digits = 2, header = F)


pred_msin_avaliable <-  predict(fit_msin_avaliable, type = "response")
df_diagnosed$pred_msin_avaliable <- pred_msin_avaliable

#library(pROC)
#g <- roc(msin_avaliable ~ pred_msin_avaliable, data = df_diagnosed)
#plot(g)
#auc(df_diagnosed$msin_avaliable, df_diagnosed$pred_msin_avaliable)

df <- df %>%
  left_join(df_diagnosed %>% select(id, pred_msin_avaliable)) %>%
  mutate(weight_msin_avaliable = ifelse(is.na(pred_msin_avaliable), 1, 1/pred_msin_avaliable)) %>%
  filter(msin_avaliable == 1 | is.na(msin_avaliable)) #filter after weights

cap_weight_msin_avaliable <- quantile(df %>% filter(weight_msin_avaliable != 1) %>% pull(weight_msin_avaliable), .99)

df <- df %>% mutate(weight_msin_avaliable = case_when(weight_msin_avaliable > cap_weight_msin_avaliable ~ cap_weight_msin_avaliable,
                                                      TRUE ~ weight_msin_avaliable))

df <- df %>%
  mutate(y1 = ifelse(y == 1, 1, 0),
         y2 = ifelse(y == 2, 1, 0))

df <- df %>%
  mutate(weight = weight_msin_avaliable)

#write_csv(df, "pseudo_dataset_for_analysis_SFACE.csv")
