library(tidyverse)
library(patchwork)
library(latex2exp)
library(xtable)

Upar_2 <- get(load("Upar_log2.Rdata"))
Upar_3 <- get(load("Upar_log3.Rdata"))
Upar_4 <- get(load("Upar_log4.Rdata"))
Upar_5 <- get(load("Upar_log5.Rdata"))
Upar_6 <- get(load("Upar_log6.Rdata"))

Upar <- list(Upar_2,Upar_3, Upar_4,Upar_5, Upar_6)

real_sface1 <- tibble(Upar_param = 2:6, real = sapply(1:5,function(i) {as.numeric(str_sub(names(Upar[[i]][1]), 20))}))
real_sface2 <- tibble(Upar_param = 2:6, real = sapply(1:5,function(i) {as.numeric(str_sub(names(Upar[[i]][2]), 20))}))
real_diff <- tibble(Upar_param = 2:6, real = sapply(1:5,function(i) {as.numeric(str_sub(names(Upar[[i]][3]), 27))}))


Upar_ssace1 <-bind_rows(
  Upar_2[[1]] %>% mutate(Upar_param = 2),
  Upar_3[[1]] %>% mutate(Upar_param = 3),
  Upar_4[[1]] %>% mutate(Upar_param = 4),
  Upar_5[[1]] %>% mutate(Upar_param = 5),
  Upar_6[[1]] %>% mutate(Upar_param = 6)
)
Upar_ssace1 <- Upar_ssace1 %>%
  arrange(Upar_param) %>%
  select(Upar_param, everything(), -size)%>%
  pivot_longer(cols = -c(Upar_param, value), names_to = "estimand", values_to = "res")%>%
  pivot_wider(names_from = value, values_from = res) %>%
  left_join(real_sface1) %>%
  mutate(relative_bias = 100*bias/real) %>%
  select(Upar_param, estimand, relative_bias, sd_estimator, se_boot) %>%
  mutate(estimand = case_when(estimand == "stand" ~ "Stand",
                              estimand == "iptw" ~ "IPTW",
                              estimand == "DR" ~ "DR",
                              estimand == "stand_total_naive"~ "Total",
                              estimand == "stand_cond_naive"~ "Cond"))

Upar_ssace1$relative_bias_jitter <- jitter(Upar_ssace1$relative_bias,100)
sface1_relative_differance <- ggplot(Upar_ssace1, aes(x = Upar_param, y = relative_bias, color = estimand)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_line() +
  geom_hline(yintercept=0, linetype = "dashed") +
  scale_color_brewer(palette = "Set1", name = "Estimator") +
 xlab(TeX("$e^{\\delta_2}$")) +
  ylab("%Bias") +
  ylim(c(-10,20))+
  theme_bw() +
  ggtitle(TeX("$SF-ACE^{(1)}_{D}$"))

Upar_ssace2 <-bind_rows(
  Upar_2[[2]] %>% mutate(Upar_param = 2),
  Upar_3[[2]] %>% mutate(Upar_param = 3),
  Upar_4[[2]] %>% mutate(Upar_param = 4),
  Upar_5[[2]] %>% mutate(Upar_param = 5),
  Upar_6[[2]] %>% mutate(Upar_param = 6)
)  %>%
  arrange(Upar_param) %>%
  select(Upar_param, everything(), -size)%>%
  pivot_longer(cols = -c(Upar_param, value), names_to = "estimand", values_to = "res")%>%
  pivot_wider(names_from = value, values_from = res) %>%
  left_join(real_sface2) %>%
  mutate(relative_bias = 100*bias/real) %>%
  select(Upar_param, estimand, relative_bias, sd_estimator, se_boot) %>%
  mutate(estimand = case_when(estimand == "stand" ~ "Stand",
                              estimand == "iptw" ~ "IPTW",
                              estimand == "DR" ~ "DR",
                              estimand == "stand_total_naive"~ "Total",
                              estimand == "stand_cond_naive"~ "Cond"))


Upar_ssace2$relative_bias_jitter <- jitter(Upar_ssace2$relative_bias,100)
sface2_relative_differance <- ggplot(Upar_ssace2, aes(x = Upar_param, y = relative_bias, color = estimand)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_line() +
  geom_hline(yintercept=0, linetype = "dashed") +
  scale_color_brewer(palette = "Set1", name = "Estimator") +
 xlab(TeX("$e^{\\delta_2}$")) +
  ylab("%Bias") +
  ylim(c(-10,20))+
  theme_bw() +
  ggtitle(TeX("$SF-ACE^{(2)}_{D}$"))



#handle the diff
Upar_diff <-bind_rows(
  Upar_2[[3]] %>% mutate(Upar_param = 2),
  Upar_3[[3]] %>% mutate(Upar_param = 3),
  Upar_4[[3]] %>% mutate(Upar_param = 4),
  Upar_5[[3]] %>% mutate(Upar_param = 5),
  Upar_6[[3]] %>% mutate(Upar_param = 6)
) %>%
  arrange(Upar_param) %>%
  select(Upar_param, everything(), -size)%>%
  pivot_longer(cols = -c(Upar_param, value), names_to = "estimand", values_to = "res")%>%
  pivot_wider(names_from = value, values_from = res) %>%
  left_join(real_diff) %>%
  mutate(relative_bias = 100*bias/real) %>%
  select(Upar_param, estimand, relative_bias, bias, sd_estimator)  %>%
  mutate(estimand = case_when(estimand == "stand" ~ "Stand",
                              estimand == "iptw" ~ "IPTW",
                              estimand == "DR" ~ "DR",
                              estimand == "stand_total_naive"~ "Total",
                              estimand == "stand_cond_naive"~ "Cond"))

Upar_diff$relative_bias_jitter <- jitter(Upar_diff$relative_bias,100)
diff_relative_differance <- ggplot(Upar_diff, aes(x = Upar_param, y = relative_bias, color = estimand)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_line() +
  geom_hline(yintercept=0, linetype = "dashed") +
  scale_color_brewer(palette = "Set1", name = "Estimator") +
 xlab(TeX("$e^{\\delta_2}$")) +
  ylab("%Bias") +
  ylim(c(-10,20))+
  theme_bw() +
  ggtitle(TeX("$\\theta_{D}$"))






Upar_2 <- get(load("Upar_RR_log2.Rdata"))
Upar_3 <- get(load("Upar_RR_log3.Rdata"))
Upar_4 <- get(load("Upar_RR_log4.Rdata"))
Upar_5 <- get(load("Upar_RR_log5.Rdata"))
Upar_6 <- get(load("Upar_RR_log6.Rdata"))
Upar <- list(Upar_2,Upar_3, Upar_4,Upar_5, Upar_6)

real_sface1 <- tibble(Upar_param = 2:6, real = sapply(1:5,function(i) {as.numeric(str_sub(names(Upar[[i]][1]), 20))}))
real_sface2 <- tibble(Upar_param = 2:6, real = sapply(1:5,function(i) {as.numeric(str_sub(names(Upar[[i]][2]), 20))}))
real_diff <- tibble(Upar_param = 2:6, real = sapply(1:5,function(i) {as.numeric(str_sub(names(Upar[[i]][3]), 27))}))


Upar_ssace1 <-bind_rows(
  Upar_2[[1]] %>% mutate(Upar_param = 2),
  Upar_3[[1]] %>% mutate(Upar_param = 3),
  Upar_4[[1]] %>% mutate(Upar_param = 4),
  Upar_5[[1]] %>% mutate(Upar_param = 5),
  Upar_6[[1]] %>% mutate(Upar_param = 6)
)
Upar_ssace1 <- Upar_ssace1 %>%
  arrange(Upar_param) %>%
  select(Upar_param, everything(), -size)%>%
  pivot_longer(cols = -c(Upar_param, value), names_to = "estimand", values_to = "res")%>%
  pivot_wider(names_from = value, values_from = res) %>%
  left_join(real_sface1) %>%
  mutate(relative_bias = 100*bias/real) %>%
  select(Upar_param, estimand, relative_bias, sd_estimator, se_boot) %>%
  mutate(estimand = case_when(estimand == "stand" ~ "Stand",
                              estimand == "iptw" ~ "IPTW",
                              estimand == "DR" ~ "DR",
                              estimand == "stand_total_naive"~ "Total",
                              estimand == "stand_cond_naive"~ "Cond"))

Upar_ssace1$relative_bias_jitter <- jitter(Upar_ssace1$relative_bias,20)
sface1_relative <- ggplot(Upar_ssace1, aes(x = Upar_param, y = relative_bias, color = estimand)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_line() +
  geom_hline(yintercept=0, linetype = "dashed") +
  scale_color_brewer(palette = "Set1", name = "Estimator") +
 xlab(TeX("$e^{\\delta_2}$")) +
  ylab("%Bias") +
  ylim(c(-1,10))+
  theme_bw() +
  ggtitle(TeX("$SF-ACE^{(1)}_{RR}$"))

Upar_ssace2 <-bind_rows(
  Upar_2[[2]] %>% mutate(Upar_param = 2),
  Upar_3[[2]] %>% mutate(Upar_param = 3),
  Upar_4[[2]] %>% mutate(Upar_param = 4),
  Upar_5[[2]] %>% mutate(Upar_param = 5),
  Upar_6[[2]] %>% mutate(Upar_param = 6)
)  %>%
  arrange(Upar_param) %>%
  select(Upar_param, everything(), -size)%>%
  pivot_longer(cols = -c(Upar_param, value), names_to = "estimand", values_to = "res")%>%
  pivot_wider(names_from = value, values_from = res) %>%
  left_join(real_sface2) %>%
  mutate(relative_bias = 100*bias/real) %>%
  select(Upar_param, estimand, relative_bias, sd_estimator, se_boot) %>%
  mutate(estimand = case_when(estimand == "stand" ~ "Stand",
                              estimand == "iptw" ~ "IPTW",
                              estimand == "DR" ~ "DR",
                              estimand == "stand_total_naive"~ "Total",
                              estimand == "stand_cond_naive"~ "Cond"))


Upar_ssace2$relative_bias_jitter <- jitter(Upar_ssace2$relative_bias,150)
sface2_relative <- ggplot(Upar_ssace2, aes(x = Upar_param, y = relative_bias, color = estimand)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_line() +
  geom_hline(yintercept=0, linetype = "dashed") +
  scale_color_brewer(palette = "Set1", name = "Estimator") +
 xlab(TeX("$e^{\\delta_2}$")) +
  ylab("%Bias") +
  ylim(c(-1,10))+
  theme_bw() +
  ggtitle(TeX("$SF-ACE^{(2)}_{RR}$"))


#handle the diff
Upar_diff <-bind_rows(
  Upar_2[[3]] %>% mutate(Upar_param = 2),
  Upar_3[[3]] %>% mutate(Upar_param = 3),
  Upar_4[[3]] %>% mutate(Upar_param = 4),
  Upar_5[[3]] %>% mutate(Upar_param = 5),
  Upar_6[[3]] %>% mutate(Upar_param = 6)
) %>%
  arrange(Upar_param) %>%
  select(Upar_param, everything(), -size)%>%
  pivot_longer(cols = -c(Upar_param, value), names_to = "estimand", values_to = "res")%>%
  pivot_wider(names_from = value, values_from = res) %>%
  left_join(real_diff) %>%
  mutate(relative_bias = 100*bias/real) %>%
  select(Upar_param, estimand, relative_bias, bias, sd_estimator)  %>%
  mutate(estimand = case_when(estimand == "stand" ~ "Stand",
                              estimand == "iptw" ~ "IPTW",
                              estimand == "DR" ~ "DR",
                              estimand == "stand_total_naive"~ "Total",
                              estimand == "stand_cond_naive"~ "Cond"))

Upar_diff$relative_bias_jitter <- jitter(Upar_diff$relative_bias,100)
diff_relative <- ggplot(Upar_diff, aes(x = Upar_param, y = relative_bias, color = estimand)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_line() +
  geom_hline(yintercept=0, linetype = "dashed") +
  scale_color_brewer(palette = "Set1", name = "Estimator") +
 xlab(TeX("$e^{\\delta_2}$")) +
  ylab("%Bias") +
  theme_bw() +
  ggtitle(TeX("$\\theta_{RR}$"))


Upar_RR <- sface1_relative + sface2_relative + diff_relative + plot_layout(guides = 'collect')



(sface1_relative_differance + sface2_relative_differance + diff_relative_differance) /
(sface1_relative + sface2_relative + diff_relative) +
  plot_layout(guides = 'collect') & theme(legend.position = 'bottom',
                                         title = element_text(size = 14),
                                         legend.title = element_text(size = 16),
                                         legend.text = element_text(size = 12),
                                         axis.title.x = element_text(size = 16),
                                         axis.text.x = element_text(size = 10),
                                         axis.title.y = element_text(size = 14))
