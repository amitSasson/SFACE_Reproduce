library(tidyverse)
library(patchwork)
library(latex2exp)


sample_size_5000 <- get(load("sample_sizes_RR_5000.RData"))
sample_size_10000 <- get(load("sample_sizes_RR_10000.RData"))
sample_size_25000 <- get(load("sample_sizes_RR_25000.RData"))
sample_size_50000 <- get(load("sample_sizes_RR_50000.RData"))

real_sface1 <- as.numeric(str_sub(names(sample_size_5000[1]), 20))
real_sface2 <- as.numeric(str_sub(names(sample_size_5000[2]), 20))
real_diff <- as.numeric(str_sub(names(sample_size_5000[3]), 27))

sample_size <- list(
  ssace1 = bind_rows(
    sample_size_5000[[1]],
    sample_size_10000[[1]],
    sample_size_25000[[1]],
    sample_size_50000[[1]]
  ),
  ssace2 = bind_rows(
    sample_size_5000[[2]],
    sample_size_10000[[2]],
    sample_size_25000[[2]],
    sample_size_50000[[2]]
  )
)
#funs to create plots for n sample
table_to_plot_sample_sizes <- function(df)
{
  ggplot(df, aes(
    x = as.factor(size),
    y = relative_bias,
    color = estimand
  )) +
    geom_point(size = 3, alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_brewer(palette = "Set1", name = "Estimand") +
    labs(x = "N", y = "%Bias") +
    theme_bw()

}




sample_size_plot(sample_size)


df1 <- sample_size[[1]] %>%
                        arrange(size) %>%
                        pivot_longer(cols = -c(size, value), names_to = "estimand", values_to = "res") %>%
                        pivot_wider(id_cols = c(size, estimand),names_from = value, values_from = res) %>%
                        mutate(relative_bias = 100*bias/ real_sface1) %>%
                        mutate(estimand = case_when(estimand == "stand" ~ "Stand",
                                                         estimand == "iptw" ~ "IPTW",
                                                         estimand == "DR" ~ "DR",
                                                         estimand == "stand_total_naive"~ "Total",
                                                         estimand == "stand_cond_naive"~ "Cond"))




df2 <- sample_size[[2]] %>%
  arrange(size) %>%
  pivot_longer(cols = -c(size, value), names_to = "estimand", values_to = "res")%>%
  pivot_wider(id_cols = c(size, estimand),names_from = value, values_from = res) %>%
  mutate(relative_bias = 100*bias/real_sface2) %>%
  mutate(estimand = case_when(estimand == "stand" ~ "Stand",
                              estimand == "iptw" ~ "IPTW",
                              estimand == "DR" ~ "DR",
                              estimand == "stand_total_naive"~ "Total",
                              estimand == "stand_cond_naive"~ "Cond"))



diff <- bind_rows(
  sample_size_5000[[3]],
  sample_size_10000[[3]],
  sample_size_25000[[3]],
  sample_size_50000[[3]]) %>%
  arrange(size) %>%
  pivot_longer(cols = -c(size, value), names_to = "estimand", values_to = "res")%>%
  pivot_wider(id_cols = c(size, estimand),names_from = value, values_from = res) %>%
  mutate(relative_bias = 100*bias/real_diff) %>%
  mutate(estimand = case_when(estimand == "stand" ~ "Stand",
                              estimand == "iptw" ~ "IPTW",
                              estimand == "DR" ~ "DR",
                              estimand == "stand_total_naive"~ "Total",
                              estimand == "stand_cond_naive"~ "Cond"))



plot1 <- table_to_plot_sample_sizes(df1) +
  ggtitle(expression(SF - ACE ^(1)))
plot2 <- table_to_plot_sample_sizes(df2) +
  ggtitle(expression(SF - ACE ^ (2)))

diff <- diff %>% select(size, estimand, relative_bias, sd_estimator)
print(xtable(diff, digits = 3), include.rownames = FALSE)


plot3 <- ggplot(diff, aes(
  x = as.factor(size),
  y = relative_bias,
  color = estimand)) +
  geom_point(size = 3, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_brewer(palette = "Set1", name = "Estimand") +
  labs(x = "N", y = "%Bias") +
  theme_bw() +
  ggtitle(TeX("$\\theta_d$"))

plot1 + plot2 + plot3 + plot_layout(guides = 'collect')




## Table
df1 <- df1 %>% dplyr::select(size, estimand, bias, relative_bias,cov_rate, sd_estimator,se_boot)
print(xtable(df1, digits = 2), include.rownames = FALSE)

df2 <- df2 %>%  dplyr::select(size, estimand, bias, relative_bias,cov_rate, sd_estimator,se_boot)
print(xtable(df2, digits = 2), include.rownames = FALSE)

df3 <- diff %>%  dplyr::select(size, estimand, bias, relative_bias,cov_rate, sd_estimator,se_boot)
print(xtable(df3, digits = 2), include.rownames = FALSE)


