library(xtable)

clean_ans <- function(ans)
{
  sface <- str_sub(names(ans)[1], 19, 30) %>% as.numeric()
  a <- ans[[1]] %>%
    select(-size) %>%
    rename(total_effect = stand_total_naive) %>%
    pivot_longer(cols=c(stand , iptw ,total_effect ,DR), names_to = "estimand", values_to = "v") %>%
    pivot_wider(names_from = value, values_from = v) %>%
    mutate(relative_bias = bias/sface) %>%
    mutate(real_val = sface) %>%
    select(estimand, mean_estimator,real_val, bias, relative_bias, cov_rate, sd_estimator, se_boot)
  return(a)
}

a <- list()

load("sim_no_break.RData")
a[[1]]<- clean_ans(ans) %>% mutate(scenario = "no_break")
load("sim_breakA.RData")
a[[2]] <- clean_ans(ans) %>% mutate(scenario = "break_A")
load("sim_breakY.RData")
a[[3]] <- clean_ans(ans) %>% mutate(scenario = "break_y")
load("sim_break_both.RData")
a[[4]] <- clean_ans(ans) %>% mutate(scenario = "break_both")

a <- bind_rows(a)

a <- a %>%
  select(estimand, mean_estimator, bias, relative_bias,cov_rate, sd_estimator,se_boot)%>%
  mutate(relative_bias = relative_bias*100)
print(xtable(a, digits = 2), include.rownames = FALSE)



