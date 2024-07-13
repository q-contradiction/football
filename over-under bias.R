# Betting on over and under and perform hypothesis testing to check whether
# the returns from under are higher than over

# data from football-data.co.uk
suppressPackageStartupMessages({
  library(broom)
  library(tidyverse)
})

rm(list = ls())

football_odds <- readRDS(file = "./data-raw/football-data.rds")

# Market/avg,max & Betbrain/avg,max are the same variables
# Substitute where avg/max market odds are NA with Betbrain odds
cols <- list(c("Max>2.5", "BbMx>2.5"), c("Max<2.5", "BbMx<2.5"),
             c("Avg>2.5", "BbAv>2.5"), c("Avg<2.5", "BbAv<2.5"))

for (v in cols)  {
  replace_indexes <- is.na(football_odds[[v[1]]])
  football_odds[[v[1]]][replace_indexes] <- football_odds[[v[2]]][replace_indexes]
}

# Remove BetBrain odds
football_odds <- football_odds %>% .[setdiff(names(.), unlist(map(cols, 2)))]

odds_ou <- list(Bet365 = c("B365>2.5", "B365<2.5"), Bet365C = c("B365C>2.5", "B365C<2.5"), 
                Pinnacle = c("P>2.5", "P<2.5"), PinnacleC = c("PC>2.5", "PC<2.5"),
                MarketMax = c("Max>2.5", "Max<2.5"), MarketMaxC = c("MaxC>2.5", "MaxC<2.5"),
                MarketAvg = c("Avg>2.5", "Avg<2.5"), MarketAvgC = c("AvgC>2.5", "AvgC<2.5"))

football_odds <- football_odds %>% 
  mutate_at(vars(unlist(odds_ou)), as.numeric) %>%
  mutate(FTL = ifelse(FTHG + FTAG > 2, ">2.5", "<2.5"))

bet_ou <- function(bookmaker, nbets, nsims = 1000) {
  
  params <- as.list(environment())
  
  over_under <- football_odds %>%
    select(all_of(odds_ou[[bookmaker]]), FTL) %>%
    setNames(c("over", "under", "FTL")) %>%
    drop_na() %>% 
    pivot_longer(cols = c(over, under), names_to = "odds_type", values_to = "odds_value") %>%
    mutate(pnl = case_when(
      odds_type == "over"  & FTL == ">2.5" ~ odds_value-1,
      odds_type == "over"  & FTL == "<2.5" ~ -1,
      odds_type == "under" & FTL == ">2.5" ~ -1,
      odds_type == "under" & FTL == "<2.5" ~ odds_value-1,
      TRUE ~ NA
    ))
  
  sampling <- function(df) {
    df %>% 
      group_by(odds_type) %>%
      slice_sample(n = nbets, replace = T) %>%
      summarise(odds = mean(odds_value),
                yield = sum(pnl)/n()) 
    }
  
  # sampling nsims times
  sim_res <- replicate(nsims, sampling(over_under), simplify = F) %>% 
    bind_rows() %>% bind_cols(params) 
}

sim_results <- names(odds_ou) %>% map_dfr(bet_ou, nbets = 5000) 

test <- sim_results %>%
  group_by(bookmaker) %>%
  summarise(tidy(
    t.test(
      yield[odds_type == "under"],
      yield[odds_type == "over"], 
      alternative = "greater"
    )
  ), .groups = 'drop')

# Market Avg & Max display over-under bias but not in their closing lines
# Pinnacle, Bet365 do not.
print(test)

