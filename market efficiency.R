# Odds Comparison between different bookmakers in Football, identifying inefficiencies in the market

# Reproducing results & methodology from: football-data.co.uk/blog/wisdom_of_the_crowd.php (Joseph Buchdahl) in R
# and extending the analysis

# data from football-data.co.uk

suppressPackageStartupMessages({
  library(viridis) #graphics
  library(scales)
  library(implied) #https://opisthokonta.net/?p=1797
  library(tidyverse)
})

rm(list = ls())
options(scipen = 999)
set.seed(1)

# 1 - Data ----------------------------------------------------------------

football_odds <- readRDS(file = "./data-raw/football-data.rds")

# Market/avg,max & Betbrain/avg,max are the same variables
# Substitute where avg/max market odds are NA with Betbrain odds
cols <- list(c("MaxH", "BbMxH"), c("MaxD", "BbMxD"), c("MaxA", "BbMxA"),
             c("AvgH", "BbAvH"), c("AvgD", "BbAvD"), c("AvgA", "BbAvA"),
             c("Max>2.5", "BbMx>2.5"), c("Max<2.5", "BbMx<2.5"),
             c("Avg>2.5", "BbAv>2.5"), c("Avg<2.5", "BbAv<2.5"))

for (v in cols)  {
  replace_indexes <- is.na(football_odds[[v[1]]])
  football_odds[[v[1]]][replace_indexes] <- football_odds[[v[2]]][replace_indexes]
}

# remove BetBrain odds
football_odds <- football_odds %>% .[setdiff(names(.), unlist(map(cols, 2)))]

# odds abbreviations for 1x2 in Football-data.co.uk. C at the end of bookmaker's name indicates closing odds
odds_1x2 <- list(Bet365 = c("B365H", "B365D", "B365A"), Bet365C = c("B365CH", "B365CD", "B365CA"),
                 Bwin = c("BWH", "BWD", "BWA"), BwinC = c("BWCH", "BWCD", "BWCA"), 
                 Interwetten = c("IWH", "IWD", "IWA"), InterwettenC = c("IWCH", "IWCD", "IWCA"),
                 Pinnacle = c("PSH", "PSD", "PSA"), PinnacleC = c("PSCH", "PSCD", "PSCA"),
                 WilliamHill = c("WHH", "WHD", "WHA"), WilliamHillC = c("WHCH", "WHCD", "WHCA"),
                 VC = c("VCH", "VCD", "VCA"), VCC = c("VCCH", "VCCD", "VCCA"),
                 MarketMax = c("MaxH", "MaxD", "MaxA"), MarketMaxC = c("MaxCH", "MaxCD", "MaxCA"),
                 MarketAvg = c("AvgH", "AvgD", "AvgA"), MarketAvgC = c("AvgCH", "AvgCD", "AvgCA"),
                 Sportingbet = c("SBH", "SBD", "SBA"),
                 Ladbrokes = c("LBH", "LBD", "LBA")) 

odds_ou <- list(Bet365 = c("B365>2.5", "B365<2.5"), Bet365C = c("B365C>2.5", "B365C<2.5"), 
                Pinnacle = c("P>2.5", "P<2.5"), PinnacleC = c("PC>2.5", "PC<2.5"),
                MarketMax = c("Max>2.5", "Max<2.5"), MarketMaxC = c("MaxC>2.5", "MaxC<2.5"),
                MarketAvg = c("Avg>2.5", "Avg<2.5"), MarketAvgC = c("AvgC>2.5", "AvgC<2.5"))

odds_abbs <- list('1X2' = odds_1x2, OU = odds_ou)

football_odds <- suppressWarnings(mutate_at(football_odds, vars(unlist(odds_abbs)), as.numeric))

# percentage of non NA odds of different bookmakers
n_obs <- colMeans(!is.na(football_odds[, unlist(c(map(odds_1x2, 1),map(odds_ou, 1)))])) %>% 
  round(2) %>% sort(decreasing = T)

# Unique match identifier & OU result (Full time Line)
football_odds <- football_odds %>%
  mutate(ID = seq(1, nrow(.), 1), 
         FTL = ifelse(FTHG + FTAG > 2, ">2.5", "<2.5"))

# 2 - Build Odds Comparison Table -----------------------------------------

fair_probabilities <- function(data, bookie_ref, market) {
  
  odds <- odds_abbs[[market]][[bookie_ref]]
  
  # filter NA or invalid odds (overround > 1, all odds > 1)
  valid_overround <- rowSums(1/data[, odds], na.rm = T) > 1
  
  data <- data[valid_overround, ] %>% filter(if_all(all_of(odds), ~. > 1))
  
  fair_props <- implied_probabilities(data[, odds], method = "power") 
  is_ok <- which(!fair_props[["problematic"]])  
  
  col_names <- if (market == "OU") c("PR>2.5", "PR<2.5") else c("PRH", "PRD", "PRA")
  
  fair_props <- as.data.frame(fair_props[["probabilities"]]) %>%
    setNames(col_names)
  
  return (bind_cols(data[is_ok, "ID", drop = F], fair_props[is_ok, ]))
}

bookies_ref_1x2 <- setdiff(names(odds_1x2), grep("Max", names(odds_1x2), value = T))
bookies_ref_ou <- setdiff(names(odds_ou), grep("Max", names(odds_ou), value = T))

props_ref_1x2 <- lapply(bookies_ref_1x2, fair_probabilities, data = football_odds, market = "1X2")
names(props_ref_1x2) <- bookies_ref_1x2

props_ref_ou <- lapply(bookies_ref_ou, fair_probabilities, data = football_odds, market = "OU")
names(props_ref_ou) <- bookies_ref_ou

props_ref <- list("1X2" = props_ref_1x2, OU = props_ref_ou)

# free memory 
rm(list = c("props_ref_1x2", "props_ref_ou"))

# long table format to compare prices between bookie_ref & bookie_bet
# reference bookie is used as an estimator of true probability of an event
build_odds_table <- function(data, bookie_ref, bookie_bet, market) {
  
  res <- if (market == "1X2") "FTR" else "FTL"
  event <- if (market == "1X2") c("H", "D", "A") else c(">2.5", "<2.5")
  
  odds_ref <- odds_abbs[[market]][[bookie_ref]]; odds_bet <- odds_abbs[[market]][[bookie_bet]]
  
  df <- data[, c("ID", odds_ref, odds_bet, res)] %>% 
    filter(if_all(all_of(odds_bet), ~. > 1)) %>%
    drop_na()
  
  if (nrow(df) == 0) return (NULL)
  
  df <- inner_join(df, props_ref[[market]][[bookie_ref]], by = "ID")
  
  # custom melt
  dfList <- list()
  
  for (i in seq_along(event)) {
    dfList[[i]] <- df %>%
      select(ID, ends_with(event[i]), !!sym(res)) %>% 
      mutate(event = event[i])
  }
  
  bets <- dfList %>% 
    lapply(setNames, c("ID", "odds_ref", "odds_bet", "prop_ref", "res", "event")) %>%
    bind_rows() %>% 
    mutate(EV = prop_ref*odds_bet - 1,
           PnL = -1 + (res == event)*odds_bet) %>%
    mutate(bookie_ref = bookie_ref,
           bookie_bet = bookie_bet, .after = ID) %>%
    arrange(ID)
  
  return (bets)
}

# combinations to compare bookmakers against each other. bookmakers are compared also to their closing line
odds_table <- expand.grid(bookie_ref = names(odds_1x2), 
                          bookie_bet = names(odds_1x2), 
                          market = c("1X2", "OU"), stringsAsFactors = FALSE) %>%
  filter(bookie_bet != bookie_ref) %>%
  filter(!(market == "OU" & (!(bookie_bet %in% names(odds_ou)) | !(bookie_ref %in% names(odds_ou))))) %>%
  filter(!grepl("Max", bookie_ref)) %>%
  
  pmap_dfr(build_odds_table, data = football_odds, .progress = T) 

# 3 - Expected Value ------------------------------------------------------

# PnL = f(EV) for all bets using bookie ref to estimate true probability of an event
ev_table <- odds_table %>% 
  group_by(bookie_ref, bookie_bet) %>%
  mutate(EV_int = ntile(EV, 40)) %>%
  group_by(EV_int, .add = T) %>%
  summarize(n = n(),
            EV_mean = mean(EV), 
            PnL_mean = mean(PnL), .groups = "drop")

plot_ev <- function(bookies_ref, bookies_bet) {
  
  pl <- ev_table %>%
    filter(bookie_ref %in% bookies_ref, bookie_bet %in% bookies_bet) %>%
    filter(between(EV_mean, -0.25, 0.25)) %>%
    ggplot(aes(x = EV_mean, y = PnL_mean, color = bookie_bet)) + geom_point(alpha = 0.5) + 
    geom_smooth(method = "lm", se = F) + 
    geom_abline(slope = 1, intercept = 0) + 
    annotate("segment", x = -Inf, xend = 0, y = 0, yend = 0, linetype = "dashed") +
    annotate("segment", x = 0, xend = 0, y = -Inf, yend = 0, linetype = "dashed") +
    facet_wrap(bookie_ref ~.) + 
    theme_bw() +
    labs(title = "Pnl[bet] ~ f(EV[bet, ref])", 
         x = "Expected Value", y = "Profit & Loss")
  print(pl)
}

lines <- c(names(odds_1x2)[!grepl("C$", names(odds_1x2))], "VC")
closing_lines <- setdiff(names(odds_1x2), lines)

# pre-closing to pre-closing
plot_ev(lines, lines)
# closing to 'closing' 
plot_ev(closing_lines, closing_lines)
# closing to 'pre-closing' 
plot_ev(closing_lines, lines)

# correlation between PnL and EV
cor_table <- ev_table %>% 
  summarize(obs = sum(n),
            cor = cor(PnL_mean, EV_mean), .by = c(bookie_ref, bookie_bet)) 

cor_table %>%
  filter(bookie_bet %in% lines) %>%
  summarise(cor_mean = mean(cor), .by = bookie_ref) %>%
  arrange(cor_mean) %>%
  mutate(bookie_ref = factor(bookie_ref, levels = bookie_ref)) %>%
  ggplot(aes(x = bookie_ref, y = cor_mean, fill = bookie_ref)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis(discrete = T, option = "B") +
  coord_flip() + guides(fill="none") + scale_y_continuous(labels = percent_format()) +  theme_bw() +
  labs(x = NULL, y = "mean correlation", title = "Mean Correlation of Bookmakers against pre-closing lines")  

cor_table %>% 
  filter(bookie_bet %in% setdiff(closing_lines, c("MarketMaxC")), 
         bookie_ref %in% setdiff(closing_lines, c("MarketMaxC"))) %>%
  ggplot(aes(x = bookie_bet, y = bookie_ref, fill = cor)) +
  geom_tile() +
  scale_fill_viridis(option = "F", direction = -1) +
  labs(title = "EV and PnL Correlation between Bookie Ref and Bet",
       subtitle = "Non closing Lines",
       x = "Bookie Bet", y = "Bookie Ref", fill = "Corr") +
  theme_bw()
# The efficiency of Pinnacle's Closing price is unparalleled and can be treated
# as a true estimator of the true probability of an event

# 4 - Odds Movement -------------------------------------------------------

# Distribution reshaping - pre-closing & closing EV compared to Pinnacle's closing price
# distr is not normal, it's negatively skewed (Gumbel or Dagum for fitting)
odds_table %>% 
  filter(bookie_ref == "PinnacleC") %>%
  filter(bookie_bet %in% names(odds_1x2)[1:6]) %>% 
  ggplot(aes(x = EV, y = after_stat(density), color = bookie_bet)) +
  geom_histogram(binwidth = 0.005) +
  geom_density() + 
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  facet_wrap(bookie_bet ~. , ncol = 2) + theme_bw() +
  coord_cartesian(xlim = c(-0.25, 0.25)) + guides(color = "none") +
  labs(title = "Expected Value distribution using Pinnacle's closing Prices",
       subtitle = "Before and at Closing odds for different bookmakers")

# Nice summary of KL divergence (https://www.perfectlynormal.co.uk/blog-kl-divergence)
# «the KL divergence represents the amount you can win from the casino by exploiting the difference between the true probabilities 
# P and the house's false beliefs Q. The closer and Q are, the harder it is to profit from your extra knowledge.»

odds_table %>% 
  filter(bookie_ref == "PinnacleC") %>%
  filter(bookie_bet %in% names(odds_1x2)[1:6]) %>% 
  mutate(DKLPQ = prop_ref*log(prop_ref*odds_bet)) %>%
  ggplot(aes(x = DKLPQ, y = after_stat(density), color = bookie_bet)) +
  geom_histogram(binwidth = 0.0025) +
  geom_density() + 
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  facet_wrap(bookie_bet ~. , ncol = 2) + theme_bw() +
  coord_cartesian(xlim = c(-0.1, 0.1)) + guides(color = "none") +
  labs(title = "KL divergence distribution using Pinnacle's closing Prices", 
       subtitle = "Before and at Closing odds for different bookmakers")

# 5 - Value Bets ----------------------------------------------------------

# percentage of value bets (EV > 0) using PinnacleC - Are some markets more error-prone than others ?
# Using Deep Learning we could identify any systematic modelling biases in pre-closing odds

# group_var can be {bookie_bet, event, Div, Season}
plot_vb_pct <- function(bookies, group_var) {
  
  var <- sym(group_var)
  
  pl <- odds_table %>%
    filter(bookie_ref == "PinnacleC", bookie_bet %in% bookies) %>%
    left_join(football_odds[, c("ID", "Div", "Season")], by = "ID") %>%
    group_by(!!var) %>%
    summarise(EV_pos_pct = (sum(EV > 0) / n()) * 100, .groups = "drop") %>%
    arrange(EV_pos_pct) %>%
    mutate(var_f = factor(!!var, levels = !!var)) %>%
    ggplot(aes(x = var_f, y = EV_pos_pct, fill = var_f)) +
    geom_bar(stat = "identity") + coord_flip() +
    scale_fill_viridis(discrete = T, option = "B") + guides(fill = "none") + theme_bw() +
    labs(x = NULL, y = "Percentage of bets with positive EV")
  print(pl)
}

plot_vb_pct(c(lines, closing_lines), "bookie_bet")
plot_vb_pct("Bet365", "event")

# 6 - Real time Value betting ---------------------------------------------

value_bet_series <- odds_table %>% 
  filter(EV > 0) %>%
  group_by(bookie_ref, bookie_bet) %>%
  arrange(ID) %>%
  mutate(NVBets = seq(1, n()),
         EV_total = cumsum(EV),
         PnL_total = cumsum(PnL)) %>% 
  ungroup() %>%
  mutate(Yield = PnL_total/NVBets)

plot_vb_series <- function(bookies_ref, bookies_bet) {
  
  pl_returns <- value_bet_series %>%
    filter(bookie_ref %in% bookies_ref, bookie_bet %in% bookies_bet) %>%
    ggplot(aes(color = bookie_bet)) + 
    geom_line(aes(x = NVBets, y = PnL_total), alpha = 0.5) + 
    geom_smooth(aes(x = NVBets, y = EV_total, linetype = "dashed"), method = "lm", formula = y ~ 0 + x) +
    facet_wrap(bookie_ref ~., scales = "free") +
    theme_bw() +
    scale_linetype_manual(name = NULL, values = "dashed", labels = "EV_total") +
    labs(title = "Cumulative ProfitnLoss and EV ~ f(Number of Value Bets)", 
         x = "Value Bets", y = "PnL") + theme_bw() 
  print(pl_returns)
}

plot_vb_series("Pinnacle", setdiff(lines, "MarketMax"))

# Inefficiencies between pre-closing lines or between closing lines
# Identifying which bookmaker can be used to profit against another
real_time_value_bets <- value_bet_series %>%
  filter( (bookie_ref %in% lines & bookie_bet %in% lines) | 
            (bookie_ref %in% closing_lines & bookie_bet %in% closing_lines)) %>%
  group_by(bookie_ref, bookie_bet) %>%
  slice_max(NVBets) %>%
  ungroup() %>%
  filter(PnL_total > 0) %>%
  select(bookie_ref, bookie_bet, NVBets:Yield)

# bootstrapping to calculate how likely is to achieve profit as calculated by value betting strategy
bootstrap <- function(bookie, nbets, nsims = 1000) {
  
  params <- as.list(environment())[1:2]
  odds <- odds_1x2[[bookie]]
  
  data <- na.omit(football_odds[, c(odds, "FTR")]) %>%
    filter(if_all(all_of(odds), ~. > 1)) %>%
    setNames(c("H", "D", "A", "FTR")) %>%
    pivot_longer(cols = -FTR, names_to = "selection", values_to = "odds") %>%
    mutate(success = (FTR == selection), 
           PnL = -1 + success*odds)   
  
  # random sample, n bets from every bin with replacement & summarize
  sampling <- function(df) {
    df %>% 
      slice_sample(n = nbets, replace = T) %>%
      summarise(yield = sum(PnL) / nbets) 
    }
  
  # sampling nsims times
  sim_res <- replicate(nsims, sampling(data), simplify = F) %>% 
    bind_rows() %>%
    summarise(yield_avg = mean(yield),
              yield_sd = sd(yield)) %>%
    bind_cols(params)

  return (sim_res)
}

bootstraping <- real_time_value_bets %>% 
  select(bookie_bet, NVBets) %>%
  distinct(.keep_all = T) %>%
  rename(nbets = NVBets, bookie = bookie_bet) %>%
  pmap_dfr(bootstrap, .progress = T)

real_time_value_bets <- left_join(real_time_value_bets, bootstraping, by = c("bookie_bet" = "bookie", 
                                                                             "NVBets" = "nbets")) %>%
  mutate(p_value = pmap_dbl(list(q = Yield, mean = yield_avg, sd = yield_sd), pnorm, lower.tail = F), 
         logp = log10(p_value)) %>% 
  select(-c(yield_avg, yield_sd)) %>%
  mutate(p_value = format(p_value, digits = 1, scientific = T)) %>% 
  arrange(desc(PnL_total)) 
