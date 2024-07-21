# Evaluating different methods to remove bookmaker's overround in football odds

# A reproduction in R of the methodology from: 
# Wisdom of Crowd - Joseph Buchdahl (www.football-data.co.uk/blog/wisdom_of_the_crowd.php)

# Additional methods included in the implied package,
# applying same logic to different bookmakers and over-under odds.

# Find the optimal method to arrive at a set of fair probabilities
# using rank probability score, linear regression and returns from fair odds

# Data from football-data.co.uk

suppressPackageStartupMessages({
  library(implied) #https://opisthokonta.net/?p=1797
  library(microbenchmark)
  library(scales)
  library(broom)
  library(tidyverse)
})

rm(list = ls())

# Data --------------------------------------------------------------------

football_odds <- readRDS(file = "./data-raw/football-data.rds")

# odds abbreviations for 1x2 in football-data.co.uk
odds_1x2 <- list(Bet365 = c("B365H", "B365D", "B365A"), Bwin = c("BWH", "BWD", "BWA"),
                 Interwetten =  c("IWH", "IWD", "IWA"), Pinnacle = c("PSH", "PSD", "PSA"), 
                 WilliamHill = c("WHH", "WHD", "WHA"),  VC = c("VCH", "VCD", "VCA"), 
                 Sportingbet =  c("SBH", "SBD", "SBA"), Blue_Square =  c("BSH", "BSD", "BSA")) 

# odds Abbreviations for over-under
odds_ou <- list(Bet365 = c("B365>2.5", "B365<2.5"), Pinnacle = c("P>2.5", "P<2.5"))

odds_abbs <- list('1X2'= odds_1x2, OU = odds_ou)

football_odds <- suppressWarnings(mutate_at(football_odds, vars(unlist(odds_abbs)), as.numeric))

# OU result (Full time 2.5 Line)
football_odds <- football_odds %>%
  mutate(ID = as.numeric(rownames(.)),
         FTL = ifelse(FTHG + FTAG > 2, "O", "U"))

n_obs <- colMeans(!is.na(football_odds[, unlist(c(map(odds_1x2, 1),map(odds_ou, 1)))])) %>% 
  round(2) %>% sort(decreasing = T)

max_odds <- football_odds %>%
  select(all_of(unlist(odds_1x2))) %>% lapply(max, na.rm = T) 

# Methods -----------------------------------------------------------------

# methods supported in 'implied_probabilities' from implied Package
methods <- c("basic", "wpo", "bb", "or", "power", "shin", "jsd", "additive")

## Complexity ----

x <- c(1.5, 3.5, 8)

microbenchmark(list = lapply(methods, function(method) {bquote(implied_probabilities(x, .(method)))})) %>%
  mutate(method = str_extract(expr, '(?<=\").*(?=\")')) %>%
  # order by mean time
  mutate(method = factor(method, 
                         levels = aggregate(time ~ method, data = ., FUN = median) %>% arrange(time) %>% pull(method))) %>%
  ggplot(aes(x = method, y = log2(time), fill = method)) + 
  geom_boxplot() + theme_bw() + guides(fill = "none")

## Applied margin ----

# mapping fair probabilities to overround probabilities, based on & method & number of events & margin 
odds_mapper <- function(method, events, margin = 0.08) {
  
  params <- as.list(environment())
  
  # Generate grid of fair probabilities - all possible combination with step 0.1
  # can't allocate enough memory if events > 3 - increase step
  fair_props <- do.call(expand.grid, replicate(events, seq(0.02, 0.98, 0.01), simplify = FALSE)) %>% 
    filter(rowSums(.) == 1) %>% setNames(paste0("P", 1:events))
  
  bookmakers_odds  <- implied_odds(fair_props, method, margin) %>% .[["odds"]] %>%
    as.data.frame() %>% setNames(paste0("O", 1:events))
  
  odds <- bind_cols(fair_props, bookmakers_odds) %>% 
    select(P1, O1) %>% rename(fair_props = P1, odds = O1) %>%
    mutate(fair_odds = 1/fair_props) %>% bind_cols(params) %>%
    mutate_at(vars("events", "margin"), as.character)
  
  return (odds)
}

# jsd is not supported with 'implied_odds', additive same transformation as wpo, exclude both
odds <- data.frame(method = rep(methods[1:6], each = 2),
                   events = c(2, 3)) %>%
  pmap_dfr(odds_mapper)

# wpo and shin identical mapping in 2 events
odds %>% 
  ggplot(aes(x = fair_props, y = 1/odds - fair_props, color = method)) +
  geom_point(aes(alpha = events)) +
  facet_wrap(events ~., labeller = label_both) +
  scale_alpha_manual(name = NULL, values = c(1, 0.1)) + 
  theme_bw() +
  labs(x = "Fair Probability", y = "Applied Margin") +
  guides(alpha = "none") 

# Fair probabilities ------------------------------------------------------

fair_probabilities <- function(method, bookmaker, market) {
  
  params <- as.list(environment())
  
  odds <- odds_abbs[[market]][[bookmaker]] 
  
  odds_cols <- if (market == "OU") c("OO", "OU") else c("OH", "OD", "OA")
  prop_cols <- if (market == "OU") c("PO", "PU") else c("PH", "PD", "PA")
  res <- if (market == "OU") "FTL" else "FTR"
  
  # filter NA or invalid odds (overround > 1, all odds > 1)
  valid_overround <- rowSums(1/football_odds[, odds], na.rm = T) > 1
  
  df <- football_odds[valid_overround, c("ID", res, odds)] %>%
    setNames(c("ID", "Res", odds_cols)) %>%
    drop_na() %>%
    filter(if_all(all_of(odds_cols), ~. > 1)) 
  
  fair_props <- implied_probabilities(df[, odds_cols], method) 
  
  is_ok <- which(!fair_props[["problematic"]])  
  fair_props <- as.data.frame(fair_props[["probabilities"]]) %>%
    setNames(prop_cols)
  
  return (bind_cols(params, df[is_ok, ], fair_props[is_ok, ]))
}

props_table_1x2 <- expand.grid(bookmaker = c("Pinnacle", "Bet365", "Bwin", "WilliamHill"),
                               method = methods[1:6], 
                               stringsAsFactors = F) %>%
  pmap_dfr(fair_probabilities, market = "1X2", .progress = T)

props_table_ou <- expand.grid(bookmaker = c("Pinnacle", "Bet365"),
                              method = methods[1:5], 
                              stringsAsFactors = F) %>%
  pmap_dfr(fair_probabilities, market = "OU", .progress = T)

props_table <- list("1X2" = props_table_1x2, "OU" = props_table_ou)

# free memory
rm(list = c("props_table_1x2", "props_table_ou"))

# Evaluation --------------------------------------------------------------

## A. Rank probability score (1x2), Brier Score (OU) ----
cum_props <- t(apply(data.matrix(select(props_table[["1X2"]], PH, PD, PA)), 1, cumsum))

res <- tribble(
  ~Res, ~r1, ~r2, ~r3,
  "H" ,  1 ,  1 ,  1 ,
  "D" ,  0 ,  1 ,  1 ,
  "A" ,  0 ,  0 ,  1 
)

res_props <- left_join(props_table[["1X2"]], res, by = "Res") %>%
  select(r1, r2, r3) %>% data.matrix()

props_table[["1X2"]]$RPS <- 0.5*rowSums((cum_props - res_props)^2)

rank_probability_score <-  props_table[["1X2"]] %>%
  group_by(bookmaker, method, market) %>%
  summarize(rps_avg = mean(RPS), .groups = "drop") %>%
  arrange(rps_avg)

brier_score <- props_table[["OU"]] %>%
  mutate(BS = case_when(
    Res == "O" ~ (1 - PO)^2, 
    Res == "U" ~ (1 - PU)^2,
    NA ~ NA)) %>%
  group_by(bookmaker, method, market) %>%
  summarize(bs_avg = mean(BS), .groups = "drop") %>%
  arrange(bs_avg)

## B. Actual & Implied prob (regression) ----
props_long <- props_table[["1X2"]] %>%
  pivot_longer(cols = c(PH, PD, PA), names_to = "selection", values_to = "implied_probability") %>%
  mutate(evaluation = (Res == substr(selection, 2, 2))) %>%
  group_by(bookmaker, method, market) %>%
  mutate(implied_prop_int = cut(implied_probability, breaks = seq(0, 1, 0.01))) %>%
  group_by(implied_prop_int, .add = T) %>%
  summarise(obs = n(), 
            implied_prop = mean(implied_probability),
            actual_prop = mean(evaluation), .groups = "drop") 

reg <- props_long %>%
  group_nest(bookmaker, method, market) %>%
  mutate(model = map(data, ~ lm(actual_prop ~ implied_prop - 1, weights = obs, data = .x)), 
         coef = map(model, tidy), 
         adj_r_squared = map_dbl(model, ~ glance(.x) %>% pull(adj.r.squared))) %>%
  select(bookmaker, method, market, coef, adj_r_squared) %>%
  unnest(cols = c(coef)) %>%
  # deviation from optimal line y = x
  mutate(dev = abs(estimate - 1)) %>%
  arrange(dev)

ggplot(props_long, aes(x = implied_prop, y = actual_prop, color = method)) +
  geom_point(alpha = 0.2) + geom_smooth(method = "lm", formula = y ~ x-1,  mapping = aes(weight = obs), se = F) + 
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(bookmaker ~.) + theme_bw() +
  labs(x = "Odds implied Probability", y = "Actual Probability")

## C. Fair PnL & Yield ----
fair_returns <- function(props) {
  
  nevents = ifelse(unique(props$market) == "1X2", 3, 2)
  
  # fair_PnL = Profit n Loss if bet on fair odds of all outcomes 
  if (nevents == 3) 
    props <- mutate(props, fair_PnL = -3 + ifelse(Res == "H", 1/PH, ifelse(Res == "D", 1/PD, 1/PA)))
  else 
    props <- mutate(props, fair_PnL = -2 + ifelse(Res == "O", 1/PO, 1/PU))
  
  fair_bets <- props %>%
    group_by(bookmaker, method, market) %>%
    mutate(NBets = seq(nevents, by = nevents, length.out = n()), 
           fair_PnL_r = cumsum(fair_PnL), 
           yield = fair_PnL_r/NBets, 
           .keep = "used") %>%
    ungroup()
  
  return (fair_bets)
}

fair_bets <- map_dfr(props_table, fair_returns) 

# power and then wpo more suitable for 1x2 - bb for OU
fair_bets %>% 
  ggplot(aes(x = NBets, y = fair_PnL_r, color = method)) + 
  geom_line() + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  facet_wrap(bookmaker ~ market, scales = "free") +
  theme_bw() + scale_x_continuous(labels = label_number(scale = 1e-3, suffix = "k")) +
  labs(title = "Profit n Loss with fair odds on all outcomes",
       subtitle = "Different overround removal methods")

fair_bets %>%
  filter(market == "1X2") %>%
  group_by(bookmaker, method) %>%
  slice_max(NBets) %>%
  ungroup() %>%
  ggplot(aes(x = method, y = yield, fill = bookmaker)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") + 
  theme_bw() + scale_y_continuous(labels = percent_format()) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 12))
