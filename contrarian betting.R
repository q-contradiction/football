# Slight modification of Contrarian Betting Strategy proposed by Joseph Buchdahl: 

# https://www.football-data.co.uk/profitable_betting_system.php 

# Contrarian betting, back the team with lower scoring rating, 
# the "colder" team that has performed relatively worse than expected

# Original Scoring System: if a team wins with probability p, 1-p, if loses -p

# New Scoring System: Draw is a positive result for "bad" teams and negative result for "good" teams, 
# so it is better not to be treated always as loss.  
# In case of a draw award: probability_draw - p. So now: {-p + [if(W, 1) if(D, PD) if(L, 0)]}

# Parametrization includes last n games in ratings & a threshold rating when backing teams

# data from football-data.co.uk

suppressPackageStartupMessages({
  library(plyr)
  library(tidyverse)
})

rm(list = ls())
options(scipen = 999)

# 1. Data ----
football_odds <- readRDS(file = "./data-raw/football-data.rds")

# Substitute where avg/max market odds are NA with Betbrain odds
cols <- list(c("MaxH", "BbMxH"), c("MaxD", "BbMxD"), c("MaxA", "BbMxA"),
             c("AvgH", "BbAvH"), c("AvgD", "BbAvD"), c("AvgA", "BbAvA"))

for (v in cols)  {
  replace_indexes <- is.na(football_odds[[v[1]]])
  football_odds[[v[1]]][replace_indexes] <- football_odds[[v[2]]][replace_indexes]
}

# Remove BetBrain odds
football_odds <- football_odds %>% .[setdiff(names(.), unlist(map(cols, 2)))]

# Odds Abbreviations for 1x2 in Football-data.co.uk
odds_abbs <- list(Bet365 = c("B365H", "B365D", "B365A"), Pinnacle = c("PSH", "PSD", "PSA"),
                  PinnacleC = c("PSCH", "PSCD", "PSCA"), WilliamHill = c("WHH", "WHD", "WHA"), 
                  Sportingbet =  c("SBH", "SBD", "SBA"), Ladbrokes = c("LBH", "LBD", "LBA"),
                  MarketMax = c("MaxH", "MaxD", "MaxA"), MarketAvg = c("AvgH", "AvgD", "AvgA"), 
                  Interwetten =  c("IWH", "IWD", "IWA")) 

# 2 - Probabilities & Ratings ----

# remove margin from odds, power method finding the value of k, such as: sum(prop^k) = 1
remove_overround <- function(bookmaker_odds) {
  
  bookmaker_odds <- as.numeric(bookmaker_odds)
  
  n <- length(bookmaker_odds)
  overround_props <- 1/bookmaker_odds

  # invalid odds
  if (n < 2 || any(is.na(overround_props)) || any(overround_props > 1) || sum(overround_props) < 1) 
    return (rep(NA, n))
  
  sum_prop <- function(k) {
    props <- overround_props^(1/k)
    return (sum(props) - 1)
  }
  
  k <- tryCatch(uniroot(sum_prop, c(0.5, 1))$root, 
                    error = function(e) { warning(conditionMessage(e)); NA })
  
  fair_props <- overround_props^(1/k) %>% {./sum(.)}
  
  return (fair_props)
}

# If we use a specific bookmaker to estimate true probability some odds might be missing, 
# so estimating overround prop as the mean of non NA reversed odds and then remove margin,
# useful for some leagues in past seasons, where we have fewer odds data
fair_probabilities <- function(df) {
  
  odds <- odds_abbs[setdiff(names(odds_abbs), c("MarketMax"))]
  
  for (i in 1:3) {
    dh <- df[, unlist(odds)] %>%
      select(ends_with(c("H", "D", "A")[i])) %>%
      data.matrix()
    
    prop <- rowMeans(1/dh, na.rm = T)
    col <- c("OH", "OD", "OA")[i]
    
    df[[col]] <- 1/prop
  }
  
  fair_props <- mapply(c, df[["OH"]], df[["OD"]], df[["OA"]], SIMPLIFY = F) %>%
    map(remove_overround) 
  
  df <- df %>% 
    mutate(PH = sapply(fair_props, `[[`, 1), 
           PD = sapply(fair_props, `[[`, 2), 
           PA = sapply(fair_props, `[[`, 3)) %>%
    filter(!is.na(PH)) %>%
    relocate(PH, PD, PA, .before = "B365H")
  
  return (df)
}

# A slight modification to the existing strategy is to take into account only the last n games
team_ratings <- function(df, last_n = Inf) {
  
  update_rating <- function(rating, p, last_n) {
    
    if (length(rating) >= last_n) rating <- rating[-1]
    
    rating <- c(rating, p)
    return(rating)
  }
  
  # Initialize team's rating  
  rating <- unique(df$HomeTeam) %>%
    lapply(function(team) numeric(0)) %>% setNames(unique(df$HomeTeam))
  
  df <- df %>% 
    arrange(Date) %>% ddply(.(ID), function(ID) {
      HT <- ID$HomeTeam
      AT <- ID$AwayTeam
      
      # pre game ratings
      ID$home_rating <- sum(rating[[HT]])
      ID$away_rating <- sum(rating[[AT]])
      
      # update team ratings 
      p_home <- -ID$PH + ifelse(ID$FTR == "H", 1, ifelse(ID$FTR == "D", ID$PD, 0))
      p_away <- -ID$PA + ifelse(ID$FTR == "A", 1, ifelse(ID$FTR == "D", ID$PD, 0))
      
      rating[[HT]] <<- update_rating(rating[[HT]], p_home, last_n)
      rating[[AT]] <<- update_rating(rating[[AT]], p_away, last_n)
      
      return(ID)
    }) %>% 
    mutate(match_rating = home_rating - away_rating) %>%
    filter(match_rating != 0) %>%
    relocate(home_rating, away_rating, match_rating, .after = AwayTeam)
  
  return (df)
}

football_odds <- football_odds %>% 
  mutate(ID = as.numeric(rownames(.))) %>%
  fair_probabilities() 

ratings <- ddply(football_odds, .(Season, Div), team_ratings, last_n = 6) 

# normally distributed match ratings
ratings %>% 
  ggplot(aes(x = match_rating)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.05, color = "red") +
  stat_function(fun = dnorm, args = list(mean = mean(ratings$match_rating), sd = sd(ratings$match_rating))) +
  ggtitle("Distribution of match ratings") 

# 3 - Contrarian Betting ----

# Contrarian betting is to back the team with the lowest rating
contrarian_betting <- function(data, bookmaker, threshold = 0) {
  
  # Odds for Home and Away - draw removed
  odds <- c(odds_abbs[[bookmaker]][1], odds_abbs[[bookmaker]][3]) 
  
  bets_contr <- data[, c("ID", "Season", "Date", "match_rating", "FTR", odds)] %>% 
    mutate_at(vars(unlist(odds)), as.numeric) %>% drop_na() %>%
    arrange(ID) %>%
    filter(abs(match_rating) > threshold) %>%
    mutate(selection_cold = ifelse(match_rating < -threshold, "H", "A"),
           odds_cold = ifelse(match_rating < -threshold, get(odds[1]), get(odds[2])),
           
           selection_hot = ifelse(match_rating > threshold, "H", "A"),
           hot_odds = ifelse(match_rating > threshold, get(odds[1]), get(odds[2])),
           
           PnL_cold = ifelse(FTR == selection_cold, odds_cold - 1, -1), 
           PnL_hot = ifelse(FTR == selection_hot, hot_odds - 1, -1))
             
  bets_contr <- bets_contr %>%
    mutate(NBets = seq(1, by = 1, length.out = nrow(.)), 
           cold = cumsum(PnL_cold), 
           hot = cumsum(PnL_hot), 
           bookmaker = bookmaker, 
           threshold = as.character(threshold)) 
  
  return (bets_contr)
}

params <- data.frame(bookmaker = rep(c("Bet365", "PinnacleC", "Pinnacle", "MarketMax"), each = 4),
                     threshold = c(0, 0.75, 1.5, 2.25))

returns <- params %>%
  pmap_dfr(contrarian_betting, data = ratings) %>% 
  pivot_longer(
    cols = c("cold", "hot"),
    names_to = "strategy",
    values_to = "PnL"
  )

# In all bookmakers returns when backing "colder"teams are higher
# And using the best price in the market we can achieve a positive yield
# We can also achieve minor profit in Pinnacle if we set a high threshold
ggplot(returns, aes(x = NBets, y = PnL, color = bookmaker, linetype = strategy)) +
  geom_line() + 
  scale_linetype_manual(values = c("cold" = "solid", "hot" = "dotted")) +
  facet_wrap(threshold ~., labeller = label_both,  scales = "free") +
  theme_bw() + 
  labs(title = "PnL backing colder and hotter teams")
