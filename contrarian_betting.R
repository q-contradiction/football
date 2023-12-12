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

rm(list = ls())

library(plyr)
library(rvest) 
library(reshape2) 
library(lubridate)
library(tidyverse)

options(scipen = 999)

# 1 - Load Football data ----

load_league_data <- function(league) {
  
  print(league)
  url_main <- "https://www.football-data.co.uk/"
  url_league <- paste0(url_main, league, "m.php")
  
  links <-  read_html(url_league) %>% html_nodes("a") %>% html_attr("href")
  
  # main leagues end in (1|0).csv
  csv_links <- paste0(url_main, links[grep("(1|0)\\.csv$", links)])  
  data <- vector("list", length(csv_links))
  
  # Add a season column 
  for (i in seq_along(csv_links)) {
    
    year <- unlist(regmatches(csv_links[i], gregexpr("\\d+", csv_links[i])))[2]  %>%
      substr(1, 2) %>% as.numeric()
    
    season <- ifelse(between(year, 0, 49), 2000+year, 1900+year)
    
    read_url <- function(URL) { 
      tryCatch(read.csv(URL, stringsAsFactors = F, check.names = F, na.strings = c("NA", "")),
               error = function(e) { warning(conditionMessage(e)); NULL })
    }
    
    data[[i]] <-  read_url(csv_links[i]) 
    
    # remove any NA rows or columns with very few observations
    if (!is.null(data[[i]])) {
      data[[i]] <- data[[i]] %>% .[rowSums(is.na(.)) != ncol(.), colSums(is.na(.)) < nrow(.) - 5] %>%
        mutate(Season = season) %>%
        relocate(Season, .after = Date)
    }
    
  }
  
  return (data)
}

# download data from Football-data.co.uk - main leagues
leagues <- c("england", "scotland", "germany", "italy", "spain", 
             "france", "netherlands", "belgium", "portugal")

# Load all CSV files into a list
football_odds <- lapply(leagues, load_league_data) 

# Merge into a single dataframe 
football_odds <- football_odds %>% do.call(c, .) %>% rbind.fill() 


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
                  Market_Max = c("MaxH", "MaxD", "MaxA"), Market_Avg = c("AvgH", "AvgD", "AvgA"), 
                  Interwetten =  c("IWH", "IWD", "IWA")) 

# 2 - Contrarian Betting ----

# remove margin from odds, power method
remove_overround <- function(bookmaker_odds) {
  
  bookmaker_odds <- as.numeric(bookmaker_odds)
  
  n <- length(bookmaker_odds)
  overround_props <- 1/bookmaker_odds

  # If NA, or any odds < 1, or margin < 0
  if (n < 2 || any(is.na(overround_props)) || any(overround_props > 1) || sum(overround_props) < 1) 
    return (rep(NA, n))
  
  sum_prop <- function(bias) {
    props <- overround_props^(1/bias)
    return (sum(props) - 1)
  }
  
  bias_ <- tryCatch(uniroot(sum_prop, c(0.5, 1))$root, 
                    error = function(e) { warning(conditionMessage(e)); NA })
  
  fair_props <- overround_props^(1/bias_) %>% {./sum(.)}
  
  
  return (fair_props)
  
}

# If we use a specific bookmaker to estimate true probability some odds might be missing, 
# so estimating overround prop as the mean of non NA reversed odds and then remove margin,
# useful for some leagues in past seasons, where we have fewer odds data
fair_probabilities <- function(df) {
  
  odds <- odds_abbs[setdiff(names(odds_abbs), c("Market_Max"))]
  
  for (i in 1:3) {
    dh <- df[, unlist(odds)] %>%
      dplyr::select(ends_with(c("H", "D", "A")[i])) %>%
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
  
  df <- df %>% arrange(Date) %>% ddply(.(ID), function(ID) {
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
  filter(between(Season, 2010, 2022)) %>%
  mutate(Date = dmy(Date)) %>% arrange(Date) %>%
  mutate(ID = as.numeric(rownames(.))) %>%
  fair_probabilities() 


ratings <- ddply(football_odds, .(Season, Div), team_ratings, last_n = 6) 

# distribution of match ratings - normally distributed
ratings %>% 
  ggplot(aes(x = match_rating)) +
  geom_histogram(aes(y = ..density..), binwidth = 0.05, color = "red") +
  stat_function(fun = dnorm, args = list(mean = mean(ratings$match_rating), sd = sd(ratings$match_rating))) +
  ggtitle("Distribution of match ratings") 


# Contrarian betting is to back the team with the lowest rating
contrarian_betting <- function(data, bookmaker, threshold = 0) {
  
  # Odds for Home and Away - draw removed
  odds <- c(odds_abbs[[bookmaker]][1], odds_abbs[[bookmaker]][3]) 
  
  data <- data[, c("ID", "Season", "match_rating", "FTR", odds)]
  
  data <- data %>% 
    arrange(ID) %>%
    filter(abs(match_rating) > threshold) %>%
    mutate(cold_selection = ifelse(match_rating < -threshold, "H", "A"),
           cold_odds = ifelse(match_rating < -threshold, get(odds[1]), get(odds[2])),
           
           hot_selection = ifelse(match_rating > threshold, "H", "A"),
           hot_odds = ifelse(match_rating > threshold, get(odds[1]), get(odds[2])),
           
           # FTR in an appropriate format
           FTR_ = ifelse(FTR == cold_selection, "N", ifelse(FTR == hot_selection, "P", "D")),
           
           cold_PnL = ifelse(FTR_ == "N", cold_odds - 1, -1), 
           hot_PnL = ifelse(FTR_ == "P", hot_odds - 1, -1)) %>%
    na.omit() 
  
  bets <- data %>%
    mutate(NBets = seq(1, by = 1, length.out = nrow(.)), 
           cold_PnL = cumsum(cold_PnL), 
           hot_PnL = cumsum(hot_PnL), 
           cold_yield = cold_PnL/NBets,
           hot_yield = hot_PnL/NBets,
           bookmaker = bookmaker, 
           threshold = as.character(threshold)) 
  
  return (bets)
}

params <- data.frame(bookmaker = rep(c("Bet365", "PinnacleC", "Pinnacle", "Market_Max"), each = 4),
                     threshold = c(0, 0.75, 1.5, 2.25))

returns <- params %>%
  pmap(contrarian_betting, data = ratings) %>% bind_rows() %>%
  melt(id.vars = c("Season", "NBets", "bookmaker", "threshold"),
       measure.vars = c("cold_PnL", "hot_PnL"))

# In all bookmakers returns when backing "colder"teams are higher
# And using the best price in the market we can achieve a positive yield
# We can also achieve minor yield in Pinnacle if we set a high threshold
pl <- ggplot(returns, aes(x = NBets, y = value, color = bookmaker, linetype = variable)) +
  geom_line() + 
  scale_linetype_manual(values = c("cold_PnL" = "solid", "hot_PnL" = "dotted")) +
  facet_wrap(threshold ~., labeller = label_both,  scales = "free") +
  labs(title = "Returns using Contrarian betting system", 
       subtitle = "Backing colder and hotter teams", 
       y = "Profit n Loss")
print(pl)


