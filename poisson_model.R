# The basic Poisson Model to estimate goal expectancy - Maher, Modelling Association Football Scores (1982)

# Applying the Poisson model to derive probabilities for the 2nd round of every season since 1920 and evaluating
# the efficiency of the model over the decades... 

# Has the Poisson approximation become better or worse to describe basic football dynamics? 

# engsoccerdata - contains results before 1900

rm(list = ls())

library(plyr)
library(reshape2) 
library(lubridate)
library(skellam)
library(engsoccerdata)
library(tidyverse)

options(scipen = 999)

# 1 - Load data ----

football_data <- england %>%
  mutate(ID = seq(1, nrow(.), 1), 
         Date = ymd(Date), 
         decade = sapply(Season, function(yr){ (yr - yr %% 10)})) %>%
  arrange(Date, division)

football_data %>% filter(between(decade, 1920, 2010)) %>%
  ggplot(aes(x = goaldif, color = factor(decade))) +
  geom_histogram(aes(y = ..density..), binwidth = 1) +
  coord_cartesian(xlim = c(-4, 6)) +
  facet_wrap(decade ~.) +
  ggtitle("Goal difference distr over the years")

# The goal difference in favor of the home is decreasing, the home effect diminishes over the years
# Goals before 1975 showed an inconsistent behavior oscillating, 
# in recent years the amount of goals scored is less variable.
football_data %>% 
  group_by(Season) %>%
  summarize(n = n(),
            hgoalm = mean(hgoal), 
            vgoalm = mean(vgoal),
            goaldiffm = mean(goaldif),
            totgoalm = mean(totgoal)) %>%
  melt(measure.vars = c("hgoalm", "vgoalm", "goaldiffm", "totgoalm")) %>%
  filter(Season >= 1920) %>%
  ggplot(aes(x = Season, y = value, color = variable)) +
  geom_point() + geom_line() +
  ggtitle("Goal scoring trends over the season in English football") +
  theme_minimal()


# 2 - Poisson model ----

add_round <- function(df) {
  
  # Initialize team's game 
  teams <- unique(c(df$home, df$visitor)) 
  games <- as.list(setNames(rep(0, length(teams)), teams))
  
  df <- df  %>% arrange(Date) %>% 
    ddply(.(ID), function(ID) {
      HT <- ID$home
      AT <- ID$visitor
      
      games[[HT]] <<- games[[HT]] + 1
      games[[AT]] <<- games[[AT]] + 1
      
      ID$Round <- max(games[[HT]], games[[AT]])
      
      return(ID)
    }) %>% arrange(tier, Round, Date) %>%
    select(division, tier, Round, Date, everything())
  
  return (df)
}

model_expectancies <- function(model, home, visitor){
  
  home_goals_mean <- tryCatch(predict(model, data.frame(Home=1, TeamA = home, TeamB = visitor), type="response"), 
                              error = function(e) {warning(conditionMessage(e)); NA})
  
  away_goals_mean <- tryCatch(predict(model, data.frame(Home=0, TeamA = visitor, TeamB = home), type="response"), 
                              error = function(e) {warning(conditionMessage(e)); NA})
  
  return (c(home_goals_mean, away_goals_mean))
}

process_second_round <- function(df) {
  
  print(paste("Processing Season:", unique(df$Season), "Division", unique(df$division)))
  
  df <- add_round(df)
  
  second_round <- seq(max(df$Round)%/%2 + 1, max(df$Round) , 1)
  results <- vector("list", length(second_round))

  for (i in seq_along(second_round)) {
    
    # train poisson model with season results before current week
    x <- filter(df, Round < second_round[i])
    
    goals_data <-rbind(
      data.frame(TeamA = x$home,    TeamB = x$visitor, Goals = x$hgoal, Home = 1), 
      data.frame(TeamA = x$visitor, TeamB = x$home,    Goals = x$vgoal, Home = 0))
    
    poisson_model <- glm(Goals ~ Home + TeamA + TeamB, family=poisson(link="log"), data = goals_data)
    
    fixtures <- filter(df, Round == second_round[i]) 
    
    lambdas <- fixtures %>% dplyr::select(home, visitor) %>%
      pmap(model_expectancies, model = poisson_model)
    
    # Goal expectancy & Win Probability
    # skellam distribution (difference of independent poisson)
    # PPH, PPD, PPA : Probability Home/Draw/Away
    results[[i]] <- fixtures %>% 
      mutate(hgoal_exp = sapply(lambdas, `[[`, 1), 
             vgoal_exp = sapply(lambdas, `[[`, 2), 
             
             PPH = pskellam(0, hgoal_exp, vgoal_exp, lower.tail = F), 
             PPD = dskellam(0, hgoal_exp, vgoal_exp),
             PPA = 1 - PPH - PPD)
  }
  
  results <- bind_rows(results)
  
  return (results)
}

results_round2 <- football_data %>% 
  filter(between(Season, 1920, 2019)) %>%
  group_split(Season, division) %>%
  lapply(process_second_round) %>%
  bind_rows()

# 3 Model Evaluation ----

results_round2$decade <- as.factor(results_round2$decade)
color_scale <- colorRampPalette(c("yellow", "red"))(length(unique(results_round2$decade)))

actual_implied_probability <- function(df) {
  
  props <- df[ , c("decade", "tier", "result", "PPH", "PPD", "PPA")] %>% na.omit() %>% 
    setNames(c("decade", "tier", "result", "H", "D", "A")) %>%
    reshape2::melt(measure.vars = c("H", "D", "A"), variable.name = "Selection", value.name = "implied_probability") %>%
    mutate(evaluation = ifelse(result == Selection, 1, 0))
  
  props_summary <- props %>%
    mutate(implied_prop_int = cut(implied_probability, breaks = seq(0, 1, 0.025))) %>%
    group_by(decade, implied_prop_int) %>%
    summarise(obs = n(), 
              implied_prop = mean(implied_probability),
              correct = sum(evaluation), 
              actual_prop = correct/obs)
  
  pl <- ggplot(props_summary, aes(x = implied_prop, y = actual_prop, color = decade)) +
    geom_point() + geom_smooth(method = "lm", se = F, mapping = aes(weight = obs)) + 
    scale_color_manual(values = color_scale) +
    geom_abline(slope = 1, intercept = 0) + 
    ggtitle("Actual Probability - Implied Probability from basic Poisson Model")
  print(pl)
  
  return(props_summary)
}

prop_table_poisson <- actual_implied_probability(results_round2)

# We underestimate low probability events and overestimate high probability events

# plot fair returns  
results_round2 %>% 
  filter(!is.na(PPH)) %>%
  split(f = .$decade) %>%
  lapply(function(df) {
    df %>% mutate(fair_PnL = -3 + ifelse(result == "H", 1/PPH, ifelse(result == "D", 1/PPD, 1/PPA)),
                  NBets = seq(3, by = 3, length.out = nrow(.)), 
                  fair_PnL_sum = cumsum(fair_PnL), 
                  yield = fair_PnL_sum/NBets)}) %>%
  bind_rows(.id = "decade") %>%
  filter(NBets > 5000) %>% 
  ggplot(aes(x = NBets, y = yield, color = decade)) + 
  geom_line() + 
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_color_manual(values = color_scale) +
  labs(title = "Yield using fair probabilities of Poisson Model over the decades")


# supremacy & actual goal difference
results_round2 %>% filter(abs(hgoal_exp - vgoal_exp) < 3) %>%
  ggplot(aes(x = hgoal_exp - vgoal_exp, y = goaldif, color = decade)) +
  geom_point(alpha = 0.02) + geom_smooth(method = "lm", se = F) +
  scale_color_manual(values = color_scale) +
  geom_abline(slope = 1, intercept = 0) + 
  coord_cartesian(xlim = c(-1.5, 2.5), ylim = c(-3, 4)) +
  ggtitle("Goal difference ~ Poisson Supremacy") 

# square distance between expectancy & actual goals
results_round2 %>%
  filter(!is.na(hgoal_exp) & !is.na(vgoal_exp)) %>%
  group_by(Season) %>% 
  summarise(matches = n(), 
            dist = median((hgoal - hgoal_exp)^2 + (vgoal - vgoal_exp)^2)) %>%
  arrange(desc(dist)) %>%
  ggplot(aes(x = Season, y = dist)) + geom_point() + 
  geom_smooth(method = "lm", formula = y ~ poly(x, 8), mapping = aes(weight = matches), se = F) +
  ggtitle("Square distance of the Poisson Aproximation over the years")

