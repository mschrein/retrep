
# Setup ----------------------------------------------------------------

# load packages
library(lme4)
library(lmerTest)
library(truncnorm)
library(simr)
library(tidyverse)
library(combinat)
library(SimDesign)
library(parallel)

# set seed
seed <- 42
set.seed(seed)


# Previous Study Estimates ------------------------------------------------

# load data
load('data_bupse/final_data.RData') # available at https://doi.org/10.17605/OSF.IO/6Z5GE
df <- filter(data_long, QUESTNNR == "simul") # only select simultaneous presentation condition
df$CASE <- factor(as.numeric(factor(df$CASE))) # make CASE numbers continuously ascending and convert to factor

# MERing
mdl_bupse <- lmer(posterior_belief ~ 0 + offset(prior_belief) +
                    (0 + dist_evidence_prior | CASE) + (0 + dist_evidence_prior | hypothesis) +
                    dist_evidence_prior +
                    dist_evidence_prior:expertise +
                    dist_evidence_prior:TP +
                    dist_evidence_prior:SI +
                    dist_evidence_prior:prior_belief + 
                    dist_evidence_prior:evidence_pres_varyingVSincreasing + 
                    dist_evidence_prior:evidence_pres_decreasingVSvarying + 
                    dist_evidence_prior:prior_belief:evidence_pres_varyingVSincreasing + 
                    dist_evidence_prior:prior_belief:evidence_pres_decreasingVSvarying
                  , data = df
                  , control = lmerControl(optimizer = 'bobyqa')
)
summary(mdl_bupse)
FE_bupse <- fixef(mdl_bupse)
RE_bupse <- VarCorr(mdl_bupse)

# predicted response
df$posterior_belief_predicted <- predict(mdl_bupse)
hist(df$prior_belief, breaks=100)
hist(df$posterior_belief, breaks=100)
hist(df$posterior_belief_predicted, breaks=100)
sum(df$posterior_belief_predicted < 0 | df$posterior_belief_predicted > 1) / nrow(df) #share of unexplainable behavior

## truncation (despite sum-to-1 constraining necessary if there are further main effects)
df[df$posterior_belief_predicted < 0,]$posterior_belief_predicted <- 0
df[df$posterior_belief_predicted > 1,]$posterior_belief_predicted <- 1
hist(df$posterior_belief_predicted, breaks=100)

## use predicted response instead of actual (from a different design)
df$posterior_belief <- df$posterior_belief_predicted


# Data Preparation --------------------------------------------------------
# remove columns
df <- df %>% 
  select(-c(STARTED, starts_with('DE'), starts_with('SI'), starts_with('SX'), 
            starts_with('X'), starts_with('TP01'), evidence_ratio)) #remove unnecessary columns

# parameters
N = length(unique(df$CASE)) #total number of participants
M = 24 #total number of items
K = 1 #total number of studies per trial
WS = 4 #number of within-subjects conditions
BS = 1 #number of between-subjects conditions
N_er = 3 # number of each  evidence ratio per condition

### add six hypotheses (randomly duplicate entries of six hypotheses)
df <- df %>%
  group_by(CASE) %>%
  mutate(row_id = row_number()) %>%
  summarise(original = list(pick(everything())), .groups = "drop") %>%
  rowwise() %>%
  mutate(extended = list({
    original_df <- original
    sampled <- original_df %>% slice_sample(n = 6, replace = TRUE)
    
    # Assign new hypothesis names for sampled rows
    new_hypotheses <- paste0("h", 19:24)
    sampled$hypothesis <- new_hypotheses
    
    bind_rows(original_df, sampled)
  })) %>%
  ungroup() %>%
  select(CASE, extended) %>%
  unnest(extended)


N = length(unique(df$CASE))

## within-subjects
df <- df %>% 
  slice(rep(1:n(), each = K)) %>% 
  group_by(CASE, hypothesis) %>% 
  mutate(k = row_number())

### evidence presentation
evidence_ratio_set <- c(0,1)

df <- df %>% 
  group_by(CASE) %>%
  mutate(condition = rep(sample(c(rep('rep1_noretr', M/WS), rep('rep2_noretr', M/WS), 
                                  rep('rep1_retr', M/WS), rep('rep2_retr', M/WS))), each=K)) %>%
  mutate(repetition_condition = case_when( # repetition condition recoded to start with 0
    condition %in% c('rep1_noretr', 'rep1_retr') ~ 0,
    condition %in% c('rep2_noretr', 'rep2_retr') ~ 1),
    retraction_condition = case_when( # effect coded
      condition %in% c('rep1_noretr', 'rep2_noretr') ~ .5,
      condition %in% c('rep1_retr', 'rep2_retr') ~ -.5
    )) %>%
  group_by(CASE, condition) %>%
  mutate(evidence_ratio = rep(sample(rep(evidence_ratio_set, each=N_er)), each=K)) %>% 
  group_by(CASE, hypothesis) %>% 
  mutate(evidence_ratio_code = which(unique(evidence_ratio)==evidence_ratio_set)) %>% 
  ungroup() %>% 
  select(-c(evidence_ratio_code)) #remove helper variables

# balanced conditions?
df %>% 
  group_by(condition) %>% 
  summarise(n()/(M*K))

# randomly assign data quality values to df
df <- df %>%
  group_by(CASE) %>%
  mutate(quality_condition = sample(c(0, 1, 2, 3), size = M, replace = TRUE))

## covariates/predictors
#distance between evidence ratio and prior belief
df$dist_evidence_prior <- df$evidence_ratio - df$prior_belief

## add need for cognitive closure (same latent scale as trust)
df <- df %>%
  group_by(CASE) %>%
  mutate(NCC = sample(df$TP, size = 1))

# Intermediate judgments and (corresponding) null model
## null model
formula_null <- posterior_belief ~ 0 + offset(prior_belief) + 
  dist_evidence_prior +
  dist_evidence_prior:expertise +
  dist_evidence_prior:TP +
  dist_evidence_prior:NCC +
  (0 + dist_evidence_prior | CASE) + (0 + dist_evidence_prior | hypothesis)

mdl_null <- lmer(formula_null, data=df)
summary(mdl_null)
FE_null <- fixef(mdl_null)
RE_null <- VarCorr(mdl_null)

RE <- RE_bupse
re_hypothesis_orig <- RE$hypothesis[1,1]
RE$hypothesis[1,1] <- RE$hypothesis[1,1]*1.3 # increase random effect of hypothesis, as hypotheses will be more heterogeneous
attributes(RE$hypothesis)$stddev <- sqrt(RE$hypothesis[1,1])
attr(RE, "sc") <- attr(RE, "sc") - sqrt(RE$hypothesis[1,1] - re_hypothesis_orig) # reduce residual variance correspondingly

# Power Simulation --------------------------------------------------------

formula <- str2lang(paste0(deparse1(formula_null), 
                            '+ 
                            dist_evidence_prior:repetition_condition + 
                            dist_evidence_prior:retraction_condition + 
                            dist_evidence_prior:quality_condition + 
                            dist_evidence_prior:repetition_condition:retraction_condition +
                            dist_evidence_prior:retraction_condition:quality_condition +
                            dist_evidence_prior:repetition_condition:quality_condition +
                            dist_evidence_prior:repetition_condition:retraction_condition:quality_condition'))

test <- list('repetition' = 'dist_evidence_prior:repetition_condition',
             'retraction' = 'dist_evidence_prior:retraction_condition',
             'quality' = 'dist_evidence_prior:quality_condition',
             'repetitionxretraction' = 'dist_evidence_prior:repetition_condition:retraction_condition',
             'retractionxquality' = 'dist_evidence_prior:retraction_condition:quality_condition',
             'repetitionxquality' = 'dist_evidence_prior:repetition_condition:quality_condition')
test_expr <- list('repetition' = simr::fixed('dist_evidence_prior:repetition_condition', method='t'),
                  'retraction' = simr::fixed('dist_evidence_prior:repetition_condition', method='t'),
                  'quality' = simr::fixed('dist_evidence_prior:quality_condition', method='t'),
                  'repetitionxretraction' = simr::fixed('dist_evidence_prior:repetition_condition:retraction_condition', method='t'),
                  'retractionxquality' = simr::fixed('dist_evidence_prior:retraction_condition:quality_condition', method='t'),
                  'repetitionxquality' = simr::fixed('dist_evidence_prior:repetition_condition:quality_condition', method='t'))
             
eff_size_repetition = .06/2
eff_size_retraction = .06
eff_size_quality = .06/3
eff_size_repetitionxretraction = .06/2
eff_size_retractionxquality = .06/3
eff_size_repetitionxquality = .06/6
regression_to_mean_factor = 0.9
FE <- FE_bupse[names(FE_bupse) %in% names(FE_null)]*regression_to_mean_factor
FE["dist_evidence_prior:TP"] <- 0.04
FE["dist_evidence_prior:expertise"] <- -0.05
FE["dist_evidence_prior:NCC"] <- -0.04
FE_ <- c('dist_evidence_prior:repetition_condition' = eff_size_repetition,
         'dist_evidence_prior:retraction_condition' = eff_size_retraction,
         'dist_evidence_prior:quality_condition' = eff_size_quality,
         'dist_evidence_prior:repetition_condition:retraction_condition' = eff_size_repetitionxretraction,
         'dist_evidence_prior:retraction_condition:quality_condition' = eff_size_retractionxquality,
         'dist_evidence_prior:repetition_condition:quality_condition' = eff_size_repetitionxquality,
         'dist_evidence_prior:repetition_condition:retraction_condition:quality_condition' = 0
)

nsim <- 1000
N_ <- seq(2000, 2060, 10)
alpha <- .05
Design <- createDesign(N = N_, hypothesis = c('repetition', 'retraction', 'quality', 'repetitionxretraction', 'retractionxquality', 'repetitionxquality'))

fixed <- list(formula = formula, test = test, test_expr = test_expr, FE = FE, FE_ = FE_, RE = RE, df = df, alpha = alpha, M = M)

## Generate
Generate <- function(condition, fixed_objects) {
  
  # randomly assign data quality values to df
  df <- df %>%
    group_by(CASE) %>%
    mutate(quality_condition = sample(c(0, 1, 2, 3), size = M, replace = TRUE))
  
  mdl <- makeLmer(formula
           , fixef = c(FE, FE_)
           , VarCorr = RE
           , sigma = attr(RE, 'sc')
           , data = df
  )
  
  mdl_extended <- extend(mdl, along='CASE', n=condition$N)
  
  resp <- doSim(mdl_extended)
  
  return(list(resp = resp, mdl_extended = mdl_extended))
  
}

## Analyse
Analyse <- function(condition, dat, fixed_objects) {
  
  mod_fit <- doFit(dat$resp, dat$mdl_extended) # fit extended model to new responses
  mod_summary <- summary(mod_fit)
  est <- mod_summary$coefficients[test[[condition$hypothesis]], 'Estimate']
  p_ts <- as.numeric(doTest(mod_fit, test_expr[[condition$hypothesis]]))
  p_os <- ifelse(sign(est) == sign(c(FE, FE_)[test[[condition$hypothesis]]]), p_ts/2, 1-p_ts/2) # convert to one-tailed p value
  
  sign_ts <- ifelse(p_ts < alpha, 1, 0)
  sign_os <- ifelse(p_os < alpha, 1, 0)
  
  return(c(p_ts = p_ts, p_os = p_os, sign_ts = sign_ts, sign_os = sign_os))
  
}

Summarise <- function(condition, results, fixed_objects) {
  
  results <- data.frame(results)
  
  power_ts = mean(results$sign_ts)
  power_os = mean(results$sign_os)
  
  return(c(power_ts = power_ts, power_os = power_os))
  
}

cl <- makeCluster(detectCores()-1)
clusterExport(cl, c('formula', 'test', 'test_expr', 'FE', 'FE_', 'RE', 'df', 'alpha', 'M'))

sim_results <- runSimulation(design = Design, replications = nsim, generate = Generate, analyse = Analyse,
                             summarise = Summarise, fixed_objects = fixed, debug = "none", seed = 1:nrow(Design),
                             save = TRUE, filename = "power_simulation_exp2_save", save_results = TRUE,
                             parallel = TRUE, ncores = detectCores()-1, cl = cl,
                             packages = c("SimDesign", "simr", "parallel"),
                             control = list(warnings_as_errors = TRUE))

stopCluster(cl)

save.image(paste0('power_exp2-effect_', eff_size_repetition, "_", eff_size_retraction, "_", eff_size_quality, "_", eff_size_repetitionxretraction, "_", eff_size_retractionxquality, "_", eff_size_repetitionxquality, '-nsim_', nsim, '-breaks_', paste(c(min(N_), max(N_), (max(N_)-min(N_))/(length(N_)-1)), collapse = '_'), '.RData'))
