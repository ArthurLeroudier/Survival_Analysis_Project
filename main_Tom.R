###IMPORTS

library(KMsurv)
library(survival)
library(ggplot2)
library(survminer)
library(survival)
library(data.table)
library(joineRML)
library(tidyr)
library(dplyr)
library(purrr)
library(mice)
library(zoo)
library(MASS) 


#setup
set.seed(0)
data(pbc2)

# Setting up facors
pbc2$histologic <- factor(pbc2$histologic)
pbc2$edema <- factor(pbc2$edema)
pbc2$drug <- factor(pbc2$drug)
pbc2$sex <- factor(pbc2$sex)
pbc2$ascites <- factor(pbc2$ascites)
 
#--------
# PART A
#--------


##Missing Values
colSums(is.na(pbc2))
missing_row <- data.frame(id= pbc2$id, missing = rowSums(is.na(pbc2)))
missing_row <-missing_row[missing_row$missing>0,]
missing_row

#Which variables are missing for each moment and patient
# Build a data.frame row-by-row
missing_info <- data.frame(
  id      = pbc2$id,
  year    = pbc2$year,
  missing = rowSums(is.na(pbc2)),
  vars    = apply(pbc2, 1, function(x) {
    vars <- names(pbc2)[is.na(x)]
    if (length(vars)>0) paste(vars, collapse = ", ") else NA
  }),
  stringsAsFactors = FALSE
)

# Keep only the visits with at least one missing value
missing_info <- missing_info[missing_info$missing > 0, ]

# Inspect
print(missing_info)

# It appeears at the last observation there are more missing values than at other measurements
#Adding a comn the idicate if it is the last observation of a given patient
pbc2$last_obs <- with(pbc2, year == ave(year, id, FUN = max))

missing_info <- data.frame(
  id       = pbc2$id,
  year     = pbc2$year,
  status   = pbc2$status,              # add the patient’s status here
  missing  = rowSums(is.na(pbc2)),
  last_obs = pbc2$last_obs,
  vars     = apply(pbc2, 1, function(x) {
    m <- names(pbc2)[is.na(x)]
    if (length(m)) paste(m, collapse = ", ") else NA
  }),
  stringsAsFactors = FALSE
)

# Keep only the visits with at least one missing value
missing_info <- missing_info[missing_info$missing > 0, ]

# Inspect
print(missing_info)

# Calculating missing values for last visit 
total_visits <- nrow(missing_info)

# 2) number of those visits that are last_obs = TRUE
last_visits   <- sum(missing_info$last_obs)

# 3) number of those visits that are last_obs = FALSE
other_visits  <- total_visits - last_visits

# 4) total missing-fields across those groups
missing_by_group <- tapply(
  missing_info$missing,
  missing_info$last_obs,
  sum
)

# 5) average # of missing measurements per visit in each group
avg_missing_by_group <- tapply(
  missing_info$missing,
  missing_info$last_obs,
  mean
)

# Print a neat summary
cat("Visits with missing data:\n")
cat("  Last visit:    ", last_visits,    " rows,", 
    missing_by_group["TRUE"], "total missing fields, avg =", 
    round(avg_missing_by_group["TRUE"],2), "\n")
cat("  Other visits:  ", other_visits,   " rows,", 
    missing_by_group["FALSE"],"total missing fields, avg =", 
    round(avg_missing_by_group["FALSE"],2), "\n")

# among the visits with missing data, keep only the last observations
last_missing <- missing_info[missing_info$last_obs, ]

# table of patient status at those last visits
table(last_missing$status)
tab <- table(last_missing$status)
pr+6412+-op <- round(100 * tab / sum(tab), 1)
data.frame(
  status      = names(tab),
  n_missing   = as.integer(tab),
  pct_missing = as.numeric(prop)
)




##Constant variables
dif_values <- data.table(pbc2)[, lapply(.SD, function(x) length(unique(x))), by=id]
#^for each individual we check how many different values for each variable
check_constants <- colSums(dif_values[,!'id'])
check_constants #if the value is 312, the variable is constant for a given individual
#constant variables are: years, status, drug, age, sex, status2
pbc2_constants <- unique(pbc2[,c('years', 'status', 'drug', 'age', 'sex', 'status2')])
hist(pbc2$years)

#creating a dataset at a baseline (years==0)
pbc2_baseline = pbc2 %>% 
  filter(year == 0) %>%
  select(-year) 



#Barplot with counts of status of patients 
ggplot(pbc2_baseline, aes(x = status, fill = status)) +
  geom_bar() +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Distribution of Patient Status", x = "Status", y = "Count") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_minimal()

# Barplot of the distribution between placebo and treatment
ggplot(pbc2_baseline, aes(x = drug, fill = drug)) +
  geom_bar() +
  scale_fill_brewer(palette = "Paired") +
  labs(title = "Treatment Group Distribution", x = "Treatment", y = "Count") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_minimal()

ggplot(pbc2_baseline, aes(x = drug, fill = status)) +
  geom_bar() +
  scale_fill_brewer(palette = "Paired") +
  labs(
    title = "Treatment Group by Patient Status",
    x     = "Treatment",
    y     = "Count",
    fill  = "Final Status"
  ) +
  geom_text(
    stat    = "count",
    aes(label = ..count..),
    position = position_stack(vjust = 0.5),
    color    = "white",
    size     = 3
  ) +
  theme_minimal()

# Histogram of the distriburion of age
ggplot(pbc2_constants, aes(x = age)) +
  geom_histogram(binwidth = 2, fill = "skyblue", color = "black") +
  labs(title = "Age Distribution at Baseline", x = "Age (years)", y = "Number of Patients") +
  theme_minimal()

# Distribution of sex
ggplot(pbc2_constants, aes(x = sex, fill = sex)) +
  geom_bar() +
  scale_fill_manual(values = c("female" = "#E69F00", "male" = "#56B4E9")) +
  labs(title = "Sex Distribution", x = "Sex", y = "Count") +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_minimal()

# SPIDERS
plot_status_by_variable <- function(varname, data) {
  ggplot(data, aes_string(x = varname, fill = "factor(status2)")) +
    geom_bar(position = "fill") +
    scale_fill_manual(
      values = c("0" = "steelblue", "1" = "firebrick"),
      labels = c("Alive or Transplanted", "Dead")
    ) +
    labs(
      title = paste("Proportion of Deaths by", varname),
      x = varname,
      y = "Proportion",
      fill = "Final Outcome"
    ) +
    theme_minimal()
}

plot_status_by_variable("sex", pbc2_baseline)
plot_status_by_variable("drug", pbc2_baseline)
plot_status_by_variable("ascites", pbc2_baseline)
plot_status_by_variable("hepatomegaly", pbc2_baseline)
plot_status_by_variable("spiders", pbc2_baseline)
plot_status_by_variable("edema", pbc2_baseline)
plot_status_by_variable("histologic", pbc2_baseline)

# Plotting all variables
plots_part_a <- list(
  sex          = plot_status_by_variable("sex", pbc2_baseline),
  drug         = plot_status_by_variable("drug", pbc2_baseline),
  ascites      = plot_status_by_variable("ascites", pbc2_baseline),
  hepatomegaly = plot_status_by_variable("hepatomegaly", pbc2_baseline),
  spiders      = plot_status_by_variable("spiders", pbc2_baseline),
  edema        = plot_status_by_variable("edema", pbc2_baseline),
  histologic   = plot_status_by_variable("histologic", pbc2_baseline)
)


purrr::iwalk(
  plots_part_a,
  ~ ggsave(
    filename = paste0(.y, ".png"),
    plot     = .x,
    width    = 6,
    height   = 4,
    dpi      = 300
  )
)

## Quantitative variables
#Mean curve
plot_biomarker_by_status <- function(data, variable, y_label = NULL, title = NULL) {
  var_sym <- rlang::sym(variable)
  
  df_filtered <- data %>%
    filter(!is.na(!!var_sym))
  
  label_map <- df_filtered %>%
    distinct(id, status) %>%
    count(status) %>%
    mutate(label = paste0(status, " (n = ", n, ")")) %>%
    select(status, label)
  
  df_labeled <- left_join(data, label_map, by = "status")
  
  # Default title/y_label if not provided
  if (is.null(title))   title   <- paste(variable, "Over Time by Status")
  if (is.null(y_label)) y_label <- variable
  
  ggplot(df_labeled, aes(x = year, y = !!var_sym, group = id, color = label)) +
    geom_line(alpha = 0.2, na.rm = TRUE) +
    geom_point(alpha = 0.2, na.rm = TRUE) +
    geom_smooth(aes(group = label), method = "loess", se = TRUE, size = 1.2, na.rm = TRUE) +
    labs(
      title = title,
      x = "Years Since Enrollment",
      y = y_label,
      color = "Final Status"
    ) +
    theme_minimal()
}


plot_biomarker_by_status(pbc2, "serBilir", y_label = "Serum Bilirubin (mg/dL)", title = "Serum Bilirubin Over Time by Status")
plot_biomarker_by_status(pbc2, "serChol",  y_label = "Serum Cholesterol (mg/dL)")
plot_biomarker_by_status(pbc2, "alkaline", y_label = "Alkaline Phosphatase (U/L)")
plot_biomarker_by_status(pbc2, "SGOT",     y_label = "SGOT (U/mL)")
plot_biomarker_by_status(pbc2, "platelets", y_label = "Platelet Count (×1000/µL)")
plot_biomarker_by_status(pbc2, "prothrombin", y_label = "Prothrombin Time (s)")



#---------------------------
# Handling Missing values 
#---------------------------

# Total missing values in serChol
sum(is.na(pbc2$serChol))

# Proportion of missing serChol values
mean(is.na(pbc2$serChol))

# Examine how many of the missing values occur at the last observation
table(pbc2$last_obs[is.na(pbc2$serChol)])  # TRUE = last visit

# Proportion of serChol NAs that occur at the last visit
mean(pbc2$last_obs[is.na(pbc2$serChol)])

#Other variables

# Create table summarizing missingness across all variables (excluding serChol)
vars_to_check <- setdiff(names(pbc2), c("serChol", "id", "year", "last_obs"))

# For each variable, calculate proportion of NAs that are in last_obs == TRUE
prop_missing_last_obs <- sapply(vars_to_check, function(var) {
  is_missing <- is.na(pbc2[[var]])
  if (sum(is_missing) == 0) return(NA)  # skip if no missing
  mean(pbc2$last_obs[is_missing])
})

# Show only variables with some missing values
prop_missing_last_obs <- prop_missing_last_obs[!is.na(prop_missing_last_obs)]

# Sort descending to see which are most concentrated at the last visit
sort(prop_missing_last_obs, decreasing = TRUE)


# Step 1: Remove serChol
pbc2_locf <- pbc2 %>%
  select(-serChol)

# Step 2: Identify variables to impute (exclude id and year)
vars_to_impute <- setdiff(names(pbc2_locf), c("id", "year"))

# Step 3: Apply LOCF imputation grouped by patient ID
pbc2_locf <- pbc2_locf %>%
  arrange(id, year) %>%
  group_by(id) %>%
  mutate(across(all_of(vars_to_impute), ~ na.locf(.x, na.rm = FALSE))) %>%
  ungroup()

colSums(is.na(pbc2_locf))

sum(is.na(pbc2_locf$platelets)) 
pbc2_locf %>% filter(is.na(platelets)) %>% select(id, year)

# For these 4 observation a median imputation grou[ed by sex and drug group 
baseline <- pbc2_locf %>% filter(year == 0)

# Step 2: Compute median platelets per drug × sex group (excluding NAs)
median_by_group <- baseline %>%
  filter(!is.na(platelets)) %>%
  group_by(drug, sex) %>%
  summarise(median_platelets = median(platelets), .groups = "drop")

# Step 3: Join and impute missing baseline values
baseline_imputed <- baseline %>%
  left_join(median_by_group, by = c("drug", "sex")) %>%
  mutate(platelets = if_else(is.na(platelets), median_platelets, platelets)) %>%
  select(-median_platelets)

# Step 4: Replace the original baseline data in the full dataset
pbc2_filled <- pbc2_locf %>%
  filter(year > 0) %>%
  bind_rows(baseline_imputed) %>%
  arrange(id, year)

sum(is.na(pbc2_filled$platelets))  # should be 0


# ————————————————————————
# PART B 
# ———————————————————-----

# histologic
# 1. Kaplan–Meier survival curves
surv_histologic <- survfit(Surv(years, status2) ~ histologic, data = pbc2_filled, conf.type = "log-log")

# 2. Plot survival curves
ggsurvplot(
  surv_histologic,
  data = pbc2_filled,
  xlab = "Years since baseline",
  ylab = "Survival probability",
  legend.title = "Histologic Stage",
  palette = "Dark2",
  ggtheme = theme_minimal()
)

# 3. Median and quartile survival times
quantile(surv_histologic, probs = c(0.10, 0.25, 0.50, 0.75, 0.90))

# 4. Log-rank test
survdiff(Surv(years, status2) ~ histologic, data = pbc2_filled)


# --- PART D (Cox model) for histologic ---

# 5. Fit Cox proportional hazards model
cox_model_hist <- coxph(Surv(years, status2) ~ factor(histologic), data = pbc2_filled)

# 6. Model summary
summary(cox_model_hist)

# 7. Predict survival curves from the Cox model
newdata_hist <- data.frame(histologic = levels(factor(pbc2_filled$histologic)))  # ensure levels match Cox input
survfit_cox_hist <- survfit(cox_model_hist, newdata = newdata_hist)

# 8. Plot Cox-estimated survival curves
ggsurvplot(
  survfit_cox_hist,
  data         = pbc2_filled,
  legend.title = "Histologic Stage (Cox Model)",
  palette      = "Dark2",
  xlab         = "Years since baseline",
  ylab         = "Estimated Survival (Cox)",
  ggtheme      = theme_minimal()
)


# Drug
# a
surv_treatment <- survfit(Surv(years, status2) ~ drug, data = pbc2_filled, conf.type='log-log')
summary(surv_treatment)
ggsurvplot(surv_treatment,
           data         = pbc2_filled,
           xlab         = "Years since baseline",
           ylab         = "Survival probability",
           legend.title = "Treatment group",
           palette      = "Dark2",
           ggtheme      = theme_minimal())

print(surv_treatment, print.rmean =TRUE)

#b 
quantile(surv_treatment,probs =c(0.10,0.25,0.50,0.75,0.90))

#c
print(survdiff(Surv(years, status2) ~ drug, data = pbc2_filled))
 
#d
# Fit Cox model with 'drug' as the only covariate
cox_model_drug <- coxph(Surv(years, status2) ~ drug, data = pbc2_filled)

# Summary to see coefficients and hazard ratios
summary(cox_model_drug)

# Predict survival curves from the Cox model for each level of 'drug'
newdata_drug <- data.frame(drug = levels(pbc2_filled$drug))  # placebo, D-penicil
survfit_cox_drug <- survfit(cox_model_drug, newdata = newdata_drug)

# Plot estimated survival functions
ggsurvplot(
  survfit_cox_drug,
  data         = pbc2_filled,
  legend.title = "Treatment Group (Cox Model)",
  palette      = "Dark2",
  xlab         = "Years since baseline",
  ylab         = "Estimated Survival (Cox)",
  ggtheme      = theme_minimal()
)


#----------------------
# Part D 
#----------------------

# Ensure categorical variables are properly encoded
pbc2_filled$drug         <- factor(pbc2_filled$drug)
pbc2_filled$sex          <- factor(pbc2_filled$sex)
pbc2_filled$ascites      <- factor(pbc2_filled$ascites)
pbc2_filled$hepatomegaly <- factor(pbc2_filled$hepatomegaly)
pbc2_filled$spiders      <- factor(pbc2_filled$spiders)
pbc2_filled$edema        <- factor(pbc2_filled$edema)
pbc2_filled$histologic   <- factor(pbc2_filled$histologic)

# Fit the initial full Cox model with all covariates exept Chol
cox_full <- coxph(Surv(years, status2) ~ 
                    age + sex + drug + ascites + hepatomegaly + spiders +
                    edema + serBilir + albumin + alkaline + SGOT +
                    platelets + prothrombin + histologic, 
                  data = pbc2_filled)

cox_selected <- stepAIC(cox_full, direction = "backward", trace = TRUE)
# Summary of the final selected model
summary(cox_selected)


# Test proportional hazards assumption
ph_test <- cox.zph(cox_selected)
print(ph_test)

# Plot for visual inspection
plot(ph_test, var = "edemaedema despite diuretics ")  # example for one variable



# Fit parametric models with different distributions
aft_weibull <- survreg(Surv(years, status2) ~ age + sex + drug + ascites +
                         hepatomegaly + spiders + edema + serBilir +
                         albumin + alkaline + SGOT, 
                       data = pbc2_filled, dist = "weibull")

aft_exponential <- survreg(Surv(years, status2) ~ age + sex + drug + ascites +
                             hepatomegaly + spiders + edema + serBilir +
                             albumin + alkaline + SGOT, 
                           data = pbc2_filled, dist = "exponential")

aft_loglogistic <- survreg(Surv(years, status2) ~ age + sex + drug + ascites +
                             hepatomegaly + spiders + edema + serBilir +
                             albumin + alkaline + SGOT, 
                           data = pbc2_filled, dist = "loglogistic")

aft_lognormal <- survreg(Surv(years, status2) ~ age + sex + drug + ascites +
                           hepatomegaly + spiders + edema + serBilir +
                           albumin + alkaline + SGOT, 
                         data = pbc2_filled, dist = "lognormal")

AIC(aft_weibull, aft_exponential, aft_loglogistic, aft_lognormal,cox_selected)

summary(aft_loglogistic)


























