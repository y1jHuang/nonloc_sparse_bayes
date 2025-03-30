library(tidyverse)
library(ggplot2)
load("../results/nutAnaly_trueK.Rdata")
colnames(eta_hat) <- paste0("Fac", 1:ncol(eta_hat))
eta_hat[, -4] <- -eta_hat[, -4]
Lambda_hat[, -4] <- -Lambda_hat[, -4]
book <- readxl::read_excel("../data/Nutrients_diff.xlsx", sheet = 3)

theta <- apply(eta_hat != 0, 2, sum) / nrow(eta_hat)
varimax.info <- varimax(Lambda_hat[, ])
Lambda.varimax <- varimax.info$loadings
# Lambda.varimax <- Lambda_hat
idx_fac <- apply(abs(Lambda.varimax), 1, 
                 order, decreasing=T) %>% t

idx <- data.frame(id=1:nrow(idx_fac), idx_fac) %>%
  arrange(across(X1:X6), decreasing=T)
Lambda.sort <- Lambda.varimax[idx$id, ]

nutrient.name <- book$Label[idx$id]
# fac.name <- c("Plant-based product", 
#               "Animal \n & \n vegetarian food", 
#               "Seafood", "Dairy food", 
#               "Animal-based product", "Antixoxidant")
fac.name <- c("Plant-based \n & \n whole grain", 
              "Ultra_processed \n & \n industrialized food", 
              "Seafood", "Dairy \n & \n processed food", 
              "Animal product", "Plant Antixoxidant")
df_Lambda <- data.frame(nutrients=as.factor(1:nrow(Lambda.sort)), 
                        -Lambda.sort) %>% reshape2::melt()

### Plot loading matrix
plt <- ggplot(df_Lambda, aes(x = variable, y = nutrients, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0, limits = c(-1, 1), name = "Loading",
                       oob = scales::squish) + # Custom gradient
  labs(x = "", y = "", fill = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom", 
        legend.direction = "horizontal") +
  scale_y_discrete(labels = nutrient.name) +
  scale_x_discrete(labels = fac.name)
ggsave("../results/plot/loadings_k6.svg", width = 10, height = 10)

### Regression analysis
library(rstanarm)
covar <- read.csv("../data/part_derv_lad1.csv") %>%
  filter(MED_ANTIHYPERT == 1 | MED_ANTIDIAB == 1) %>%
  left_join(read.csv("../data/anta_lad1.csv"), by = "PID") %>%
  mutate(BMI = ANTA4 / (HEIGHT/100)^2,
         OBESITY = if_else(BMI >= 30, 1, 0))

factor <- paste0("Fac", 1:6)
var <- c("PAG2008YN", "CESD10", "EMPLOYED", "YRSUS",
         "MARITAL_STATUS", "INCOME_C3", "EDUCATION_C3",
         "DTIA20", "BKGRD1_C7", "GENDERNUM", "CENTER", factor)
coef <- data.frame()
for (outcome in c("HYPERTENSION", "DIABETES_SELF", "HIGH_TOTAL_CHOL", "OBESITY")) {
  formula <- reformulate(var, response = outcome)
  df <- read.csv("../../data/nutrients/dtia.csv") %>%
    mutate(SCSFA = DTIA69,
           MCSFA = rowSums(across(DTIA70:DTIA73), na.rm = TRUE),
           LCSFA = rowSums(across(DTIA74:DTIA78), na.rm = TRUE),
           LCMFA = rowSums(across(DTIA80:DTIA83), na.rm = TRUE)) %>%
    filter(if_all(book$Variable, ~ . >= 0 | is.na(.))) %>%
    inner_join(covar, by = "PID") %>%
    # cbind(data.frame(eta_hat %*% varimax.info$rotmat) %>%
    #         rename_with(~ paste0("Fac", seq_along(.)))) %>%
    cbind(eta_hat) %>%
    left_join(read.csv("../data/pa_derv_lad1.csv"), by = "PID") %>%
    select(all_of(c(var, outcome)))
  
    set.seed(1)
    t_prior <- student_t(df = 7, location = 0, scale = 2.5)
    fit <- stan_glm(formula, 
                    prior = t_prior,
                    data = df,
                    family = binomial(link = "logit"))
                    # family = gaussian())
    ci <- posterior_interval(fit, prob = 0.95)
    coef <- data.frame(OR = exp(coef(fit)), exp(ci)) %>%
      filter(str_detect(row.names(.), "Fac")) %>%
      mutate(outcome = outcome,
             sig = if_else((X2.5. > 1 & X97.5. > 1) | (X2.5. < 1 & X97.5. < 1), 1, 0),
             factor = 1:n()) %>%
      rbind(coef)
}

plt <- ggplot(coef, aes(x = OR, y = factor(factor), color = outcome)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = X2.5., xmax = X97.5.),
                 position = position_dodge(width = 0.5)) +
  labs(x = "Odds Ratio (95% CI)", y = NULL) +
  scale_y_discrete(labels = fac.name) +
  scale_color_discrete(labels = c("DIABETES_SELF" = "Diabetes",
                                  "HIGH_TOTAL_CHOL" = "High cholesterol",
                                  "HYPERTENSION" = "Hypertension",
                                  "OBESITY" = "Obesity")) +
  geom_vline(xintercept = 1, 
             linetype = "dashed", 
             linewidth = 1,
             color = "gray") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom", 
        legend.direction = "horizontal")
ggsave("../results/plot/OR_4diseases.svg", width = 8, height = 4)
saveRDS(coef, "../results/plot/OR_4diseases.RDS")