library(tidyverse)
library(ggplot2)
file_path <- list.files("../results", 
                        full.names=TRUE,
                        pattern = "estK")

results = theta = k_estMat <- data.frame()
# Combine metrics from different scenarios
for (file in file_path) {
  load(file)
  results <- data.frame(RV_eta, RV_Lambda) %>%
    pivot_longer(cols = everything(),
                 names_to = "metric",
                 values_to = "val") %>%
    mutate(scen = str_extract(file, "scenario[^_]*"),
           # step = str_extract(file, "(?<=step)\\d+(\\.\\d+)?"),
           model = if_else(grepl("BFMAN", file, ignore.case = TRUE),
                           "BFMAN", "MGPS")) %>%
    rbind(results)
  if (grepl("BFMAN", file, ignore.case = TRUE)) {
    num_rep <- length(RV_eta)
    theta <- data.frame(theta_est,
                        fac = rep(1:(length(theta_est) / num_rep), num_rep),
                        # step = str_extract(file, "(?<=step)\\d+(\\.\\d+)?"),
                        scen = str_extract(file, "scenario[^_]*")) %>%
      rbind(theta)
  } else {
    k_estMat <- data.frame(n_fac = k_est,
                           scen = str_extract(file, "scenario[^_]*")) %>%
      rbind(k_estMat)
  }
}

# Boxplot for RV coefficients under different scenarios by BFMAN & MGPS
plt <- results %>%
  mutate(scen = scen %>%
           str_replace_all("scenario([0-9])",
                           function(x) paste0("Scenario ", str_extract(x, "[0-9]")
                           ))) %>%
  # filter(step==1|is.na(step)) %>%
  ggplot() +
    geom_boxplot(aes(x = metric, y = val, fill = model, alpha=metric),
                 linewidth = 1) +
    facet_wrap(~ scen, nrow = 1) + # Facet by the `scenario` variable
    scale_x_discrete(labels = c("RV_eta" = expression(eta~eta^{~T}), 
                                "RV_Lambda" = expression(Lambda~Lambda^{~T})
                                )) +
    scale_alpha_discrete(range = c(1, 0.4)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12),
          plot.title = element_text(hjust = 0.5),
          legend.position = "bottom",
          panel.grid.major = element_line(color = "grey", linewidth = 0.5), # Set major grid lines to grey
          panel.grid.minor = element_blank(), # Set minor grid lines to lighter grey
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background = element_rect(fill = "transparent", color = NA)) + # Adjust x-axis labels
    labs(x = "",
         y = "RV",
         fill = "Model") +
    guides(alpha = "none")
ggsave("../results/plot/RV_scens_estK.svg", plot = plt, width = 6, height = 4)

# Histogram for the number of factors estimated by MGPS
plt <- k_estMat %>%
  mutate(scen = scen %>%
           str_replace_all("scenario([0-9])",
                           function(x) paste0("Scenario ", str_extract(x, "[0-9]")
                           ))) %>%
  ggplot() +
  geom_histogram(aes(x = n_fac),
                 binwidth = 0.5) +
  facet_wrap(~ scen, nrow = 1) + # Facet by the `scenario` variable
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_line(color = "grey", linewidth = 0.5), # Set major grid lines to grey
        panel.grid.minor = element_blank(), # Set minor grid lines to lighter grey
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA)) + # Adjust x-axis labels
  labs(x = "Number of Factors Estimated by MGPS",
       y = "Frequency")
ggsave("../results/simul/MGPS_estK_start.svg", plot = plt, width = 6, height = 3)

# Boxplot for theta under different scenarios by BFMAN
theta <- theta %>%
  mutate(fac,
         scen = scen %>%
           str_replace_all("scenario([0-9])",
                           function(x) paste0("Scenario ", str_extract(x, "[0-9]")
                                              )))
plt <- ggplot() +
  geom_boxplot(data = theta,
               aes(x = factor(fac), y = theta_est),
               linewidth = 1) +
  facet_wrap(vars(scen), nrow = 1, scales = "free_x") + 
  geom_segment(data = data.frame(scen="Scenario 1"), 
               aes(x = 0.5, y = 0.4, xend = 3.5, yend = 0.4), 
               color = "#FF4040", linetype = "solid", linewidth = 1.5) +
  geom_segment(data = data.frame(scen="Scenario 2"), 
               aes(x = 0.5, y = 0.8, xend = 1.5, yend = 0.8), 
               color = "#FF4040", linetype = "solid", linewidth = 1.5) +
  geom_segment(data = data.frame(scen="Scenario 2"), 
               aes(x = 1.5, y = 0.6, xend = 2.5, yend = 0.6), 
               color = "#FF4040", linetype = "solid", linewidth = 1.5) +
  geom_segment(data = data.frame(scen="Scenario 2"), 
               aes(x = 2.5, y = 0.4, xend = 3.5, yend = 0.4), 
               color = "#FF4040", linetype = "solid", linewidth = 1.5) +
  geom_segment(data = data.frame(scen="Scenario 3"), 
               aes(x = 0.5, y = 0.9, xend = 1.5, yend = 0.9), 
               color = "#FF4040", linetype = "solid", linewidth = 1.5) +
  geom_segment(data = data.frame(scen="Scenario 3"), 
               aes(x = 1.5, y = 0.8, xend = 2.5, yend = 0.8), 
               color = "#FF4040", linetype = "solid", linewidth = 1.5) +
  geom_segment(data = data.frame(scen="Scenario 3"), 
               aes(x = 2.5, y = 0.7, xend = 3.5, yend = 0.7), 
               color = "#FF4040", linetype = "solid", linewidth = 1.5) +
  geom_segment(data = data.frame(scen="Scenario 3"), 
               aes(x = 3.5, y = 0.6, xend = 4.5, yend = 0.6), 
               color = "#FF4040", linetype = "solid", linewidth = 1.5) +
  geom_segment(data = data.frame(scen="Scenario 3"), 
               aes(x = 4.5, y = 0.5, xend = 5.5, yend = 0.5), 
               color = "#FF4040", linetype = "solid", linewidth = 1.5) +
  labs(x = "Factors", y = expression(theta)) +
  scale_y_continuous(breaks = seq(0.8, 0, -0.2)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_line(color = "grey", linewidth = 0.5), # Set major grid lines to grey
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA))
ggsave("../results/plot/theta_scens_estK.svg", plot = plt, width = 6, height = 3)