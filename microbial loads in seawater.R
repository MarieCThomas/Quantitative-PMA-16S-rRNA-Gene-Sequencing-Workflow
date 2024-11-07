
###################  Comparison of total microbial loads in seawater estimated by flow cytometry and droplet digital PCR  ###################################
# Marie Thomas

# Load libraries --------------------------------------------------------
library(phyloseq)
library(ggplot2)
library(dplyr)
library(vegan)
library(RVAideMemoire)
library(speedyseq)
library(wesanderson)
library(grid)
library(scales)
library(ggpubr)

# Data preparation ------------------------------------------------------
new_order <- c("Flow cytometry", "ddPCR")
new_order2 <- c("100", "80", "60", "40", "20", "0")

# Filter datasets for each method
ddpcr_data <- subset(count_data, method == "ddPCR")
fc_data <- subset(count_data, method == "Flow cytometry")
qubit_data <- subset(count_data, method == "Qubit")

# Flow cytometry plot ---------------------------------------------------
fc_data$Sample <- as.character(fc_data$Sample)

fc_plot <- ggplot(fc_data, aes(x = Sample, y = average, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = average - sd, ymax = average + sd), width = 0.2, position = position_dodge(0.9)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(margin = unit(c(6, 0, 0, 0), "mm")),
    axis.title.y = element_text(margin = unit(c(0, 6, 0, 0), "mm")),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey90"),
    panel.background = element_blank(),
    panel.border = element_rect(color = "grey80", size = 0.3),
    axis.line = element_line(colour = "grey80", size = 0.3),
    axis.ticks = element_line(color = "grey80", size = 0.3)
  ) +
  labs(x = expression("% natural seawater"), y = expression("Cells mL"^"-1")) +
  scale_y_continuous(labels = scientific_format()) +
  scale_fill_manual(values = c("Flow cytometry" = "#73a69c"))

fc_plot$data$Sample <- factor(fc_plot$data$Sample, levels = new_order2)

print(fc_plot)

# Qubit plot ------------------------------------------------------------
qubit_plot <- ggplot(qubit_data, aes(x = Sample, y = average, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = average - sd, ymax = average + sd), width = 0.2, position = position_dodge(0.9)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18),
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(margin = unit(c(6, 0, 0, 0), "mm")),
    axis.title.y = element_text(margin = unit(c(0, 6, 0, 0), "mm")),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey90"),
    panel.background = element_blank(),
    panel.border = element_rect(color = "grey80", size = 0.3),
    axis.line = element_line(colour = "grey80", size = 0.3),
    axis.ticks = element_line(color = "grey80", size = 0.3)
  ) +
  labs(x = expression("% natural seawater"), y = expression("DNA yields (ng μL"^{-1})) +
  scale_fill_manual(values = c("Qubit" = "#4892a8"))

qubit_plot$data$Sample <- factor(qubit_plot$data$Sample, levels = new_order2)

print(qubit_plot)

# ddPCR plot ------------------------------------------------------------
ddpcr_plot <- ggplot(ddpcr_data, aes(x = Sample, y = average, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = average - sd, ymax = average + sd), width = 0.2, position = position_dodge(0.9)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(margin = unit(c(6, 0, 0, 0), "mm")),
    axis.title.y = element_text(margin = unit(c(0, 6, 0, 0), "mm")),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey90"),
    panel.background = element_blank(),
    panel.border = element_rect(color = "grey80", size = 0.3),
    axis.line = element_line(colour = "grey80", size = 0.3),
    axis.ticks = element_line(color = "grey80", size = 0.3)
  ) +
  labs(x = expression("% natural seawater"), y = expression("16S rRNA gene copies mL"^"-1")) +
  scale_y_continuous(labels = scientific_format()) +
  scale_fill_manual(values = c("ddPCR" = "#a7cad4"))

ddpcr_plot$data$Sample <- factor(ddpcr_plot$data$Sample, levels = new_order2)

print(ddpcr_plot)

# Correlation plot ------------------------------------------------------
cor <- correlation_data %>%
  mutate(cells_log = log10(cells), copies_log = log10(copies), yields_log = log10(yields))

cor_plot <- ggplot(cor, aes(x = cells, y = copies)) + 
  geom_point(color = "#003366", size = 3, alpha = 0.6) +
  theme_bw() +
  theme(axis.title = element_text(size = 18), axis.text = element_text(size = 16)) +
  geom_smooth(method = "lm", se = TRUE, color = "#0099CC", linewidth = 0.5) +
  labs(x = "FC counts", y = "ddPCR counts") +
  scale_x_continuous(labels = scientific_format()) +
  scale_y_continuous(labels = scientific_format()) +
  stat_regline_equation(label.x = 10000, label.y = 639000) +
  stat_cor(aes(label = ..rr.label..), label.x = 10000, label.y = 620000)

model <- lm(copies ~ cells, data = cor)
summary(model)
confint(model, level = 0.95)

print(cor_plot)

# CV plot ---------------------------------------------------------------
cv_plot <- ggplot(CV_data, aes(x = method, y = CV, fill = method)) +
  stat_boxplot(geom = 'errorbar', width = 0.2) +
  geom_boxplot(width = 0.5) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 20),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16)
  ) +
  ylab("Coefficient of variation (%)") +
  xlab("Method") +
  scale_fill_manual(values = c("Flow cytometry" = "#73a69c", "ddPCR" = "#a7cad4"))

print(cv_plot)

# Save plots -----------------------------------------------------
plots.combined <- ggarrange(fc_plot, ddpcr_plot, cv_plot, cor_plot,
                             ncol = 2,
                             nrow = 2,   # Changed nrows to nrow
                             labels = c("A)", "B)", "C)", "D)"), 
                             label.x = 0, 
                             label.y = 1, 
                             hjust = 0, 
                             vjust = 1.5,  
                             # width = c(1, 1), # Uncomment if you want to adjust widths
                             font.label = list(size = 18, color = "black", face = "bold"))


ggsave(file = "Figure 2.jpg", plots.combined, width = 13, height = 11) 

# Perform regression analysis -----------------------------------------------------
ddPCR <- subset(preds_obs_data, method %in% c("ddPCR"))
FC <- subset(preds_obs_data, method %in% c("Flow cytometry"))

model <- lm(FC$predicted ~ FC$observed)
summary(model)
confint(model, level = 0.95)

plot_list <- list()

methods <- unique(data$method)

for (method in methods) {
  
  subset <- data %>% filter(method == !!method)
  
  p <- ggplot(subset, aes(x = predicted, y = observed)) +
    geom_point(size = 3, color = "#003366") +
    geom_smooth(method = 'lm', color = '#0099CC') +
    labs(title = paste(method, ": Predicted vs Observed"),
         x = "Predicted Values",
         y = "Observed Values") +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 18),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(margin = unit(c(6, 0, 0, 0), "mm")),
      axis.title.y = element_text(margin = unit(c(0, 6, 0, 0), "mm")),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey90"),
      panel.background = element_blank(),
      panel.border = element_rect(color = "grey80", size = 0.3),
      axis.line = element_line(colour = "grey80", size = 0.3),
      axis.ticks = element_line(color = "grey80", size = 0.3))
  
  plot_list[[method]] <- p
}


print(plot_list[["method_name"]])


for (method in methods) {
  print(plot_list[[method]])
}

FC <- plot_list$`Flow cytometry`
ddPCR <- plot_list$ddPCR

FC <- FC + stat_regline_equation(label.x=2, label.y=97) +
  stat_cor(aes(label=..rr.label..), label.x= 2, label.y=93)

ddPCR <- ddPCR + stat_regline_equation(label.x=2, label.y=99) +
  stat_cor(aes(label=..rr.label..), label.x= 2, label.y=95)

plots.combined2 <-  ggarrange(FC, ddPCR, labels = c("A)", "B)"),
                             ncol = 2, nrow= 1, font.label = list(size = 18, color = "black", face = "bold"))

ggsave(file = "Figure S2.jpg", plots.combined2, width = 12, height = 5) 


# Calculate Pearson correlation coefficient for each method  ------------------------------------------------
correlation_results <- sapply(methods, function(method) {
  subset <- data %>% filter(method == !!method)
  cor(subset$predicted, subset$observed)
})

correlation_results

# Calculate error metrics ------------------------------------------------
error_metrics <- sapply(methods, function(method) {
  subset <- data %>% filter(method == !!method)
  mae <- mean(abs(subset$observed - subset$predicted))
  rmse <- sqrt(mean((subset$observed - subset$predicted)^2))
  c(MAE = mae, RMSE = rmse)
})

error_metrics
