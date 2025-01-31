# Load necessary libraries
library(dplyr)
library(ggplot2)
library(latex2exp)
library(egg)
library(cowplot)    # For get_legend() and plot_grid()
library(patchwork)



# Read and group data
df_np <- read.csv("df_1.1_n_P_results.csv", stringsAsFactors = TRUE)
df_np <- df_np %>%
  group_by(n_P, Method) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    .groups = "drop" # Ensure the data is ungrouped after summarising
  )
df_np$Method <- factor(df_np$Method)
levels(df_np$Method) <- gsub("TLDP", "bottom", levels(df_np$Method))

df_np



df_nq0 <- read.csv("df_1.1_nq_0_results.csv", stringsAsFactors = TRUE)
df_nq0 <- df_nq0 %>%
  group_by(n_Q, Method) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    .groups = "drop" # Ensure the data is ungrouped after summarising
  )
df_nq0$Method <- factor(df_nq0$Method, labels = gsub("TLDP", "Target-Only TLDP", levels(df_nq0$Method)))

df_nq0


df_nq2 <- read.csv("df_1.1_nq_2_results.csv", stringsAsFactors = TRUE)
# Perform the summarisation as before
df_nq2 <- df_nq2 %>%
  group_by(n_Q, Method) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    .groups = "drop" # Ensure the data is ungrouped after summarising
  )

df_nq2

df_nq_com <- read.csv("df_1.1_2_results.csv", stringsAsFactors = TRUE)
df_nq_com <- df_nq_com %>%
  group_by(n_Q, Method) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    .groups = "drop" # Ensure the data is ungrouped after summarising
  )

# Combine the data frames
df_nq <- dplyr::bind_rows(df_nq0, df_nq2, df_nq_com_updated)



df_gamma <- read.csv("df_1.1_gamma_results.csv", stringsAsFactors = TRUE)
df_gamma <- df_gamma %>%
  group_by(gamma, Method) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    .groups = "drop" # Ensure the data is ungrouped after summarising
  ) 
df_gamma$Method <- factor(df_gamma$Method, labels = gsub("TLDP", "top", levels(df_gamma$Method)))

df_gamma

df_kappa <- read.csv("df_1.1_kappa_results.csv", stringsAsFactors = TRUE)
df_kappa <- df_kappa  %>%
  group_by(kappa, Method) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    .groups = "drop" # Ensure the data is ungrouped after summarising
  ) 
df_kappa$Method <- factor(df_kappa$Method, labels = gsub("TLDP", "bottom", levels(df_kappa$Method)))

df_kappa

df_C_I <- read.csv("df_1.1_C_I_results.csv", stringsAsFactors = TRUE)
df_C_I <- df_C_I %>%
  group_by(C_I, Method) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    .groups = "drop" # Ensure the data is ungrouped after summarising
  ) 
df_C_I$Method <- factor(df_C_I$Method, labels = gsub("TLDP", "top", levels(df_C_I$Method)))

df_C_I


df_C_r <- read.csv("df_1.1_C_r_results.csv", stringsAsFactors = TRUE)
df_C_r <- df_C_r %>%
  group_by(C_r, Method) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    .groups = "drop" # Ensure the data is ungrouped after summarising
  ) 
df_C_r$Method <- factor(df_C_r$Method, labels = gsub("TLDP", "bottom", levels(df_C_r$Method)))

df_C_r



all_levels <- c("bottom", "top", 
                "TLDP", "Target-Only TLDP", 
                "ABE", "ExUCB")




# Data Preparation

base_width <- 2000

# Plot 1: n_P vs mean
plot_1 <- ggplot(df_np, aes(x = n_P, y = mean, color = "bottom")) +
  geom_line() + 
  geom_point() + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = base_width) +
  labs(y = "Regret", x = TeX("$n_P$")) +
  scale_color_manual(values = c("bottom" = "#8da0cb")) +
  theme_classic()

x_range_1 <- max(df_np$n_P) - min(df_np$n_P)
x_range_2 <- max(df_nq$n_Q) - min(df_nq$n_Q)
adjusted_width <- base_width * (x_range_2 / x_range_1)

df_nq$Method <- factor(df_nq$Method, levels =  c("TLDP","Target-Only TLDP", "ABE", "ExUCB"))


color_list <- setNames(c("#66c2a5", "#e29578", "#A5A58D", "#598392"), levels(df_nq$Method))

plot_2 <- ggplot(df_nq, mapping = aes(x = n_Q, y = mean, colour = Method)) +
  geom_line() + 
  geom_point() + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = adjusted_width) +
  labs(y = "Regret", x = TeX("$n_Q$")) +
  scale_color_manual(
    values = color_list,
    labels = c("TLDP","Target-Only TLDP", "ABE", "ExUCB")
  ) + 
  theme_classic()

plot_2

x_range_2 <- max(df_kappa$kappa) - min(df_kappa$kappa)
adjusted_width <- base_width * (x_range_2 / x_range_1)


df_gamma2 <- df_gamma %>%
  mutate(kappa_mapped = 0.2 + 0.4*(gamma - 0.5))

plot_3 <-ggplot() +
  geom_line(data = df_kappa, aes(x = kappa, y = mean, color = "bottom" )) +
  geom_point(data = df_kappa, aes(x = kappa, y = mean, color = "bottom")) +
  geom_errorbar(
    data = df_kappa,
    aes(x = kappa, ymin = mean - sd, ymax = mean + sd, color = "bottom"),
    width = adjusted_width
  ) +
  # Plot the gamma data (mapped onto the kappa scale)
  geom_line(data = df_gamma2, aes(x = kappa_mapped, y = mean, color = "top")) +
  geom_point(data = df_gamma2, aes(x = kappa_mapped, y = mean, color = "top")) +
  geom_errorbar(
    data = df_gamma2,
    aes(x = kappa_mapped, ymin = mean - sd, ymax = mean + sd, color ="top"),
    width = adjusted_width
  ) +
  # Main (bottom) axis: kappa
  scale_color_manual(values = c("bottom" = "#8da0cb", "top" = "#DDB892")) +
  scale_x_continuous(
    name   = TeX("$\\kappa$"),
    limits = c(0.13, 1.07),
    # We want ticks at the bottom for kappa
    breaks = c(0.2, 0.4, 0.6, 0.8, 1),
    # Add a small "expand" so 0.2 doesn't get clipped
    expand = c(0, 0),
    
    # Secondary (top) axis: gamma
    sec.axis = sec_axis(
      ~ . , 
      name = TeX("$\\gamma$"),
      
      # We place ticks at those bottom‐axis (kappa) points
      breaks = c(0.2, 0.4, 0.6, 0.8, 1),
      
      # Then label them with the corresponding gamma values
      labels = c(0.5, 1.0, 1.5, 2.0, 2.5)
    )
  ) +
  ylab("Regret") + 
  theme_classic()



x_range_3 <- max(df_C_r$C_r) - min(df_C_r$C_r)
adjusted_width_1 <- base_width * (x_range_3 / x_range_1)


df_C_I2 <- df_C_I %>%
  mutate(C_r_mapped = 0.15 + 0.2*(C_I - 0.25))

plot_4 <-ggplot() +
  geom_line(data = df_C_r, aes(x = C_r, y = mean, color = "bottom")) +
  geom_point(data = df_C_r, aes(x = C_r, y = mean, color = "bottom")) +
  geom_errorbar(
    data = df_C_r,
    aes(x = C_r, ymin = mean - sd, ymax = mean + sd, color = "bottom"),
    width = adjusted_width_1
  ) +
  # Plot the gamma data (mapped onto the kappa scale)
  geom_line(data = df_C_I2, aes(x = C_r_mapped, y = mean, color = "top")) +
  geom_point(data = df_C_I2, aes(x = C_r_mapped, y = mean, color = "top")) +
  geom_errorbar(
    data = df_C_I2,
    aes(x = C_r_mapped, ymin = mean - sd, ymax = mean + sd, color = "top"),
    width = adjusted_width_1
  ) +
  scale_color_manual(values = c("bottom" = "#8da0cb", "top" = "#DDB892")) +
  scale_x_continuous(
    name   = TeX("$C_r$"),
    limits = c(0.135, 0.365),
    # We want ticks at the bottom for kappa
    breaks = c(0.15, 0.2, 0.25, 0.3, 0.35),
    # Add a small "expand" so 0.2 doesn't get clipped
    expand = c(0, 0),
    
    # Secondary (top) axis: gamma
    sec.axis = sec_axis(
      ~ . , 
      name = TeX("$C_I$"),
      
      # We place ticks at those bottom‐axis (kappa) points
      breaks = c(0.15, 0.2, 0.25, 0.3, 0.35),
      
      # Then label them with the corresponding gamma values
      labels = c(0.25, 0.5, 0.75, 1, 1.25)
    )
  ) +
  ylab("Regret") + 
  theme_classic()





plot_1_no_legend <- plot_1 + theme(legend.position = "none")
plot_2_no_legend <- plot_2 + theme(legend.position = "none")
plot_3_no_legend <- plot_3 + theme(legend.position = "none")
plot_4_no_legend <- plot_4 + theme(legend.position = "none")



# Define separate color mappings for upper and lower sections
upper_methods <- c("TLDP", "Target-Only TLDP", "ABE", "ExUCB")
lower_methods <- c("TLDP (varies w/ bottom axis)", "TLDP (varies w/ top axis)")

upper_colors <- c(
  "TLDP" = "#66c2a5",
  "Target-Only TLDP" = "#e29578",
  "ABE" = "#A5A58D",
  "ExUCB" = "#598392"
)



lower_colors <- c(
  "TLDP (varies w/ bottom axis)" = "#8da0cb",
  "TLDP (varies w/ top axis)" = "#DDB892"
)

# Create dummy data for upper and lower legends
dummy_df_upper <- data.frame(
  x = seq_along(upper_methods),
  y = seq_along(upper_methods),
  Method = factor(upper_methods, levels = upper_methods)
)

dummy_df_lower <- data.frame(
  x = seq_along(lower_methods),
  y = seq_along(lower_methods),
  Method = factor(lower_methods, levels = lower_methods)
)

# Generate upper legend plot (focus on color guide only)
upper_legend_plot <- ggplot(dummy_df_upper, aes(x = x, y = y, color = Method)) +
  geom_point() +
  scale_color_manual(values = upper_colors) +
  theme_minimal() +
  theme(legend.position = "right", legend.title = element_text(size = 10)) +
  guides(color = guide_legend(title = "Method (Panel A)"), size = "none", shape = "none")

# Generate lower legend plot (focus on color guide only)
lower_legend_plot <- ggplot(dummy_df_lower, aes(x = x, y = y, color = Method)) +
  geom_point() +
  scale_color_manual(values = lower_colors) +
  theme_minimal() +
  theme(legend.position = "right", legend.title = element_text(size = 10)) +
  guides(color = guide_legend(title = "Method (Panels B, C & D)"), size = "none", shape = "none")

# Extract the legends
upper_legend_grob <- cowplot::get_legend(upper_legend_plot)
lower_legend_grob <- cowplot::get_legend(lower_legend_plot)

# Arrange legends vertically
right_legends <- plot_grid(
  upper_legend_grob, 
  lower_legend_grob, 
  ncol = 1, # Stack legends vertically
  align = "v",
  rel_heights = c(1, 1) # Equal space for both legends
)

# Combine plots without legends
combined_plots_no_legend <- plot_grid(
  plot_grid(plot_2_no_legend, plot_1_no_legend, ncol = 2, labels = c("A", "B")),
  plot_grid(plot_3_no_legend, plot_4_no_legend, ncol = 2, labels = c("C", "D")),
  nrow = 2,
  rel_heights = c(1, 1.2)
)

# Combine the plots with the legends aligned on the right
final_combined_plot <- plot_grid(
  combined_plots_no_legend, 
  right_legends, 
  ncol = 2, # Plots on the left, legends on the right
  rel_widths = c(3, 1) # Adjust the width ratio
)

# Display the final plot
final_combined_plot
