# ===============================================================
# COMPLETE R SCRIPT: TARGETED SCFA BOX PLOTS WITH STATISTICS
# ===============================================================

# -----------------------------
# 1. Load packages
# -----------------------------
install.packages(c("tidyverse", "rstatix", "multcompView"))  # run once if needed
library(tidyverse)
library(rstatix)
library(multcompView)

# -----------------------------
# 2. Load data
# -----------------------------
file_path <- "your_path.csv"
data <- read.csv(file_path)

# -----------------------------
# 3. Select metabolite
# -----------------------------
metabolite_name <- "Acetate"   # change to any metabolite column
plot_title <- "Acetate"

df <- data %>%
  select(SampleID = 1, Group = 2, Value = all_of(metabolite_name)) %>%
  mutate(logValue = log10(Value + 1e-6))

# -----------------------------
# 4. Define correct group order
# -----------------------------
order_levels <- c(
  "Germ-free",
  "Mono-colonized",
  "OMM12",
  "SPF")

df$Group <- factor(df$Group, levels = order_levels)

# -----------------------------
# 5. Kruskal-Wallis test
# -----------------------------
kw <- kruskal_test(df, logValue ~ Group)
print(kw)

# -----------------------------
# 6. Dunns test
# -----------------------------
pairwise <- df %>%
  dunn_test(logValue ~ Group, p.adjust.method = "none") %>%
  mutate(
    p.adj.BH = p.adjust(p, method = "BH"),
    p.adj.BY = p.adjust(p, method = "BY")
  )

# -----------------------------
# 7. CLD (FIXED, internal hyphen handling only)
# -----------------------------
pairwise <- pairwise %>%
  mutate(
    group1_safe = gsub("-", "~", group1),
    group2_safe = gsub("-", "~", group2)
  )

pvals <- pairwise$p.adj.BH
names(pvals) <- paste(pairwise$group1_safe, pairwise$group2_safe, sep = "-")

cld <- multcompLetters(pvals)

cld_letters <- cld$Letters
names(cld_letters) <- gsub("~", "-", names(cld_letters))

cld_df <- data.frame(
  Group = names(cld_letters),
  Letters = cld_letters
)

# -----------------------------
# 8. Merge CLD into data
# -----------------------------
plot_data <- df %>%
  left_join(cld_df, by = "Group")

plot_data$Group <- factor(plot_data$Group, levels = order_levels)

# -----------------------------
# 9. Plot
# -----------------------------
ggplot(plot_data, aes(x = Group, y = logValue, fill = Group)) +
  
  geom_boxplot(outlier.shape = NA, linewidth = 1.2) +
  
  geom_jitter(color = "black", width = 0.15, size = 5, alpha = 0.8) +
  
  geom_text(
    data = plot_data %>%
      group_by(Group) %>%
      summarise(y = max(logValue, na.rm = TRUE) + 0.5, Letters = unique(Letters)),
    aes(x = Group, y = y, label = Letters),
    size = 8,
    fontface = "bold",
    inherit.aes = FALSE
  ) +
  
  scale_fill_manual(values = c(
    "Germ-free" = "#E69F00",
    "Mono-colonized)" = "#009E73",
    "OMM12" = "#0072B2",
    "SPF" = "#D55E00")) +
  
  labs(
    title = plot_title,
    y = expression(log[10]~"(Concentration (μmol/g))"),
    x = NULL
  ) +
  
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 1.2, color = "black"),
    axis.text.x = element_text(size = 14, angle = 30, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.position = "none"
  )

