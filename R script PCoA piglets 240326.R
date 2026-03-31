# 0. Load libraries
# ------------------------------
library(vegan)
library(ggplot2)
library(ggrepel)
library(grid)

# ------------------------------
# 1. Load your data
# ------------------------------
data_path <- "C:/Users/gaute/OneDrive/Dokumenter/Postdoc 2026-2029/Targeted SCFA/Piglets.csv"
df <- read.csv(data_path, header = TRUE, check.names = FALSE, row.names = 1)

# Separate metadata and metabolites
groups <- df[,1]
metab_df <- df[, -1]

# ------------------------------
# Optional: custom group labels
# ------------------------------
use_custom_labels <- FALSE  # TRUE or FALSE

custom_labels <- c(
  "Germ free" = "Germ-free",
  "Mono colonized" = "Mono-colonized")

# ------------------------------
# 2. Filter metabolites (optional)
# ------------------------------
keep <- colMeans(is.na(metab_df)) < 0.5
metab_df <- metab_df[, keep]

# ------------------------------
# 3. Half-min imputation
# ------------------------------
metab_imputed <- metab_df
for (i in 1:ncol(metab_imputed)) {
  min_val <- min(metab_imputed[metab_imputed[,i] > 0, i], na.rm = TRUE)
  metab_imputed[is.na(metab_imputed[,i]) | metab_imputed[,i] == 0, i] <- min_val / 2
}

# ------------------------------
# 4. Log transform
# ------------------------------
metab_log <- log10(metab_imputed)

# ------------------------------
# 5. Distance + PCoA
# ------------------------------
distm <- vegdist(metab_log, method = "canberra")
pcoa_res <- cmdscale(distm, k = 2, eig = TRUE)

pcoa_coords <- as.data.frame(pcoa_res$points)
colnames(pcoa_coords) <- c("PCoA1", "PCoA2")
variance <- round(pcoa_res$eig[1:2] * 100 / sum(pcoa_res$eig), 1)

pcoa_coords$Group <- groups

# ------------------------------
# 6. PERMANOVA
# ------------------------------
adon_res <- adonis2(distm ~ Group, data = pcoa_coords)
r2 <- round(adon_res$R2[1] * 100, 1)
pval <- signif(adon_res$`Pr(>F)`[1], 3)

# ------------------------------
# 7. envfit vectors
# ------------------------------
fit <- envfit(pcoa_res, metab_log, permutations = 999)

vecs_scaled <- as.data.frame(fit$vectors$arrows * sqrt(fit$vectors$r))
colnames(vecs_scaled) <- c("PCoA1", "PCoA2")

vecs_scaled$Metabolite <- rownames(fit$vectors$arrows)
vecs_scaled$pval <- fit$vectors$pvals
vecs_scaled$signif <- vecs_scaled$pval < 0.05

# ------------------------------
# 8. Manual arrow scaling
# ------------------------------
arrow_manual_scale <- 0.4  # <-- adjust manually to change arrow length
vecs_scaled$PCoA1 <- vecs_scaled$PCoA1 * arrow_manual_scale
vecs_scaled$PCoA2 <- vecs_scaled$PCoA2 * arrow_manual_scale

# ------------------------------
# 9. Plot
# ------------------------------
p <- ggplot(pcoa_coords, aes(x = PCoA1, y = PCoA2, color = Group, fill = Group)) +
  
  # Points
  geom_point(size = 4, alpha = 0.9) +
  
  # Ellipses
  stat_ellipse(geom = "polygon", alpha = 0.2, color = NA) +
  
  # Biplot arrows
  geom_segment(data = vecs_scaled,
               aes(x = 0, y = 0, xend = PCoA1, yend = PCoA2),
               arrow = arrow(length = unit(0.3, "cm")),
               linetype = ifelse(vecs_scaled$signif, "solid", "dashed"),
               size = 0.8,
               color = "black",
               inherit.aes = FALSE,
               show.legend = FALSE) +
  
  # Metabolite labels (ggrepel)
  geom_text_repel(data = vecs_scaled,
                  aes(x = PCoA1, y = PCoA2, label = Metabolite,
                      fontface = ifelse(signif, "bold", "plain")),
                  size = 4,
                  segment.color = "black",
                  segment.size = 0.5,
                  box.padding = 0.5,
                  point.padding = 0.5,
                  max.overlaps = Inf,
                  inherit.aes = FALSE,
                  nudge_x = vecs_scaled$PCoA1 * 0.05,
                  nudge_y = vecs_scaled$PCoA2 * 0.1) +
  
  # Colors + optional custom labels
  {
    if (use_custom_labels) {
      list(
        scale_color_manual(values = c("#E69F00", "#009E73"),
                           labels = custom_labels),
        scale_fill_manual(values = c("#E69F00", "#009E73"),
                          labels = custom_labels)
      )
    } else {
      list(
        scale_color_manual(values = c("#E69F00", "#009E73")),
        scale_fill_manual(values = c("#E69F00", "#009E73"))
      )
    }
  } +
  
  # Axis labels
  xlab(bquote("PCoA1 ("*.(variance[1])*"%)")) +
  ylab(bquote("PCoA2 ("*.(variance[2])*"%)")) +
  
  # Theme (white background)
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.line = element_line(color = "black"),
    legend.text = element_text(size = 12),
    legend.title = element_blank()
  ) +
  
  # Bigger legend points
  guides(color = guide_legend(override.aes = list(size = 6)),
         fill  = guide_legend(override.aes = list(size = 6)))

# ------------------------------
# PERMANOVA annotation - top right
# ------------------------------
x_max <- max(pcoa_coords$PCoA1)
y_max <- max(pcoa_coords$PCoA2)

p <- p + annotate("text",
                  x = x_max * 2.5,       # adjust horizontally
                  y = y_max * 2.5,       # adjust vertically
                  label = paste0("\nR² = ", r2, "%\np = ", pval),
                  size = 5,               # font size
                  hjust = 1,
                  vjust = 1)

# ------------------------------
# Show plot
# ------------------------------
print(p)

# ------------------------------
# 10. Save plot
# ------------------------------
ggsave("PCoA_piglets_final.png", p, width = 10, height = 8, dpi = 300)