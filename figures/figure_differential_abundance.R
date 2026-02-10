################################################################################
# Publication Figure: Differential Abundance Visualization
# Description: Creates publication-ready figures for ALDEx2 differential
#              abundance results with functional annotations
################################################################################

# Load required libraries
library(ggplot2)
library(dplyr)
library(cowplot)
library(grid)

################################################################################
# LOAD DATA
################################################################################

# Load ALDEx2 results with taxonomy and annotations
# This assumes you have added functional annotations to the significant taxa
# (e.g., anaerobic, biofilm former, antimicrobial tolerant)

# For demonstration, loading from workspace
# In practice, load from saved CSV file:
# enrich <- read.csv("significant_da_taxa_with_annotations.csv")

# Example structure expected:
# - Species: Species name
# - effect: ALDEx2 effect size
# - Anaerobic: "Yes"/"No"
# - Biofilm.former: "Yes"/"No"
# - Antimicrobial.tolerant: "Yes"/"No"

################################################################################
# DATA PREPARATION
################################################################################

# Create functional annotation summary
df <- enrich %>%
  mutate(
    # Create functional group categories
    Functional_group = case_when(
      Anaerobic == "Yes" & Biofilm.former == "Yes" ~ "Anaerobe + Biofilm",
      Anaerobic == "Yes" ~ "Anaerobe",
      Biofilm.former == "Yes" ~ "Biofilm former",
      TRUE ~ "Other"
    ),
    # Enrichment direction
    Direction = ifelse(effect < 0, "Sink A", "Sink B"),
    # Absolute effect for plotting
    abs_effect = abs(effect)
  )

# Sort by effect size
df <- df %>%
  arrange(effect) %>%
  mutate(Species = factor(Species, levels = Species))

# Define colors matching other figures
sink_colors <- c("Sink A" = "#4575b4", "Sink B" = "#1a9850")

# Filter to top 25 from each direction (top 50 total)
top_taxa <- df %>%
  group_by(Direction) %>%
  slice_max(abs_effect, n = 25) %>%
  ungroup() %>%
  arrange(effect) %>%
  mutate(Species = factor(Species, levels = Species))

# Add annotation position columns for markers
top_taxa <- top_taxa %>%
  mutate(
    anaerobe_x = ifelse(effect > 0, effect + 0.12, effect - 0.12),
    biofilm_x = ifelse(effect > 0, effect + 0.22, effect - 0.22),
    amr_x = ifelse(effect > 0, effect + 0.32, effect - 0.32)
  )

# Functional group colors
functional_colors <- c(
  "Anaerobe + Biofilm" = "#d73027",
  "Anaerobe" = "#fc8d59",
  "Biofilm former" = "#91bfdb",
  "Other" = "#969696"
)

top_taxa$Functional_display <- factor(
  top_taxa$Functional_group,
  levels = c("Anaerobe + Biofilm", "Anaerobe", "Biofilm former", "Other")
)

################################################################################
# MAIN FIGURE: DIVERGING BAR CHART
################################################################################

p_final <- ggplot(top_taxa, aes(x = effect, y = Species)) +
  # Background shading for direction
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
           fill = "#4575b4", alpha = 0.05) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#1a9850", alpha = 0.05) +
  # Bars colored by functional group
  geom_col(aes(fill = Functional_display), width = 0.75) +
  # Center line
  geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
  # Add annotation markers for traits
  geom_point(data = filter(top_taxa, Anaerobic == "Yes"),
             aes(x = anaerobe_x, y = Species),
             shape = 16, color = "#d73027", size = 1.5, inherit.aes = FALSE) +
  geom_point(data = filter(top_taxa, Biofilm.former == "Yes"),
             aes(x = biofilm_x, y = Species),
             shape = 17, color = "#7570b3", size = 1.5, inherit.aes = FALSE) +
  geom_point(data = filter(top_taxa, Antimicrobial.tolerant == "Yes"),
             aes(x = amr_x, y = Species),
             shape = 8, color = "#e7298a", size = 1.5, inherit.aes = FALSE) +
  # Color scale
  scale_fill_manual(values = functional_colors, name = "Functional Group") +
  scale_x_continuous(
    breaks = seq(-1.5, 2, 0.5),
    limits = c(-1.7, 2.5)
  ) +
  # Labels
  labs(
    x = "ALDEx2 Effect Size",
    y = NULL,
    title = "Differentially Abundant Taxa Between Sinks"
  ) +
  # Theme
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.3),
    axis.text.y = element_text(size = 8, face = "italic"),
    axis.text.x = element_text(size = 9),
    axis.title.x = element_text(size = 10, margin = margin(t = 10)),
    legend.position = "bottom",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    plot.title = element_text(size = 12, face = "bold"),
    plot.margin = margin(10, 15, 10, 10)
  ) +
  guides(fill = guide_legend(nrow = 1))

# Add direction labels at top
p_final <- p_final +
  annotate("text", x = -0.85, y = Inf, label = "Enriched in Sink A",
           color = "#4575b4", size = 3.5, fontface = "bold", vjust = -0.5) +
  annotate("text", x = 1.1, y = Inf, label = "Enriched in Sink B",
           color = "#1a9850", size = 3.5, fontface = "bold", vjust = -0.5) +
  coord_cartesian(clip = "off")

# Save final version
ggsave("figure_aldex2_differential_abundance.pdf", p_final, 
       width = 9, height = 13)

################################################################################
# ALTERNATIVE: SPLIT LAYOUT (SIDE-BY-SIDE)
################################################################################

# Split into two panels: Sink A enriched (left), Sink B enriched (right)
sink_a_taxa <- top_taxa %>% 
  filter(Direction == "Sink A") %>%
  arrange(desc(abs_effect)) %>%
  mutate(Species = factor(Species, levels = Species))

sink_b_taxa <- top_taxa %>% 
  filter(Direction == "Sink B") %>%
  arrange(desc(abs_effect)) %>%
  mutate(Species = factor(Species, levels = Species))

# Panel A (Sink A enriched)
p_a <- ggplot(sink_a_taxa, aes(x = abs_effect, y = Species, fill = Functional_display)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = functional_colors, guide = "none") +
  scale_x_reverse(limits = c(1.6, 0)) +
  labs(x = NULL, y = NULL, title = "Enriched in Sink A") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 8, face = "italic", hjust = 1),
    plot.title = element_text(size = 11, face = "bold", color = "#4575b4", hjust = 0.5),
    plot.margin = margin(10, 5, 10, 10)
  )

# Panel B (Sink B enriched)
p_b <- ggplot(sink_b_taxa, aes(x = abs_effect, y = Species, fill = Functional_display)) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = functional_colors, name = "Functional Group") +
  scale_x_continuous(limits = c(0, 2.2)) +
  labs(x = NULL, y = NULL, title = "Enriched in Sink B") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 8, face = "italic"),
    plot.title = element_text(size = 11, face = "bold", color = "#1a9850", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = 9, face = "bold"),
    plot.margin = margin(10, 10, 10, 5)
  ) +
  guides(fill = guide_legend(nrow = 2))

# Combine side by side
p_split <- plot_grid(p_a, p_b, ncol = 2, rel_widths = c(0.45, 0.55), align = "h")

# Add shared x-axis label
p_split <- ggdraw(p_split) +
  draw_label("ALDEx2 Effect Size (absolute)", x = 0.5, y = 0.02, size = 10)

ggsave("figure_aldex2_split_layout.pdf", p_split, width = 12, height = 10)

cat("\n=== Publication Figures Saved ===\n")
cat("  - figure_aldex2_differential_abundance.pdf (diverging bars)\n")
cat("  - figure_aldex2_split_layout.pdf (side-by-side layout)\n")

################################################################################
# NOTES:
# - Top 25 taxa from each direction (50 total) based on absolute effect size
# - Functional annotations: Anaerobe, Biofilm former, AMR tolerant
# - Color scheme: Blue (Sink A), Green (Sink B)
# - Markers indicate phenotypic traits (circle = anaerobe, triangle = biofilm, 
#   asterisk = AMR tolerant)
# - Effect size interpretation: Negative = Sink A enriched; Positive = Sink B enriched
################################################################################
