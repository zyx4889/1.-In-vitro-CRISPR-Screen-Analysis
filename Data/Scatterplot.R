library(ggplot2)
library(plotly)
library(gapminder)
library(ggrepel)
data <- read.csv ("drugz_results.csv")

p <- ggplot(data, aes(x = rank_synth, y = normZ, color=normZ, label=GENE)) +
  geom_point(show.legend = FALSE) + 
  scale_color_gradient(low = "red", high = "blue") +
  labs(title="MyC-CaP + MYCi-975",
       x="Gene Rank", y = "NormZ Score") +
  # Label top genes in red
  geom_text_repel(
    data = subset(data, rank_synth < 26),
    aes(label = GENE),
    size = 4,
    family = "sans",
    colour = "red",
    segment.color = "grey",
    max.overlaps = 50,
    direction = "both"
  ) +
  # Label bottom genes in blue
  geom_text_repel(
    data = subset(data, rank_synth > 19455),
    aes(label = GENE),
    size = 4,
    family = "sans",
    colour = "blue",
    segment.color = "grey",
    max.overlaps = 50
  ) +
  theme_classic() +
  expand_limits(x = c(min(data$rank_synth) - 20, max(data$rank_synth) + 20)) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1)))

# Print the plot
p

