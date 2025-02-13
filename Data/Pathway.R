install.packages("pathfindR")
install.packages("pak")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("ReactomePA")
install.packages("hrbrthemes")


library(ReactomePA)
library(pathfindR)
library(ggplot2)
library(biomaRt)
library(plyr)
library(dplyr)
library(stringr)
library(ggrepel)
library(hrbrthemes)


input_df <- read.csv ("RA_input_synergy.csv")
output_df <- run_pathfindR(input_df, gene_sets = "KEGG")
output_df <- run_pathfindR(input_df, gene_sets = "Reactome")
output_df <- run_pathfindR(input_df, gene_sets = "GO-MF")
output_df <- run_pathfindR(input_df, gene_sets = "GO-BP")
output_df <- run_pathfindR(input_df, gene_sets = "GO-CC")
output_df <- run_pathfindR(input_df, gene_sets = "GO-All")
output_df <- run_pathfindR(input_df, gene_sets = "BioCarta")

write.csv(output_df, "output_GO_ALL.csv", row.names=TRUE)

output_df <- read.csv ("output_GO_ALL.csv")
output_df$Term_Description <- factor(output_df$Term_Description, levels = output_df$Term_Description[order(output_df$Fold_Enrichment, decreasing=FALSE)])
output_df$Term_Description # notice the changed order of factor levels

output_df$freq<-str_count(output_df$Down_regulated,'\\w+')

output_df %>% count(Down_regulated)

#Top 15
top15<- output_df %>% 
  top_n(15, Fold_Enrichment)

# Calculate -log10(p-value)
top15$log_p <- -log10(top15$lowest_p)


ggplot(data = top15, aes(x = Fold_Enrichment, 
                         y = reorder(Term_Description, Fold_Enrichment), 
                         color = log_p, 
                         size = freq)) + 
  geom_point() +
  scale_color_gradient2(low = "darkgreen", mid = "seagreen3", high = "indianred2", 
                        midpoint = median(top15$log_p)) +
  scale_size_continuous(range = c(2, 10)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(family = "sans", size = 12),
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    axis.ticks.x = element_line(colour = "black"),
    axis.ticks.y = element_blank(),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_blank(),
    legend.position = "right",
    plot.margin = margin(10, 10, 10, 10),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    aspect.ratio = 3  # Adjust this value to change the plot's aspect ratio
  ) +
  labs(
    y = NULL,
    x = "Fold Enrichment",
    color = "-log10(p-value)",
    size = "Number of Genes"
  )

