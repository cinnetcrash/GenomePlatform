#!/usr/bin/env Rscript
# Annotates a phylogenetic tree with sample metadata.
# Args: treefile metadata.csv output.svg

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(tidyverse)
  library(RColorBrewer)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: annotate_tree.R <treefile> <metadata.csv> <output.svg>")

treefile   <- args[1]
meta_file  <- args[2]
output_svg <- args[3]

tree <- read.tree(treefile)
meta <- read.csv(meta_file, stringsAsFactors = FALSE, na.strings = c("", "NA"))

if (!"sample_id" %in% colnames(meta)) stop("metadata.csv must have a 'sample_id' column")

countries <- unique(na.omit(meta$country))
n_col     <- max(length(countries), 1)
pal       <- colorRampPalette(brewer.pal(min(n_col, 12), "Set3"))(n_col)
names(pal) <- countries

source_shapes <- c(hospital=16, food=17, environment=15,
                   animal=18, human=19, reference=8, unknown=4)

p <- ggtree(tree, layout = "rectangular") %<+% meta +
  geom_tippoint(aes(color = country, shape = source), size = 2.5, na.rm = TRUE) +
  geom_tiplab(aes(label = paste0(label, ifelse(!is.na(mlst_st) & mlst_st != "",
                                               paste0(" [", mlst_st, "]"), ""))),
              size = 2.2, offset = 0.001, na.rm = TRUE) +
  scale_color_manual(values = pal, na.value = "grey70", name = "Country") +
  scale_shape_manual(values = source_shapes, na.value = 4, name = "Source") +
  theme_tree2() +
  theme(legend.position = "right", legend.text = element_text(size = 7)) +
  ggtitle("Global Phylogenetic Tree — Core Genome (IQ-TREE, GTR+G)")

ggsave(output_svg, plot = p,
       width = 14, height = max(6, length(tree$tip.label) * 0.25),
       units = "in", device = "svg")

cat("Tree SVG written to:", output_svg, "\n")
