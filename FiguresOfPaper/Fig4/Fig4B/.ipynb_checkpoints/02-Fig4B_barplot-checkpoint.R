# Load required packages
library(dplyr)
library(svglite)
library(ggplot2)
library(RColorBrewer)
library(extrafont)

# Import custom fonts
font_import(paths = "your_fonts_path/fonts/")
fonts()

### Read genus data
genus = read.table("./data/Response_Non-Response_genus.txt", header = TRUE, sep = "\t", row.names = NULL)

# Sort response data and create ordered factors
sorted_response <- genus %>% filter(category == "Response") %>% arrange(value)
genus$subcategory <- factor(genus$subcategory, levels = unique(genus$subcategory[genus$category == "Response"]))

### Read species data
species = read.table("./data/Response_Non-Response_species.txt", header = TRUE, sep = "\t", row.names = NULL)
species$subcategory <- gsub("_", " ", species$subcategory)
sorted_response <- species %>% filter(category == "Response") %>% arrange(value)
species$subcategory <- factor(species$subcategory, levels = unique(species$subcategory[species$category == "Response"]))

# Define color palette
colors <- c("#9d001f","#d64863", "#540051","#ff82ab","#d06bc6",
            "#ff70ab", "#0083cd", "#acecff", "#00beff","#002f50", "#2a749a")

# Generate genus plot
p1 <- ggplot(genus, aes(x = category, y = value, fill = subcategory)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colors[1:(nrow(genus)/2)]) +
  guides(fill = guide_legend(title = "", ncol = 1)) +
  ylab("Relative abundances of microbes") +
  theme(
    axis.title.y = element_text(size = 9, face = "bold", family = "Helvetica"),
    axis.text.x = element_text(angle = 25, hjust = 1, size = 8, face = "bold", color = "black", family = "Helvetica"),
    axis.text.y = element_text(size = 8, face = "bold", color = "black", family = "Helvetica"),
    legend.key.size = unit(0.6, "lines"),
    legend.text = element_text(size = 10, face = "bold.italic", color = "black", family = "Helvetica"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  xlab("") +
  scale_x_discrete(limits = c("Response", "Non_Response"),
                   labels = c("Response", "Non-Response"))

# Save genus plot
ggsave("./figures/Fig4B_barplot_genus.pdf", p1, width = 3, height = 3, device = "pdf")

# Generate species plot
p2 <- ggplot(species, aes(x = category, y = value, fill = subcategory)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = colors[1:(nrow(genus)/2)]) +
  guides(fill = guide_legend(title = "", ncol = 1)) +
  ylab("Relative abundances of microbes") +
  theme(
    axis.title.y = element_text(size = 9, face = "bold", family = "Helvetica"),
    axis.text.x = element_text(angle = 25, hjust = 1, size = 8, face = "bold", color = "black", family = "Helvetica"),
    axis.text.y = element_text(size = 8, face = "bold", color = "black", family = "Helvetica"),
    legend.key.size = unit(0.6, "lines"),
    legend.text = element_text(size = 10, face = "bold.italic", color = "black", family = "Helvetica"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  xlab("") +
  scale_x_discrete(limits = c("Response", "Non_Response"),
                   labels = c("Response", "Non-Response"))

# Save species plot
ggsave("./figures/Fig4B_barplot_species.pdf", p2, width = 4, height = 3, device = "pdf")