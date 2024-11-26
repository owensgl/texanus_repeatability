library(tidyverse)

climate <- read_tsv("texanus_repeatability/meta/common_garden_climate.txt",
         col_names = c("name","long","lat",paste0("bio",rep(1:19))),
         skip=1)
# Specify the desired order of samples
sample_order <- c("LBJ", "BFL", "KPC", "HCC")

# Filter and reorder the climate data
climate <- climate %>%
  filter(name %in% sample_order) %>%
  arrange(factor(name, levels = sample_order))
# Euclidean distance
dist <- dist(climate[ , c(4:22)] , diag=TRUE)

library(reshape2)
  

dist_matrix <- as.matrix(dist)
# Melt the matrix for ggplot2
rownames(dist_matrix) <- climate$name
colnames(dist_matrix) <- climate$name

# Melt the matrix for ggplot2
melted_dist <- melt(dist_matrix)
colnames(melted_dist) <- c("Var1", "Var2", "Distance")

# Create the heatmap
pdf("2021/plots/environmental_distance.v1.pdf")
melted_dist %>%
  filter(as.numeric(Var1) <= as.numeric(Var2)) %>%
ggplot(aes(x = Var1, y = Var2, fill = Distance)) +
  geom_tile() +
  scale_fill_viridis_c(name="Environmental\ndistance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank()) +
  coord_fixed()
dev.off()
