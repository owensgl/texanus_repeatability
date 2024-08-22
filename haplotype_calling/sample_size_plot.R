library(tidyverse)
library(cowplot)
sample_info <- read_tsv("/media/drive_5_usb/Childs/texanus_2021/meta/texanus.sampleinfo.sequenced.txt",
                        col_names = c("sample","year","gen","loc_type","type","location",
                                      "folder","seq_1","seq_2","seq_type"))%>%
  filter(seq_type == "both" | seq_type == "GBS2") %>%
  filter(type == "BC")

pdf("2021/plots/samplesize.v1.pdf")
sample_info %>%
  filter(gen %% 1 == 0) %>% 
  group_by(gen, location) %>%
  summarize(n=n()) %>%
  ggplot(.,aes(x=gen,y=n)) +
  geom_bar(stat="identity",position="dodge") +
  theme_cowplot() +
  facet_wrap(~location) +
  ylab("Sample size") +
  xlab("Generation") +
  ggtitle("Location")
dev.off()
