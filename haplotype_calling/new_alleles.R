library(tidyverse)
library(cowplot)
library(patchwork)
new_alleles <- read_tsv("../../../Childs/texanus_2021/vcf/texanus.freebayes.newalleles.txt")

sample_info <- read_tsv("../../../Childs/texanus_2021/meta/texanus.sampleinfo.sequenced.nonqtl.txt",
                        col_names = c("sample","year","gen","location_class","class",
                                      "location","group","seq_1","seq_2","seq_3"),
                        col_types = "cddcccccc")

n1 <- new_alleles %>%
  inner_join(sample_info) %>%
  filter(gen != 1 & gen != 2.5 & gen != 7.5) %>%
  ggplot(.,aes(x=gen,y=new_alleles)) +
  #geom_boxplot(aes(group=gen),outlier.shape = NA) +
  geom_jitter(alpha=0.2,width=0.2) +
  facet_wrap(~location) +
  theme_cowplot() +
  ylab("Number of new alleles") +
  xlab("Generation")

n3 <- new_alleles %>%
  inner_join(sample_info) %>%
  filter(location == "KPC") %>%
  filter(gen != 1 & gen != 2.5 & gen != 7.5) %>%
  mutate(gen = paste0("Generation: ",gen)) %>%
  ggplot(.,aes(new_alleles)) +
  #geom_boxplot(aes(group=gen),outlier.shape = NA) +
  geom_histogram() +
  facet_wrap(~gen,ncol=1) +
  theme_cowplot() +
  ylab("Sample count") +
  xlab("Number of new alleles") +
  ggtitle("KPC")

n2 <- new_alleles %>%
  inner_join(sample_info) %>%
  filter(gen != 1 & gen != 2.5 & gen != 7.5) %>%
  mutate(novel = case_when(new_alleles > 1000 ~ "Contaminated",
                           TRUE ~ "Pure")) %>%
  group_by(gen,location,novel ) %>%
  summarize(count=n()) %>%
  ggplot(.,aes(x=gen,y=count,fill=novel)) +
  geom_bar(position="fill", stat="identity") +
  facet_wrap(~location) +
  theme_cowplot() +
  ylab("Proportion of samples") +
  xlab("Generation") +
  scale_fill_manual(values = c("#368a8d", "#d23d27"),
                    name="Non-parent alleles",
                    labels=c("Present","Absent")) +
  theme(legend.position='top')
pdf("2021/plots/new_alleles.v2.pdf",height=9,width=10)
((n1/ n2) | n3) + plot_annotation(tag_levels = 'A')
dev.off()
