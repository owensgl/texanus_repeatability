#Plotting the amount of the genome with extra or fewer haplotypes from outside gene flow
library(tidyverse)
library(cowplot)
library(broom)
library(colorspace)
library(patchwork)
sample_info <- read_tsv("/media/drive_5_usb/Childs/texanus_2021/meta/texanus.sampleinfo.sequenced.txt",
                        col_names = c("sample","year","gen","loc_type","type","location",
                                      "folder","seq_1","seq_2","seq_type"))%>%
  filter(seq_type == "both" | seq_type == "GBS2") %>%
  filter(type == "BC")

state_fixer <- tibble(state=c("00","01","11",NA),
                      genotype=c(0,1,2,NA))



parents <- c("deb","ann1","ann2","ann3")
directory <- "/media/drive_5_usb/Childs/texanus_2021/ancestry_hmm/"
all_ancestry <- read_tsv(paste0(directory,"allmarkers.txt.gz"),
                         col_names = c("chr","pos","state","sample","parent"),
                         col_types=c("cdccc"))

all_ancestry <- all_ancestry %>%
  inner_join(state_fixer) %>%
  inner_join(sample_info) 

#Copy first generation to all locations
first_gen <- tibble()
for (loc in c("HCC","BFL","LBJ","KPC")){
  first_gen <- rbind(first_gen, all_ancestry %>%
                       filter(gen == 1) %>%
                       mutate(location = loc))
}
all_ancestry <- all_ancestry %>%
  filter(gen != 1) %>%
  rbind(first_gen)

summed_ancestry <- all_ancestry %>%
  filter(gen != 2.5) %>%
  group_by(chr, pos, sample, location, gen) %>%
  summarize(total_genotypes = sum(genotype)) 

summed_ancestry <- summed_ancestry %>%
group_by(sample, location, gen, total_genotypes) %>%
summarize(genome_copies = n())




sample_order <- summed_ancestry %>%
  filter(total_genotypes == 2) %>%
  group_by(location, gen, sample) %>%
  summarize(genome_copies_sum = sum(genome_copies)) %>%
  arrange(location, gen, desc(genome_copies_sum)) %>%  # Sort by genome_copies descending
  group_by(location, gen) %>%
  mutate(sample_number = row_number()) %>%
  ungroup()
  


pdf("2021/plots/extra_genome_copies.v1.pdf",height=6,width=10)
summed_ancestry %>%
  filter(gen != "7.5") %>%
  inner_join(sample_order) %>%
  ggplot(.,aes(x=sample_number,y=genome_copies)) +
  geom_bar(stat="identity",position="stack",aes(fill=as.factor(total_genotypes))) +
  scale_fill_brewer(palette = "Set1",name="Haplotype copies") +
  facet_grid(as.factor(location)~as.factor(gen),scales="free_x") +
  ylab("Genomic regions") +
  xlab("Sample") +
  theme_cowplot() +
  theme(axis.text.x = element_blank()) 
dev.off()

#3058 windows

p1 <- summed_ancestry %>%
  filter(gen != 7.5) %>%
  filter(total_genotypes == 2) %>%
  mutate(percent_correct = 100*genome_copies/3058) %>%
  ggplot(.,aes(x=location,y=percent_correct,color=as.factor(gen))) +
  geom_boxplot() +
  scale_color_viridis_d(name="Generation") +
  theme_cowplot() +
  ylab("Percent of genome with\ntwo haplotype copies") +
  xlab("Location")
  

p2 <- summed_ancestry %>%
  filter(gen != 7.5) %>%
  group_by(sample, location, gen) %>%
  mutate(copies = (total_genotypes*genome_copies)) %>%
  summarize(percent_local = sum(copies,na.rm=T)/(3058*2)) %>%
  ggplot(.,aes(x=location,y=percent_local*100,color=as.factor(gen))) +
  geom_boxplot() +
  scale_color_viridis_d(name="Generation") +
  theme_cowplot() +
  ylab("Percent local haplotypes") +
  xlab("Location")
  

pdf("2021/plots/percent_geneflow.v1.pdf")
p1 / p2 + plot_annotation(tag_levels = 'a') +
  plot_layout(guides = 'collect')
dev.off()
