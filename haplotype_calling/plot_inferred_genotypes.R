library(tidyverse)
library(cowplot)
library(patchwork)
sample_info <- read_tsv("/media/drive_5_usb/Childs/texanus_2021/meta/texanus.sampleinfo.sequenced.txt",
                        col_names = c("sample","year","gen","loc_type","type","location",
                                      "folder","seq_1","seq_2","seq_type"))%>%
  filter(seq_type == "both" | seq_type == "GBS2") %>%
  filter(type == "BC")
parents <- c("deb","ann1","ann2","ann3")

#Load in actual genotypes
all_genotypes <- tibble()
for (parent in parents){
  data <- read_tsv(paste0("/media/drive_5_usb/Childs/texanus_2021/vcf/texanus.freebayes.parentagesnps.",parent,".tidy.txt.gz"),
                   col_types = c("cdcc")) %>%
    mutate(parent=parent) 
  all_genotypes <- rbind(all_genotypes, data)
}


sample_chosen <- "2008.696"
directory <- "/media/drive_5_usb/Childs/texanus_2021/ancestry_hmm/"
all_ancestry <- tibble()
for (parent in parents){
  data <- read_tsv(paste0(directory,parent,"/",sample_chosen,".blocks.txt")) %>% 
    mutate(parent=parent) 
  all_ancestry <- rbind(all_ancestry, data)
}
tmp <- all_genotypes %>%
  filter(sample == sample_chosen)
p1 <- tmp %>%
  filter(chr == "Ha412HOChr03") %>%
  mutate(parent = fct_relevel(parent, c("ann3","ann2","ann1","deb"))) %>%
  ggplot(aes(x=pos/1000000,y=parent,color=genotype_count)) +
  geom_point(shape=124) +
 # facet_wrap(~chr) +
  scale_color_brewer(palette = "Set1",name="Genotype state") + theme_cowplot() +
  xlab("Mbp")

p2 <- all_ancestry %>%
  filter(chr == "Ha412HOChr03") %>%
  mutate(parent = fct_relevel(parent, c("ann3","ann2","ann1","deb"))) %>%
  ggplot(.,aes(x=start/1000000,xend=end/1000000,y=parent,yend=parent,color=as.factor(state),group=parent)) +
  geom_segment(size=3) +
  #facet_wrap(~chr) +
  scale_color_brewer(palette = "Set1",name="HMM state") + theme_cowplot() +
  xlab("Mbp")


sample_chosen <- "2009.002"
directory <- "/media/drive_5_usb/Childs/texanus_2021/ancestry_hmm/"
all_ancestry <- tibble()
for (parent in parents){
  data <- read_tsv(paste0(directory,parent,"/",sample_chosen,".blocks.txt")) %>% 
    mutate(parent=parent) 
  all_ancestry <- rbind(all_ancestry, data)
}
tmp <- all_genotypes %>%
  filter(sample == sample_chosen)
p3 <- tmp %>%
  filter(chr == "Ha412HOChr03") %>%
  mutate(parent = fct_relevel(parent, c("ann3","ann2","ann1","deb"))) %>%
  ggplot(aes(x=pos/1000000,y=parent,color=genotype_count)) +
  geom_point(shape=124) +
  # facet_wrap(~chr) +
  scale_color_brewer(palette = "Set1",name="Genotype state") + theme_cowplot() +
  xlab("Mbp")

p4 <- all_ancestry %>%
  filter(chr == "Ha412HOChr03") %>%
  mutate(parent = fct_relevel(parent, c("ann3","ann2","ann1","deb"))) %>%
  ggplot(.,aes(x=start/1000000,xend=end/1000000,y=parent,yend=parent,color=as.factor(state),group=parent)) +
  geom_segment(size=3) +
  #facet_wrap(~chr) +
  scale_color_brewer(palette = "Set1",name="HMM state") + theme_cowplot() +
  xlab("Mbp")


sample_chosen <- "2009.808"
directory <- "/media/drive_5_usb/Childs/texanus_2021/ancestry_hmm/"
all_ancestry <- tibble()
for (parent in parents){
  data <- read_tsv(paste0(directory,parent,"/",sample_chosen,".blocks.txt")) %>% 
    mutate(parent=parent) 
  all_ancestry <- rbind(all_ancestry, data)
}
tmp <- all_genotypes %>%
  filter(sample == sample_chosen)
p5 <- tmp %>%
  filter(chr == "Ha412HOChr03") %>%
  mutate(parent = fct_relevel(parent, c("ann3","ann2","ann1","deb"))) %>%
  ggplot(aes(x=pos/1000000,y=parent,color=genotype_count)) +
  geom_point(shape=124) +
  # facet_wrap(~chr) +
  scale_color_brewer(palette = "Set1",name="Genotype state") + theme_cowplot() +
  xlab("Mbp")

p6 <- all_ancestry %>%
  filter(chr == "Ha412HOChr03") %>%
  mutate(parent = fct_relevel(parent, c("ann3","ann2","ann1","deb"))) %>%
  ggplot(.,aes(x=start/1000000,xend=end/1000000,y=parent,yend=parent,color=as.factor(state),group=parent)) +
  geom_segment(size=3) +
  #facet_wrap(~chr) +
  scale_color_brewer(palette = "Set1",name="HMM state") + theme_cowplot() +
  xlab("Mbp")

all_markers <- tibble()
for (parent in parents){
  markers <- read_tsv(paste0("/media/drive_5_usb/Childs/texanus_2021/ancestry_hmm/",parent,"/",sample_chosen,".markers.txt"),
                      col_names = c("chr","pos","state")) %>%
    mutate(parent = parent)
  all_markers <- rbind(all_markers,markers)
}
all_markers %>%
  filter(chr == "Ha412HOChr03") %>%
  mutate(parent = fct_relevel(parent, c("ann3","ann2","ann1","deb"))) %>%
  ggplot(.,aes(x=pos/1000000,y=parent,color=as.factor(state),group=parent)) +
  geom_point(shape=124,size=4) +
  #facet_wrap(~chr) +
  scale_color_brewer(palette = "Set1",name="HMM state") + theme_cowplot() +
  xlab("Mbp")


p1 / p2 / p3 / p4 / p5 / p6
#Minimum covered range of the genome
all_ancestry %>% 
  group_by(chr,parent) %>% 
  summarize(min_pos = min(pos), max_pos = max(pos)) %>%
  group_by(chr) %>%
  summarize(start = max(min_pos),end=min(max_pos)) %>%
  write_tsv(.,"/media/drive_5_usb/Childs/texanus_2021/meta/texanus.parentagecovered.bed",
            col_names = F)
