library(tidyverse)

ld <- read_tsv("/media/drive_5_usb/Childs/texanus_2021/ld/texanus.freebayes.BC1qtl.parentagesnps.deb.geno.ld") %>%
  rename(r2 = `R^2`)

ld %>%
  filter(CHR == "Ha412HOChr01") %>%
  ggplot(.,aes(x=as.factor(POS1),y=as.factor(POS2),fill=r2)) + 
  geom_tile() +
  scale_fill_viridis_c()

ld %>%
  rename(CHR1 = CHR) %>%
  mutate(CHR2 = CHR1) -> ld_for_inter

inter_ld <- read_tsv("/media/drive_5_usb/Childs/texanus_2021/ld/texanus.freebayes.BC1qtl.parentagesnps.deb.interchrom.geno.ld") %>%
  rename(r2 = `R^2`)
inter_ld$CHR1 <- gsub("Ha412HO","",inter_ld$CHR1)
inter_ld$CHR2 <- gsub("Ha412HO","",inter_ld$CHR2)

pdf("2021/BC1_interchromosome_LD.pdf",height=5,width=6)

inter_ld %>%
  group_by(CHR1,CHR2) %>%
  summarise(mean_ld = mean(r2,na.rm=T)) %>%
  ggplot(.,aes(x=CHR1,y=CHR2,fill=mean_ld)) +
  geom_tile() +
  scale_fill_viridis_c(name="Mean linkage") +
  theme(axis.text.x = element_text(angle = 60, hjust=1)) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


inter_ld %>%
  filter(CHR1 == "Chr06", CHR2 == "Chr15") %>%
  ggplot(.,aes(x=as.factor(POS1),y=as.factor(POS2),fill=r2)) +
  geom_tile() +
  scale_fill_viridis_c(name="r2") +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank()) +
  xlab("Chr06") + ylab("Chr15")

inter_ld %>%
  filter(CHR1 == "Chr12", CHR2 == "Chr16") %>%
  ggplot(.,aes(x=as.factor(POS1),y=as.factor(POS2),fill=r2)) +
  geom_tile() +
  scale_fill_viridis_c(name="r2") +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank()) +
  xlab("Chr12") + ylab("Chr16") 

inter_ld %>%
  filter(CHR1 == "Chr12", CHR2 == "Chr17") %>%
  ggplot(.,aes(x=as.factor(POS1),y=as.factor(POS2),fill=r2)) +
  geom_tile() +
  scale_fill_viridis_c(name="r2") +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank()) +
  xlab("Chr12") + ylab("Chr17")

inter_ld %>%
  filter(CHR1 == "Chr16", CHR2 == "Chr17") %>%
  ggplot(.,aes(x=as.factor(POS1),y=as.factor(POS2),fill=r2)) +
  geom_tile() +
  scale_fill_viridis_c(name="r2") +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank()) +
  xlab("Chr16") + ylab("Chr17") 

inter_ld %>%
  filter(CHR1 == "Chr04", CHR2 == "Chr07") %>%
  ggplot(.,aes(x=as.factor(POS1),y=as.factor(POS2),fill=r2)) +
  geom_tile() +
  scale_fill_viridis_c(name="r2") +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank()) +
  xlab("Chr04") + ylab("Chr07") 


inter_ld %>%
  ggplot(.,aes(x=as.factor(POS1),y=as.factor(POS2),fill=r2)) +
  geom_tile() +
  scale_fill_viridis_c(name="r2") +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank()) +
  facet_wrap(CHR1 ~ CHR2)

inter_ld %>%
  mutate(high_ld = case_when(r2 > 0.3 ~ "high",
                             TRUE ~ "low")) %>%
  group_by(CHR1,CHR2,high_ld) %>%
  summarise(n = n()) %>%
  group_by(CHR1,CHR2) %>%
  mutate(total_n = sum(n)) %>%
  filter(high_ld == "low") %>%
  mutate(percent_high = 1-(n/total_n)) %>%
  filter(CHR1 == "Ha412HOChr04" | CHR2 == "Ha412HOChr04") %>%
  ggplot(.,aes(x=CHR1,y=CHR2,fill=percent_high)) +
  geom_tile() +
  scale_fill_viridis_c(name="high linkage") +
  theme(axis.text.x = element_text(angle = 60, hjust=1))

inter_ld %>%

  group_by(CHR1,CHR2) %>%
  summarise(max_ld = max(r2)) %>%
  ggplot(.,aes(x=CHR1,y=CHR2,fill=max_ld)) +
  geom_tile() +
  scale_fill_viridis_c(name="max linkage") +
  theme(axis.text.x = element_text(angle = 60, hjust=1))
