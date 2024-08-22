library(tidyverse)
library(zoo)
library(cowplot)
deb <- read_tsv("/media/drive_5_usb/Childs/texanus_2021/vcf/texanus.freebayes.5dp.parentagesnps.vcf.recode.debfreq.txt")

first_gen <- tibble()
for (loc in c("HCC","BFL","LBJ","KPC")){
  first_gen <- rbind(first_gen, deb %>%
                       filter(generation == 1) %>%
                       group_by(chr,pos,generation) %>%
                       summarize(percent = sum(percent*total_genotyped)/sum(total_genotyped),
                                 total_genotyped = sum(total_genotyped)) %>%
                       mutate(location = loc))
}
deb <- deb %>%
  filter(generation != 1) %>%
  rbind(first_gen)

pdf("figures/texanus.freebayes.5dp.parentagesnps.allele_freq_change.v1.pdf",height=6,width=10)
deb %>%
  filter(chr == "Ha412HOChr01") %>%
  filter(generation%%1==0) %>%
  #filter(pos >= 65150699, pos <= 65151183 ) %>%
  ggplot(.,aes(x=pos/1000000,y=percent)) +
  geom_point(alpha=0.05) +
  #geom_smooth() +
  facet_grid(location~generation) +
  coord_cartesian(ylim=c(0,1)) +
  theme_cowplot() +
  xlab("MBp") + ylab("Allele frequency") +
  ggtitle("Deb alleles on Chr01") +
  theme(axis.text.x = element_text(size=7))
dev.off()
deb %>%
  dplyr::select(chr,pos) %>%
  distinct() %>% View()

deb %>%
  filter(chr == "Ha412HOChr02") %>%
  #inner_join(deviations %>% filter(mean_deviation < 0.1)) %>%
  group_by(generation,location) %>%
  mutate(median_percent=rollapply(percent,10,median,align='right',fill=NA)) %>% 
  ggplot(.) +
  geom_point(aes(x=pos,y=percent)) +
  geom_line(aes(x=pos,y=median_percent),color="blue") +
  facet_grid(location~generation) +
  scale_color_viridis_c() +
  coord_cartesian(ylim=c(0,1))

deviations <- deb %>%
  group_by(chr,location, generation) %>%
  mutate(median_percent=rollapply(percent,10,median,align='right',partial=T)) %>% 
  mutate(deviation = percent - median_percent) %>%
  group_by(chr,pos, location) %>%
  summarize(mean_deviation = mean(deviation)) %>%
  mutate(id = paste0(chr,".",pos))
  

deviations %>%
  group_by(id) %>%
  mutate(mean_mean_deviation = mean(mean_deviation)) %>%
  ggplot(.) +
  geom_jitter(aes(x=fct_reorder(id,mean_mean_deviation),y=mean_deviation,color=location)) +
  geom_point(aes(x=fct_reorder(id,mean_mean_deviation),y=mean_mean_deviation),color="red") +
  facet_wrap(~chr)

deb %>%
  filter(chr == "Ha412HOChr04",location == "HCC") %>%
  group_by(chr,location, generation) %>%
  mutate(median_percent=rollapply(percent,10,median,align='right',partial=T)) %>% 
  mutate(deviation = percent - median_percent) %>%
  group_by(chr,pos, location) %>%
  mutate(mean_deviation = mean(deviation)) %>%
  mutate(id = paste0(chr,".",pos))
  ggplot(.) +
  geom_jitter(aes(x=fct_reorder(id,mean_deviation),y=deviation)) +
  geom_point(aes(x=fct_reorder(id,mean_deviation),y=mean_deviation),color="blue")
