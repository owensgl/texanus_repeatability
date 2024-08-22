library(poolSeq)

library(psych) 
library(tidyverse)
library(cowplot)
library(broom)
library(colorspace)
library(patchwork)
library(geosphere)
sample_info <- read_tsv("/media/drive_5_usb/Childs/texanus_2021/meta/texanus.sampleinfo.sequenced.txt",
                        col_names = c("sample","year","gen","loc_type","type","location",
                                      "folder","seq_1","seq_2","seq_type"))%>%
  filter(seq_type == "both" | seq_type == "GBS2") %>%
  filter(type == "BC")


population_list <- c("HCC","BFL","LBJ","KPC")
parent_list <- c("ann1","ann2","ann3","deb")

selection_output <- tibble()
for (population_chosen in population_list){
  for (parent_chosen in parent_list){
    tmp <- read_tsv(paste0("/media/drive_5_usb/Childs/texanus_2021/vcf/poolseq_selection_output.",population_chosen,
                           ".",parent_chosen,".txt.gz"))
    
    selection_output <- rbind(selection_output, tmp)
  }
}


#Get genetic recombination rate per position
genetic_map <- read_tsv("/media/drive_5_usb/Fuchs/ref/Ha412HOv2.0-20181130.Nov22k22.geneticmap.extradivisions.txt") %>% 
  mutate(rate = lead(cM)-cM) %>% 
  rename(start = pos) %>% mutate(end = start + 999999) %>%
  mutate(rate = case_when(is.na(rate) ~ cM - lag(cM),
                          rate < 0 ~ cM - lag(cM),
                          TRUE ~ rate))
chr_lengths <- read_tsv("../wild_gwas_2018/Ha412HO.chrlengths.txt")

windows <- selection_output %>%
  ungroup() %>%
  select(chr,pos) %>%
  distinct()


windows$cm_rate <- NA
for (i in 1:nrow(windows)){
  print(i)
  chosen_start <- windows$pos[i]
  chosen_chr <- windows$chr[i]
  rate <- genetic_map %>% filter(chr == chosen_chr) %>% 
    filter(start <= chosen_start, end >= chosen_start) %>%
    pull(rate)
  windows$cm_rate[i] <- rate
}
select_output <- inner_join(selection_output,windows)

selection_output %>%
  #filter(chr == "Ha412HOChr17") %>%
  filter(parent == "ann1") %>%
  #filter(location == "LBJ") %>%
  ggplot(.) +
  geom_point(aes(x=pos/1000000,y=s,color=location)) +
  #geom_linerange(aes(x=pos/1000000,ymin=s_min,ymax=s_max)) +
  facet_grid(~chr) +
  scale_color_brewer(palette = "Set1") +
  theme_cowplot()  +
  geom_hline(yintercept = 0,linetype="dotted") +
  xlab("MBp") +
  facet_grid(~chr)


selection_output %>%
  ggplot(.,aes(x=s,fill=parent)) +
  geom_density(alpha=0.8) +
  facet_wrap(~location) +
  theme_cowplot() +
  geom_vline(xintercept=0,linetype="dashed") +
  xlab("Selection coefficient") +
  scale_fill_brewer(palette = "Set1",name="Location")



###


select_output_cum <- chr_lengths %>%
  select(chr,end) %>%
  rename(chr_len = end) %>%
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(select_output, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate( poscum=pos+tot)

axisdf = select_output_cum %>% group_by(chr) %>% summarize(center=( max(poscum) + min(poscum) ) / 2 )

good_seed_qtl_locations <- read_tsv("2021/good_seed_deb_qtl.locations.txt")
good_seed_qtl_locations_cum <- chr_lengths %>%
  select(chr,end) %>%
  rename(chr_len = end) %>%
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(good_seed_qtl_locations, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, pos) %>%
  mutate( poscum=pos+tot)
good_seed_qtl_locations_cum <- good_seed_qtl_locations_cum %>%
  select(chr,type,poscum) %>%
  pivot_wider(names_from = type,values_from=poscum)
pdf("2021/selection_estimate_debilis.v1.pdf",height=8,width=10)
p1 <- select_output_cum %>%
  filter(parent == "deb") %>%
  ggplot(.) +
  geom_line( aes(x=poscum, y=s,color=as.factor(chr)), alpha=0.8, size=1) +
  # Show all points
  scale_color_manual(values = rep(c("#24492e", "#2c6184"), 17 )) +
  
  # custom X axis:
  scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$chr), breaks= axisdf$center ) +
  theme_cowplot() +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  xlab("Chromosome") +
  facet_wrap(~location,nrow=4,scales="free_y") +
  geom_hline(yintercept=0,linetype="dashed") +
  geom_linerange(data=good_seed_qtl_locations_cum,
                 aes(y=0.5,xmin=start,xmax=end)) +
  geom_point(data=good_seed_qtl_locations_cum,
                 aes(y=0.5,x=max)) 

p2 <- select_output_cum %>%
  filter(parent == "deb") %>%
  ggplot(.) +
  geom_line( aes(x=poscum, y=cm_rate,color=as.factor(chr)), alpha=0.8, size=1) +
  # Show all points
  scale_color_manual(values = rep(c("#24492e", "#2c6184"), 17 )) +
  
  # custom X axis:
  scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$chr), breaks= axisdf$center ) +
  theme_cowplot() +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  xlab("Chromosome") +
  ylab("Recombination\nrate") 
p1 / p2 +
  plot_layout(heights = c(4, 1))
dev.off()

####
library(lme4)
library(car)
select_output %>%
  filter(parent == "deb") %>%
  ggplot(.,aes(x=log10(cm_rate),y=s,color=location)) +
  geom_point() +
  geom_smooth(method="loess") +
  facet_wrap(~chr)

a <- lm(formula = s ~ location +chr + cm_rate, 
        data = select_output %>% filter(parent == "deb"))

summary(a)
anova(a)
af <- anova(a)
afss <- af$"Sum Sq"
print(cbind(af,PctExp=afss/sum(afss)*100))
Anova(a)

a <- lm(formula = s ~ location + chr + cm_rate, 
        data = select_output %>% filter(parent == "deb") %>%
          filter(chr == "Ha412HOChr06" |
                   chr == "Ha412HOChr12" |
                   chr == "Ha412HOChr15" |
                   chr == "Ha412HOChr16" |
                   chr == "Ha412HOChr17" ))
summary(a)
anova(a)
af <- anova(a)
afss <- af$"Sum Sq"
print(cbind(af,PctExp=afss/sum(afss)*100))
Anova(a)

a <- lm(formula = s ~ location+ cm_rate, 
        data = select_output %>% filter(parent == "deb") %>%
          filter(chr == "Ha412HOChr02"  ))
summary(a)
anova(a)
af <- anova(a)
afss <- af$"Sum Sq"
print(cbind(af,PctExp=afss/sum(afss)*100))
Anova(a)



#####
p1 <- select_output_cum %>%
  filter(parent == "deb") %>%
  filter(chr == "Ha412HOChr02" ) %>%
  group_by(chr) %>%
  mutate( quintile_rank = ntile(cm_rate,5)) %>%
  group_by(location, quintile_rank,chr) %>%
  summarize(mean_s = mean(s),sd = sd(s),
            se = sd(s)/sqrt(n())) %>%
  ggplot(.) +
  geom_point(aes(x=quintile_rank,y=mean_s,color=location,group=location),
             position = position_dodge(width=0.6)) +
  geom_linerange(aes(x=quintile_rank,ymin=mean_s-2*se,ymax=mean_s+2*se,
                     color=location,group=location),
                 position = position_dodge(width=0.6)) +
  facet_wrap(~chr) +
  theme_cowplot()  +
  ylab("Selection coefficient") +
  xlab("Recombination quintile") +
  scale_color_brewer(palette = "Set1",name="Location")
p2 <- select_output_cum %>%
  filter(parent == "deb") %>%
  filter(chr == "Ha412HOChr12" ) %>%
  group_by(chr) %>%
  mutate( quintile_rank = ntile(cm_rate,5)) %>%
  group_by(location, quintile_rank,chr) %>%
  summarize(mean_s = mean(s),sd = sd(s),
            se = sd(s)/sqrt(n())) %>%
  ggplot(.) +
  geom_point(aes(x=quintile_rank,y=mean_s,color=location,group=location),
             position = position_dodge(width=0.6)) +
  geom_linerange(aes(x=quintile_rank,ymin=mean_s-2*se,ymax=mean_s+2*se,
                     color=location,group=location),
                 position = position_dodge(width=0.6)) +
  facet_wrap(~chr) +
  theme_cowplot()  +
  ylab("Selection coefficient") +
  xlab("Recombination quintile") +
  scale_color_brewer(palette = "Set1",name="Location")
 
pdf("2021/selection_and_recombination.v1.pdf",width=8,height=4)
p1 + p2 
dev.off()

########Plotting the other parental genome selection

pdf("2021/selection_coefficients_ann.v1.pdf",height=10,width=8)
select_output_cum %>%
  filter(parent == "ann1") %>%
  ggplot(.) +
  geom_line( aes(x=poscum, y=s,color=as.factor(chr)), alpha=0.8, size=1) +
  # Show all points
  scale_color_manual(values = rep(c("#24492e", "#2c6184"), 17 )) +
  
  # custom X axis:
  scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$chr), breaks= axisdf$center ) +
  theme_cowplot() +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  xlab("Chromosome") +
  facet_wrap(~location,nrow=4,scales="free_y") +
  geom_hline(yintercept=0,linetype="dashed") +
  ggtitle("ann1")

select_output_cum %>%
  filter(parent == "ann2") %>%
  ggplot(.) +
  geom_line( aes(x=poscum, y=s,color=as.factor(chr)), alpha=0.8, size=1) +
  # Show all points
  scale_color_manual(values = rep(c("#24492e", "#2c6184"), 17 )) +
  
  # custom X axis:
  scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$chr), breaks= axisdf$center ) +
  theme_cowplot() +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  xlab("Chromosome") +
  facet_wrap(~location,nrow=4,scales="free_y") +
  geom_hline(yintercept=0,linetype="dashed") +
  ggtitle("ann2")

select_output_cum %>%
  filter(parent == "ann3") %>%
  ggplot(.) +
  geom_line( aes(x=poscum, y=s,color=as.factor(chr)), alpha=0.8, size=1) +
  # Show all points
  scale_color_manual(values = rep(c("#24492e", "#2c6184"), 17 )) +
  
  # custom X axis:
  scale_x_continuous( label = gsub("Ha412HOChr","",axisdf$chr), breaks= axisdf$center ) +
  theme_cowplot() +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  xlab("Chromosome") +
  facet_wrap(~location,nrow=4,scales="free_y") +
  geom_hline(yintercept=0,linetype="dashed") +
  ggtitle("ann3")
dev.off()
