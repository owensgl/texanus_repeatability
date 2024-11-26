#Plotting non-parental ancestry from the number of haplotypes found in each position of the genome after HMM
library(tidyverse)
library(cowplot)
library(zoo)
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
  

tmp <- all_ancestry %>%
  group_by(chr,pos, sample,location,gen) %>%
  summarize(total_chrs = sum(genotype))
tmp <- tmp %>%
  group_by(sample,location,gen,total_chrs) %>%
  summarize(count=n())

diploid_per <- tmp %>%
  group_by(sample) %>%
  mutate(total= sum(count),
         percent_diploid = count/total) %>%
  filter(total_chrs == 2) %>%
  select(sample,percent_diploid)

first_gen <- tibble()
for (loc in c("HCC","BFL","LBJ","KPC")){
  first_gen <- rbind(first_gen, tmp %>%
    filter(gen == 1) %>%
    mutate(location = loc))
}

tmp %>%
  filter(gen != 1) %>%
  rbind(first_gen) %>%
  inner_join(diploid_per) %>%
  filter(location == "BFL") %>%
  ggplot(.,aes(x=fct_reorder(sample,-percent_diploid),y=count,fill=as.factor(total_chrs))) +
  geom_bar(position="stack",stat="identity") +
  facet_wrap(~gen,scales="free_x") +
  theme_cowplot() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 
  
  
#Percent of non-parental genome.
tmp %>%
  filter(gen != 1) %>%
  rbind(first_gen) %>%
  mutate(total_chrs = case_when(total_chrs >= 2 ~ 2,
                                TRUE ~ total_chrs)) %>%
  group_by(sample,location,gen,total_chrs) %>%
  summarize(count = sum(count)) %>%
  pivot_wider(names_from = total_chrs,values_from = count) %>%
  mutate(total_windows = (`0` + `1` + `2`)*2) %>%
  mutate(number_unaccounted = (`0` * 2) + `1`) %>%
  mutate(percent_unaccounted = number_unaccounted/total_windows) %>%
  ggplot(.,aes(x=as.numeric(gen),y=percent_unaccounted)) +
  geom_jitter() +
  facet_wrap(~location) +
  ylab("Proportion non-parental") +
  xlab("Generation") +
  theme_cowplot()

#Outside ancestry across the genome.
first_gen <- tibble()
for (loc in c("HCC","BFL","LBJ","KPC")){
  first_gen <- rbind(first_gen, all_ancestry %>%
                       filter(gen == 1) %>%
                       mutate(location = loc))
}


outside_ancestry <- all_ancestry %>%
  filter(gen != 1) %>%
  rbind(first_gen) %>%
  group_by(chr,pos,location,gen,sample) %>%
  summarize(total_chrs = sum(genotype,na.rm=T)) %>%
  mutate(total_chrs = case_when(total_chrs > 2 ~ 2,
                                TRUE ~ total_chrs)) %>%
  mutate(outside_chrs = 2 - total_chrs)


summarized_outside <- outside_ancestry %>%
  group_by(chr,pos,location,gen) %>%
  summarize(total_outside_chrs = sum(outside_chrs),
            total_samples = n(),
            percent_outside=total_outside_chrs/total_samples)

slopes_outside <- summarized_outside %>% 
  filter(gen != 7.5, gen != 2.5) %>%
  ungroup() %>%
  nest_by(chr,pos,location) %>% 
  mutate(model = list(lm(I(percent_outside-0) ~ 0 + gen, data = data))) %>% 
  summarise(tidy(model))

slopes_outside %>%
  #filter(chr == "Ha412HOChr15") %>%
  ggplot(.,aes(x=pos/1000000,y=estimate)) +
  geom_point() +
  facet_grid(location~chr)

slopes_outside %>%
  group_by(location) %>%
  mutate(rank = rank(estimate, ties.method = "first")) %>%
  group_by(chr,pos) %>%
  summarize(mean_rank = mean(rank)) %>%
  mutate(rolling_rank=rollapply(mean_rank,10,mean,align='right',fill=NA,partial=T)) %>%
  ggplot(.,aes(x=pos,y=mean_rank)) +
  geom_point(alpha=0.4) +
  facet_wrap(~chr) +
  geom_line(aes(x=pos,y=rolling_rank),color="blue") +
  theme_cowplot()

outside_correlations <- tibble()
outside_correlation_plots <- list()
counter=1
for (i in 1:4){
  for (j in 2:4){
    if (i >=j){next}
    location1 <- locations[i]
    location2 <- locations[j]
    tmp <- slopes_outside %>%
      filter(location == location1 | location == location2) %>%
      select(chr,pos,estimate) %>%
      pivot_wider(names_from = location,values_from = estimate)  %>%
      rename(loc1 = 3, loc2 = 4)
    correlation <- cor.test(tmp$loc1,tmp$loc2) 
    outside_correlations <- rbind(outside_correlations,
                              tibble(location_1 = location1, location_2 = location2,
                                     parent="outside",
                                     r = correlation$estimate[[1]], 
                                     pvalue = correlation$p.value[[1]]))
    plot <- ggplot(tmp, aes(x=loc1,y=loc2)) +
      geom_point() + 
      geom_smooth(method="lm") +
      theme_cowplot() +
      xlab(location1) + ylab(location2) +
      ggtitle("Outside")
    outside_correlation_plots[[counter]] <- plot
    counter=counter+1
    
  }
}

outside_correlation_plots[[1]] +
  outside_correlation_plots[[2]] + outside_correlation_plots[[3]] + outside_correlation_plots[[4]] +
  outside_correlation_plots[[5]] + outside_correlation_plots[[6]] 
