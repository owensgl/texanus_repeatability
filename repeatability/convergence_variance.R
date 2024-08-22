#Change in allele frequency covariance from https://www.pnas.org/doi/10.1073/pnas.1919039117
library(tidyverse)
library(tidymodels)
library(cowplot)
sample_info <- read_tsv("/media/drive_5_usb/Childs/texanus_2021/meta/texanus.sampleinfo.sequenced.txt",
                        col_names = c("sample","year","gen","loc_type","type","location",
                                      "folder","seq_1","seq_2","seq_type"))%>%
  filter(seq_type == "both" | seq_type == "GBS2") %>%
  filter(type == "BC")
state_fixer <- tibble(state=c("00","01","11",NA),
                      genotype=c(0,1,2,NA))
parents <- c("deb","ann1","ann2","ann3")
locations <- c("BFL","HCC","KPC","LBJ")
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
  filter(gen != 2.5, gen != 7.5) %>%
  rbind(first_gen)


summed_ancestry <- all_ancestry %>%
  group_by(chr,pos,gen,location,parent) %>%
  summarize(summed_genotype = sum(genotype,na.rm=T),
            n_samples = n(),
            percent_overall = summed_genotype/n_samples) %>%
  group_by(chr,pos,gen,location) %>%
  mutate(total_chrs = sum(summed_genotype, na.rm=T)) %>%
  mutate(percent_parent = summed_genotype/total_chrs)

formatted_1_6 <- summed_ancestry %>%
  filter(gen == 1 | gen == 6) %>%
  filter(parent == "deb") %>%
  select(chr,pos,gen,location,total_chrs, percent_parent) %>%
  rename(pr = percent_parent,reads=total_chrs) %>%
  pivot_wider(names_from = c(gen,location),values_from = c(pr,reads))

write_tsv(formatted_1_6, "2021/deb_1_6_cvtk.v1.tsv")

formatted_1_6 <- summed_ancestry %>%
  filter(gen == 1 | gen == 6) %>%
  filter(parent == "ann1") %>%
  select(chr,pos,gen,location,total_chrs, percent_parent) %>%
  rename(pr = percent_parent,reads=total_chrs) %>%
  pivot_wider(names_from = c(gen,location),values_from = c(pr,reads))

write_tsv(formatted_1_6, "2021/ann1_1_6_cvtk.v1.tsv")

formatted_1_6 <- summed_ancestry %>%
  filter(gen == 1 | gen == 6) %>%
  filter(parent == "ann2") %>%
  select(chr,pos,gen,location,total_chrs, percent_parent) %>%
  rename(pr = percent_parent,reads=total_chrs) %>%
  pivot_wider(names_from = c(gen,location),values_from = c(pr,reads))

write_tsv(formatted_1_6, "2021/ann2_1_6_cvtk.v1.tsv")

formatted_1_6 <- summed_ancestry %>%
  filter(gen == 1 | gen == 6) %>%
  filter(parent == "ann3") %>%
  select(chr,pos,gen,location,total_chrs, percent_parent) %>%
  rename(pr = percent_parent,reads=total_chrs) %>%
  pivot_wider(names_from = c(gen,location),values_from = c(pr,reads))

write_tsv(formatted_1_6, "2021/ann3_1_6_cvtk.v1.tsv")
parent_colors <- c("#368a8d","#1f6285","#114069", "#d23d27")

covariance_results_cvtk <- tibble(parent=c("deb","ann1","ann2","ann3"),
                                  percent_convergence_mid = c(0.88118155,0.58804163,0.27483379,0.48093516),
                                  percent_convergence_top = c(0.90962051,0.60925361,0.33638771,0.52584025),
                                  percent_convergence_bottom = c(0.84758051,0.4961866,0.22024338, 0.38380524))
pdf("2021/convergence_variance_plot.v2.pdf",height=4,width=4)

covariance_results_cvtk %>%
  ggplot(.,aes(color=parent)) +
  geom_linerange(aes(x=parent, ymin=percent_convergence_bottom*100,ymax=percent_convergence_top*100),size=1.1) +
  geom_point(aes(x=parent,y=percent_convergence_mid*100),size=3) +
  xlab("Parent") + ylab("Percent variance explained\nby shared selection") +
  theme_cowplot() +
  scale_color_manual(values=parent_colors,name="Parent") +
  coord_cartesian(ylim=c(0,100))
dev.off()

##################
#Code not used, replaced with cvtk code.
##################
covariance_results_1_6 <- tibble()
for (i in 1:4){
  for (j in 2:4){
    loc_1 <- locations[i]
    loc_2 <- locations[j]
    if (i >= j){next}
    if (loc_1 == loc_2){next}
    for (parent_chosen in parents){
      x <- 1
      y <- 6
      
      tmp_1 <- summed_ancestry %>%
        filter(parent == parent_chosen, gen == x | gen == y, location == loc_1) %>%
        select(chr,pos,gen,percent_parent) %>%
        pivot_wider(names_from = gen,values_from = percent_parent) %>%
        rename(time_1 = 4,
               time_2 = 5) %>%
        mutate(dif = time_1 - time_2) %>% pull(dif)
      tmp_2 <- summed_ancestry %>%
        filter(parent == parent_chosen, gen == x | gen == y, location == loc_2) %>%
        select(chr,pos,gen,percent_parent) %>%
        pivot_wider(names_from = gen,values_from = percent_parent) %>%
        rename(time_1 = 4,
               time_2 = 5) %>%
        mutate(dif = time_1 - time_2) %>% pull(dif)
      covariance <- cov(tmp_1,tmp_2)
      variance_1 <- var(tmp_1)
      variance_2 <- var(tmp_2)
      sd = sqrt(variance_1 * variance_2)
      result_tmp <- tibble(parent = parent_chosen,
                           location_1 = loc_1,
                           location_2 = loc_2,
                           time = paste0(x,"-",y),
                           covariance = covariance,
                           var_1 = variance_1,
                           var_2 = variance_2,
                           sd = sd)
      covariance_results_1_6 <- rbind(covariance_results_1_6,result_tmp)
      print(paste(x,loc_1,loc_2,parent_chosen))
      
      
    }
    
  }
}

covariance_results_1_6 %>%
  mutate(pair_location = paste0(location_1,"-",location_2)) %>%
  group_by(parent,pair_location) %>%
  summarize(total_variance = mean(var_1 + var_2),
            covariance = mean(covariance),
            percent_convergent = covariance/total_variance) %>%
  ggplot(.,aes(x=parent,y=percent_convergent,color=parent)) +
  geom_jitter(width=0.1)

mean_variance <- covariance_results_1_6 %>%
  pivot_longer(c(var_1,var_2),names_to = "locations",values_to = "variance") %>%
  select(parent, variance) %>%
  distinct() %>%
  group_by(parent) %>%
  summarize(mean_variance = mean(variance))
covariance_results_1_6 %>%
  inner_join(mean_variance) %>%
  group_by(parent) %>%
  summarize(average_covariance = mean(covariance),
            average_sd = mean(sd),
            mean_variance = mean(mean_variance),
            convergent_correlation = average_covariance/average_sd,
            percent_convergent = average_covariance/mean_variance)
covariance_results_1_6 %>%
  mutate(percent_convergent = covariance/sd) %>%
  mutate(id = paste0(location_1,"-",location_2)) %>%
  ggplot(.,aes(x=parent,y=percent_convergent,color=id)) +
  geom_point()

#Bootstrapping somehow
subset_ancestry <- all_ancestry %>%
  filter(gen == 1 |  gen == 6)

#Doing it the stupid way :(
gens <- c(1,6)
parents <- c("ann1","ann2","ann3","ann4")



bootstrapped_1 <- read_tsv("/media/drive_5_usb/Childs/texanus_2021/ancestry_hmm/allmarkers.freqboot.gen1.txt.gz",
                           col_names = c("chr","pos","gen","location","parent",
                                         "percent_parent","rep"))
first_boot <- tibble()
for (loc in c("HCC","BFL","LBJ","KPC")){
  first_boot <- rbind(first_boot, bootstrapped_1 %>%
                       mutate(location = loc))
}

bootstrapped_6 <- read_tsv("/media/drive_5_usb/Childs/texanus_2021/ancestry_hmm/allmarkers.freqboot.gen6.txt.gz",
                           col_names = c("chr","pos","gen","location","parent",
                                         "percent_parent","rep"))

bootstrapped <- rbind(first_boot,bootstrapped_6)


covariance_results_1_6_boot <- tibble()
for (rep_chosen in 1:250){
  print(paste(rep_chosen))
  bootstrap_rep <- bootstrapped %>%
    filter(rep == rep_chosen)
  for (i in 1:4){
    for (j in 2:4){
      loc_1 <- locations[i]
      loc_2 <- locations[j]
      if (i >= j){next}
      if (loc_1 == loc_2){next}
      for (parent_chosen in parents){
        
        
        tmp_1 <- bootstrap_rep %>%
          group_by(chr,pos,parent) %>%
          filter(location == loc_1,parent==parent_chosen,) %>%
          select(chr,pos,gen,percent_parent) %>%
          pivot_wider(names_from = gen,values_from = percent_parent) %>%
          rename(time_1 = 4,
                 time_2 = 5) %>%
          mutate(dif = time_1 - time_2) %>% pull(dif)
        tmp_2 <- bootstrap_rep %>%
          group_by(chr,pos,parent) %>%
          filter(parent == parent_chosen, location == loc_2) %>%
          select(chr,pos,gen,percent_parent) %>%
          pivot_wider(names_from = gen,values_from = percent_parent) %>%
          rename(time_1 = 4,
                 time_2 = 5) %>%
          mutate(dif = time_1 - time_2) %>% pull(dif)
        covariance <- cov(tmp_1,tmp_2)
        variance_1 <- var(tmp_1)
        variance_2 <- var(tmp_2)
        sd = sqrt(variance_1 * variance_2)
        result_tmp <- tibble(parent = parent_chosen,
                             location_1 = loc_1,
                             location_2 = loc_2,
                             time = paste0(1,"-",6),
                             covariance = covariance,
                             var_1 = variance_1,
                             var_2 = variance_2,
                             sd = sd,
                             rep=rep_chosen)
        covariance_results_1_6_boot <- rbind(covariance_results_1_6_boot,result_tmp)
        
        
      }
    }
  }
}

#write_tsv(covariance_results_1_6_boot, "/media/drive_5_usb/Childs/texanus_2021/ancestry_hmm/allmarkers.freqboot.covariance.txt")
covariance_results_1_6_boot <- read_tsv("/media/drive_5_usb/Childs/texanus_2021/ancestry_hmm/allmarkers.freqboot.covariance.txt")

#Plotting without outside convergence

test_values <- covariance_results_1_6 %>%
  inner_join(mean_variance) %>%
  group_by(parent) %>%
  summarize(average_covariance = mean(covariance),
            average_sd = mean(sd),
            mean_variance = mean(mean_variance),
            convergent_correlation = average_covariance/average_sd,
            percent_convergent = average_covariance/mean_variance)

mean_variance <- covariance_results_1_6_boot %>%
  group_by(rep) %>%
  mutate(id = paste0(location_1,"-",location_2)) %>%
  pivot_longer(c(var_1,var_2),names_to = "locations",values_to = "variance") %>%
  select(parent, variance,id,rep) %>%
  distinct() %>%
  group_by(parent,rep,id) %>%
  summarize(mean_variance = mean(variance))

pdf("2021/convergence_variance_plot.pdf",height=4,width=4)
covariance_results_1_6_boot %>%
  mutate(id = paste0(location_1,"-",location_2)) %>%
  inner_join(mean_variance) %>%
  group_by(parent,rep) %>% 
  summarize(average_covariance = mean(covariance),
            average_sd = mean(sd),
            mean_variance = mean(mean_variance),
            convergent_correlation = average_covariance/average_sd,
            percent_convergent = average_covariance/mean_variance) %>%
  group_by(parent) %>%
  summarise(enframe(quantile(percent_convergent, c(0.025, 0.975)), "quantile", "percent_convergent")) %>%
  pivot_wider(names_from = quantile,values_from = percent_convergent) %>%
  rename(min5 = 2,
         max5 = 3) %>%
  ggplot(.,aes(color=parent)) +
  geom_linerange(aes(x=parent, ymin=min5*100,ymax=max5*100),position="dodge") +
  geom_point(data=test_values,aes(x=parent,y=percent_convergent*100)) +
  xlab("Parent") + ylab("Percent variance explained\nby shared selection") +
  theme_cowplot() +
  scale_color_manual(values=parent_colors,name="Parent") 
dev.off()
#Running test and bootstrap for outside geneflow.
covariance_results_1_6out <- tibble()
for (i in 1:4){
  for (j in 2:4){
    loc_1 <- locations[i]
    loc_2 <- locations[j]
    if (i >= j){next}
    if (loc_1 == loc_2){next}
    x <- 1
    y <- 6
    
    tmp_1  <- summed_ancestry %>%
      group_by(chr,pos,gen,location) %>%
      mutate(total_genotypes = sum(summed_genotype)) %>%
      summarize(outside_freq = 1-total_genotypes/(n_samples*2)) %>% 
      distinct() %>%
      filter(gen == x | gen == y, location == loc_1) %>%
      select(chr,pos,gen,outside_freq) %>%
      pivot_wider(names_from = gen,values_from = outside_freq) %>%
      rename(time_1 = 4,
             time_2 = 5) %>%
      mutate(dif = time_1 - time_2) %>% pull(dif)
    tmp_2 <- summed_ancestry %>%
      group_by(chr,pos,gen,location) %>%
      mutate(total_genotypes = sum(summed_genotype)) %>%
      summarize(outside_freq = 1-total_genotypes/(n_samples*2)) %>% 
      distinct() %>%
      filter(gen == x | gen == y, location == loc_2) %>%
      select(chr,pos,gen,outside_freq) %>%
      pivot_wider(names_from = gen,values_from = outside_freq) %>%
      rename(time_1 = 4,
             time_2 = 5) %>%
      mutate(dif = time_1 - time_2) %>% pull(dif)
    covariance <- cov(tmp_1,tmp_2)
    variance_1 <- var(tmp_1)
    variance_2 <- var(tmp_2)
    sd = sqrt(variance_1 * variance_2)
    result_tmp <- tibble(parent = "outside",
                         location_1 = loc_1,
                         location_2 = loc_2,
                         time = paste0(x,"-",y),
                         covariance = covariance,
                         var_1 = variance_1,
                         var_2 = variance_2,
                         sd = sd)
    covariance_results_1_6out <- rbind(covariance_results_1_6out,result_tmp)
    print(paste(x,loc_1,loc_2,"outside"))
    
    
    
  }
}
mean_variance <- covariance_results_1_6out %>%
  pivot_longer(c(var_1,var_2),names_to = "locations",values_to = "variance") %>%
  select(parent, variance) %>%
  distinct() %>%
  group_by(parent) %>%
  summarize(mean_variance = mean(variance))
covariance_results_1_6out %>%
  inner_join(mean_variance) %>%
  group_by(parent) %>%
  summarize(average_covariance = mean(covariance),
            average_sd = mean(sd),
            mean_variance = mean(mean_variance),
            convergent_correlation = average_covariance/average_sd,
            percent_convergent = average_covariance/mean_variance)

#Bootstrapped outside convergence

bootstrapped_1out <- read_tsv("/media/drive_5_usb/Childs/texanus_2021/ancestry_hmm/allmarkers.freqboot.gen1out.txt.gz",
                           col_names = c("chr","pos","gen","location","parent",
                                         "percent_parent","rep"))
first_boot <- tibble()
for (loc in c("HCC","BFL","LBJ","KPC")){
  first_boot <- rbind(first_boot, bootstrapped_1out %>%
                        mutate(location = loc))
}

bootstrapped_6out <- read_tsv("/media/drive_5_usb/Childs/texanus_2021/ancestry_hmm/allmarkers.freqboot.gen6out.txt.gz",
                           col_names = c("chr","pos","gen","location","parent",
                                         "percent_parent","rep"))

bootstrappedout <- rbind(first_boot,bootstrapped_6out)


covariance_results_1_6out_boot <- tibble()
for (rep_chosen in 1:250){
  print(paste(rep_chosen))
  bootstrap_rep <- bootstrappedout %>%
    filter(rep == rep_chosen)
  for (i in 1:4){
    for (j in 2:4){
      loc_1 <- locations[i]
      loc_2 <- locations[j]
      if (i >= j){next}
      if (loc_1 == loc_2){next}
      
      
      
      tmp_1 <- bootstrap_rep %>%
        group_by(chr,pos,parent) %>%
        filter(location == loc_1,parent=="outside") %>%
        select(chr,pos,gen,percent_parent) %>%
        pivot_wider(names_from = gen,values_from = percent_parent) %>%
        rename(time_1 = 4,
               time_2 = 5) %>%
        mutate(dif = time_1 - time_2) %>% pull(dif)
      tmp_2 <- bootstrap_rep %>%
        group_by(chr,pos,parent) %>%
        filter(parent == "outside", location == loc_2) %>%
        select(chr,pos,gen,percent_parent) %>%
        pivot_wider(names_from = gen,values_from = percent_parent) %>%
        rename(time_1 = 4,
               time_2 = 5) %>%
        mutate(dif = time_1 - time_2) %>% pull(dif)
      covariance <- cov(tmp_1,tmp_2)
      variance_1 <- var(tmp_1)
      variance_2 <- var(tmp_2)
      sd = sqrt(variance_1 * variance_2)
      result_tmp <- tibble(parent = "outside",
                           location_1 = loc_1,
                           location_2 = loc_2,
                           time = paste0(1,"-",6),
                           covariance = covariance,
                           var_1 = variance_1,
                           var_2 = variance_2,
                           sd = sd,
                           rep=rep_chosen)
      covariance_results_1_6out_boot <- rbind(covariance_results_1_6out_boot,result_tmp)
      
      
    }
  }
}
write_tsv(covariance_results_1_6out_boot, "/media/drive_5_usb/Childs/texanus_2021/ancestry_hmm/allmarkers.freqboot.out.covariance.txt")



#Final plotting
test_values <- rbind(covariance_results_1_6out %>%
                       inner_join(mean_variance) %>%
                       group_by(parent) %>%
                       summarize(average_covariance = mean(covariance),
                                 average_sd = mean(sd),
                                 mean_variance = mean(mean_variance),
                                 convergent_correlation = average_covariance/average_sd,
                                 percent_convergent = average_covariance/mean_variance),
                     covariance_results_1_6 %>%
                       inner_join(mean_variance) %>%
                       group_by(parent) %>%
                       summarize(average_covariance = mean(covariance),
                                 average_sd = mean(sd),
                                 mean_variance = mean(mean_variance),
                                 convergent_correlation = average_covariance/average_sd,
                                 percent_convergent = average_covariance/mean_variance))
  


mean_variance <- covariance_results_1_6_boot %>%
  rbind(covariance_results_1_6out_boot) %>%
  group_by(rep) %>%
  mutate(id = paste0(location_1,"-",location_2)) %>%
  pivot_longer(c(var_1,var_2),names_to = "locations",values_to = "variance") %>%
  select(parent, variance,id,rep) %>%
  distinct() %>%
  group_by(parent,rep,id) %>%
  summarize(mean_variance = mean(variance))
covariance_results_1_6_boot %>%
  rbind(covariance_results_1_6out_boot) %>%
  mutate(id = paste0(location_1,"-",location_2)) %>%
  inner_join(mean_variance) %>%
  group_by(parent,rep) %>% 
  summarize(average_covariance = mean(covariance),
            average_sd = mean(sd),
            mean_variance = mean(mean_variance),
            convergent_correlation = average_covariance/average_sd,
            percent_convergent = average_covariance/mean_variance) %>%
  group_by(parent) %>%
  summarise(enframe(quantile(percent_convergent, c(0.025, 0.975)), "quantile", "percent_convergent")) %>%
  pivot_wider(names_from = quantile,values_from = percent_convergent) %>%
  rename(min5 = 2,
         max5 = 3) %>%
  ggplot(.,aes(color=parent)) +
  geom_errorbar(aes(x=parent, ymin=min5*100,ymax=max5*100),position="dodge",
                width=0.2) +
  geom_point(data=test_values,aes(x=parent,y=percent_convergent*100)) +
  xlab("Parent") + ylab("Percent variance explained by shared selection") +
  theme_cowplot()



