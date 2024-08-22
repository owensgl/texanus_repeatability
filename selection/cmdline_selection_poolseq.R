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

state_fixer <- tibble(state=c("00","01","11",NA),
                      genotype=c(0,1,2,NA))

locations <- read_tsv("/media/drive_5_usb/Childs/texanus_2021/meta/location.gps.txt")
location_distance <- tibble()
for (i in 1:4){
  for (j in 1:4){
    pop1 <- locations$location[i]
    pop2 <- locations$location[j]
    dist <- distGeo(c(locations$long[i],locations$lat[i]),
                    c(locations$long[j],locations$lat[j]))/1000
    tmp <- tibble(location_1 = pop1, location_2 = pop2,
                  geodist = dist)
    location_distance <- rbind(location_distance, tmp)
  }
}


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
  group_by(chr,pos,gen,location,parent) %>%
  summarize(summed_genotype = sum(genotype,na.rm=T),
            n_samples = n(),
            percent_overall = summed_genotype/n_samples) %>%
  filter(gen != 7.5, gen != 2.5) %>%
  group_by(chr,pos,gen,location) %>%
  mutate(total_chrs = sum(summed_genotype, na.rm=T)) %>%
  mutate(percent_parent = summed_genotype/total_chrs)

#Population sizes
#Use the harmonic mean
pop_sizes <- read_tsv("/media/drive_5_usb/owens/bin/texanus/2021/texanus_popsizes.txt")
harmonic_pop_size <- pop_sizes %>%
  group_by(Plot) %>%
  filter(PopSize > 0) %>%
  summarize(mean_size = round(harmonic.mean(PopSize))) 

args = commandArgs(trailingOnly=TRUE)
#input_string <- "BFL-ann2"
input_string = args[1]

split_string = str_split(input_string,"-")
location_chosen <- split_string[[1]][1]
parent_chosen <- split_string[[1]][2]
#Calculate selection using PoolSeq
selection_output <- tibble()

popsize <- harmonic_pop_size %>%
  filter(Plot == location_chosen) %>%
  pull(mean_size)
print(paste(location_chosen,parent_chosen))
tmp <- summed_ancestry %>%
  filter(location == location_chosen) %>%
  filter(parent == parent_chosen) %>%
  select(chr,pos,gen,percent_parent) %>%
  pivot_wider(names_from = gen,values_from=percent_parent)
tmp.matrix <- tmp %>%
  ungroup() %>%
  select(-location,-chr,-pos) %>%
  as.matrix()
output <- tmp %>%
  select(location,chr,pos) %>%
  mutate(parent = parent_chosen, s = NA, s_min = NA, s_max = NA)
for (i in 1:nrow(tmp)){
  print(i)
  sh_out <- estimateSH(tmp.matrix[i,], Ne=popsize, t=colnames(tmp.matrix),method="LLS",h=0.5)
  confidence_interval <- confint.estsh(sh_out)
  output$s[i] <- sh_out$s
  output$s_min[i] <- confidence_interval[1,1]
  output$s_max[i] <- confidence_interval[1,2]
}
selection_output <- rbind(selection_output, output)



write_tsv(selection_output, 
          paste0("/media/drive_5_usb/Childs/texanus_2021/vcf/poolseq_selection_output.",
                 location_chosen,".",parent_chosen,".txt.gz"))