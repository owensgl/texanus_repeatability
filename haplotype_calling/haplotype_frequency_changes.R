library(tidyverse)
library(tidymodels)
library(cowplot)
sample_info <- read_tsv("/media/drive_5_usb/Childs/texanus_2021/meta/texanus.sampleinfo.sequenced.txt",
                        col_names = c("sample","year","gen","loc_type","type","location",
                                      "folder","seq_1","seq_2","seq_type"))%>%
  filter(seq_type == "both" | seq_type == "GBS2") %>%
  filter(type == "BC")

allele_freqs <- tibble()

for (parent_chosen in c("deb","ann1","ann2","ann3")){
  
  alleles <- read_tsv(paste0("/media/drive_5_usb/Childs/texanus_2021/vcf/texanus.freebayes.parentagesnps.",
                             parent_chosen,".tidy.txt.gz"),
                      col_types = "cdcd") %>%
    inner_join(sample_info) %>% mutate(parent = parent_chosen)
  
  
  #Copy first generation to all locations
  first_gen <- tibble()
  for (loc in c("HCC","BFL","LBJ","KPC")){
    first_gen <- rbind(first_gen, alleles %>%
                         filter(gen == 1) %>%
                         mutate(location = loc))
  }
  alleles <- alleles %>%
    filter(gen != 1) %>%
    filter(gen != 2.5, gen != 7.5) %>%
    rbind(first_gen)
  
  freq <- alleles %>%
    group_by(chr,pos,gen,location,parent) %>%
    summarize(allele_freq = mean(genotype_count)/2) 
  
  allele_freqs <- rbind(allele_freqs, freq)
}

deb_freq %>%
  filter(chr == "Ha412HOChr02") %>%
  ggplot(.,aes(x=pos/1000000,y=allele_freq)) +
  geom_point(alpha=0.2) +
  facet_grid(location~gen) +
  theme_classic() +
  xlab("MBp") +
  ylab("Allele frequency") +
  geom_hline(yintercept=0.25,linetype="dashed")

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


chr_n <- "Ha412HOChr01"
parent_chosen <- "deb"

pdf("2021/plots/frequency_changes.all.v1.pdf",height=12,width=10)
for (n in c(1:17)){
  chr_n = paste0("Ha412HOChr",sprintf("%02d", n))
  n_full <- sprintf("%02d", n)
  for (parent_chosen in c("deb","ann1","ann2","ann3")){
    print(paste("printing ",chr_n," ",parent_chosen))
    
    
    n1 <- allele_freqs %>%
      filter(parent == parent_chosen) %>%
      filter(chr == chr_n) %>%
      ggplot(.,aes(x=pos/1000000,y=allele_freq)) +
      geom_point(alpha=0.2) +
      facet_grid(location~gen) +
      theme_classic() +
      xlab("MBp") +
      ylab("Allele frequency") +
      geom_hline(yintercept=0.25,linetype="dashed") +
      geom_hline(yintercept=0,linetype="solid") +
      ggtitle(paste0("Chr",n_full,":",parent_chosen," parental allele frequency")) +
      theme(axis.text.x = element_text(size=8))
    
    n2 <- summed_ancestry %>%
      filter(chr == chr_n) %>%
      filter(parent == parent_chosen) %>%
      ggplot(.,aes(x=pos/1000000,y=percent_parent)) +
      geom_point(alpha=0.2) +
      facet_grid(location~gen) +
      theme_classic() +
      xlab("MBp") +
      ylab("Allele frequency") +
      geom_hline(yintercept=0.25,linetype="dashed") +
      geom_hline(yintercept=0,linetype="solid") +
      ggtitle(paste0("Chr",n_full,":",parent_chosen," parental haplotype frequency")) +
      theme(axis.text.x = element_text(size=8))
    
    print(
      n1 / n2
    )
  }
}
dev.off()


ancestry_change <- summed_ancestry %>%
  select(chr, pos, gen, parent, percent_parent) %>%
  group_by(location, chr, pos, parent) %>%  # Group by location, chr, pos, and parent
  arrange(gen) %>%                          # Ensure that the data is ordered by generation
  mutate(diff_percent_parent = percent_parent - lag(percent_parent),
         prev_perc = lag(percent_parent),
         gen_dif = gen - lag(gen),
         change_rate = abs(diff_percent_parent)/gen_dif) %>%  # Calculate the difference
  ungroup()   

ancestry_change %>%
  ggplot(.,aes(x=as.factor(gen),y=change_rate,fill=parent)) +
  geom_boxplot() +
  facet_wrap(~location)
