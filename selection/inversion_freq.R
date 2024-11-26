#Plotting fst outliers, but actual calls. for pet and neg F1 maps

library(tidyverse)
library(ggthemes)
library(forcats)
library(ggrastr)
library(gridExtra)
library(scatterpie)
library(viridis)

#Load in all inversion genotypes
folder <- "MDS_outliers"
chosen_species <- "annuus"
chosen_species_file <- "annuus"
chosen_species_abbreviation <- c("Ann")
filtering <- "pcasites"
prefix <- "Ha412HO_inv.v3"
base_directory <- "/media/drive_5_usb/speedy/working/texanus"
labels <- read_tsv("/media/drive_5_usb/speedy/working/texanus/texanus.sampleinfo.txt",col_names = T) %>%
  rename(sample=ID)
inversion_list <- read_tsv("/media/drive_5_usb/owens/bin/wild_gwas_2018/MDS_outliers/Ha412HO/all/Ha412HO_inv.v3.txt") %>%
  rename(chr_n = chr) %>%
  mutate(chr = paste0("Ha412HOChr",sprintf("%02d", chr_n)),
         mds = paste0(direction,mds))
inversions <- inversion_list %>% filter(spe == chosen_species) %>%
  select(chr,mds) %>% unique()

#pdf("texanus.v3.inversions.pdf",height=5,width=9)

all_inv_calls <- data.frame()
for (n in (c(4,5,6,7,8))){
  chosen_chr <- inversions[n,1] %>% pull()
  chosen_mds <- inversions[n,2] %>% pull()
  
  
  
  
  
  fst_calls <- read_tsv(paste(base_directory,"/texanus",".", chosen_chr,".",chosen_mds,
                              ".txt.gz",sep=""),guess_max= 100000,
                        col_names = c("XRQchr","XRQpos","chr","pos","sample","genotype","depth"),
                        col_types = "cdcdccd") %>%
    filter(genotype == "00" | genotype == "01" | genotype == "11")
  
  min_depth_genotyping <- 5
  min_depth_visualizing <- 2
  fst_calls %>%
    filter(depth >= min_depth_genotyping) %>%
    group_by(sample,genotype) %>%
    dplyr::count() %>%
    filter(!is.na(genotype)) %>%
    group_by(genotype) %>%
    dplyr::count()-> check
  if (nrow(check) < 3){next}
  
  
  fst_calls %>%
    filter(depth >= min_depth_genotyping) %>%
    group_by(sample,genotype) %>%
    dplyr::count() %>%
    filter(!is.na(genotype)) %>%
    spread(genotype, n,fill=0) %>%
    mutate(total =  (`00` + `01` + `11`)*2,
           perc_1 = ((`11`*2) + `01`)/total,
           het = (`01`*2)/total,
           left_formula = ((-(2/3)*perc_1) + (2/3)),
           right_formula = (((2/3)*perc_1))) %>%
    mutate(fst_call = case_when(perc_1 < 0.5  & left_formula >= het ~ 0,
                                perc_1 < 0.5  & left_formula < het ~ 1,
                                perc_1 >= 0.5  & right_formula >= het ~ 2,
                                perc_1 >= 0.5  & right_formula < het ~ 1)) %>%
    mutate(inv = paste0(chosen_chr,".",chosen_mds)) -> fst_samples
  
  all_inv_calls <- rbind(all_inv_calls, fst_samples)
  
  
}
sample_types <- all_ancestry %>%
  filter(chr == "Ha412HOChr01", pos == 2000000) %>% filter(parent == "ann1") %>%
  select(sample, gen, location, seq_type)

all_inv_calls %>%
  inner_join(sample_types) %>%
  filter(inv != "Ha412HOChr14.neg4") %>%
  filter(seq_type == "both" | seq_type == "GBS2") %>%
  ggplot(.,aes(x=total)) +
  geom_histogram() +
  facet_wrap(~inv)


all_inv_calls %>%
  inner_join(sample_types) %>%
  filter(inv != "Ha412HOChr14.neg4") %>%
  filter(seq_type == "both" | seq_type == "GBS2") %>%
  group_by(inv) %>%
  summarize(mean_markers = mean(total), sd = sd(total))
#Identifying which haplotype has which inversion genotype.

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

#chr13.1

chr13.1 <- all_ancestry %>%
  filter(chr == "Ha412HOChr13", pos == 50000000) %>%
  filter(gen == 1) 

chr13.1 %>%
  filter(parent == "ann1") %>%
  inner_join(all_inv_calls %>% filter(inv == "Ha412HOChr13.pos1")) %>%
  ggplot(.,aes(x=genotype,y=fst_call)) +
  geom_count()

#chr13.2

chr13.2 <- all_ancestry %>%
  filter(chr == "Ha412HOChr13", pos == 150000000) %>%
  filter(gen == 1) 

chr13.2 %>%
  filter(parent == "ann2") %>%
  inner_join(all_inv_calls %>% filter(inv == "Ha412HOChr13.pos2")) %>%
  ggplot(.,aes(x=genotype,y=fst_call)) +
  geom_count()


#chr14.1

chr14.1 <- all_ancestry %>%
  filter(chr == "Ha412HOChr14", pos == 111000000) %>%
  filter(gen == 1) 

chr14.1 %>%
  filter(parent == "ann2") %>%
  inner_join(all_inv_calls %>% filter(inv == "Ha412HOChr14.neg1")) %>%
  ggplot(.,aes(x=genotype,y=fst_call)) +
  geom_count()


#chr14.2
#HOMOGENOUS

#chr15.1
chr15.1 <- all_ancestry %>%
  filter(chr == "Ha412HOChr15", pos == 120000000) %>%
  filter(gen == 1) 

chr15.1 %>%
  filter(parent == "ann1") %>%
  inner_join(all_inv_calls %>% filter(inv == "Ha412HOChr15.pos1")) %>%
  ggplot(.,aes(x=genotype,y=fst_call)) +
  geom_count()


####Plotting selection coefficients for the texanus type individual
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

#Chr13.1, ann1 (mismatch), ann2, ann3 (match) chr13 10-109
selection_output %>%
  filter(chr == "Ha412HOChr13") %>%
  filter(pos > 10000000, pos < 109000000) %>%
  filter(parent != "deb") %>%
  group_by(parent) %>%
  summarize(
    mean = mean(s),
    sd = sd(s),
    n = n(),
    se = sd/sqrt(n)  # Calculate standard error
  ) %>%
  # Create the plot
  ggplot(aes(x = parent, y = mean)) +
  geom_point(size = 3) +  # Mean points
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2) +  # Error bars
  labs(
    x = "Parent",
    y = "Mean Selection Coefficient",
    title = "Mean Selection Coefficient by Parent"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Chr13.2, ann1 (match), ann2, ann3 (mismatch) chr13 140-157
selection_output %>%
  filter(chr == "Ha412HOChr13") %>%
  filter(pos > 140000000, pos < 157000000) %>%
  filter(parent != "deb") %>%
  group_by(parent) %>%
  summarize(
    mean = mean(s),
    sd = sd(s),
    n = n(),
    se = sd/sqrt(n)  # Calculate standard error
  ) %>%
  # Create the plot
  ggplot(aes(x = parent, y = mean)) +
  geom_point(size = 3) +  # Mean points
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2) +  # Error bars
  labs(
    x = "Parent",
    y = "Mean Selection Coefficient",
    title = "Mean Selection Coefficient by Parent"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Chr14.1, ann3 (match), ann2, ann1 (mismatch) chr14 102-128
selection_output %>%
  filter(chr == "Ha412HOChr14") %>%
  filter(pos > 102000000, pos < 128000000) %>%
  filter(parent != "deb") %>%
  group_by(parent) %>%
  summarize(
    mean = mean(s),
    sd = sd(s),
    n = n(),
    se = sd/sqrt(n)  # Calculate standard error
  ) %>%
  # Create the plot
  ggplot(aes(x = parent, y = mean)) +
  geom_point(size = 3) +  # Mean points
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2) +  # Error bars
  labs(
    x = "Parent",
    y = "Mean Selection Coefficient",
    title = "Mean Selection Coefficient by Parent"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Chr15.1, ann3 (match), ann2, ann1 (mismatch) chr14 105-175
selection_output %>%
  filter(chr == "Ha412HOChr15") %>%
  filter(pos > 105000000, pos < 175000000) %>%
  filter(parent != "deb") %>%
  group_by(parent) %>%
  summarize(
    mean = mean(s),
    sd = sd(s),
    n = n(),
    se = sd/sqrt(n)  # Calculate standard error
  ) %>%
  # Create the plot
  ggplot(aes(x = parent, y = mean)) +
  geom_point(size = 3) +  # Mean points
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2) +  # Error bars
  labs(
    x = "Parent",
    y = "Mean Selection Coefficient",
    title = "Mean Selection Coefficient by Parent"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#Merged figure
tmp1 <- selection_output %>%
  filter(chr == "Ha412HOChr13") %>%
  filter(pos > 10000000, pos < 109000000) %>%
  filter(parent != "deb") %>%
  group_by(parent) %>%
  summarize(
    mean = mean(s),
    sd = sd(s),
    n = n(),
    se = sd/sqrt(n)  # Calculate standard error
  ) %>%
  mutate(inv = "ann13.01", match = case_when(parent == "ann1" ~ "mismatch",
                                             TRUE ~ "match"))

tmp2 <- selection_output %>%
  filter(chr == "Ha412HOChr13") %>%
  filter(pos > 140000000, pos < 157000000) %>%
  filter(parent != "deb") %>%
  group_by(parent) %>%
  summarize(
    mean = mean(s),
    sd = sd(s),
    n = n(),
    se = sd/sqrt(n)  # Calculate standard error
  )%>%
  mutate(inv = "ann13.02", match = case_when(parent == "ann1" ~ "match",
                                             TRUE ~ "mismatch"))

tmp3 <- selection_output %>%
  filter(chr == "Ha412HOChr14") %>%
  filter(pos > 102000000, pos < 128000000) %>%
  filter(parent != "deb") %>%
  group_by(parent) %>%
  summarize(
    mean = mean(s),
    sd = sd(s),
    n = n(),
    se = sd/sqrt(n)  # Calculate standard error
  )%>%
  mutate(inv = "ann14.01", match = case_when(parent == "ann3" ~ "match",
                                             TRUE ~ "mismatch"))

tmp4 <- selection_output %>%
  filter(chr == "Ha412HOChr15") %>%
  filter(pos > 105000000, pos < 175000000) %>%
  filter(parent != "deb") %>%
  group_by(parent) %>%
  summarize(
    mean = mean(s),
    sd = sd(s),
    n = n(),
    se = sd/sqrt(n)  # Calculate standard error
  ) %>%
  mutate(inv = "ann15.01", match = case_when(parent == "ann3" ~ "match",
                                             TRUE ~ "mismatch"))

pdf("2021/inversion_selection.v1.pdf",height=5,width=8)
rbind(tmp1, tmp2, tmp3, tmp4) %>%
  ggplot(aes(x = parent, y = mean,color=match)) +
  geom_point(size = 3) +  # Mean points
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2) +  # Error bars
  labs(
    x = "Parent",
    y = "Mean Selection Coefficient",
    title = expression(paste("Selection coefficient in ", italic("H. annuus"), " inversions"))
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~inv,scales="free") +
  scale_color_brewer(palette = "Set1", name="Texanus genotype match?",
                     labels=c("Matched","Mismatched"))
dev.off()
