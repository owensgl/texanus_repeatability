library(tidyverse)
library(cowplot)
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
all_genotypes <- all_genotypes %>%
  inner_join(sample_info)


first_gen <- tibble()
for (loc in c("HCC","BFL","LBJ","KPC")){
  first_gen <- rbind(first_gen, all_genotypes %>%
                       filter(gen == 1) %>%
                       mutate(location = loc))
}
all_genotypes <- all_genotypes %>%
  filter(gen != 1) %>%
  rbind(first_gen)

all_genotypes <- all_genotypes %>%
  filter(gen != 2.5, gen != 7.5)

summarized_genotypes <- all_genotypes %>%
  group_by(chr,pos,gen,location,parent) %>%
  mutate(genotyped = case_when(!is.na(genotype_count) ~ 1,
                               TRUE ~ 0))
summarized_genotypes <- summarized_genotypes %>% 
  summarize(summed_genotype = sum(as.numeric(genotype_count),na.rm=T),
            n_samples = sum(genotyped,na.rm=T)*2,
            allele_freq = summed_genotype/n_samples)

formatted_1_6 <- summarized_genotypes %>%
  filter(gen == 1 | gen == 6) %>%
  filter(parent == "deb") %>%
  select(chr,pos,gen,location,n_samples, allele_freq) %>%
  rename(pr = allele_freq,reads=n_samples) %>%
  pivot_wider(names_from = c(gen,location),values_from = c(pr,reads)) %>%
  filter(reads_6_BFL > 0, reads_6_HCC > 0, reads_6_LBJ > 0, reads_6_KPC > 0)

write_tsv(formatted_1_6, "2021/deb_1_6_snps_cvtk.v1.tsv")

formatted_1_6 <- summarized_genotypes %>%
  filter(gen == 1 | gen == 6) %>%
  filter(parent == "ann1") %>%
  select(chr,pos,gen,location,n_samples, allele_freq) %>%
  rename(pr = allele_freq,reads=n_samples) %>%
  pivot_wider(names_from = c(gen,location),values_from = c(pr,reads)) %>%
  filter(reads_6_BFL > 0, reads_6_HCC > 0, reads_6_LBJ > 0, reads_6_KPC > 0)


write_tsv(formatted_1_6, "2021/ann1_1_6_snps_cvtk.v1.tsv")

formatted_1_6 <- summarized_genotypes %>%
  filter(gen == 1 | gen == 6) %>%
  filter(parent == "ann2") %>%
  select(chr,pos,gen,location,n_samples, allele_freq) %>%
  rename(pr = allele_freq,reads=n_samples) %>%
  pivot_wider(names_from = c(gen,location),values_from = c(pr,reads)) %>%
  filter(reads_6_BFL > 0, reads_6_HCC > 0, reads_6_LBJ > 0, reads_6_KPC > 0)


write_tsv(formatted_1_6, "2021/ann2_1_6_snps_cvtk.v1.tsv")

formatted_1_6 <- summarized_genotypes %>%
  filter(gen == 1 | gen == 6) %>%
  filter(parent == "ann3") %>%
  select(chr,pos,gen,location,n_samples, allele_freq) %>%
  rename(pr = allele_freq,reads=n_samples) %>%
  pivot_wider(names_from = c(gen,location),values_from = c(pr,reads)) %>%
  filter(reads_6_BFL > 0, reads_6_HCC > 0, reads_6_LBJ > 0, reads_6_KPC > 0)


write_tsv(formatted_1_6, "2021/ann3_1_6_snps_cvtk.v1.tsv")

parent_colors <- c("#368a8d","#1f6285","#114069", "#d23d27")

covariance_results_cvtk <- tibble(parent=c("deb","ann1","ann2","ann3"),
                                  percent_convergence_bottom = c(0.8504246,0.64521931,0.56979981,0.58322454),
                                  percent_convergence_mid = c(0.88985396,0.73290524,0.65419923,0.6627292),
                                  percent_convergence_top = c(0.95448491, 0.76056189,0.6875503,0.72787547))
                                  
pdf("2021/convergence_variance_plot_snps.v2.pdf",height=4,width=4)

covariance_results_cvtk %>%
  ggplot(.,aes(color=parent)) +
  geom_linerange(aes(x=parent, ymin=percent_convergence_bottom*100,ymax=percent_convergence_top*100),size=1.1) +
  geom_point(aes(x=parent,y=percent_convergence_mid*100),size=3) +
  xlab("Parent") + ylab("Percent variance explained\nby shared selection") +
  theme_cowplot() +
  scale_color_manual(values=parent_colors,name="Parent")  +
  coord_cartesian(ylim=c(0,100))
dev.off()