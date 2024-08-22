library(tidyverse)
library(gdsfmt)
library(SNPRelate)
library(cowplot)

sample_info <- read_tsv("/media/drive_5_usb/Childs/texanus_2021/meta/texanus.sampleinfo.sequenced.txt",
                        col_names = c("sample","year","gen","loc_type","type","location",
                                      "folder","seq_1","seq_2","seq_type"))%>%
  filter(seq_type == "both" | seq_type == "GBS2") %>%
  filter(type == "BC")

# snpgdsVCF2GDS("/media/drive_5_usb/Childs/texanus_2021/vcf/texanus.freebayes.parentagesnps.noQTL.vcf", 
#               "/media/drive_5_usb/Childs/texanus_2021/vcf/texanus.freebayes.parentagesnps.noQTL.gds", 
#               method="biallelic.only")

snpgdsVCF2GDS("/media/drive_5_usb/Childs/texanus_2021/ancestry_hmm/allmarkers.vcf.gz", 
              "/media/drive_5_usb/Childs/texanus_2021/ancestry_hmm/allmarkers.gds", 
              method="biallelic.only")

genofile <- snpgdsOpen("/media/drive_5_usb/Childs/texanus_2021/ancestry_hmm/allmarkers.gds")
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,autosome.only = F)
snpset.id <- unlist(unname(snpset))

#NO LD filtering!
pca <- snpgdsPCA(genofile, num.thread=2, autosome.only = F,eigen.cnt=50)

tab <- data.frame(sample = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
pc.percent <- pca$varprop*100

first_gen <- tibble()
for (location_chosen in c("BFL","HCC","KPC","LBJ")){
  tmp <- inner_join(tab,sample_info) %>%
    filter(gen == 1) %>%
    mutate(location = location_chosen)
  first_gen <- rbind(first_gen, tmp)
}

inner_join(tab,sample_info) %>%
  filter(gen != 1) %>%
  rbind(first_gen) %>%
  group_by(location,gen) %>%
  ggplot(.,aes(x=EV1, y=EV2,color=as.factor(gen))) +
  geom_point() +
  scale_color_viridis_d() +
  facet_wrap(~location) +
  stat_ellipse()

pdf("2021/texanus.ancestry_hmm.markers.pca.v1.pdf",
    height=6,width=6)
inner_join(tab,sample_info) %>%
  filter(gen != 1) %>%
  filter(gen != 2.5, gen != 7.5) %>%
  rbind(first_gen) %>%
  group_by(location,gen) %>%
  summarize(mean_PC1 = mean(EV1),mean_PC2 = mean(EV2)) %>%
  ggplot(.,aes(x=mean_PC1, y=mean_PC2,color=as.factor(gen))) +
  geom_point(data=inner_join(tab,sample_info) %>%
               filter(gen != 2.5, gen != 7.5),
             aes(x=EV1,y=EV2),alpha=0.1) +
  geom_point(size=0.2) +
  scale_color_viridis_d(option="magma",name="Generation") +
  geom_path(aes(group=location),
            arrow = arrow(length = unit(0.3, "cm")),
            lineend = "round",size=1) +
  theme_cowplot() +
  xlab("PC1 (PVE=6.4%)") + ylab("PC2 (PVE=2.8%)") +
  theme(axis.title = element_text(size = 14))
dev.off()
inner_join(tab,sample_info) %>%
  filter(gen != 1) %>%
  filter(gen != 2.5, gen != 7.5) %>%
  rbind(first_gen) %>%
  group_by(location,gen) %>%
  summarize(mean_PC1 = mean(EV1),mean_PC2 = mean(EV2)) %>%
  ggplot(.,aes(x=mean_PC1, y=mean_PC2,color=as.factor(gen))) +
  geom_point(data=inner_join(tab,sample_info) %>%
               filter(gen != 2.5, gen != 7.5),
             aes(x=EV1,y=EV2),alpha=0.1) +
  geom_point(size=0.2) +
  scale_color_viridis_d(option="magma") +
  geom_path(aes(group=location),
            arrow = arrow(length = unit(0.3, "cm")),
            lineend = "round",size=1) +
  theme_cowplot() +
  facet_wrap(~location)

#Distance between timepoints.
usable_dimensions <- sum(pca$varprop > 0.01,na.rm=T)

tab <- as_data_frame(cbind(pca$sample.id, pca$eigenvect)) %>%
  rename(sample = V1)

first_gen <- tibble()
for (location_chosen in c("BFL","HCC","KPC","LBJ")){
  tmp <- inner_join(tab,sample_info) %>%
    filter(gen == 1) %>%
    mutate(location = location_chosen)
  first_gen <- rbind(first_gen, tmp)
}

tab <-inner_join(tab, sample_info) %>%
  filter(gen != 1) %>%
  rbind(first_gen) 

test <- tab %>%
  select(-sample)
test <- sapply( test, as.numeric ) 

pcs <- select(as.tibble(test), V2:V34)

pc_data <- cbind(tab$sample, pcs, tab$location, tab$gen ) %>%
  rename(sample=`tab$sample`,
         gen=`tab$gen`,
         location=`tab$location`) 

pc_data_sum <- pc_data %>%
  filter(gen != 2.5, gen != 7.5) %>%
  group_by(location,gen) %>%
  summarise(across(V2:V34, ~ mean(.x, na.rm = TRUE)))

euc_dist <- function(x1, x2){
  return(sqrt(sum((x1 - x2)^2)))
}
pca_distances <- tibble()
for (i in 1:(nrow(pc_data_sum)-1)){
  j <- i + 1
  location_1 <- pc_data_sum$location[i]
  location_2 <- pc_data_sum$location[j]
  if (location_1 != location_2){next}
  gen_1 <- pc_data_sum$gen[i]
  gen_2 <- pc_data_sum$gen[j]
  tmp_1 <- pc_data_sum[i,c(3:5)] %>%
    unlist(., use.names=FALSE)
  tmp_2 <- pc_data_sum[j,c(3:5)] %>%
    unlist(., use.names=FALSE)
  dist <- euc_dist(tmp_1,tmp_2)
  result <- tibble(location=location_1,gen_1 = gen_1, gen_2 = gen_2, pca_dist = dist)
  pca_distances <- rbind(pca_distances, result)
}

pdf("2021/PCA_distances.v1.pdf",height=5,width=5)
pca_distances %>%
  mutate(length = gen_2 - gen_1,
         norm_dist = pca_dist / length) %>%
  ggplot(.,aes(x=gen_1,y=norm_dist,color=location)) +
  geom_line(size=2) +
  theme_cowplot() +
  ylab("PCA distance") +
  xlab("Starting generation") +
  scale_color_brewer(palette = "Set1",name="Location") 
dev.off()

            