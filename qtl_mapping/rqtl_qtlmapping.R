library(qtl)
library(tidyverse)
library(cowplot)
#This uses the scanone and scantwo functions to get qtl size and effect.
deb <- read.cross(format="csvr", "/media/drive_5_usb/Childs/texanus_2021/vcf/","texanus.freebayes.BC1qtl.parentagesnps.dp3.deb.combined.csvr",estimate.map=FALSE)
deb.maf <- read_tsv("/media/drive_5_usb/Childs/texanus_2021/vcf/texanus.freebayes.BC1qtl.parentagesnps.dp3.deb.combined.maf.txt")
deb.maf %>% 
  filter(percent_A > 0.2, percent_A < 0.8) %>%
  filter(genotyped > 40) %>% pull(marker)-> keeping_samples

todrop <- deb.maf %>% filter(!marker %in% keeping_samples) %>% pull(marker)
deb <- drop.markers(deb, todrop)
deb.sim <- sim.geno(deb, step=0, n.draws=1000, err=0.05)

ann <- read.cross(format="csvr", "/media/drive_5_usb/Childs/texanus_2021/vcf/","texanus.freebayes.BC1qtl.parentagesnps.dp3.ann2.combined.csvr",estimate.map=FALSE)
ann.maf <- read_tsv("/media/drive_5_usb/Childs/texanus_2021/vcf/texanus.freebayes.BC1qtl.parentagesnps.dp3.ann2.combined.maf.txt")
ann.maf %>% 
  filter(percent_A > 0.2, percent_A < 0.8) %>%
  filter(genotyped > 40) %>% pull(marker)-> keeping_samples

todrop <- ann.maf %>% filter(!marker %in% keeping_samples) %>% pull(marker)
ann <- drop.markers(ann, todrop)
ann.sim <- sim.geno(ann, step=0, n.draws=1000, err=0.05)




nperm <- 200
phenotypes <- colnames(ann.sim$pheno)
#pdf("2021/texanus.freebayes.BC1qtl.parentagesnps.dp3.deb.combined.pdf", height=10,width=10)
par(mfrow=c(2,2))
for (i in 1:25){
  deb.em <- scanone(deb.sim, pheno.col=i)
  deb.perm <- scanone(deb.sim, pheno.col=i,n.perm=nperm)
  
  plot(deb.em)
  title(paste("Debilis -",phenotypes[i]))
  abline(h=summary(deb.perm)[[1]],col="green",lty=2,lwd=2)
  abline(h=summary(deb.perm)[[2]],col="orange",lty=2,lwd=2)
  
  effect_result_deb <- effectscan(deb.sim,pheno.col=i,get.se=T,add.legend = F)
  title(paste("Debilis -",phenotypes[i],"Effect Size"))
  
  ann.em <- scanone(ann.sim, pheno.col=i)
  ann.perm <- scanone(ann.sim, pheno.col=i,n.perm=nperm)
  
  plot(ann.em)
  title(paste("Annuus BC parent -",phenotypes[i]))
  abline(h=summary(ann.perm)[[1]],col="green",lty=2,lwd=2)
  abline(h=summary(ann.perm)[[2]],col="orange",lty=2,lwd=2)
  
  effect_result_ann <- effectscan(ann.sim,pheno.col=i,get.se=T,add.legend = F)
  title(paste("Annuus BC parent -",phenotypes[i],"Effect Size"))
} 
dev.off()

ft_qtl <- scanone(deb.sim, pheno.col=4)
as_tibble(ft_qtl) %>%
 cbind(rownames(ft_qtl)) %>%
  rename(marker = "rownames(ft_qtl)") %>%
  filter(chr == "06") %>%
  separate(marker,c("chr","pos"),"_") %>%
  ggplot(.,aes(x=as.numeric(pos)/1000000,y=lod)) + geom_line() +
  theme_cowplot() +
  geom_vline(xintercept = 132,linetype="dashed") +
  xlab("MBp") + ylab("DaysToFlower LOD") +
  ggtitle("Chr06")


#Finding the max effect size
max_effects <- tibble()
for (i in 1:25){
  print(i)
  deb.em <- scanone(deb.sim, pheno.col=i)
  effect_result_deb <- effectscan(deb.sim,pheno.col=i,get.se=T,add.legend = F)
  max_effect_deb <- max(abs(effect_result_deb$a))
  tmp <- tibble(parent="F1",phenotype=phenotypes[i],effect_size = max_effect_deb)
  max_effects <- rbind(max_effects,tmp )
  
  ann.em <- scanone(ann.sim, pheno.col=i)
  effect_result_ann <- effectscan(ann.sim,pheno.col=i,get.se=T,add.legend = F)
  max_effect_ann <- max(abs(effect_result_ann$a))
  tmp <- tibble(parent="BC",phenotype=phenotypes[i],effect_size = max_effect_ann)
  max_effects <- rbind(max_effects,tmp )
  
} 
max_effects %>%
  group_by(phenotype) %>%
  mutate(min_size = min(effect_size),
         normalized_size = effect_size/min_size)) %>%
  ggplot(.,aes(x=parent,y=normalized_size)) + 
  geom_point() + 
  geom_line(aes(group=phenotype))
  
phenotype_categories <- read_tsv("/media/drive_5_usb/Childs/texanus_2021/meta/phenotype_categories.txt") %>%
  rename(phenotype=trait)
max_effects <- read_tsv("2021/max_qtl_effect_sizes.v1.txt")
max_effects_sum <- max_effects %>%
  group_by(phenotype) %>%
  pivot_wider(names_from = parent,values_from = effect_size) %>%
  mutate(deb_effect = F1/BC) %>%
  inner_join(phenotype_categories) 

#write_tsv(max_effects, "2021/max_qtl_effect_sizes.v1.txt")
pdf("2021/QTL_effect_size_comparison.v1.pdf",
    height=4,width=6)
qtl_plot <- max_effects_sum %>%
  ungroup() %>%
  summarize(mean_effect = mean(deb_effect),
            sd = sd(deb_effect,na.rm=T),
            se = sd/sqrt(n())) %>% 
  rename(deb_effect = mean_effect) %>%
  mutate(type = "all") %>%
  ggplot() +
  geom_point(aes(x=type,y=deb_effect),size=3) +
  geom_linerange(aes(x=type,ymin=deb_effect-2*se,
                     ymax=deb_effect+2*se),size=1.1) +
  geom_jitter(data=max_effects_sum,
              aes(x=type,y=deb_effect),
              width=0.1) +
  theme_cowplot() +
  geom_hline(yintercept=1,linetype="dashed") +
  ylab("F1 effect size / BC effect size") +
  xlab("Phenotype category")
dev.off()
  
t.test(x=max_effects_sum$deb_effect,mu=1)
  
###
#Finding intervals for QTL
i <- 7
deb.em <- scanone(deb.sim, pheno.col=i)
deb.perm <- scanone(deb.sim, pheno.col=i,n.perm=nperm)

good_seed_qtl_locations <- tibble()
qtl_chrs <- c("06","12","15","16","17")
for (qtl_chr in qtl_chrs){
  tmp <- bayesint(deb.em, qtl_chr, 0.95)
  tmp_strings <- str_split(rownames(tmp),"_")
  tmp_tibble <- tibble(chr=c(tmp_strings[[1]][1],tmp_strings[[1]][1],tmp_strings[[1]][1]),
                       type=c("start","max","end"),
                       pos=as.numeric(c(tmp_strings[[1]][2],tmp_strings[[2]][2],tmp_strings[[3]][2])))
  good_seed_qtl_locations <- rbind(good_seed_qtl_locations,tmp_tibble)
}
write_tsv(good_seed_qtl_locations,"2021/good_seed_deb_qtl.locations.txt")


####Finder all qtl intervals

deb.ehk<- scanone(deb.sim, pheno.col=1:25)
deb.perm <- scanone(deb.sim, pheno.col=1:25,n.perm=nperm)

ann.ehk<- scanone(ann.sim, pheno.col=1:25)
ann.perm <- scanone(ann.sim, pheno.col=1:25,n.perm=nperm)


allpeaks_deb <- summary(deb.ehk, format="allpeaks", perms=deb.perm, alpha=0.05, pvalues=TRUE)
allpeaks_ann <- summary(ann.ehk, format="allpeaks", perms=ann.perm, alpha=0.05, pvalues=TRUE)

tidy_peaks <- data.frame(trait=as.character(),chr=as.character(),
                         left_edge=as.numeric(),right_edge=as.numeric(),
                         peak=as.numeric())
for (i in 1:25){
  pcol <- (i-1)*3+4
  chrcol <- 1
  loccol <- (i-1)*3+2
  trait <- colnames(allpeaks_deb)[(i-1)*3+3]
  for (j in 1:nrow(allpeaks_deb)){
    pvalue <- allpeaks_deb[j, pcol]
    chr <- allpeaks_deb[j,chrcol]
    loc <- allpeaks_deb[j, loccol]
    if (pvalue < 0.05){
      interval <- bayesint(deb.ehk, chr, 0.95, lodcolumn=i)
      interval[3,2]
      tmp_peak <- data.frame(trait=as.character(trait),chr=as.character(chr),
                             left_edge=as.numeric(interval[1,2]),right_edge=as.numeric(interval[3,2]),
                             peak=as.numeric(loc),parent="deb")
      tidy_peaks <- rbind(tidy_peaks,tmp_peak)
      
    }
  }
}

for (i in 1:25){
  pcol <- (i-1)*3+4
  chrcol <- 1
  loccol <- (i-1)*3+2
  trait <- colnames(allpeaks_ann)[(i-1)*3+3]
  for (j in 1:nrow(allpeaks_ann)){
    pvalue <- allpeaks_ann[j, pcol]
    chr <- allpeaks_ann[j,chrcol]
    loc <- allpeaks_ann[j, loccol]
    if (pvalue < 0.05){
      interval <- bayesint(ann.ehk, chr, 0.95, lodcolumn=i)
      interval[3,2]
      tmp_peak <- data.frame(trait=as.character(trait),chr=as.character(chr),
                             left_edge=as.numeric(interval[1,2]),right_edge=as.numeric(interval[3,2]),
                             peak=as.numeric(loc),parent="ann")
      tidy_peaks <- rbind(tidy_peaks,tmp_peak)
      
    }
  }
}
tidy_peaks %>%
  inner_join(phenotypes %>% select(trait, full_name, class)) %>%
  write_tsv("2021/all_qtl_ranges.txt")

phenotypes <- read_tsv("/media/drive_5_usb/Childs/texanus_2021/meta/phenotype_full_names.txt")
tidy_peaks$chr <- factor(tidy_peaks$chr, levels= tidy_peaks$chr[order(as.numeric(as.character(tidy_peaks$chr)))])
pdf("2021/all_qtl_peaks.v1.pdf",height=8,width=16)
tidy_peaks %>% 
  inner_join(phenotypes) %>%
  ggplot(.) + 
  geom_segment(aes(x=left_edge,xend=right_edge,y=fct_reorder(as.factor(full_name),-order),
                   yend=fct_reorder(as.factor(full_name),-order),color=as.factor(class))) +
  geom_point(aes(x=peak,y=fct_reorder(as.factor(full_name),-order),color=as.factor(class))) +
  facet_grid(parent~chr) +
  scale_y_discrete(name="", limits = (levels(tidy_peaks$full_name))) +
  xlab("cM") + 
  theme_cowplot() +
  theme(panel.grid.major.x=element_line(linetype = "solid",color="light grey"),
        panel.grid.major.y=element_line(linetype = "dashed",color="light grey"),
        axis.text.x=element_text(size=8)) +
  scale_color_brewer(palette = "Set1",name="Trait\nCategory")
dev.off()
