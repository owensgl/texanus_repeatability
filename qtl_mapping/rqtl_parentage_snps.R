library(tidyverse)
library(qtl)

chrlist <- paste0(sprintf("%02d", seq(17)))
marker_sets <- read_tsv("2021/texanus.BC1.linkagegroups.txt") 
linkagegroups <- read_tsv("/media/drive_5_usb/Childs/texanus_2021/meta/texanus.lgmanualcuration.txt") %>%
  mutate(chr = paste0("Ha412HOChr",sprintf("%02d",Chr))) %>%
  select(-Chr)
marker_info <- inner_join(marker_sets, linkagegroups) %>%
  mutate(marker= paste0(chr,".",pos))
pdf("2021/texanus.BC1.linkageplots_testgroups.v1.pdf")
for (i in 1:17){
  for (j in 1:4){
  
    chr <- chrlist[i]
    mapthis <- read.cross(format="csvr", "/media/drive_5_usb/Childs/texanus_2021/vcf/rqtl",paste("texanus.freebayes.BC1.10AF.50missing.5dp.bi.15-35MAF.chr",chr,".csvr",sep=""),estimate.map=FALSE,sep="\t")
    nmarkers <- summary.map(mapthis)$n.mar[1]
    min.markers <- round(nmarkers * 0.5)
    mapthis <- subset(mapthis, ind=(ntyped(mapthis)>min.markers))
    full_chrname <- paste0("Ha412HOChr",chr)
    marker_info %>%
      filter(chr == full_chrname) %>%
      filter(!is.na(haplotype)) %>%
      filter(haplotype == j) %>%
      arrange(haplotype,pos)  -> marker_order
    kept_markers <- marker_order %>% pull(marker)
    nt.bymar <- ntyped(mapthis, "mar")
    all_markers <- names(nt.bymar)
    todrop <- setdiff(all_markers, kept_markers)
    mapthis <- drop.markers(mapthis, todrop)
    mapthis <- est.rf(mapthis)
    # current_order <- tibble(marker = names(ntyped(mapthis, "mar"))) %>%
    #   mutate(n=row_number())
    # inner_join(marker_order,current_order) %>%
    #   pull(n) -> new_order
    print(
    plotRF(mapthis, chr=c(1:10),alternate.chrid=TRUE, main=paste0("Chr",chr," - Haplotype_",j))
    )
  }
}
dev.off()

snp_polarity <- read_tsv("/media/owens/Childs/texanus_2021/vcf/texanus.freebayes.BC1.10AF.50missing.5dp.bi.15-35MAF.SNPified_indels.txt") %>%
  select(-pos) %>% rename(pos = original_pos)


marker_info %>%
  inner_join(snp_polarity) %>%
  mutate(haplotype = case_when(haplotype == 1 ~ "deb",
                               haplotype == 2 ~ "ann1",
                               haplotype == 3 ~ "ann2",
                               haplotype == 4 ~ "ann3")) %>%
  filter(!is.na(haplotype)) %>%
  select(chr,pos,haplotype,ref,alt,minor) -> marker_together

write_tsv(marker_together, 
          "/media/owens/Childs/texanus_2021/vcf/texanus.freebayes.BC1.10AF.50missing.5dp.bi.15-35MAF.snpparentage.txt")

