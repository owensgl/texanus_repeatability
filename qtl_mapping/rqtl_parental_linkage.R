library(qtl)
library(tidyverse)

deb_alleles <- read_tsv("/media/drive_5_usb/Childs/texanus_2021/vcf/texanus.freebayes.BC1.10AF.50missing.5dp.bi.15-35MAF.SNPified.deb.txt") %>%
  filter(deb_af > 0.5) %>% filter(minor == "alt") %>% mutate(sign = "deb") %>%
  select(chr,pos,sign)

chrlist <- paste0("chr",sprintf("%02d", seq(17)))
marker_identity <- data.frame(chr=as.character(),pos=as.numeric(),parent=as.character())
full_marker_set <- data.frame(chr=as.character(),pos=as.numeric(),lg=as.numeric())

pdf("2021/texanus.BC1.linkageplots.pdf")
for (i in 1:length(chrlist)){
  chr <- chrlist[i]
  mapthis <- read.cross(format="csvr", "/media/drive_5_usb/Childs/texanus_2021/vcf/rqtl",paste("texanus.freebayes.BC1.10AF.50missing.5dp.bi.15-35MAF.",chr,".csvr",sep=""),estimate.map=FALSE,sep="\t")
  nmarkers <- summary.map(mapthis)$n.mar[1]
  min.markers <- round(nmarkers * 0.5)
  mapthis <- subset(mapthis, ind=(ntyped(mapthis)>min.markers))
  nt.bymar <- ntyped(mapthis, "mar")
  todrop <- names(nt.bymar[nt.bymar < 45])
  mapthis <- drop.markers(mapthis, todrop)
  mapthis <- est.rf(mapthis)
  min_lod <- 16
  if (i == 17){
    min_lod <- 18
  }
  mapthis <- formLinkageGroups(mapthis, max.rf=0.10, min.lod=min_lod, reorgMarkers=TRUE)
  plotRF(mapthis, chr=c(1:10),alternate.chrid=TRUE, main=paste(chr))
  lg <- tibble(chr=as.character(),pos=as.numeric(),lg=as.numeric())
  for (j in 1:10){
    tmp <- markernames(mapthis, chr=j)
    tmp <- as.data.frame(tmp)
    tmp %>% separate(tmp, c("chr","pos")) -> tmp
    tmp$lg <- j
    tmp$pos <- as.numeric(tmp$pos)
    lg <- rbind(lg,tmp)
  }
  print(
    full_join(lg,deb_alleles) %>%
      mutate(sign=  case_when(is.na(sign) ~ "ann",
                              TRUE ~"deb")) %>%
      filter(!is.na(lg))%>%
      ggplot(.) + geom_bar(aes(x=as.factor(lg),fill=sign),position = "fill") +
        ylab("Proportion Debilis sites") +
        xlab("Linkage group") +
        ggtitle(paste(chr)) 
  )
  all_markers <- lg %>% select(chr, pos, lg)
  full_marker_set <- rbind(full_marker_set,all_markers)
}
dev.off()


#write_tsv(full_marker_set,"2021/texanus.BC1.linkagegroups.txt")



