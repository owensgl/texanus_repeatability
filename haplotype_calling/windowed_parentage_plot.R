library(tidyverse)
window <- read_tsv("/media/drive_5_usb/Childs/texanus_2021/vcf/texanus.freebayes.5dp.5mbwindow.parentage.txt",
                   col_types = "cccnncc",
                   col_names=c("parent","sample","chr","start","end","counts","call"),
                   skip=1)


info <- read_tsv("/media/drive_1/working/texanus/texanus.sampleinfo.txt",guess_max=10000)
info %>% select(ID, gen, location) -> info
colnames(info) <- c("sample","gen","location")
window <- inner_join(window, info)
bw_pallet <- c("#cccccc","#808080","#0c0c0c","#ffffff")
poplist <- window %>% select(location) %>% unique() %>% head(4) %>% pull()
chrlist <- window %>% select(chr) %>% unique()  %>% pull()
pdf("2021/texanus.freebayes.5dp.5mbwindow.bw.pdf",height=4,width=6)
for (i in 1:4){
  print(poplist[i])
  window %>% filter(location == poplist[i]) %>% select(gen) %>% unique() %>% pull() -> genlist
  for (j in 1:length(genlist)){
    print(genlist[j])
    for (k in 1:length(chrlist)){
      print(
        window %>%
          filter(location == poplist[i],gen == genlist[j]) %>% 
          filter(chr == chrlist[k]) %>%
          ggplot(.) +
          geom_segment(aes(x=start, xend=end,y=as.factor(sample),yend=as.factor(sample),color=call),size=2) +
          facet_grid(~parent) +
          scale_colour_manual(values=bw_pallet, 
                              name="Allelic state",
                              breaks=c("0", "1", "2","Draw"),
                              labels=c("0", "1", "2","Unknown")) +
          theme_minimal() +
          theme(axis.text.y=element_blank(),axis.ticks=element_blank(),
                axis.title.x=element_blank(),
                axis.title.y=element_blank()) +
          ggtitle(paste(poplist[i], genlist[j], chrlist[k],"5MB window", sep=" "))
        
      )
    }
  }
}
dev.off()