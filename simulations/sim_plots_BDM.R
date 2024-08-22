library(tidyverse)

directory <- "/media/drive_5_usb/owens/bin/texanus/2021/slim_sims"

convergence.s <- tibble()
for (i in c(1:100)){
  i <- sprintf("%03d",i)
  for (j in c(0:9)){
    tmp <- read_delim(paste0(directory,"/","conv.sel.",i,".",j,".txt"),
                    col_names=c("convergence"),delim = "\t") %>%
      mutate(selection = j/10,rep=i)
    convergence.s <- rbind(convergence.s,tmp)
  }
}

convergence.s %>%
  ggplot(.,aes(x=as.factor(selection),y=convergence)) + geom_boxplot()
    
pdf("2021/selection_simulations.v2.pdf",height=5,width=5)
selection_sim_plot <- convergence.s %>%
  group_by(selection) %>%
  summarize(bottom = quantile(convergence, probs = 0.0275),
            mid = quantile(convergence, probs = 0.5),
            top = quantile(convergence,probs = 0.975)) %>%
  ggplot() + geom_point(aes(x=selection,y=mid),size=3) +
  geom_linerange(aes(x=selection,ymin=bottom,ymax=top),size=1.1) +
  theme_cowplot() + geom_hline(aes(yintercept=0),linetype="dotted") +
  ylab("Porportion variance due to selection") +
  xlab("DMI Selection")     
selection_sim_plot
dev.off()

convergence.r <- tibble()
for (i in c(1:100)){
  i <- sprintf("%03d",i)
  for (j in c("1e-3","1e-4","5e-4","1e-5","5e-5","1e-6","5e-6","5e-7")){
    tmp <- read_delim(paste0(directory,"/","conv.re.",i,".",j,".txt"),
                      col_names=c("convergence"),delim="\t") %>%
      mutate(recombination = as.numeric(j),rep=i)
    convergence.r <- rbind(convergence.r,tmp)
  }
}

pdf("2021/recombination_simulations.v2.pdf",height=5,width=5)

recomb_sim_plot <- convergence.r %>%
  group_by(recombination) %>%
  summarize(bottom = quantile(convergence, probs = 0.0275),
            mid = quantile(convergence, probs = 0.5),
            top = quantile(convergence,probs = 0.975)) %>%
  mutate(cm_length = recombination * 1000000) %>%
  ggplot() + geom_point(aes(x=cm_length,y=mid),size=3) +
  geom_linerange(aes(x=cm_length,ymin=bottom,ymax=top),size=1.1) +
  theme_cowplot() + geom_hline(aes(yintercept=0),linetype="dotted") +
  scale_x_continuous(trans='log10') +
  ylab("Porportion variance due to selection") +
  xlab("Chromosome length (cM)")
recomb_sim_plot
dev.off()