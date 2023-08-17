library(UpSetR)
library(ggpubr)
library(tidyverse)
library(showtext)

font_path <-'/System/Library/Fonts/Supplemental/' %>% font_paths()
fonts <- font_files() %>% as_tibble()
font_add('Arial','/System/Library/Fonts/Supplemental/Arial.ttf')
font_families()

df = read.csv("upset_input_host_phage.txt",sep='\t',header=TRUE)
pdf(file='upset_input_host_phage.pdf',width = 9, height = 6)
par(mar = c(5, 6, 2, 2) + 0.5) 
plot_theme = theme_gray() %+replace%  
  theme(text = element_text(family = "sans",face="bold", size=12))
theme_set(plot_theme)
p<-upset(df,sets=c('CRISPR','tRNA','Prophage'),mb.ratio=c(0.7,0.3),order.by='freq',
         nsets=3,number.angles=0,point.size=3,line.size=1,mainbar.y.label='Number',
         sets.bar.color = ggpubr::get_palette('npg',3),
         main.bar.color = "#A19BA9",matrix.color="#A19BA9",
         sets.x.label='Number')
p
dev.off()
