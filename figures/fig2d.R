library(ggplot2)
library(gcookbook)
library(dplyr)
library(ggpubr)

df0= read.csv("inext_input_eco.txt",sep='\t',header=TRUE)
df0$t <- as.numeric(as.character(df0$t))
pdf(file='inext_out_eco.pdf',width = 6, height = 9)
plot_theme =theme(panel.background = element_rect(fill = "white"),text = element_text(family = "sans",face="bold", size=12),legend.title=element_blank(),legend.text=element_text(size=12)) #,legend.position = 'none'
theme_set(plot_theme)

p<-ggplot() +
  geom_line(data=subset(df0,method="Extrapolation"),aes(x = t, y = qD,group=order,color=order),linetype=2,linewidth=0.5)+
  geom_line(data=subset(df0,method!="Extrapolation"),aes(x = t, y = qD,group=order,color=order),linetype=1,linewidth=0.5)+
  geom_hline(yintercept=254, lwd=0.5,lty=6,col='grey30')+
  geom_point(aes(2228,1574),color="#ED0000",size=3,shape=21,stroke = 0.3)+
  geom_point(aes(5790,3224),color="#42B540",size=3,shape=21,stroke = 0.3)+
  geom_point(aes(16442,19864),color="#00468B",size=3,shape=21,stroke = 0.3)+
  geom_point(aes(1490,1295),color="#925E9F",size=3,shape=21,stroke = 0.3)+
  geom_point(aes(670,1023),color="#0099B4",size=3,shape=21,stroke = 0.3)+
  geom_ribbon(data=df0,aes(ymin = qD_LCL,
                           ymax = qD_UCL,group=order,fill=order,x=t),
              alpha = 0.2)+
  #scale_linetype_manual(df0$method,values=c("twodash","solid","solid"))+
  scale_color_manual(values=c("#ED0000","#42B540","#00468B","#925E9F","#0099B4"))+
  scale_fill_manual(values=c("#ED0000","#42B540","#00468B","#925E9F","#0099B4"))+
  xlab('Number of genomes')+
  ylab('Number of RiPP families')




p
dev.off()
