setwd("~/PosDoc/UNICAMP/Codes/CompleteCode/Cluster/")

method = 2
snail_pop = 1000

n_pop_ga = 1000
n_gen = 30

folder  ="fixed_seed/mu_h_bad/"

fitness = matrix(0,n_pop_ga,n_gen)

for(i in 1:n_gen){
    dat = read.table(paste(folder,"/GA_result_fitness_",i,"_",method,"_",snail_pop,".dat",sep=""),h=F,stringsAsFactors = F)$V1
    fitness[,i] = dat
}

boxplot(fitness)

m = order(fitness[,i],decreasing = T)[1]
sort(fitness[,i],decreasing=T)
fitness[m,i]


dat = read.table(paste(folder,"/GA_result_pop_",i,"_",method,"_",snail_pop,".dat",sep=""),h=F,stringsAsFactors = F)
dat[,m]

a = c(0.00493120,0.00017434,0.01065208,7.49600000,0.29320000)
round(a,digits = 6)
for(i in 1:30){
  dat = read.table(paste(folder,"/GA_result_pop_",i,"_",method,"_",snail_pop,".dat",sep=""),h=F,stringsAsFactors = F)
  for(j in 1:ncol(dat)){
    s = sum(round(dat[,j],digits = 6) == round(a,digits = 6))
    if(s==5){
      print(c(i,j))
    }
  }
}


 dim(dados_pop) 
 dados_pop[,dados$V1==max(dados$V1)]
 
 
 ##########################################################################################################
 #################################################
 
 sizes = c(500,1000,2000)
 
 method = 2
 snail_pop = 500
 
 n_pop_ga = 1000
 n_gen = 17
 
 folder  ="fixed_seed/bkp1010"
 
 fitness = matrix(0,n_pop_ga,n_gen)
 
 for(i in 1:n_gen){
   dat = read.table(paste(folder,"/GA_result_fitness_",i,"_",method,"_",snail_pop,".dat",sep=""),h=F,stringsAsFactors = F)$V1
   fitness[,i] = dat
 }
 
 aux = stack(as.data.frame(fitness))
 
 bl = unique(as.vector(aux$ind))
 i = seq(1,17,2)
 bl[i]
 
 mp1 = ggplot()+geom_boxplot(data = aux,aes(x = ind,y = values),color = "purple",fill="purple1",size = 1.5,alpha = 0.5)+
   scale_x_discrete(name = TeX('Generations'),labels = as.character(i),breaks = as.factor(bl[i]))+
   scale_y_continuous(name = TeX("Fitness score"))+
   labs(tag = paste("A",sep=""))+
   theme_bw()+
   #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
   #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
   theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
           panel.grid.minor = element_blank(),
           #panel.border = element_blank(),
           panel.background = element_blank(),
           #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
           #text=element_text(family = "Tahoma"),
           axis.title = element_text(face="bold"),
           axis.text.x = element_text(colour="black", size = 35),
           axis.text.y = element_text(colour="black", size = 35,angle=90,hjust = 0.5),
           axis.title.y = element_text(colour="black", size = 40),
           axis.title.x = element_text(colour="black", size = 40),
           # plot.tag = element_text(size = 40),
           legend.text = element_text(size = 25),
           legend.title = element_text(size = 16),
           legend.position = "top",
           plot.tag.position = "topleft",
           plot.tag = element_text(size = 40),
           #axis.line = element_line(size=0.5, colour = "black")
   )
 mp1
 
 method = 2
 snail_pop = 1000
 
 n_pop_ga = 1000
 n_gen = 17
 
 folder  ="fixed_seed/mu_h_bad/"
 
 fitness = matrix(0,n_pop_ga,n_gen)
 
 for(i in 1:n_gen){
   dat = read.table(paste(folder,"/GA_result_fitness_",i,"_",method,"_",snail_pop,".dat",sep=""),h=F,stringsAsFactors = F)$V1
   fitness[,i] = dat
 }
 
 aux = stack(as.data.frame(fitness))
 
 bl = unique(as.vector(aux$ind))
 i = seq(1,17,2)
 bl[i]
 
 mp2 = ggplot()+geom_boxplot(data = aux,aes(x = ind,y = values),color = "red",fill="red1",size = 1.5,alpha = 0.5)+
   scale_x_discrete(name = TeX('Generations'),labels = as.character(seq(1,17,2)),breaks = as.factor(bl[i]))+
   scale_y_continuous(name = TeX("Fitness score"))+
   labs(tag = paste("B",sep=""))+
   theme_bw()+
   #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
   #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
   theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
           panel.grid.minor = element_blank(),
           #panel.border = element_blank(),
           panel.background = element_blank(),
           #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
           #text=element_text(family = "Tahoma"),
           axis.title = element_text(face="bold"),
           axis.text.x = element_text(colour="black", size = 35),
           axis.text.y = element_text(colour="black", size = 35,angle=90,hjust = 0.5),
           axis.title.y = element_text(colour="black", size = 40),
           axis.title.x = element_text(colour="black", size = 40),
           # plot.tag = element_text(size = 40),
           legend.text = element_text(size = 25),
           legend.title = element_text(size = 16),
           legend.position = "top",
           plot.tag.position = "topleft",
           plot.tag = element_text(size = 40),
           #axis.line = element_line(size=0.5, colour = "black")
   )
 mp2
 
 
 method = 2
 snail_pop = 2000
 
 n_pop_ga = 1000
 n_gen = 19
 
 folder  ="fixed_seed/mu_h_bad/"
 
 fitness = matrix(0,n_pop_ga,n_gen)
 
 for(i in 1:n_gen){
   dat = read.table(paste(folder,"/GA_result_fitness_",i,"_",method,"_",snail_pop,".dat",sep=""),h=F,stringsAsFactors = F)$V1
   fitness[,i] = dat
 }
 
 aux = stack(as.data.frame(fitness))
 
 bl = unique(as.vector(aux$ind))
 i = seq(1,19,2)
 bl[i]
 
 
 mp3 = ggplot()+geom_boxplot(data = aux,aes(x = ind,y = values),color = "blue",fill="blue1",size = 1.5,alpha = 0.5)+
   scale_x_discrete(name = TeX('Generations'),labels = as.character(i),breaks = as.factor(bl[i]))+
   scale_y_continuous(name = TeX("Fitness score"))+
   labs(tag = paste("C",sep=""))+
   theme_bw()+
   #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
   #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
   theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
           panel.grid.minor = element_blank(),
           #panel.border = element_blank(),
           panel.background = element_blank(),
           #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
           #text=element_text(family = "Tahoma"),
           axis.title = element_text(face="bold"),
           axis.text.x = element_text(colour="black", size = 35),
           axis.text.y = element_text(colour="black", size = 35,angle=90,hjust = 0.5),
           axis.title.y = element_text(colour="black", size = 40),
           axis.title.x = element_text(colour="black", size = 40),
           # plot.tag = element_text(size = 40),
           legend.text = element_text(size = 25),
           legend.title = element_text(size = 16),
           legend.position = "top",
           plot.tag.position = "topleft",
           plot.tag = element_text(size = 40),
           #axis.line = element_line(size=0.5, colour = "black")
   )
 mp3
 
 t1 = theme(axis.text.x = element_text(colour="black", size = 32),
            axis.text.y = element_text(colour="black", size = 35,angle=90,hjust = 0.5),
            axis.title.y = element_text(colour="black", size = 45),
            axis.title.x = element_text(colour="black", size = 40))
 
 mpt = ggarrange(mp1+t1,mp2+t1+rremove("y.title"),mp3+t1+rremove("y.title"),ncol=3,nrow=1,widths = c(0.83,0.8,0.8))
 mpt
 
 ggsave(paste("Fitness.png",sep=""),plot=mpt,device="png",width = 50,height = 25,units = "cm")
 
 ##############################################################################################################
 ##############################33 Plot prevalence ##############################################################
 ##########################################################################################################
 
 getwd()
 setwd("~/PosDoc/UNICAMP/Codes/CompleteCode/")
 
 library(ggplot2)
 library(latex2exp)
 library(ggpubr)
 
 method = 2

 
 snail_pop = 500
 tag_ = "A"
 file = 0
 n_pop_ga = 500
 n_gen = 7
 n_boots = 150
 
 folder = paste("result_",snail_pop,"_method_",method,sep="")
 
 matrix = read.table(paste(folder,"/matrix_",file,"_",snail_pop,".dat",sep = ""),h=F,stringsAsFactors = F)
 matrix = as.matrix(matrix)
 prevalence = read.table("prevalence_field.dat",h=F,stringsAsFactors = F)

 m = matrix(0,nrow(matrix),3)
for(i in 1:nrow(matrix)){
  m[i,1] = quantile(matrix[i,],0.025,names=F)
  m[i,2] = quantile(matrix[i,],0.975,names=F)
  m[i,3] = quantile(matrix[i,],0.5,names=F)
} 
 
df1 = data.frame(group = prevalence$V1, CI1=m[,1], CI2=m[,2], mediana = m[,3])

colnames(prevalence) = c("group","KK","HTX")
prevalence = data.frame(prevalence)


mp=ggplot()+geom_ribbon(data=df1,aes(x=1:length(group),ymin=CI1,ymax=CI2),color = "blue",fill = "blue",alpha = 0.5,size = 0.5)+
  geom_line(data=df1,aes(x=1:length(group),y = mediana),color = "blue",size = 1.5)+
  geom_line(data=prevalence,aes(x=1:length(group),y = HTX),color = "red",size = 1.5)+
  geom_point(data=prevalence,aes(x=1:length(group),y = HTX),color = "red",fill="red",size = 2.5)+
  scale_x_continuous(name = 'Age groups',labels = prevalence$group,breaks = 1:length(prevalence$group))+
  scale_y_continuous(name = "Prevalence",limits = c(0,1))+
  labs(tag = paste(tag_,sep=""))+
  theme_bw()+
  #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25,angle=45,hjust = 1.0),
          axis.text.y = element_text(colour="black", size = 25,angle=90,hjust = 0.5),
          axis.title.y = element_text(colour="black", size = 25),
          axis.title.x = element_text(colour="black", size = 25),
          # plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 16),
          legend.position = "top",
          plot.tag.position = "topleft",
          plot.tag = element_text(size = 40),
          #axis.line = element_line(size=0.5, colour = "black")
  )
mp

snail_pop = 1000
tag_ = "B"
file = 1
n_pop_ga = 500
n_gen = 7
n_boots = 150

folder = paste("result_",snail_pop,"_method_",method,sep="")

matrix = read.table(paste(folder,"/matrix_",file,"_",snail_pop,".dat",sep = ""),h=F,stringsAsFactors = F)
matrix = as.matrix(matrix)
prevalence = read.table("prevalence_field.dat",h=F,stringsAsFactors = F)

m = matrix(0,nrow(matrix),3)
for(i in 1:nrow(matrix)){
  m[i,1] = quantile(matrix[i,],0.025,names=F)
  m[i,2] = quantile(matrix[i,],0.975,names=F)
  m[i,3] = quantile(matrix[i,],0.5,names=F)
} 

df1 = data.frame(group = prevalence$V1, CI1=m[,1], CI2=m[,2], mediana = m[,3])

colnames(prevalence) = c("group","KK","HTX")
prevalence = data.frame(prevalence)


mp1=ggplot()+geom_ribbon(data=df1,aes(x=1:length(group),ymin=CI1,ymax=CI2),color = "blue",fill = "blue",alpha = 0.5,size = 0.5)+
  geom_line(data=df1,aes(x=1:length(group),y = mediana),color = "blue",size = 1.5)+
  geom_line(data=prevalence,aes(x=1:length(group),y = HTX),color = "red",size = 1.5)+
  geom_point(data=prevalence,aes(x=1:length(group),y = HTX),color = "red",fill="red",size = 2.5)+
  scale_x_continuous(name = TeX('Age groups'),labels = prevalence$group,breaks = 1:length(prevalence$group))+
  scale_y_continuous(name = TeX("Prevalence"),limits = c(0,1))+
  labs(tag = paste(tag_,sep=""))+
  theme_bw()+
  #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25,angle=45,hjust = 1.0),
          axis.text.y = element_text(colour="black", size = 25,angle=90,hjust = 0.5),
          axis.title.y = element_text(colour="black", size = 25),
          axis.title.x = element_text(colour="black", size = 25),
          # plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 16),
          legend.position = "top",
          plot.tag.position = "topleft",
          plot.tag = element_text(size = 40),
         # axis.line = element_line(size=0.5, colour = "black")
  )
mp1

snail_pop = 2000
tag_ = "C"
file = 0
n_pop_ga = 500
n_gen = 7
n_boots = 150

folder = paste("result_",snail_pop,"_method_",method,sep="")

matrix = read.table(paste(folder,"/matrix_",file,"_",snail_pop,".dat",sep = ""),h=F,stringsAsFactors = F)
matrix = as.matrix(matrix)
prevalence = read.table("prevalence_field.dat",h=F,stringsAsFactors = F)

m = matrix(0,nrow(matrix),3)
for(i in 1:nrow(matrix)){
  m[i,1] = quantile(matrix[i,],0.025,names=F)
  m[i,2] = quantile(matrix[i,],0.975,names=F)
  m[i,3] = quantile(matrix[i,],0.5,names=F)
} 

df1 = data.frame(group = prevalence$V1, CI1=m[,1], CI2=m[,2], mediana = m[,3])

colnames(prevalence) = c("group","KK","HTX")
prevalence = data.frame(prevalence)


mp2=ggplot()+geom_ribbon(data=df1,aes(x=1:length(group),ymin=CI1,ymax=CI2),color = "blue",fill = "blue",alpha = 0.5,size = 0.5)+
  geom_line(data=df1,aes(x=1:length(group),y = mediana),color = "blue",size = 1.5)+
  geom_line(data=prevalence,aes(x=1:length(group),y = HTX),color = "red",size = 1.5)+
  geom_point(data=prevalence,aes(x=1:length(group),y = HTX),color = "red",fill="red",size = 2.5)+
  scale_x_continuous(name = TeX('Age groups'),labels = prevalence$group,breaks = 1:length(prevalence$group))+
  scale_y_continuous(name = TeX("Prevalence"),limits = c(0,1))+
  labs(tag = paste(tag_,sep=""))+
  theme_bw()+
  #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25,angle=45,hjust = 1.0),
          axis.text.y = element_text(colour="black", size = 25,angle=90,hjust = 0.5),
          axis.title.y = element_text(colour="black", size = 25),
          axis.title.x = element_text(colour="black", size = 25),
          # plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 16),
          legend.position = "top",
          plot.tag.position = "topleft",
          plot.tag = element_text(size = 40),
          #axis.line = element_line(size=0.5, colour = "black")
  )
mp2

mpt = ggarrange(mp,mp1+rremove("y.text")+rremove("y.title"),mp2+rremove("y.text")+rremove("y.title"),ncol=3,nrow=1,widths = c(1,0.8,0.8))

ggsave(paste("prevalence.png",sep=""),plot=mpt,device="png",width = 45,height = 25,units = "cm")


ggsave(paste("prevalence_",snail_pop,"_",method,".png",sep=""),plot=mp,device="png",width = 45,height = 25,units = "cm")

ggsave(paste("prevalence_",snail_pop,"_",method,".eps",sep=""),plot=mp,device="eps",width = 25,height = 25,units = "cm")




##############################################################################################################
##############################33 Plot time ##############################################################
##########################################################################################################

setwd("/data/thomas/schisto_treat/")

library(ggplot2)
library(latex2exp)
library(ggpubr)
method = 2
snail_pop = 500
tag_ = "A"
file = 1
n_pop_ga = 500
n_gen = 7
n_boots = 150

folder = paste("result_",snail_pop,"_method_",method,sep="")

matrix = as.matrix(read.table(paste(folder,"/matrix_time_",file,"_",snail_pop,".dat",sep = ""),h=F,stringsAsFactors = F))


m = matrix(0,nrow(matrix),4)
for(i in 1:nrow(matrix)){
  m[i,1] = quantile(matrix[i,],0.025,names=F)
  m[i,2] = quantile(matrix[i,],0.975,names=F)
  m[i,3] = quantile(matrix[i,],0.5,names=F)
  m[i,4] = mean(matrix[i,])
} 

df1 = data.frame(time_d= seq(1,nrow(m)), CI1=m[,1], CI2=m[,2], mediana = m[,3],mean = m[,4])


mp1=ggplot()+geom_ribbon(data=df1,aes(x=time_d/365,ymin=CI1,ymax=CI2),color = "lightgrey",fill = "lightgrey",alpha = 0.5,size = 0.5)+
  geom_line(data=df1,aes(x=time_d/365,y = mean),color = "black",size = 1.5)+
  scale_x_continuous(name = 'Time (years)',limits = c(0,60))+
  scale_y_continuous(name = "Number of infected individuals")+
  #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  labs(tag = paste(tag_,sep=""))+
  theme_bw()+
  #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25,angle=90,hjust = 0.5),
          axis.title.y = element_text(colour="black", size = 25),
          axis.title.x = element_text(colour="black", size = 25),
          # plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 16),
          legend.position = "top",
          plot.tag.position = "topleft",
          plot.tag = element_text(size = 40),
          #axis.line = element_line(size=0.5, colour = "black")
  )
mp1


method = 2
snail_pop = 1000
tag_ = "B"
file = 1
n_pop_ga = 500
n_gen = 7
n_boots = 150

folder = paste("result_",snail_pop,"_method_",method,sep="")

matrix = as.matrix(read.table(paste(folder,"/matrix_time_",file,"_",snail_pop,".dat",sep = ""),h=F,stringsAsFactors = F))


m = matrix(0,nrow(matrix),4)
for(i in 1:nrow(matrix)){
  m[i,1] = quantile(matrix[i,],0.025,names=F)
  m[i,2] = quantile(matrix[i,],0.975,names=F)
  m[i,3] = quantile(matrix[i,],0.5,names=F)
  m[i,4] = mean(matrix[i,])
} 

df1 = data.frame(time_d= seq(1,nrow(m)), CI1=m[,1], CI2=m[,2], mediana = m[,3],mean = m[,4])


mp2=ggplot()+geom_ribbon(data=df1,aes(x=time_d/365,ymin=CI1,ymax=CI2),color = "lightgrey",fill = "lightgrey",alpha = 0.5,size = 0.5)+
  geom_line(data=df1,aes(x=time_d/365,y = mean),color = "black",size = 1.5)+
  scale_x_continuous(name = TeX('Time (years)'),limits = c(0,60))+
  scale_y_continuous(name = TeX("Number of infected individuals"))+
  #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  labs(tag = paste(tag_,sep=""))+
  theme_bw()+
  #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25,angle=90,hjust = 0.5),
          axis.title.y = element_text(colour="black", size = 25),
          axis.title.x = element_text(colour="black", size = 25),
          # plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 16),
          legend.position = "top",
          plot.tag.position = "topleft",
          plot.tag = element_text(size = 40),
          #axis.line = element_line(size=0.5, colour = "black")
  )
mp2


method = 2
snail_pop = 2000
tag_ = "C"
file = 1
n_pop_ga = 500
n_gen = 7
n_boots = 150

folder = paste("result_",snail_pop,"_method_",method,sep="")

matrix = as.matrix(read.table(paste(folder,"/matrix_time_",file,"_",snail_pop,".dat",sep = ""),h=F,stringsAsFactors = F))


m = matrix(0,nrow(matrix),4)
for(i in 1:nrow(matrix)){
  m[i,1] = quantile(matrix[i,],0.025,names=F)
  m[i,2] = quantile(matrix[i,],0.975,names=F)
  m[i,3] = quantile(matrix[i,],0.5,names=F)
  m[i,4] = mean(matrix[i,])
} 

df1 = data.frame(time_d= seq(1,nrow(m)), CI1=m[,1], CI2=m[,2], mediana = m[,3],mean = m[,4])


mp3=ggplot()+geom_ribbon(data=df1,aes(x=time_d/365,ymin=CI1,ymax=CI2),color = "lightgrey",fill = "lightgrey",alpha = 0.5,size = 0.5)+
  geom_line(data=df1,aes(x=time_d/365,y = mean),color = "black",size = 1.5)+
  scale_x_continuous(name = TeX('Time (years)'),limits = c(0,60))+
  scale_y_continuous(name = TeX("Number of infected individuals"))+
  #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  labs(tag = paste(tag_,sep=""))+
  theme_bw()+
  #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25,angle=90,hjust = 0.5),
          axis.title.y = element_text(colour="black", size = 25),
          axis.title.x = element_text(colour="black", size = 25),
          # plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 16),
          legend.position = "top",
          plot.tag.position = "topleft",
          plot.tag = element_text(size = 40),
          #axis.line = element_line(size=0.5, colour = "black")
  )
mp3


#mpt = ggarrange(mp,mp1+rremove("y.text")+rremove("y.title"),mp2+rremove("y.text")+rremove("y.title"),ncol=3,nrow=1,widths = c(1,0.8,0.8))
mpt = ggarrange(mp1,mp2+rremove("y.title"),mp3+rremove("y.title"),ncol=3,nrow=1,widths = c(0.9,0.8,0.8))

ggsave(paste("time_series_",snail_pop,"_",method,".png",sep=""),plot=mp,device="png",width = 25,height = 25,units = "cm")

ggsave(paste("time_series.png",sep=""),plot=mpt,device="png",width = 45,height = 25,units = "cm")


##################################################################################################
########################3 Plot time with intervention
###################################################################################################



setwd("~/PosDoc/UNICAMP/Codes/CompleteCode/")

library(ggplot2)
library(latex2exp)
library(ggpubr)
method = 2
snail_pop = 500
tag_ = "A"
file = 1
n_pop_ga = 500
n_gen = 7
n_boots = 150
rounds = 10
interval = 1

folder = paste("result_",snail_pop,"_method_",method,"_",rounds,"_",interval,sep="")

matrix = as.matrix(read.table(paste(folder,"/matrix_time_",file,"_",snail_pop,".dat",sep = ""),h=F,stringsAsFactors = F))


m = matrix(0,nrow(matrix),4)
for(i in 1:nrow(matrix)){
  m[i,1] = quantile(matrix[i,],0.025,names=F)
  m[i,2] = quantile(matrix[i,],0.975,names=F)
  m[i,3] = quantile(matrix[i,],0.5,names=F)
  m[i,4] = mean(matrix[i,])
} 

df1 = data.frame(time_d= seq(1,nrow(m)), CI1=m[,1], CI2=m[,2], mediana = m[,3],mean = m[,4])


mp1=ggplot()+geom_ribbon(data=df1,aes(x=time_d/365,ymin=CI1,ymax=CI2),color = "lightgrey",fill = "lightgrey",alpha = 0.5,size = 0.5)+
  geom_line(data=df1,aes(x=time_d/365,y = mean),color = "black",size = 1.5)+
  scale_x_continuous(name = 'Time (years)',limits = c())+
  scale_y_continuous(name = "Number of infected individuals")+
  #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  labs(tag = paste(tag_,sep=""))+
  theme_bw()+
  #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25,angle=90,hjust = 0.5),
          axis.title.y = element_text(colour="black", size = 25),
          axis.title.x = element_text(colour="black", size = 25),
          # plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 16),
          legend.position = "top",
          plot.tag.position = "topleft",
          plot.tag = element_text(size = 40),
          #axis.line = element_line(size=0.5, colour = "black")
  )
mp1


method = 2
snail_pop = 1000
tag_ = "B"
file = 1
n_pop_ga = 500
n_gen = 7
n_boots = 150

folder = paste("result_",snail_pop,"_method_",method,sep="")

matrix = as.matrix(read.table(paste(folder,"/matrix_time_",file,"_",snail_pop,".dat",sep = ""),h=F,stringsAsFactors = F))


m = matrix(0,nrow(matrix),4)
for(i in 1:nrow(matrix)){
  m[i,1] = quantile(matrix[i,],0.025,names=F)
  m[i,2] = quantile(matrix[i,],0.975,names=F)
  m[i,3] = quantile(matrix[i,],0.5,names=F)
  m[i,4] = mean(matrix[i,])
} 

df1 = data.frame(time_d= seq(1,nrow(m)), CI1=m[,1], CI2=m[,2], mediana = m[,3],mean = m[,4])


mp2=ggplot()+geom_ribbon(data=df1,aes(x=time_d/365,ymin=CI1,ymax=CI2),color = "lightgrey",fill = "lightgrey",alpha = 0.5,size = 0.5)+
  geom_line(data=df1,aes(x=time_d/365,y = mean),color = "black",size = 1.5)+
  scale_x_continuous(name = TeX('Time (years)'),limits = c(0,60))+
  scale_y_continuous(name = TeX("Number of infected individuals"))+
  #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  labs(tag = paste(tag_,sep=""))+
  theme_bw()+
  #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25,angle=90,hjust = 0.5),
          axis.title.y = element_text(colour="black", size = 25),
          axis.title.x = element_text(colour="black", size = 25),
          # plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 16),
          legend.position = "top",
          plot.tag.position = "topleft",
          plot.tag = element_text(size = 40),
          #axis.line = element_line(size=0.5, colour = "black")
  )
mp2


method = 2
snail_pop = 2000
tag_ = "C"
file = 1
n_pop_ga = 500
n_gen = 7
n_boots = 150

folder = paste("result_",snail_pop,"_method_",method,sep="")

matrix = as.matrix(read.table(paste(folder,"/matrix_time_",file,"_",snail_pop,".dat",sep = ""),h=F,stringsAsFactors = F))


m = matrix(0,nrow(matrix),4)
for(i in 1:nrow(matrix)){
  m[i,1] = quantile(matrix[i,],0.025,names=F)
  m[i,2] = quantile(matrix[i,],0.975,names=F)
  m[i,3] = quantile(matrix[i,],0.5,names=F)
  m[i,4] = mean(matrix[i,])
} 

df1 = data.frame(time_d= seq(1,nrow(m)), CI1=m[,1], CI2=m[,2], mediana = m[,3],mean = m[,4])


mp3=ggplot()+geom_ribbon(data=df1,aes(x=time_d/365,ymin=CI1,ymax=CI2),color = "lightgrey",fill = "lightgrey",alpha = 0.5,size = 0.5)+
  geom_line(data=df1,aes(x=time_d/365,y = mean),color = "black",size = 1.5)+
  scale_x_continuous(name = TeX('Time (years)'),limits = c(0,60))+
  scale_y_continuous(name = TeX("Number of infected individuals"))+
  #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  labs(tag = paste(tag_,sep=""))+
  theme_bw()+
  #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25,angle=90,hjust = 0.5),
          axis.title.y = element_text(colour="black", size = 25),
          axis.title.x = element_text(colour="black", size = 25),
          # plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 16),
          legend.position = "top",
          plot.tag.position = "topleft",
          plot.tag = element_text(size = 40),
          #axis.line = element_line(size=0.5, colour = "black")
  )
mp3


#mpt = ggarrange(mp,mp1+rremove("y.text")+rremove("y.title"),mp2+rremove("y.text")+rremove("y.title"),ncol=3,nrow=1,widths = c(1,0.8,0.8))
mpt = ggarrange(mp1,mp2+rremove("y.title"),mp3+rremove("y.title"),ncol=3,nrow=1,widths = c(0.9,0.8,0.8))

ggsave(paste("time_series_",snail_pop,"_",method,".png",sep=""),plot=mp,device="png",width = 25,height = 25,units = "cm")

ggsave(paste("time_series.png",sep=""),plot=mpt,device="png",width = 45,height = 25,units = "cm")






##################################################################################################
########################3 Plot time with intervention found
###################################################################################################



setwd("~/PosDoc/UNICAMP/Codes/CompleteCode/")

library(ggplot2)
library(latex2exp)
library(ggpubr)
method = 2
snail_pop = 500
tag_ = "A"
file = 1
n_pop_ga = 500
n_gen = 7
n_boots = 150
rounds = 10
interval = 1

folder = paste("result_",snail_pop,"_method_",method,"_",rounds,"_",interval,sep="")

matrix = as.matrix(read.table(paste(folder,"/matrix_time_found_",file,"_",snail_pop,".dat",sep = ""),h=F,stringsAsFactors = F))


m = matrix(0,nrow(matrix),4)
for(i in 1:nrow(matrix)){
  m[i,1] = quantile(matrix[i,],0.025,names=F)
  m[i,2] = quantile(matrix[i,],0.975,names=F)
  m[i,3] = quantile(matrix[i,],0.5,names=F)
  m[i,4] = mean(matrix[i,])
} 

df1 = data.frame(time_d= seq(1,nrow(m)), CI1=m[,1], CI2=m[,2], mediana = m[,3],mean = m[,4])


mp1=ggplot()+geom_ribbon(data=df1,aes(x=time_d/365,ymin=CI1,ymax=CI2),color = "lightgrey",fill = "lightgrey",alpha = 0.5,size = 0.5)+
  geom_line(data=df1,aes(x=time_d/365,y = mean),color = "black",size = 1.5)+
  scale_x_continuous(name = TeX('Time (years)'),limits = c(60,80))+
  scale_y_continuous(name = TeX("Number of infected individuals"))+
  #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  labs(tag = paste(tag_,sep=""))+
  theme_bw()+
  #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25,angle=90,hjust = 0.5),
          axis.title.y = element_text(colour="black", size = 25),
          axis.title.x = element_text(colour="black", size = 25),
          # plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 16),
          legend.position = "top",
          plot.tag.position = "topleft",
          plot.tag = element_text(size = 40),
          #axis.line = element_line(size=0.5, colour = "black")
  )
mp1


method = 2
snail_pop = 1000
tag_ = "B"
file = 1
n_pop_ga = 500
n_gen = 7
n_boots = 150

folder = paste("result_",snail_pop,"_method_",method,sep="")

matrix = as.matrix(read.table(paste(folder,"/matrix_time_",file,"_",snail_pop,".dat",sep = ""),h=F,stringsAsFactors = F))


m = matrix(0,nrow(matrix),4)
for(i in 1:nrow(matrix)){
  m[i,1] = quantile(matrix[i,],0.025,names=F)
  m[i,2] = quantile(matrix[i,],0.975,names=F)
  m[i,3] = quantile(matrix[i,],0.5,names=F)
  m[i,4] = mean(matrix[i,])
} 

df1 = data.frame(time_d= seq(1,nrow(m)), CI1=m[,1], CI2=m[,2], mediana = m[,3],mean = m[,4])


mp2=ggplot()+geom_ribbon(data=df1,aes(x=time_d/365,ymin=CI1,ymax=CI2),color = "lightgrey",fill = "lightgrey",alpha = 0.5,size = 0.5)+
  geom_line(data=df1,aes(x=time_d/365,y = mean),color = "black",size = 1.5)+
  scale_x_continuous(name = TeX('Time (years)'),limits = c(0,60))+
  scale_y_continuous(name = TeX("Number of infected individuals"))+
  #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  labs(tag = paste(tag_,sep=""))+
  theme_bw()+
  #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25,angle=90,hjust = 0.5),
          axis.title.y = element_text(colour="black", size = 25),
          axis.title.x = element_text(colour="black", size = 25),
          # plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 16),
          legend.position = "top",
          plot.tag.position = "topleft",
          plot.tag = element_text(size = 40),
          #axis.line = element_line(size=0.5, colour = "black")
  )
mp2


method = 2
snail_pop = 2000
tag_ = "C"
file = 1
n_pop_ga = 500
n_gen = 7
n_boots = 150

folder = paste("result_",snail_pop,"_method_",method,sep="")

matrix = as.matrix(read.table(paste(folder,"/matrix_time_",file,"_",snail_pop,".dat",sep = ""),h=F,stringsAsFactors = F))


m = matrix(0,nrow(matrix),4)
for(i in 1:nrow(matrix)){
  m[i,1] = quantile(matrix[i,],0.025,names=F)
  m[i,2] = quantile(matrix[i,],0.975,names=F)
  m[i,3] = quantile(matrix[i,],0.5,names=F)
  m[i,4] = mean(matrix[i,])
} 

df1 = data.frame(time_d= seq(1,nrow(m)), CI1=m[,1], CI2=m[,2], mediana = m[,3],mean = m[,4])


mp3=ggplot()+geom_ribbon(data=df1,aes(x=time_d/365,ymin=CI1,ymax=CI2),color = "lightgrey",fill = "lightgrey",alpha = 0.5,size = 0.5)+
  geom_line(data=df1,aes(x=time_d/365,y = mean),color = "black",size = 1.5)+
  scale_x_continuous(name = TeX('Time (years)'),limits = c(0,60))+
  scale_y_continuous(name = TeX("Number of infected individuals"))+
  #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  labs(tag = paste(tag_,sep=""))+
  theme_bw()+
  #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25,angle=90,hjust = 0.5),
          axis.title.y = element_text(colour="black", size = 25),
          axis.title.x = element_text(colour="black", size = 25),
          # plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 16),
          legend.position = "top",
          plot.tag.position = "topleft",
          plot.tag = element_text(size = 40),
          #axis.line = element_line(size=0.5, colour = "black")
  )
mp3


#mpt = ggarrange(mp,mp1+rremove("y.text")+rremove("y.title"),mp2+rremove("y.text")+rremove("y.title"),ncol=3,nrow=1,widths = c(1,0.8,0.8))
mpt = ggarrange(mp1,mp2+rremove("y.title"),mp3+rremove("y.title"),ncol=3,nrow=1,widths = c(0.9,0.8,0.8))

ggsave(paste("time_series_",snail_pop,"_",method,".png",sep=""),plot=mp,device="png",width = 25,height = 25,units = "cm")

ggsave(paste("time_series.png",sep=""),plot=mpt,device="png",width = 45,height = 25,units = "cm")

################################################################################################3
### time series several data ##########################################3

library(ggplot2)
library(latex2exp)
library(ggpubr)
library(viridis)
library(RColorBrewer)
setwd("/data/thomas/schisto_treat/")

method = 2
snail_pop = 500
file = 4
ly = c(0,250)

for(snail_pop in c(500,1000,2000)){
  for(file in 4:5){
  
    print(c(snail_pop,file))
r = "1"
ii = "0.0"


folder = paste("result_",snail_pop,"_method_",method,"_",r,"_",ii,sep="")

matrix = as.matrix(read.table(paste(folder,"/matrix_time_",file,".dat",sep = ""),h=F,stringsAsFactors = F))

m = matrix(0,nrow(matrix),4)
for(i in 1:nrow(matrix)){
  m[i,1] = quantile(matrix[i,],0.025,names=F)
  m[i,2] = quantile(matrix[i,],0.975,names=F)
  m[i,3] = quantile(matrix[i,],0.5,names=F)
  m[i,4] = mean(matrix[i,])
} 

df1 = data.frame(time_d= seq(1,nrow(m)), CI1=m[,1], CI2=m[,2], mediana = m[,3],mean = m[,4],round=rep(r,nrow(m)),interval=rep(ii,nrow(m)),idx=rep(6,nrow(m)))


intr = c("0.5","1.0","2.0")

rnd = c("2","4","6","8","10")

for(kk in 1:length(rnd)){
  for(ii in intr){
    r = rev(rnd)[kk]
    print(c(r,ii))
    folder = paste("result_",snail_pop,"_method_",method,"_",r,"_",ii,sep="")
    
    matrix = as.matrix(read.table(paste(folder,"/matrix_time_",file,".dat",sep = ""),h=F,stringsAsFactors = F))
    
    m = matrix(0,nrow(matrix),4)
    for(i in 1:nrow(matrix)){
      m[i,1] = quantile(matrix[i,],0.025,names=F)
      m[i,2] = quantile(matrix[i,],0.975,names=F)
      m[i,3] = quantile(matrix[i,],0.5,names=F)
      m[i,4] = mean(matrix[i,])
    } 
    
    df1 = rbind(df1,data.frame(time_d= seq(1,nrow(m)), CI1=m[,1], CI2=m[,2], mediana = m[,3],mean = m[,4],round=rep(r,nrow(m)),interval=rep(ii,nrow(m)),idx=rep(kk,nrow(m))))
  }
}


cv = viridis(6)
br1 = c(1,2,3,4,5,6)
br = c("1","2","4","6","8","10")

library(vctrs)


intf = "0.5"
df2 = df1[df1$interval %in% c("0.0",intf),]

tag_ = "A"
mp1=ggplot()+#geom_ribbon(data=df2,aes(x=time_d/365,ymin=CI1,ymax=CI2,color = round,fill = round),alpha = 0.5,size = 0.5)+
  geom_line(data=df2,aes(x=time_d/365,y = mean,color = as.factor(idx)),size = 1.5)+
  scale_x_continuous(name = TeX('Time (years)'),limits = c(55,90))+
  scale_y_continuous(name = TeX("Number of infected individuals"),limits = ly)+
  scale_color_manual(values = rev(cv),name = "Rounds",breaks = rev(br1),labels = br)+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  labs(tag = paste(tag_,sep=""))+
  theme_bw()+
  #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25,angle=90,hjust = 0.5),
          axis.title.y = element_text(colour="black", size = 25),
          axis.title.x = element_text(colour="black", size = 25),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 16),
          legend.position = "top",
          plot.tag.position = "topleft",
          plot.tag = element_text(size = 40),
          #axis.line = element_line(size=0.5, colour = "black")
  )
mp1
ggsave(paste("plots/time_series_",snail_pop,"_",method,"_",intf,"_",file,".png",sep=""),plot=mp1,device="png",width = 25,height = 25,units = "cm")


intf = "1.0"
df2 = df1[df1$interval %in% c("0.0",intf),]
tag_ = "B"
mp2=ggplot()+#geom_ribbon(data=df2,aes(x=time_d/365,ymin=CI1,ymax=CI2,color = round,fill = round),alpha = 0.5,size = 0.5)+
  geom_line(data=df2,aes(x=time_d/365,y = mean,color = as.factor(idx)),size = 1.5)+
  scale_x_continuous(name = TeX('Time (years)'),limits = c(55,90))+
  scale_y_continuous(name = TeX("Number of infected individuals"),limits = ly)+
  scale_color_manual(values = rev(cv),name = "Rounds",breaks = rev(br1),labels = br)+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  labs(tag = paste(tag_,sep=""))+
  theme_bw()+
  #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25,angle=90,hjust = 0.5),
          axis.title.y = element_text(colour="black", size = 25),
          axis.title.x = element_text(colour="black", size = 25),
          # plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 16),
          legend.position = "top",
          plot.tag.position = "topleft",
          plot.tag = element_text(size = 40),
          #axis.line = element_line(size=0.5, colour = "black")
  )
mp2
ggsave(paste("plots/time_series_",snail_pop,"_",method,"_",intf,"_",file,".png",sep=""),plot=mp2,device="png",width = 25,height = 25,units = "cm")


intf = "2.0"
df2 = df1[df1$interval %in% c("0.0",intf),]
tag_ = "C"
mp3=ggplot()+#geom_ribbon(data=df2,aes(x=time_d/365,ymin=CI1,ymax=CI2,color = round,fill = round),alpha = 0.5,size = 0.5)+
  geom_line(data=df2,aes(x=time_d/365,y = mean,color = as.factor(idx)),size = 1.5)+
  scale_x_continuous(name = TeX('Time (years)'),limits = c(55,90))+
  scale_y_continuous(name = TeX("Number of infected individuals"),limits = ly)+
  scale_color_manual(values = rev(cv),name = "Rounds",breaks = rev(br1),labels = br)+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  labs(tag = paste(tag_,sep=""))+
  theme_bw()+
  #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25,angle=90,hjust = 0.5),
          axis.title.y = element_text(colour="black", size = 25),
          axis.title.x = element_text(colour="black", size = 25),
          # plot.tag = element_text(size = 40),
          legend.text = element_text(size = 25),
          legend.title = element_text(size = 16),
          legend.position = "top",
          plot.tag.position = "topleft",
          plot.tag = element_text(size = 40),
          #axis.line = element_line(size=0.5, colour = "black")
  )
mp3

ggsave(paste("time_series_",snail_pop,"_",method,"_",intf,"_",file,".png",sep=""),plot=mp3,device="png",width = 25,height = 25,units = "cm")

legend <- cowplot::get_legend(mp1+guides(colour = guide_legend(override.aes = list(size=10)))+ theme(legend.position = "bottom",legend.text = element_text(size = 40),legend.title = element_text(size = 40)))

mpt = ggarrange(mp1+theme(legend.position = "none"),NULL,mp2+theme(legend.position = "none")+rremove("y.title")+rremove("y.text"),NULL,mp3+theme(legend.position = "none")+rremove("y.title")+rremove("y.text"),NULL,NULL,legend,NULL,NULL,ncol = 5,nrow = 2,widths = c(1.1,0.1,0.95,0.1,0.95,0.05,0.05,3,0.05,0.05),heights = c(1,0.2))
mpt


ggsave(paste("plots/time_series_treatment_",snail_pop,"_",method,"_",file,".png",sep=""),plot=mpt,device="png",width = 60,height = 25,units = "cm")

ggsave(paste("plots/time_series_treatment_",snail_pop,"_",method,"_",file,".eps",sep=""),plot=mpt,device="eps",width = 60,height = 25,units = "cm")
  }
}



#######################################################################################################################
############################## proportion of zeros######################################################################
#####################################################################################################################

setwd("~/PosDoc/UNICAMP/Codes/CompleteCode/")

library(ggplot2)
library(latex2exp)
library(ggpubr)
method = 2
snail_pop = 500

file = 2
n_pop_ga = 500
n_gen = 7
n_boots = 150


int = c("0.5","1.0","2.0")
rounds = c(2,4,6,8,10)


folder = paste("Cluster/treatment/result_",snail_pop,"_method_",method,"/",sep="")


data0 = read.table(paste(folder,"inf_time_series_r_0.dat",sep=""),h = F,stringsAsFactors = F)
n0 = length(data0[nrow(data0),data0[nrow(data0),]==0])

file = 2

m1 = c()
m2 = c()
m3 = c()


for(ii in 1:length(int)){
  i = int[ii]
  aux1 = c()
  folder = paste("Cluster/treatment/result_",snail_pop,"_method_",method,"_1_0.0/",sep="")
  data = read.table(paste(folder,"inf_time_series_r_",file,".dat",sep=""),h = F,stringsAsFactors = F)
  n = length(data[nrow(data),data[nrow(data),]==0])
  aux1[1] = (n-n0)/(1000-n0)
  for(rr in 1:length(rounds)){
    
    r = rounds[rr]
    print(c(i,r))
    folder = paste("Cluster/treatment/result_",snail_pop,"_method_",method,"_",r,"_",i,"/",sep="")
    data = read.table(paste(folder,"inf_time_series_r_",file,".dat",sep=""),h = F,stringsAsFactors = F)
    n = length(data[nrow(data),data[nrow(data),]==0])
    aux1[rr+1] = (n-n0)/(1000-n0)
  }
  m1 = c(m1,aux1)
  m2 = c(m2,c(1,rounds))
  m3 = c(m3,rep(i,length(aux1)))
}

df2 = data.frame(prop = m1,rounds = m2,intervalo = m3)
write.table(df2,"proportion_red_1.dat",col.names = F)
tag_ = "A"
mp1=ggplot()+#geom_ribbon(data=df2,aes(x=time_d/365,ymin=CI1,ymax=CI2,color = round,fill = round),alpha = 0.5,size = 0.5)+
  geom_line(data=df2,aes(x=as.factor(rounds),y = prop,color = as.factor(intervalo),group = as.factor(intervalo)),size = 3)+
  geom_point(data=df2,aes(x=as.factor(rounds),y = prop,color = as.factor(intervalo)),size = 5)+
  scale_x_discrete(name = TeX('Number of treatment rounds'),breaks = c(1,rounds),labels = c(1,rounds))+
  scale_y_continuous(name = TeX("Proportion of extinction"),limits = c(0,1))+
  scale_color_manual(values = c("navy","red","purple"),name = "Interval of treatment",breaks = int,labels = c("half year","one year","two years"))+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  labs(tag = paste("A",sep=""))+
  theme_bw()+
  #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25,angle=90,hjust = 0.5),
          axis.title.y = element_text(colour="black", size = 40),
          axis.title.x = element_text(colour="black", size = 25),
          # plot.tag = element_text(size = 40),
          legend.text = element_text(size = 30),
          legend.title = element_text(size = 30),
          legend.box.background = element_rect(color="black", size=1.5),
          legend.position = c(0.2, 0.85),
          plot.tag.position = "topleft",
          plot.tag = element_text(size = 40),
          #axis.line = element_line(size=0.5, colour = "black")
  )
mp1

###
file = 4
m1 = c()
m2 = c()
m3 = c()


for(ii in 1:length(int)){
  i = int[ii]
  aux1 = c()
  folder = paste("Cluster/treatment/result_",snail_pop,"_method_",method,"_1_0.0/",sep="")
  data = read.table(paste(folder,"inf_time_series_r_",file,".dat",sep=""),h = F,stringsAsFactors = F)
  n = length(data[nrow(data),data[nrow(data),]==0])
  aux1[1] = (n-n0)/(1000-n0)
  for(rr in 1:length(rounds)){
    
    r = rounds[rr]
    print(c(i,r))
    folder = paste("Cluster/treatment/result_",snail_pop,"_method_",method,"_",r,"_",i,"/",sep="")
    data = read.table(paste(folder,"inf_time_series_r_",file,".dat",sep=""),h = F,stringsAsFactors = F)
    n = length(data[nrow(data),data[nrow(data),]==0])
    aux1[rr+1] = (n-n0)/(1000-n0)
  }
  m1 = c(m1,aux1)
  m2 = c(m2,c(1,rounds))
  m3 = c(m3,rep(i,length(aux1)))
}

df3 = data.frame(prop = m1,rounds = m2,intervalo = m3)
write.table(df3,"proportion_red_2.dat",col.names = F)
tag_ = "B"
mp2=ggplot()+#geom_ribbon(data=df2,aes(x=time_d/365,ymin=CI1,ymax=CI2,color = round,fill = round),alpha = 0.5,size = 0.5)+
  geom_line(data=df3,aes(x=as.factor(rounds),y = prop,color = as.factor(intervalo),group = as.factor(intervalo)),size = 3)+
  geom_point(data=df3,aes(x=as.factor(rounds),y = prop,color = as.factor(intervalo)),size = 5)+
  scale_x_discrete(name = TeX('Number of treatment rounds'),breaks = c(1,rounds),labels = c(1,rounds))+
  scale_y_continuous(name = TeX("Proportion of extinction"),limits = c(0,1))+
  scale_color_manual(values = c("navy","red","purple"),name = "Interval of treatment",breaks = int,labels = c("half year","one year","two years"))+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  labs(tag = paste("B",sep=""))+
  theme_bw()+
  #scale_color_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  #scale_fill_manual(values = c("red","blue"),labels = c("HCW","Residents"),name="Individual")+
  theme(  panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.minor = element_blank(),
          #panel.border = element_blank(),
          panel.background = element_blank(),
          #plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
          #text=element_text(family = "Tahoma"),
          axis.title = element_text(face="bold"),
          axis.text.x = element_text(colour="black", size = 25),
          axis.text.y = element_text(colour="black", size = 25,angle=90,hjust = 0.5),
          axis.title.y = element_text(colour="black", size = 40),
          axis.title.x = element_text(colour="black", size = 25),
          # plot.tag = element_text(size = 40),
          legend.text = element_text(size = 30),
          legend.title = element_text(size = 30),
          legend.box.background = element_rect(color="black", size=1.5),
          legend.position = c(0.2, 0.85),
          plot.tag.position = "topleft",
          plot.tag = element_text(size = 40),
          #axis.line = element_line(size=0.5, colour = "black")
  )
mp2

mpt = ggarrange(mp1+theme(legend.position = c(0.25,0.9), legend.text = element_text(size = 20),legend.title = element_text(size = 20)),NULL,mp2+rremove("y.title")+rremove("y.text")+theme(legend.position = c(0.25,0.9), legend.text = element_text(size = 20),legend.title = element_text(size = 20)),ncol = 3,nrow = 1,widths = c(1.1,0.1,1.0))
mpt


ggsave(paste("Cluster/treatment/plots/probability",snail_pop,".eps",sep=""),plot=mpt,device="eps",width = 35,height = 25,units = "cm")

ggsave(paste("Cluster/treatment/plots/probability",snail_pop,".png",sep=""),plot=mpt,device="png",width = 35,height = 25,units = "cm")

