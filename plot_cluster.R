
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
    
    ggsave(paste("plots/time_series_",snail_pop,"_",method,"_",intf,"_",file,".png",sep=""),plot=mp3,device="png",width = 25,height = 25,units = "cm")
    
    legend <- cowplot::get_legend(mp1+guides(colour = guide_legend(override.aes = list(size=10)))+ theme(legend.position = "bottom",legend.text = element_text(size = 40),legend.title = element_text(size = 40)))
    
    mpt = ggarrange(mp1+theme(legend.position = "none"),NULL,mp2+theme(legend.position = "none")+rremove("y.title")+rremove("y.text"),NULL,mp3+theme(legend.position = "none")+rremove("y.title")+rremove("y.text"),NULL,NULL,legend,NULL,NULL,ncol = 5,nrow = 2,widths = c(1.1,0.1,0.95,0.1,0.95,0.05,0.05,3,0.05,0.05),heights = c(1,0.2))
    mpt
    
    
    ggsave(paste("plots/time_series_treatment_",snail_pop,"_",method,"_",file,".png",sep=""),plot=mpt,device="png",width = 60,height = 25,units = "cm")
    
    ggsave(paste("plots/time_series_treatment_",snail_pop,"_",method,"_",file,".eps",sep=""),plot=mpt,device="eps",width = 60,height = 25,units = "cm")
  }
}



# histograms --------------------------------------------------------------


xx = read.table("result_500_method_2/inf_time_series_r_0.dat", h=F)
x1 = rowMeans(xx)


xx = read.table("result_500_method_2/inf_time_series_r_2.dat", h=F)
x6 = rowMeans(xx)

xx = read.table("old_data/result_500_method_2_8_0.5/inf_time_series_r_2.dat", h=F)
x2 = rowMeans(xx)

xx = read.table("old_data/result_500_method_2_10_0.5/inf_time_series_r_2.dat", h=F)
x3 = rowMeans(xx)


xx = read.table("result_500_method_2_8_0.5/inf_time_series_r_6.dat", h=F)
x4 = rowMeans(xx)

xx = read.table("result_500_method_2_10_0.5/inf_time_series_r_6.dat", h=F)
x5 = rowMeans(xx)


# 
# x3 = x3 %>% group_by(ind) %>% mutate(soma = sum(values)) %>%
#   filter(soma > 0)


plot(x1, xlim = c(20000, 30000), type = "l")
lines(x2, col = "red")
lines(x3, col = "blue")
lines(x4, col = "green")
lines(x5, col = "purple")
lines(x6, col = "orange")

  ggplot()+
  geom_histogram(data = x1, aes(x = values), color = "blue",fill = "blue", alpha = 0.2)+
    geom_histogram(data = x2, aes(x = values), color = "red",fill ="red",  alpha = 0.2)+
    geom_histogram(data = x3, aes(x = values), color = "green",fill ="green",  alpha = 0.2)+
  scale_y_continuous(trans = "log10")


  
  x1 %>% group_by(ind) %>% summarise(soma = sum(values)) %>%
    filter(soma == 0) %>% nrow
  
  
  x2 %>% group_by(ind) %>% summarise(soma = sum(values)) %>%
    filter(soma == 0) %>% nrow
  
  x3 %>% group_by(ind) %>% summarise(soma = sum(values)) %>%
    filter(soma == 0) %>% nrow
  
  
  
  x1 %>% 
    filter(values>0) %>% nrow/1000
  
  
  x2 %>% 
    filter(values>0) %>% nrow/1000
  
  
  x3 %>% 
    filter(values>0) %>% nrow/1000

  

# testing -----------------------------------------------------------------
  
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
  
  
  folder = paste("result_",snail_pop,"_method_",method,"/",sep="")
  
  
  data0 = read.table(paste(folder,"inf_time_series_r_0.dat",sep=""),h = F,stringsAsFactors = F)
  n0 = length(data0[nrow(data0),data0[nrow(data0),]==0])
  
  file = 2
  
  m1 = c()
  m2 = c()
  m3 = c()
  
  
  for(ii in 1:length(int)){
    i = int[ii]
    aux1 = c()
    folder = paste("old_data/result_",snail_pop,"_method_",method,"_1_0.0/",sep="")
    data = read.table(paste(folder,"inf_time_series_r_",file,".dat",sep=""),h = F,stringsAsFactors = F)
    n = length(data[nrow(data),data[nrow(data),]==0])
    aux1[1] = (n-n0)/(1000-n0)
    for(rr in 1:length(rounds)){
      
      r = rounds[rr]
      print(c(i,r))
      folder = paste("old_data/result_",snail_pop,"_method_",method,"_",r,"_",i,"/",sep="")
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
    folder = paste("old_data/result_",snail_pop,"_method_",method,"_1_0.0/",sep="")
    data = read.table(paste(folder,"inf_time_series_r_",file,".dat",sep=""),h = F,stringsAsFactors = F)
    n = length(data[nrow(data),data[nrow(data),]==0])
    aux1[1] = (n-n0)/(1000-n0)
    for(rr in 1:length(rounds)){
      
      r = rounds[rr]
      print(c(i,r))
      folder = paste("old_data/result_",snail_pop,"_method_",method,"_",r,"_",i,"/",sep="")
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
  