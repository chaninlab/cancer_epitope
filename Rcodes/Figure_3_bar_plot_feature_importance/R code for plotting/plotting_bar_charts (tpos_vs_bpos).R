#setwd("d:/R_test/plotting bar charts")
#setwd("~/Downloads")
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

#COLOR codes for epitope classes
#dual #ff66CC (purple)
#mono #66CCFF (cyan)
#tpos "blue"
#bpos "gree"

library(readxl)
library(ggplot2)
#library(grid)
library(gridExtra)
library(ggthemes)


df<-read_excel("tpos_vs_bpos_desc.xlsx")
df<-na.omit(df)


subplot_num=20

plots<-list()

for(i in 1:subplot_num)
  {
   df_row<-df[i,]
   formatted_df<-data.frame(matrix(ncol=4,nrow=2))
   x<-c("Descriptors","Epitope_type","Mean","Stdev")
   colnames(formatted_df)<-x
   
   #print(df_row)
   #print("\n")
   #saves the name of the descriptors for the 2 cells of first row
   formatted_df[,1]<-(df_row)[1]    #df_row$Descriptors
   #saves the epitope types (saved as column names)
   formatted_df[1,2]<-colnames(df_row)[2]
   formatted_df[2,2]<-colnames(df_row)[4]
   #saves the mean values of the epitope types
   formatted_df[1,3]<-(df_row)[2]
   formatted_df[2,3]<-(df_row)[4]
   #saves standard dev 
   formatted_df[1,4]<-(df_row)[3]
   formatted_df[2,4]<-(df_row)[5]
   
   
   #print(formatted_df)
   #print(" ")
   #line<-readline()
   #reordering the epitope types to the order as their row index so that ggplots plots
   #the bars the way they appear in the df
   formatted_df$Epitope_type<-as.factor(formatted_df$Epitope_type)
   formatted_df$Epitope_type<-factor(formatted_df$Epitope_type, levels=formatted_df$Epitope_type[order(as.numeric(rownames(df)))])
  
   #this shows the variable classes of the dataframe
   #str(formatted_df)
   
   dodge<-position_dodge(width=0.7)
   limits<-aes(ymax=formatted_df$Mean + (formatted_df$Mean>0)*formatted_df$Stdev, ymin=formatted_df$Mean - (formatted_df$Mean<0)*formatted_df$Stdev)
    
   #p<-ggplot(data=formatted_df,aes(x=Descriptors,y=Mean, fill=Epitope_type))
   p<-ggplot(data=formatted_df,aes(x=Descriptors,y=Mean, fill=Epitope_type))
                                                                                #adjust this for size of frame around bar
   p<-p + geom_bar(stat="identity", position = dodge, width = 0.7,colour="black",size=1)
   
   p<-p + scale_fill_manual("legend",values=c("blue","green"))
   #dual #ff66CC (purple)
   #mono #66CCFF (cyan)
   #tpos "blue"
   #bpos "green"
   
   p<-p + geom_errorbar(limits, position = dodge, width =0.3,size=1)
   
   
   p <- p + labs(x = " ", y = " ")
   p <- p + theme_few()
   p <- p + theme(legend.position="none")
   
   #adjusts the text size
   p <- p + theme(text=element_text(face="bold",size=20))
   
   #ignore this one for now
   #p <- p + theme(axis.line = element_line(size= 0.5, colour="black"))
   
   #adjusts the text distances from the axis, adjust the first numerical parameter
   p <- p + theme(axis.text.x = element_text(margin=margin(5,5,10,5,"pt")),
            axis.text.y = element_text(margin=margin(5,5,10,5,"pt")))
   
   #adjusts the frame around the picture
   p <- p + theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5))
   
   #adjest  of the tick marks
   p <- p + theme(axis.ticks.length = unit(.2, "cm"), axis.ticks = element_line(size = 1, colour="black")) 
   
  # print(p)
   
   plots[[i]]<-p

    
#this part saves the pictures to file   
   filename=paste0("tpos_vs_bpos", i, ".pdf")
   pdf(filename)
   print(p)
   dev.off()
   
  }


#do.call(grid.arrange(arrangeGrob(c(plots, ncol=5))))
#grid.arrange(arrangeGrobs(plots))


#this part attempts to save multiple plots together. Note that grid.arrange can be used with the whole plot vector instead of individual elements of it.
#the plotting together works up to 8 subplots. More will result in errornous display that are not corrigible with scaling the picture and is
#the error is retained inside the saved image file. Irrelevant of the file format (e.g. png, pdf)

#grid.arrange(plots[[1]],plots[[2]],plots[[3]],nrow=1,ncol=3) 
#pdf("plot.pdf")
#do.call("grid.arrange",c(plots,nrow=1))#,nrow=1,ncol=length(plots)))
#dev.off()