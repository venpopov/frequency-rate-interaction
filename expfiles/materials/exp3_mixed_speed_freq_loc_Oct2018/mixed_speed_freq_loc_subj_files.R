#############
#Setup
#############
rm(list=ls())
library(here)
setwd(here())
set.seed(75952)
screen_x<-1680
screen_y<-1050
number_of_subj_files<-50
all_words <- read.csv(file = "expfiles/materials/exp6_mixed_speed_freq_loc_Oct2018/words.csv",stringsAsFactors = FALSE)
config_to_order<- list(c(1,1,1,-1,-1),
                       c(1,1,-1,1,-1),
                       c(1,1,-1,-1,1),
                       c(1,-1,1,1,-1),
                       c(1,-1,1,-1,1),
                       c(1,-1,-1,1,1),
                       c(-1,1,1,1,-1),
                       c(-1,1,1,-1,1),
                       c(-1,1,-1,1,1),
                       c(-1,-1,1,1,1))

for (n in (1:number_of_subj_files))
{
all_trials<-data.frame(matrix(ncol=32,nrow=300),stringsAsFactors = FALSE)
colnames(all_trials)<-c('word1','word2','word3','word4','word5','angle1','angle2','angle3','angle4','angle5','target_word_place', 'target_word','target_angle','x1','y1','x2','y2','x3','y3','x4','y4','x5','y5','dom_freq','duration','Procedure','freq_config','freq1','freq2','freq3','freq4','freq5')
all_words <- all_words[sample(nrow(all_words)),]
HF_pool<-all_words[all_words$freq=='HF',]
LF_pool<-all_words[all_words$freq=='LF',]


all_trials[,'dom_freq']<-rep(c(rep('HF',15),rep('LF',15)),10)
all_trials[,'duration']<-rep(c(rep(500,5),rep(750,5),rep(1000,5)),20)
all_trials[,'target_word_place']<-rep(1:5,60)
all_trials[,'Procedure']<-rep('TestProc',300)
all_trials[,'freq_config']<-c(sample(0:9),sample(0:9),sample(0:9))

for (i in (1:300))
{
  bad <- TRUE
  while(bad)
  {
    new_angle_row<-sort(sample(360,5))
    diff_1 <- min(new_angle_row[2] - new_angle_row[1])
    diff_2 <- min(new_angle_row[3] - new_angle_row[2])
    diff_3 <- min(new_angle_row[4] - new_angle_row[3])
    diff_4 <- min(new_angle_row[5] - new_angle_row[4])
    diff_5 <- min(360-new_angle_row[5]+new_angle_row[1]) #this is to account for wraparound effects
    if (diff_1 > 60 & diff_2 > 60 & diff_3 > 60 & diff_4 > 60 & diff_5 > 60)
    {
      bad<-FALSE
    }
  }
  all_trials[i,c('angle1','angle2','angle3','angle4','angle5')]<-sample(new_angle_row)
  all_trials[i,'target_angle']<-all_trials[i,paste0('angle',all_trials[i,'target_word_place'])]
  
}



for (j in (1:300))
{
  word1_coords<-c((screen_x/2)+400*sin(all_trials[j,'angle1']*pi/180),(screen_y/2)-400*cos(all_trials[j,'angle1']*pi/180))
  word2_coords<-c((screen_x/2)+400*sin(all_trials[j,'angle2']*pi/180),(screen_y/2)-400*cos(all_trials[j,'angle2']*pi/180))
  word3_coords<-c((screen_x/2)+400*sin(all_trials[j,'angle3']*pi/180),(screen_y/2)-400*cos(all_trials[j,'angle3']*pi/180))
  word4_coords<-c((screen_x/2)+400*sin(all_trials[j,'angle4']*pi/180),(screen_y/2)-400*cos(all_trials[j,'angle4']*pi/180))
  word5_coords<-c((screen_x/2)+400*sin(all_trials[j,'angle5']*pi/180),(screen_y/2)-400*cos(all_trials[j,'angle5']*pi/180))
  all_trials[j,c('x1','y1')]<-word1_coords
  all_trials[j,c('x2','y2')]<-word2_coords
  all_trials[j,c('x3','y3')]<-word3_coords
  all_trials[j,c('x4','y4')]<-word4_coords
  all_trials[j,c('x5','y5')]<-word5_coords
}


for (k in (1:9))
{
  all_trials[(k*30)+(1:30),'freq_config']<-(all_trials[1:30,'freq_config']+k)%%10
}

HF_index<-1
LF_index<-1
freq_order<-c()
for (l in (1:300))
{
 freq_order<-unlist(config_to_order[all_trials[l,'freq_config']+1]) 
 if (all_trials[l,'dom_freq'] == 'LF')
 {
   freq_order <- (-1*freq_order)
 }
 
 for (u in (1:5))
 {
   if (freq_order[u]>0)
   {
     all_trials[l,paste0('word',toString(u))]<-HF_pool[HF_index,1]
     all_trials[l,paste0('freq',toString(u))]<-HF_pool[HF_index,2]
     HF_index<-HF_index+1
   }
   else
   {
     all_trials[l,paste0('word',toString(u))]<-LF_pool[LF_index,1]
     all_trials[l,paste0('freq',toString(u))]<-LF_pool[LF_index,2]
     LF_index<-LF_index+1
   }
  }
 all_trials[l,'target_word']<-all_trials[l,paste0('word',all_trials[l,'target_word_place'])]
  
 }
  
  
data<-data.frame(matrix(ncol=32,nrow=310),stringsAsFactors = FALSE)
colnames(data)<-c('word1','word2','word3','word4','word5','angle1','angle2','angle3','angle4','angle5','target_word_place', 'target_word','target_angle','x1','y1','x2','y2','x3','y3','x4','y4','x5','y5','dom_freq','duration','Procedure','freq_config','freq1','freq2','freq3','freq4','freq5')  
for (v in 0:9)
{
  new_block <- all_trials[30*v+(1:30),]
  data[31*v+(1:30),]<-new_block[sample(nrow(new_block)),]
  data[31*(v+1),] <- rep(NA,32)
  data[31*(v+1),'Procedure'] <- c('BreakProc') 
}
  
  
w_n<-data.frame(matrix(ncol=2,nrow=310))
w_n[,1]<-rep(1,310)
w_n[,2]<-rep('',310)
colnames(w_n)<-c('Weight','Nested')

data<-data[,c('Procedure','dom_freq','word1','word2','word3','word4','word5','angle1','angle2','angle3','angle4','angle5','duration','target_word_place', 'target_word','target_angle','x1','y1','x2','y2','x3','y3','x4','y4','x5','y5','freq1','freq2','freq3','freq4','freq5','freq_config')]
data<-cbind(w_n,data)

write.table(data, paste0("expfiles/materials/exp6_mixed_speed_freq_loc_Oct2018/Subject_Files/mixed_speed_freq_loc_subject_",toString(n),".txt"), sep="\t", quote= FALSE, row.names = FALSE)


}