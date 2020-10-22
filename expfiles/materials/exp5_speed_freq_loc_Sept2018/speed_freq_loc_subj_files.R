#############
#Setup
#############
rm(list=ls())
library(here)
setwd(here())
set.seed(18952)
screen_x<-1680
screen_y<-1050
number_of_subj_files<-50
all_words <- read.csv(file = "expfiles/materials/exp5_speed_freq_loc_Sept2018/words.csv",stringsAsFactors = FALSE)


for (n in (1:number_of_subj_files))
{
  all_trials<-data.frame(matrix(ncol=26,nrow=310),stringsAsFactors = FALSE)
  colnames(all_trials)<-c('word1','word2','word3','word4','word5','angle1','angle2','angle3','angle4','angle5','target_word_place', 'target_word','target_angle','x1','y1','x2','y2','x3','y3','x4','y4','x5','y5','freq','duration','Procedure')
  all_words <- all_words[sample(nrow(all_words)),]
  HF_pool<-all_words[all_words$freq=='HF',]
  LF_pool<-all_words[all_words$freq=='LF',]
  for (m in (1:10))
  {
    offset <- (m-1)*75
    #################
    #Create HF trials
    #################
    HF_trials<-data.frame(matrix(nrow=15,ncol=25),stringsAsFactors = FALSE)
    colnames(HF_trials)<-c('word1','word2','word3','word4','word5','angle1','angle2','angle3','angle4','angle5','target_word_place', 'target_word','target_angle','x1','y1','x2','y2','x3','y3','x4','y4','x5','y5','freq','duration')
    HF_trials[,1] <- HF_pool[offset+(1:15),1]
    HF_trials[,2] <- HF_pool[offset+(16:30),1]
    HF_trials[,3] <- HF_pool[offset+(31:45),1]
    HF_trials[,4] <- HF_pool[offset+(46:60),1]
    HF_trials[,5] <- HF_pool[offset+(61:75),1]
    
    ########Assign random angles, make sure no two are closer than 60 degrees
    for (j in (1:15))
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
      HF_trials[j,6:10]<-sample(new_angle_row)
      
    }
    
    ###### Generate order for target word position
    HF_word_position_order_duration<-data.frame(matrix(nrow=15,ncol=2))
    
    HF_word_position_order_duration[,1]<-c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3))
    HF_word_position_order_duration[,2]<-rep(c(500,750,1000),5)
    HF_word_position_order_duration<-HF_word_position_order_duration[sample(nrow(HF_word_position_order_duration)),]
    
    HF_trials[,c(11,25)]<-HF_word_position_order_duration
    
    for (i in (1:15))
    {
      HF_trials[i,12] <- HF_trials[i,HF_word_position_order_duration[i,1]]
      HF_trials[i,13] <- HF_trials[i,HF_word_position_order_duration[i,1]+5]
    }
    
    
    ##### Generate coordinates for each word based on the angle, adjust for different screen resolutions
    
    for (k in (1:15))
    {
      word1_coords<-c((screen_x/2)+400*sin(HF_trials[k,'angle1']*pi/180),(screen_y/2)-400*cos(HF_trials[k,'angle1']*pi/180))
      word2_coords<-c((screen_x/2)+400*sin(HF_trials[k,'angle2']*pi/180),(screen_y/2)-400*cos(HF_trials[k,'angle2']*pi/180))
      word3_coords<-c((screen_x/2)+400*sin(HF_trials[k,'angle3']*pi/180),(screen_y/2)-400*cos(HF_trials[k,'angle3']*pi/180))
      word4_coords<-c((screen_x/2)+400*sin(HF_trials[k,'angle4']*pi/180),(screen_y/2)-400*cos(HF_trials[k,'angle4']*pi/180))
      word5_coords<-c((screen_x/2)+400*sin(HF_trials[k,'angle5']*pi/180),(screen_y/2)-400*cos(HF_trials[k,'angle5']*pi/180))
      HF_trials[k,14:23]<-c(word1_coords,word2_coords,word3_coords,word4_coords,word5_coords)
    }
    
    HF_trials[,24]<-rep('HF',15)
    
    
    #################
    #Create LF trials
    #################
    LF_trials<-data.frame(matrix(nrow=15,ncol=25),stringsAsFactors = FALSE)
    colnames(LF_trials)<-c('word1','word2','word3','word4','word5','angle1','angle2','angle3','angle4','angle5','target_word_place', 'target_word','target_angle','x1','y1','x2','y2','x3','y3','x4','y4','x5','y5','freq','duration')
    LF_trials[,1] <- LF_pool[offset+(1:15),1]
    LF_trials[,2] <- LF_pool[offset+(16:30),1]
    LF_trials[,3] <- LF_pool[offset+(31:45),1]
    LF_trials[,4] <- LF_pool[offset+(46:60),1]
    LF_trials[,5] <- LF_pool[offset+(61:75),1]
    
    ########Assign random angles, make sure no two are closer than 60 degrees
    for (j in (1:15))
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
      LF_trials[j,6:10]<-sample(new_angle_row)
      
    }
    
    ###### Generate order for target word position
    LF_word_position_order_duration<-data.frame(matrix(nrow=15,ncol=2))
    
    LF_word_position_order_duration[,1]<-c(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3))
    LF_word_position_order_duration[,2]<-rep(c(500,750,1000),5)
    LF_word_position_order_duration<-LF_word_position_order_duration[sample(nrow(LF_word_position_order_duration)),]
    
    LF_trials[,c(11,25)]<-LF_word_position_order_duration
    
    for (i in (1:15))
    {
      LF_trials[i,12] <- LF_trials[i,LF_word_position_order_duration[i,1]]
      LF_trials[i,13] <- LF_trials[i,LF_word_position_order_duration[i,1]+5]
    }
    
    
    ##### Generate coordinates for each word based on the angle, adjust for different screen resolutions
    
    for (k in (1:15))
    {
      word1_coords<-c((screen_x/2)+400*sin(LF_trials[k,'angle1']*pi/180),(screen_y/2)-400*cos(LF_trials[k,'angle1']*pi/180))
      word2_coords<-c((screen_x/2)+400*sin(LF_trials[k,'angle2']*pi/180),(screen_y/2)-400*cos(LF_trials[k,'angle2']*pi/180))
      word3_coords<-c((screen_x/2)+400*sin(LF_trials[k,'angle3']*pi/180),(screen_y/2)-400*cos(LF_trials[k,'angle3']*pi/180))
      word4_coords<-c((screen_x/2)+400*sin(LF_trials[k,'angle4']*pi/180),(screen_y/2)-400*cos(LF_trials[k,'angle4']*pi/180))
      word5_coords<-c((screen_x/2)+400*sin(LF_trials[k,'angle5']*pi/180),(screen_y/2)-400*cos(LF_trials[k,'angle5']*pi/180))
      LF_trials[k,14:23]<-c(word1_coords,word2_coords,word3_coords,word4_coords,word5_coords)
    }
    
    LF_trials[,24]<-rep('LF',15)
    
    
    #####################
    #Construct block
    #####################
    full_block<-data.frame(matrix(nrow=31,ncol=26),stringsAsFactors = FALSE)
    colnames(full_block)[26]<-c('Procedure')

    all_HF_trials<-HF_trials[sample(nrow(HF_trials)),]
    all_LF_trials<-LF_trials[sample(nrow(LF_trials)),]
    
    
    #####generate frequency order, make sure no groups of 4+ consecutive trials with the same frequency
    
    ordering<-c(rep(1,15),rep(-1,15))
    is_list_bad <- TRUE
    while(is_list_bad)
    {
      attempt<-sample(ordering)
      indic_list<-c()
      for (k in 1:27)
      {
        
        consec_indic<-abs(attempt[k]+
                            attempt[k+1]+
                            attempt[k+2]+
                            attempt[k+3])
        indic_list<-c(indic_list,consec_indic)
      }
      is_list_bad <- any(indic_list>3)
    }
    
    freq_order <- attempt
    ############ Assemble block by pulling trials based on frequency order
    HF_index<-1
    LF_index<-1
    for (l in 1:length(freq_order))
    {
      if (freq_order[l]>0)
      {
        full_block[l,1:25]<-all_HF_trials[HF_index,]
        HF_index<-HF_index+1
      }
      else
      {
        full_block[l,1:25]<-all_LF_trials[LF_index,]
        LF_index<-LF_index+1
      }
    }
    
    #add break row
    full_block[1:30,26] <- rep('TestProc',30)
    full_block[31,]<-c(rep(NA,25),'BreakProc')
    
    
    offset2<- (m-1)*31
    all_trials[(offset2)+(1:31),]<-full_block
    
  }
  #generate 'weight' and 'nested' columns
  w_n<-data.frame(matrix(ncol=2,nrow=310))
  w_n[,1]<-rep(1,310)
  w_n[,2]<-rep('',310)
  colnames(w_n)<-c('Weight','Nested')
  #rearrange columns
  all_trials <- all_trials[c(26,24,1:10,25,11:23)]
  data<-cbind(w_n,all_trials)
  
  
  write.table(data, paste0("expfiles/materials/exp5_speed_freq_loc_Sept2018/Subject_Files/speed_freq_loc_subject_",toString(n),".txt"), sep="\t", quote= FALSE, row.names = FALSE)
  
}




