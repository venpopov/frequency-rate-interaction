rm(list=ls())

library(tidyverse)
library(here)

setwd(here())

set.seed(1404)


all_words<-read.csv('expfiles/materials/exp4_freq_precision_Sept2018/SUBTLEX_all_words.csv', stringsAsFactors = FALSE)

  
#filter out words
word_pool<-filter(all_words,
                Length <= 7,
                Length >= 5,
                Dom_PoS_SUBTLEX == 'Noun',
                substring(Word,nchar(Word)) != 's', #eliminate most plural words
                FREQcount >= 50, #eliminate very uncommon words and some non-words
                !(Word %in% c('asshole', 'bastard', 'doesn', 'fella', 'grandpa', 'honour', 'madame', 'nigger', 'outta', 'pussy', 'suicide','weren', 'whore','wouldn','armour', 'bijou', 'blowjob', 'bonsoir', 'brothel', 'commie', 'creme', 'dildo', 'gangsta', 'gigolo', 'grammy', 'heinie', 'hombre', 'hubba','hubby', 'imaging','inbound','mahjong','mardi','massa','mugging','nonny', 'parlour', 'perdy', 'potty', 'prima', 'quickie', 'rummy','rumba', 'seder', 'shite', 'sicko', 'sirree', 'swede', 'senora', 'bitch','kraut','pissant','auntie','betcha','boner','booby', 'buick','buffy','bugger', 'carte','centre','chico','chink','colour','condom','couldn','doggie','faggot','fatso','favour','flavour','fucker','gideon','shylock','givin','girlie','gringo','hickey','homie','homey','homeboy','hooch','hooker','hottie','humour','kleenex','latino','lesbian','lordy','mammy','merci','midget','milord','negro','nigga','orgasm','pardner','pecker','peewee','prick','rapist','rumour','sayid','semen','signor','signora','sissy','snort','telly','thong','titty','vagina','wanker','weenie','wiener','women','zorro','sperm','squaw','mommy','daddy','pickin','mitzvah','anytime','tootsie','renoir','enema','tellin','theatre','femme','gonzo', 'leone','messin','beaut','matey','fuhrer','laddie','cajun','kaiser','norse','altair','mongol','saviour','english','baldy','chemo','banzai','chewie','cabbie'))
                # manually eliminate obscenities, cut off contractions, racist words, words with multiple spellings, slang, non-nouns and proper nouns that made it through the filter, and foreign words
                )

HF_pool <- word_pool[word_pool$Lg10WF >= 2.7,] #HF words have Lg10WF greater than 2.68 
LF_pool <- word_pool[word_pool$Lg10WF <= 2.03,] #LF words have LF10WF less than 2.03
#weird numbers are to allow sufficient pool size

#randomly select words
HF_words <- HF_pool[sample(1:nrow(HF_pool), 750),]
LF_words <- LF_pool[sample(1:nrow(LF_pool), 750),] 

#randomly select practice words from remaining HF words
remaining_HF <- HF_pool[!(HF_pool$Word %in% HF_words$Word),]

#practice_words <- remaining_HF[sample(1:nrow(remaining_HF),32),1]


#prepare word lists for final output
HF_words_final <- data.frame(cbind(HF_words[,1],rep('HF',750)))
LF_words_final <- data.frame(cbind(LF_words[,1],rep('LF',750)))

stim <- rbind(HF_words_final,LF_words_final)
colnames(stim) <- c('word','freq')

#practice_words_final <- cbind(practice_words,rep('HF',32))
#colnames(practice_words_final) <- c('word','freq')

#write word list files
write.csv(stim,'expfiles/materials/exp4_freq_precision_Sept2018/words.csv', row.names = FALSE)
#write.csv(practice_words_final, 'practice_words.csv', row.names = FALSE)