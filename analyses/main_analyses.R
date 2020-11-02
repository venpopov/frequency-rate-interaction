#############################################################################
# GOAL: analyze data from the three continuous precision freq experiments, July 2018
# AUTHOR: Ven Popov
# DATE UPDATED: 22 Oct 2020
#############################################################################

rm(list=ls())
library(tidyverse)
library(stats4)
library(circular)
library(lme4)
theme_set(theme_classic(base_size=12))
source('analyses/main_mixture_functions.R')

#############################################################################!
# DATA                                                                   ####
#############################################################################!

# load and clean up data
# source('analyses/exp5_clean_data.R')
exp1 <- read.csv('data/exp1_raw.csv')
exp2 <- read.csv('data/exp2_raw.csv')
exp3 <- read.csv('data/exp3_raw.csv')

# add exp number and set the duration for exp1, set pure vs mixed
exp1$exp <- 1
exp2$exp <- 2
exp3$exp <- 3
exp1$duration <- 0.75
exp2$trial_type <- '1) Pure'
exp3$trial_type <- '2) Mixed'

# shift subject number so that they are unique when combining the datasets
exp2$subject <- exp2$subject + max(exp1$subject)
exp3$subject <- exp3$subject + max(exp2$subject)
dat <- bind_rows(exp1, exp2, exp3)
dat$exp <- as.factor(dat$exp)

# # extract SUBTLEX context diversity measure and attach it
# subtlex <- read.csv('expfiles/materials/SUBTLEX_all_words.csv')
# subtlex <- select(subtlex, c('Word','Lg10CD', 'SUBTLCD'))
# dat <- left_join(dat, subtlex, by=c('target_word'='Word'))

#############################################################################!
# OVERALL DENSITY PLOTS                                                  ####
#############################################################################!

# extend the data on the left and right of the error scale in order to
# get wrap-up effects for plots (to display circular density correctly)
dat3 <- bind_rows(dat, dat, dat)
dat3$angle_diff <- dat3$angle_diff + rep(c(-360,0,360), each=nrow(dat3)/3)
dat3$resp_angle <- dat3$resp_angle + rep(c(-360,0,360), each=nrow(dat3)/3)

# overall error histogram Looks like density for wrong answers peaks at the other locations
(f1 <- dat %>%
    ggplot() +
    geom_histogram(bins=60, color='darkgrey', fill='grey', aes(x=angle_diff*pi/180, y=..density..)) +
    # stat_density(position='identity', bw=3, fill='white', color='black', alpha=0) +
    coord_cartesian(xlim=c(-pi,pi+0.1), expand=0) +
    scale_x_continuous(name='Error (in degrees)', breaks=c(-180,-120,-60,0,60,120,180)*pi/180, labels=c(-180,-120,-60,0,60,120,180)))

# Fit 3p model to overall data and plot
pars <- dat %>% fit_mixture3(empirical=T)
pars$rad_sigma = pars$sigma * pi /180
pars$kappa = sd2k(pars$rad_sigma)

centers <- dat[,c('V1','V2','V4','V5')] %>%       # extract location of non-targets
  gather(key, value, V1:V5) %>% 
  mutate(value = round(value * pi / 180,3)) %>% 
  count(key, value) %>% 
  group_by(key) %>% 
  mutate(prop = n/sum(n))

x1 <- seq(-pi,pi,0.001)
y1 <- pars$p_correct * brms::dvon_mises(x1,median(dat$angle_diff)*pi/180,pars$kappa) +
  pars$p_other*colSums(centers$prop*t(sapply(centers$value, function(i) brms::dvon_mises(x1,i,pars$kappa))))/4 +
  pars$p_guess*dunif(x1,min=-pi,max=pi)

f1 + geom_line(data=data.frame(x=x1, y=y1), aes(x,y))
ggsave('figures/histogram.png', width=3.38, height=3, units='in')


# plot some model for illustration
pars <- data.frame(p_correct=0.3, p_other=0.2, p_guess = 0.5, sigma = 14)
pars$rad_sigma = pars$sigma * pi /180
pars$kappa= (1/pars$rad_sigma) ** 2
x1 <- seq(-pi,pi,0.001)
y1 <- pars$p_correct*brms::dvon_mises(x1,median(dat$angle_diff)*pi/180,pars$kappa) +
  pars$p_other*rowSums(sapply(colMeans(dat[,c('V1','V2','V4','V5')])/180*pi, function(i) brms::dvon_mises(x1,i,pars$kappa)))/4 +
  pars$p_guess*dunif(x1,min=-pi,max=pi)
sim <- data.frame(x=x1,y=y1)
ggplot(sim, aes(x,y)) +
  geom_line() +
  scale_x_continuous(breaks=colMeans(dat[,c('V1','V2','V3', 'V4','V5')])*pi/180, labels=c('NT','NT','T','NT','NT'))
ggsave('figures/model.png', width=4, height=2, units='in')


#############################################################################!
# MAIN PLOTS                                                             ####
#############################################################################!

# overall position effect on error
dat %>% 
  # filter(exp==3) %>% 
  group_by(subject, target_word_place, exp) %>% 
  summarise(error = sd(angle_diff)) %>%
  group_by(exp) %>%
  do(Rmisc::normDataWithin(.,'subject', 'error')) %>%
  ggplot(aes(target_word_place, errorNormed, color=exp)) +
  stat_summary(fun.data = mean_se, geom="pointrange") +
  stat_summary(fun.data = mean_se, geom="errorbar", width=0.1) +
  stat_summary(fun.data = mean_se, geom="line") +
  scale_x_continuous('Serial position of test word') +
  scale_y_continuous('Absolute angle error (in degrees)')


# overall duration effect on error
dat %>% 
  group_by(subject, duration, exp) %>% 
  filter(exp != 1) %>% 
  summarise(error = sd(angle_diff)) %>%
  group_by(exp) %>% 
  do(Rmisc::normDataWithin(., 'subject', 'error')) %>% 
  ggplot(aes(duration, errorNormed, color=exp)) +
  stat_summary(fun.data = mean_se, geom="pointrange") +
  stat_summary(fun.data = mean_se, geom="errorbar", width=0.05) +
  stat_summary(fun.data = mean_se, geom="line") +
  scale_x_continuous('Presentation rate', breaks=c(0.5,0.75,1)) +
  scale_y_continuous('Absolute angle error (in degrees)') 



# overall frequency (exp 1 & 2)
theme_set(theme_classic(base_size=8))
dat %>% 
  filter(exp != 3) %>% 
  group_by(subject, freq, exp) %>% 
  summarise(error = sd(angle_diff)) %>%
  group_by(exp) %>% 
  do(Rmisc::normDataWithin(., 'subject', 'error')) %>% 
  ggplot(aes(freq, errorNormed, color=exp, group=exp)) +
  stat_summary(fun.data = mean_se, geom="pointrange") +
  stat_summary(fun.data = mean_se, geom="errorbar", width=0.05) +
  stat_summary(fun.data = mean_se, geom="line") +
  scale_y_continuous('Absolute angle error (in degrees)') +
  scale_x_discrete('Word frequency') +
  scale_color_discrete('Experiment', labels=c('1: Pure lists','2: Pure lists','3: Mixed lists')) +
  theme(legend.position = c(1,0.05), legend.justification = c(1,0)) 
ggsave('figures/raw_freq_e12.png', width=3.38, height=3, units='in')


# overall frequency (exp 1 & 2 & 3)
theme_set(theme_classic(base_size=8))
dat %>% 
  group_by(subject, freq, exp) %>% 
  summarise(error = sd(angle_diff)) %>%
  group_by(exp) %>% 
  do(Rmisc::normDataWithin(., 'subject', 'error')) %>% 
  ggplot(aes(freq, errorNormed, color=exp, group=exp)) +
  stat_summary(fun.data = mean_se, geom="pointrange") +
  stat_summary(fun.data = mean_se, geom="errorbar", width=0.05) +
  stat_summary(fun.data = mean_se, geom="line") +
  scale_y_continuous('Absolute angle error (in degrees)') +
  scale_x_discrete('Word frequency') +
  scale_color_discrete('Experiment', labels=c('1: Pure lists','2: Pure lists','3: Mixed lists')) +
  theme(legend.position = c(1,0.05), legend.justification = c(1,0)) 
ggsave('figures/raw_freq.png', width=3.38, height=3, units='in')


# overall frequency by duration by exp
dat %>% 
  group_by(subject, duration, exp, freq) %>% 
  filter(exp != 1, !(subject %in% c(55,41, 99))) %>%
  summarise(error = sd(angle_diff)) %>%
  group_by(exp) %>% 
  do(Rmisc::normDataWithin(., 'subject', 'error')) %>% 
  ungroup() %>% 
  mutate(exp = ifelse(exp==2, "2) Pure freq exp", "3) Mixed freq exp")) %>% 
  ggplot(aes(duration, errorNormed, color=freq)) +
  stat_summary(fun.data = mean_se, geom="pointrange") +
  stat_summary(fun.data = mean_se, geom="errorbar", width=0.05) +
  stat_summary(fun.data = mean_se, geom="line") +
  scale_x_continuous('Presentation rate', breaks=c(0.5,0.75,1)) +
  scale_y_continuous('Absolute angle error (in degrees)') +
  scale_color_discrete('Word frequency') +
  facet_wrap(~exp) +
  theme(legend.position = 'bottom')
ggsave('figures/freq_by_duration_raw.png', width=3.38, height=2.5, units='in')


# overall duration effect by word position on error
dat %>% 
  filter(exp != 1) %>%
  group_by(subject, target_word_place, duration) %>% 
  summarise(error = sd(angle_diff)) %>%
  Rmisc::normDataWithin('subject', 'error') %>% 
  ggplot(aes(target_word_place, errorNormed, color=as.factor(duration))) +
  stat_summary(fun.data = mean_se, geom="pointrange") +
  stat_summary(fun.data = mean_se, geom="errorbar", width=0.1) +
  stat_summary(fun.data = mean_se, geom="line") +
  scale_x_continuous('Serial position of test word') +
  scale_y_continuous('Absolute angle error (in degrees)')


# error as a function of continuous frequency
dat %>% 
  mutate(LogWF = round(LogWF_target, 2)) %>% 
  filter(exp != 3) %>% 
  group_by(subject, LogWF, freq, exp) %>% 
  summarise(error = sd(angle_diff)) %>%
  ggplot(aes(LogWF, error, color=freq, group=freq)) +
  stat_summary(fun.data = mean_se, geom="point", size=0.5) +
  geom_smooth(method='lm', se=F, color='black', size=0.5) +
  coord_cartesian(ylim=c(0,60)) +
  scale_x_continuous('Word frequency (per million words)', 
                     breaks=c(log10(2),log10(3),log10(10), log10(100), log10(1000)), 
                     labels=c(1,2,10,100,1000)) +
  scale_y_continuous('Absolute angle error (in degrees)') +
  scale_color_discrete('') +
  theme_classic(base_size=9) +
  theme(legend.position=c(0.02,0.02), legend.justification=c(0,0))
ggsave('figures/continuous_freq.png', width=3.38, height=2.5, units='in')


# table of freq effect mean diff and se by experiment
dat %>% 
  group_by(subject, freq, exp) %>% 
  filter(serial_position <= 3) %>% 
  summarise(error = mean(abs_angle_diff)) %>% 
  spread(freq, error) %>% 
  mutate(diff=HF-LF) %>% 
  group_by(exp) %>% 
  summarise(mean = mean(diff),
            se = sd(diff)/sqrt(length(exp)))

#############################################################################!
# FIT MIXTURE MODELs                                                     ####
#############################################################################!

# fit mixture model to serial position
sp_pars <- dat %>% 
  group_by(subject, target_word_place, exp) %>% 
  do(fit_mixture3(.))

sp_pars %>% 
  do(Rmisc::normDataWithin(., 'subject', 'p_correct')) %>% 
  do(Rmisc::normDataWithin(., 'subject', 'p_other')) %>% 
  do(Rmisc::normDataWithin(., 'subject', 'p_guess')) %>% 
  do(Rmisc::normDataWithin(., 'subject', 'sigma')) %>% 
  gather(par, value, p_correctNormed, p_otherNormed, p_guessNormed, sigmaNormed) %>% 
  mutate(par = recode(par, p_correctNormed='P(correct)', p_otherNormed='P(misbinding)', p_guessNormed='P(guess)', sigmaNormed='Sigma')) %>% 
  ggplot(aes(target_word_place, value, color=exp)) +
  stat_summary(fun.data = mean_se, geom="point") +
  stat_summary(fun.data = mean_se, geom="line") +
  stat_summary(fun.data = mean_se, geom="errorbar", width=0.1) +
  scale_x_continuous('Serial position of probe word') +
  facet_wrap(~par, scales='free') +
  facet_wrap(~par, scales='free', nrow=1) +
  theme_bw(base_size=8) +
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.title.y = element_blank())
ggsave('figures/serial_position_pars_by_exp.png', width=7, height=1.9, units='in')




# overall parameters for every subject
subjpars <- dat %>% 
  group_by(subject) %>% 
  do(fit_mixture3(.))

subjpars %>% 
  gather(key, value, p_correct, p_other, p_guess, sigma) %>% 
  ggplot(aes(subject, value, color=key)) +
  stat_summary(fun.data = mean_se, geom="point") +
  stat_summary(fun.data = mean_se, geom="line") +
  facet_wrap(~key, scales='free')

subjpars %>% 
  rename(p_misbinding = p_other) %>% 
  gather(key, value, p_correct, p_misbinding, p_guess, sigma) %>% 
  ggplot(aes(value, y=..density..)) +
  geom_histogram(fill='lightgrey', color='black', size=0.2, bins=15) +
  facet_wrap(~key, scales='free', nrow=1) +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) 
ggsave('figures/par_distribution.png', width=3.38, height=1)

#correlation between parameters for every subject
subjpars %>% 
  ungroup() %>% 
  select(-negll) %>% 
  as.matrix() %>% 
  Hmisc::rcorr()




# fit mixture model to freq
freqpars12 <- dat %>% 
  filter(exp != 3) %>%
  group_by(subject, freq, exp) %>%
  do(fit_mixture3(.))

t.test(p_correct ~ freq, data=filter(freqpars12, exp==1), paired=T)
t.test(p_other ~ freq, data=filter(freqpars12, exp==1), paired=T)
t.test(p_guess ~ freq, data=filter(freqpars12, exp==1), paired=T)
t.test(sigma ~ freq, data=filter(freqpars12, exp==1), paired=T)

t.test(p_correct ~ freq, data=filter(freqpars12, exp==2), paired=T)
t.test(p_other ~ freq, data=filter(freqpars12, exp==2), paired=T)
t.test(p_guess ~ freq, data=filter(freqpars12, exp==2), paired=T)
t.test(sigma ~ freq, data=filter(freqpars12, exp==2), paired=T)

freqpars12 %>% 
  group_by(exp) %>% 
  do(Rmisc::normDataWithin(., 'subject', 'p_correct')) %>% 
  do(Rmisc::normDataWithin(., 'subject', 'p_other')) %>% 
  do(Rmisc::normDataWithin(., 'subject', 'p_guess')) %>% 
  do(Rmisc::normDataWithin(., 'subject', 'sigma')) %>% 
  gather(par, value, p_correctNormed, p_otherNormed, p_guessNormed, sigmaNormed) %>% 
  mutate(par = recode(par, p_correctNormed='P(correct)', p_otherNormed='P(misbinding)', p_guessNormed='P(guess)', sigmaNormed='Sigma')) %>% 
  ggplot(aes(freq, value, linetype=exp, group=exp, shape=exp)) +
  stat_summary(fun.data = mean_se, geom="point") +
  stat_summary(fun.data = mean_se, geom="errorbar", width=0.15, linetype=1) +
  stat_summary(fun.data = mean_se, geom="line") +
  scale_linetype_discrete('Exp.') +
  scale_shape_discrete('Exp.') +
  facet_wrap(~par, scales='free', nrow=1) +
  theme_bw(base_size=8) +
  theme(panel.grid = element_blank(),
        legend.position = 'bottom',
        axis.title = element_blank())
ggsave('figures/freq_pars_e12.png', width=7, height=2.2, units='in')



# parameters as a function of duration
durpars <- dat %>% 
  filter(exp != 1) %>%
  group_by(subject, exp, duration) %>%
  do(fit_mixture3(.))


durpars %>% 
  # filter(exp == 3) %>%
  Rmisc::normDataWithin('subject', 'p_correct') %>% 
  Rmisc::normDataWithin('subject', 'p_other') %>% 
  Rmisc::normDataWithin('subject', 'p_guess') %>% 
  Rmisc::normDataWithin('subject', 'sigma') %>% 
  gather(par, value, p_correctNormed, p_otherNormed, p_guessNormed, sigmaNormed) %>% 
  mutate(par = recode(par, p_correctNormed='P(correct)', p_otherNormed='P(misbinding)', p_guessNormed='P(guess)', sigmaNormed='Sigma')) %>% 
  ggplot(aes(duration, value)) +
  stat_summary(fun.data = mean_se, geom="point") +
  stat_summary(fun.data = mean_se, geom="errorbar", width=0.05) +
  stat_summary(fun.data = mean_se, geom="line") +
  scale_x_continuous("Presentation rate", breaks=c(0.5,0.75,1), labels=c('500ms','750ms','1000ms')) +
  facet_wrap(~par, scales='free_y', nrow=1) +
  theme_bw(base_size=8) +
  theme(panel.grid = element_blank(),
        legend.position = 'bottom',
        axis.title = element_blank())
ggsave('figures/dur_pars_e12.png', width=7, height=1.8, units='in')



# parameters by block
par_block <- dat %>% 
  # mutate(block = floor((block-3)/10)*10) %>% 
  group_by(subject, block) %>% 
  do(fit_mixture3(.))

par_block %>%
  gather(par, value, p_correct:p_guess) %>% 
  ggplot(aes(block, value)) +
  stat_summary(fun.data = mean_se, geom="pointrange") +
  stat_summary(fun.data = mean_se, geom="line") +
  # geom_point() +
  # geom_smooth(se=FALSE) +
  # geom_line() +
  facet_wrap(~par, scales='free')




#############################################################################
# CONFUSIBILITY (SOURCE OF ERRORS)
#############################################################################

LL_circ1 <- function(dat, sigma) {
  # transform x from degrees to radians
  rad = dat[c('y','angle1_s', 'angle2_s','angle3_s','angle4_s','angle5_s')] * pi / 180
  # transform the normal sd into radians kappa for circular vonmises concentration parameter
  rad_sigma = sigma * pi /180
  kappa = (1/rad_sigma) ** 2
  l_norm1 <- brms::dvon_mises(rad$y, mu=rad$angle1_s, kappa=kappa)
  l_norm2 <- brms::dvon_mises(rad$y, mu=rad$angle2_s, kappa=kappa)
  l_norm3 <- brms::dvon_mises(rad$y, mu=rad$angle3_s, kappa=kappa)
  l_norm4 <- brms::dvon_mises(rad$y, mu=rad$angle4_s, kappa=kappa)
  l_norm5 <- brms::dvon_mises(rad$y, mu=rad$angle5_s, kappa=kappa)
  l_unif <- brms::dvon_mises(rad$y, mu=0, kappa=0)
  return(data.frame(p1 = l_norm1, p2 = l_norm2, p3 = l_norm3, p4 = l_norm4, p5 = l_norm5, g = l_unif))
}

ps <- LL_circ1(dat, sigma=fit_mixture3(dat)$sigma)
ps <- ps / rowSums(ps)

dat2 <- cbind(dat, round(ps, 4)) %>% filter(subject != 19)

dat2 %>% 
  group_by(subject, serial_position, freq) %>% 
  summarise(p1 = mean(p1), p2 = mean(p2), p3 = mean(p3), p4 = mean(p4), p5 = mean(p5)) %>% 
  gather(key, value, p1, p2, p3, p4, p5) %>%
  mutate(key = as.numeric(as.factor(key))-as.numeric(serial_position)) %>% 
  # filter(key != 0) %>%
  ggplot(aes(key, value, color=as.factor(serial_position), group=as.factor(serial_position))) +
  stat_summary(geom='pointrange') +
  stat_summary(geom='point') +
  stat_summary(geom='line') +
  facet_grid(~serial_position) +
  # scale_x_discrete(labels=c(1:5)) +
  scale_color_discrete('Serial position of word probe') +
  xlab('Serial position of location') 

dat2 %>% 
  gather(key, value, p1, p2, p3, p4, p5) %>% 
  ggplot(aes(key, value, color=freq, group=freq)) +
  stat_summary(geom='pointrange') +
  stat_summary(geom='line') +
  facet_grid(~target_word_place) +
  scale_x_discrete(labels=c(1:5)) +
  scale_color_discrete('Serial position of word probe') +
  xlab('Serial position of location')


dat2 %>% 
  gather(key, value, p1, p2, p3, p4, p5) %>% 
  ggplot(aes(key, value, color=as.factor(duration), group=as.factor(duration))) +
  stat_summary(geom='pointrange') +
  stat_summary(geom='line') +
  facet_grid(~target_word_place) +
  scale_x_discrete(labels=c(1:5)) +
  scale_color_discrete('Serial position of word probe') +
  xlab('Serial position of location')


dat2 %>% 
  group_by(subject, serial_position) %>% 
  summarise(p1 = mean(p1), p2 = mean(p2), p3 = mean(p3), p4 = mean(p4), p5 = mean(p5)) %>% 
  gather(key, value, p1, p2, p3, p4, p5) %>% 
  mutate(key = as.numeric(as.factor(key))-as.numeric(serial_position)) %>% 
  ggplot(aes(key, value)) +
  stat_summary(geom='pointrange') +
  stat_summary(geom='point') +
  stat_summary(geom='line') +
  xlab('Serial position of location')


dat2 %>% 
  group_by(subject, serial_position) %>% 
  summarise(p1 = mean(p1), p2 = mean(p2), p3 = mean(p3), p4 = mean(p4), p5 = mean(p5)) %>% 
  gather(key, value, p1, p2, p3, p4, p5) %>% 
  mutate(key = as.numeric(as.factor(key))-as.numeric(serial_position)) %>% 
  filter(key != 0) %>% 
  ggplot(aes(key, value, group=as.factor(key>0))) +
  stat_summary(geom='pointrange') +
  stat_summary(geom='point') +
  stat_summary(geom='line') +
  xlab('Lag')

dat2 %>% 
  # filter(exp == 2) %>% 
  group_by(subject, serial_position, freq) %>% 
  summarise(p1 = mean(p1), p2 = mean(p2), p3 = mean(p3), p4 = mean(p4), p5 = mean(p5)) %>% 
  gather(key, value, p1, p2, p3, p4, p5) %>% 
  mutate(key = as.numeric(as.factor(key))-as.numeric(serial_position)) %>% 
  filter(key != 0) %>% 
  ggplot(aes(key, value, group=interaction(as.factor(key>0),as.factor(freq)), color=as.factor(freq))) +
  stat_summary(geom='pointrange') +
  stat_summary(geom='point') +
  stat_summary(geom='line') +
  xlab('Serial position of location')

dat2 %>% 
  group_by(subject, serial_position) %>% 
  summarise(p1 = mean(p1), p2 = mean(p2), p3 = mean(p3), p4 = mean(p4), p5 = mean(p5)) %>% 
  gather(key, value, p1, p2, p3, p4, p5) %>% 
  mutate(key = as.numeric(as.factor(key))-as.numeric(serial_position)) %>% 
  filter(key != 0) %>% 
  ggplot(aes(key, value, group=interaction(as.factor(key>0)))) +
  stat_summary(geom='pointrange') +
  stat_summary(geom='point') +
  stat_summary(geom='line') +
  xlab('Serial position of location') +
  facet_wrap(~subject)



# -------------------------------------------------------------------
# Performance change over time
# -------------------------------------------------------------------

dat <- dat %>% 
  group_by(subject, block) %>% 
  mutate(withinblocktrial = 1:length(block))


dat %>% 
  filter(block > 1) %>% 
  group_by(withinblocktrial, subject) %>% 
  summarise(error = sd(angle_diff)) %>% 
  ggplot(aes(withinblocktrial, error)) +
  stat_summary(fun.data = mean_se, geom="point") +
  coord_cartesian(ylim=c(47,59)) +
  xlab('Within-block trial number') +
  ylab('Absolute angle error') +
  geom_smooth(method='lm', se=F) +
  theme_classic()

dat %>% 
  group_by(block, subject) %>% 
  summarise(error = sd(angle_diff)) %>% 
  ggplot(aes(block, error)) +
  stat_summary(fun.data = mean_se, geom="pointrange") +
  xlab('Block number') +
  ylab('Absolute angle error') + 
  stat_summary(fun.data = mean_se, geom="line") +
  theme_classic()


ml0 <- lmer(abs_angle_diff ~ block + (1|subject) + (1|target_word), data=dat)
ml1 <- lmer(abs_angle_diff ~ withinblocktrial + block + (1|subject) + (1|target_word), data=dat)
anova(ml0, ml1)

dat <- dat %>% 
  group_by(subject, block) %>% 
  prior_item_analysis('rt','rt')

dat %>% 
  mutate(rt = as.numeric(rt)) %>% 
  filter(block > 1, !is.na(rt_prioritem)) %>% 
  group_by(subject) %>% 
  mutate(rt_prioritem = as.numeric(rt_prioritem),
         rt_priorcat = cut(rt_prioritem, quantile(rt_prioritem, probs=seq(0,1,0.1)), labels=F)) %>% 
  group_by(rt_priorcat, subject) %>% 
  summarise(error = sd(angle_diff, na.rm=T)) %>%
  ungroup() %>% 
  Rmisc::normDataWithin('subject','error', na.rm=T) %>%
  ggplot(aes(rt_priorcat, errorNormed)) +
  stat_summary(fun.data = mean_se, geom="pointrange") +
  geom_smooth(method='lm', se=F)


# if preceding trial is a HF vs LF
e12 <- dat %>%
  group_by(subject, block) %>% 
  prior_item_analysis('freq','freq','LF')


e12 %>% 
  filter(!is.na(freq_prioritem), exp==3, freq=='LF') %>% 
  group_by(freq_prioritem, subject, duration) %>% 
  # summarise(error = abs(angle_diff)) %>% 
  # ungroup() %>% 
  # Rmisc::normDataWithin('subject', 'error') %>% 
  ggplot(aes(duration, abs_angle_diff, color=freq_prioritem, group=freq_prioritem)) +
  stat_summary(fun.data = mean_se, geom="pointrange") +
  stat_summary(fun.data = mean_se, geom="line")




