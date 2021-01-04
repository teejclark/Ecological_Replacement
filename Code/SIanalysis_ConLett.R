#### Stable Isotope Overlap Analysis Code
#### T.J. Clark - 10/7/19
#### Code used for "A wolf in fox's clothing? Using stable isotopes to quantify ecological replacement"
#### accepted for Conservation Letters

library("ggplot2")
library("dplyr")
library("nicheROVER")

###########################################################################

#### 1) Niche Overlap

## Data manip.
samp <- read.csv("SamplesMaster_final.csv")
sagf <- data.frame(samp[which(samp$sp=="SAGF"),])

# falkland wolf data - only 1 sample per individual
fiw <- samp %>% filter(sp=="FIW" & HairPart != "Underfur") %>%
  select(individual,d15N,d13C) %>%
  group_by(individual) %>%
  summarize(d15N = mean(d15N), d13C = mean(d13C)) # average over individual

# adjust for SUESS effect
# by averaging over deltad13C b/t historical and current samples
kgosamp <- samp %>% filter(sp=="KGO") %>%
  group_by(provenance) %>%
  summarize(d13C = mean(d13C))
kgodiff <- diff(kgosamp$d13C) # diff = 2.243
ugosamp <- samp %>% filter(sp=="UGO") %>%
  group_by(provenance) %>%
  summarize(d13C = mean(d13C))
ugodiff <- diff(ugosamp$d13C) # diff = 4.152

avdiff <- mean(c(ugodiff,kgodiff)) # adjustment = 3.198; rate of -0.033 per year

# southam. fox data - av. over prox and distal samples
sagf <- samp %>% filter(sp=="SAGF" & HairPart != "Underfur") %>%
  select(individual,d15N,d13C) %>%
  group_by(individual) %>%
  summarize(d15N = mean(d15N), d13C = mean(d13C))

# combine dataframes
can <- bind_rows(fiw,sagf) %>%
  mutate(sp=c(rep("FIW",8),rep("SAGF",32)))
can <- data.frame(can) # fix dplyr goofiness...

# check out means by species
aggregate(can[2:3], can[4], mean)

######################################################################################

# compare temporal shifts in diet - SAGF
sagf_temp <- samp %>% filter(sp=="SAGF" & HairPart != "Underfur")
sagf_temp$HairPart <- factor(sagf_temp$HairPart)
windows();plot(sagf_temp$d15N ~ sagf_temp$HairPart,border="red")
windows();plot(sagf_temp$d13C ~ sagf_temp$HairPart,border="red")
t.test(sagf_temp$d15N ~ sagf_temp$HairPart) # t=0.39, p = .69
t.test(sagf_temp$d13C ~ sagf_temp$HairPart) # t=0.71, p = .48

# compare temporal shifts in diet - FIW
fiw_temp <- samp %>% filter(sp=="FIW" & HairPart != "Underfur" & HairPart != "Whole")
fiw_temp$HairPart <- factor(fiw_temp$HairPart)
windows();plot(fiw_temp$d15N ~ fiw_temp$HairPart,border="blue")
windows();plot(fiw_temp$d13C ~ fiw_temp$HairPart,border="blue")
t.test(fiw_temp$d15N ~ fiw_temp$HairPart) # t=0.24, p=0.82
t.test(fiw_temp$d13C ~ fiw_temp$HairPart) # t=0.01, p=0.99

######################################################################################

## Create 2-D projections of niche regions

clrs <- c("red","blue")
nsamples <- 100 #
# create priors
can.par <- tapply(1:nrow(can),
                  can$sp,
                  function(ii) niw.post(nsamples=nsamples,
                                        X=can[ii,3:2]))

# format data to plot & plot
can.data <- tapply(1:nrow(can),
                   can$sp,
                   function(ii) X=can[ii,3:2])
windows()
niche.plot(niche.par = can.par,
           niche.data = can.data,
           pfrac = 0.05, # alpha level
           iso.names = expression(delta^{13}*C, delta^{15}*N),
           col=clrs,
           xlab = "Isotope Ratio")

# Calculate and statistically compare Niche Area between Species

# calculate distance from midpoint (mu) to both ends
# easy way
ctr <- c(can.par$SAGF$mu[1,1],can.par$SAGF$mu[1,2])
a_points <- sqrt(rowSums((t(t(el_mat)-ctr))^2)) # distance from center
pi*min(a_points)*max(a_points)

# double check - nerd eigenvalue way
cov_dat <- can.par$SAGF$Sigma[,,1]
eig_dat <- eigen(cov_dat)$values # eigenvalues of covariance matrix
vec <- sqrt(5.991*eig_dat) # half length of major and minor axis for 95% CIs
pi*vec[1]*vec[2]

# for loop - SAGF - use can.par - 100 draws
el_area_sagf <- rep(NA,nrow(can.par$SAGF$mu))
for (i in 1:length(el_area_sagf)){
  cov_dat <- can.par$SAGF$Sigma[,,i]
  eig_dat <- eigen(cov_dat)$values
  vec <- sqrt(5.991*eig_dat)
  el_area_sagf[i] <- pi*vec[1]*vec[2]
}

# for loop - FIW - use can.par - 100 draws
el_area_fiw <- rep(NA,nrow(can.par$FIW$mu))
for (i in 1:length(el_area_fiw)){
  cov_dat <- can.par$FIW$Sigma[,,i]
  eig_dat <- eigen(cov_dat)$values
  vec <- sqrt(5.991*eig_dat)
  el_area_fiw[i] <- pi*vec[1]*vec[2]
}

# combine and t-test?
el_area_total <- data.frame(sagf=el_area_sagf,
                            fiw=el_area_fiw)
t.test(el_area_total$sagf,el_area_total$fiw,alternative="greater") # 2-sample t-test, p<.00001

# shows:
# 1) scatterplots for isotopes
# 2) ellipses of probablistic niche regions
# 3) density estimates

## Calculate and Display Niche Overlap
nsamples <- 10000
can.par <- tapply(1:nrow(can),
                  can$sp,
                  function(ii) niw.post(nsamples=nsamples,
                                        X=can[ii,2:3]))

# calculate overlap
over.stat <- overlap(can.par, nreps=nsamples,
                     nprob=nsamples,
                     alpha=c(0.95))
# show mean and 95% CI
over.mean <- apply(over.stat, c(1:2),mean)*100
round(over.mean,2)
over.cred <- apply(over.stat*100,
                   c(1:2),
                   quantile,
                   prob=c(0.025,0.975),na.rm=T)
round(over.cred[,,1:2],2)

#############################################################################################################################
# EXTRA: Randomization Test
# use to see if small sample size in FIW could be biasing the results...Supp Fig. S6.

# 1) Calculate overlap for different sizes of N (5 to 32)
par <- c(5,10,15,20,25,32)

blah <- array(NA,dim=c(2,2,2,100,length(par))) #95% CI
bloo <- array(NA,dim=c(2,2,100,length(par))) # mean

for (j in 1:length(par)){
for (i in 1:100){
  # RANDOM SAMPLER DATA
  sagf_rand <- sample_n(sagf, size = par[j])
  can <- bind_rows(fiw,sagf_rand)%>%
    mutate(sp=c(rep("FIW",8),rep("SAGF",par[j])))
  can <- data.frame(can)

  # CALCULATE
  nsamples <- 10000
  can.par <- tapply(1:nrow(can),
                    can$sp,
                    function(ii) niw.post(nsamples=nsamples,
                                          X=can[ii,2:3]))
  # calculate overlap
  over.stat <- overlap(can.par, nreps=nsamples,
                       nprob=nsamples,
                       alpha=c(0.95))

  bloo[,,i,j] <- apply(over.stat, c(1:2),mean)*100

  blah[,,,i,j] <- apply(over.stat*100,
                     c(1:2),
                     quantile,
                     prob=c(0.025,0.975),na.rm=T)
}}

# REVISE TO INCORPORATE THE VARIED PARAMS
bloo1 <- bloo[,2,,]; bloo1 <- bloo1[-2,,]; bloo1_mean <- apply(bloo1,2,mean) # FIW -> SAGF (mean)
bloo2 <- bloo[,1,,]; bloo2 <- bloo2[-1,,]; bloo2_mean <- apply(bloo2,2,mean) # # SAGF -> FIW (mean)

blah1 <- blah[,,1,,]; blah1 <- blah1[,-1,,];
blah1low <- blah1[1,,]; blah1low <- apply(blah1low,2,mean) # SAGF -> FIW (95% CI)
blah1high <- blah1[2,,]; blah1high <- apply(blah1high,2,mean)

blah2 <- blah[,,2,,]; blah2 <- blah2[,-2,,];
blah2low <- blah2[1,,]; blah2low <- apply(blah2low,2,mean) # FIW -> SAGF (95% CI)
blah2high <- blah2[2,,]; blah2high <- apply(blah2high,2,mean)

# OTHER CALC?
# % of simulations where the 95% CIs are greater than a threshold (10,20,30,40,50%); FIW -> SAGF
blah_pct <- blah2[1,,]

blah_pct_10 <- rep(NA,length(par))
blah_pct_20 <- rep(NA,length(par))
blah_pct_30 <- rep(NA,length(par))
blah_pct_40 <- rep(NA,length(par))
blah_pct_50 <- rep(NA,length(par))

for (i in 1:length(par)){
  bloop <- blah_pct[,i]
  blah_pct_10[i] <- length(bloop[bloop>=10])/100
  blah_pct_20[i] <- length(bloop[bloop>=20])/100
  blah_pct_30[i] <- length(bloop[bloop>=30])/100
  blah_pct_40[i] <- length(bloop[bloop>=40])/100
  blah_pct_50[i] <- length(bloop[bloop>=50])/100
}

# GRAPH
library(ggplot2)

windows()
ggplot()+
  geom_line(aes(x=par,y=bloo1_mean),size=2,color="blue")+ # FIW -> SAGF
  geom_line(aes(x=par,y=blah2low),size=1,color="blue",lty="dashed")+
  geom_line(aes(x=par,y=blah2high),size=1,color="blue",lty="dashed")+
  geom_line(aes(x=par,y=bloo2_mean),size=2,color="red")+ # SAGF -> FIW
  geom_line(aes(x=par,y=blah1low),size=1,color="red",lty="dashed")+
  geom_line(aes(x=par,y=blah1high),size=1,color="red",lty="dashed")+
  xlab("Number of SAGF Samples")+
  ylab("% Overlap")+
  ylim(0,100)+
  theme_bw() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.text = element_blank(),
    legend.title = element_blank())

windows()
ggplot()+
  geom_line(aes(x=par,y=blah_pct_10),size=2,color="#ffff3f")+
  geom_line(aes(x=par,y=blah_pct_20),size=2,color="#dddf00")+
  geom_line(aes(x=par,y=blah_pct_30),size=2,color="#bfd200")+
  geom_line(aes(x=par,y=blah_pct_40),size=2,color="#55a630")+
  geom_line(aes(x=par,y=blah_pct_50),size=2,color="#007f5f")+
  xlab("Number of SAGF Samples")+
  ylab("Proportion of Simulated Overlaps\n > Threshold for FIW in SAGF")+
  ylim(0,1)+
  theme_bw() +
  theme(
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.title = element_text(size=20),
    axis.text = element_text(size=16, color="black"),
    legend.title = element_blank())


#############################################################################################################################
# Overlap Plot

clrs <- c("red","blue")
windows();
overlap.plot(over.stat,
             col=clrs,
             mean.cred.col="black",
             equal.axis=T,
             xlab="Overlap Probability/Niche Region Size (95%)")

#############################################################################################################################

