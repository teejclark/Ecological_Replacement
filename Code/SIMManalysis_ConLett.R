#### Stable Isotope Mixing Model Analysis Code
#### T.J. Clark - 10/30/19
#### Code used for "A wolf in fox's clothing? Using stable isotopes to quantify ecological replacement"
#### accepted for Conservation Letters

library("ggplot2")
library("dplyr")
library("simmr")

###########################################################################

#### Data manip.
samp <- read.csv("SamplesMaster_final.csv")
sagf <- data.frame(samp[which(samp$sp=="SAGF"),])

# average over proximal and distal samples
sagf <- samp %>% filter(sp=="SAGF" & HairPart != "Underfur") %>%
  select(individual,d15N,d13C) %>%
  group_by(individual) %>%
  summarize(d15N = mean(d15N), d13C = mean(d13C))

# reorganize SAGF for analysis
mix <- matrix(c(sagf$d13C,sagf$d15N), ncol=2, nrow=32)
colnames(mix) <- c("d13C","d15N")

# calculate means and sds for food!
ugo <- samp %>% filter(sp=="UGO" & provenance == "field") %>%
  select(individual, d13C, d15N) %>%
  summarize(sdC = sd(d13C), sdN = sd(d15N),
    d13C = mean(d13C), d15N = mean(d15N))

kgo <- samp %>% filter(sp=="KGO" & provenance == "field") %>%
  select(individual, d13C, d15N) %>%
  group_by(individual) %>% # average over individual
  summarize(d13C = mean(d13C), d15N = mean(d15N)) %>%
  summarize(sdC = sd(d13C), sdN = sd(d15N),
            d13C = mean(d13C), d15N = mean(d15N))

penguin <- samp %>% filter(sp == "Penguin sp.") %>%
  filter(sample != "BV0303a_B") %>%
  select(individual, d13C, d15N) %>%
  group_by(individual) %>% # average over individual
  summarize(d13C = mean(d13C), d15N = mean(d15N)) %>%
  summarize(sdC = sd(d13C), sdN = sd(d15N),
            d13C = mean(d13C), d15N = mean(d15N))

cricket <- samp %>% filter(sp == "camel_cricket") %>%
  select(individual, d13C, d15N) %>%
  summarize(sdC = sd(d13C), sdN = sd(d15N),
           d13C = mean(d13C), d15N = mean(d15N))

diddle_dee <- samp %>% filter(sp == "diddle_dee") %>%
  select(individual, d13C, d15N) %>%
  summarize(sdC = sd(d13C), sdN = sd(d15N),
            d13C = mean(d13C), d15N = mean(d15N))

mtn_berry <- samp %>% filter(sp == "mountain_berry") %>%
  select(individual, d13C, d15N) %>%
  summarize(sdC = sd(d13C), sdN = sd(d15N),
            d13C = mean(d13C), d15N = mean(d15N))

limpet <- samp %>% filter(sp == "common_limpet" | sp == "keyhole_limpet") %>%
  select(individual, d13C, d15N) %>%
  summarize(sdC = sd(d13C), sdN = sd(d15N),
            d13C = mean(d13C), d15N = mean(d15N))

#################################################################################
#### EXTRA ANALYSIS - use SIDER to convert SI feather values to muscle...

library(SIDER)

SIDER_data <- scrumpSider(iso.data="all")
SIDER_trees <- scrumpSider(tree="all")

ugo_test <- recipeSider(species = "Chrysocyon_brachyurus",
                        taxonomic.class = "mammalia",
                        tissue = "hair",
                        diet.type = "omnivore",
                        habitat = "terrestrial", # could also be marine?
                        tree = SIDER_trees)

UGO_tdf_data_n <- prepareSider(data.estimate = ugo_test,
                               data.isotope = SIDER_data,
                               tree = SIDER_trees,
                               isotope="nitrogen")
formula_n <- delta15N ~ diet.type + habitat
random_terms <- ~ animal + sp.col + tissue

UGO_tdf_est <- imputeSider(mulTree.data = UGO_tdf_data_n,
                           formula=formula_n,
                           output="test_n_run")

summary(UGO_tdf_est$tdf_global)

#################################################################################
#### Re-adjust correction values for Trophic Enrichment Factors
# model adds whatever value is inputted into the mean and sd
# order is: ugo, kgo, penguin, cricket, diddle-dee, mountain berry, limpet
# fox TEF = add 2.6 for carbon, 3.2 for nitrogen
# UGO TEF = subtract 0.55 for carbon, 1.3 for nitrogen
# KGO TEF = subtract 0.53 for carbon, 1.27 for nitrogen
# Pen TEF = subtract 0.57 for carbon, 1.29 for nitrogen
# Seabird TEF = subtract 0.56 for carbon, 1.27 for nitrogen
# Sealion TEF = subtract 0.8 for carbon, 0.04 for nitrogen

c_means <-matrix(c(2.05,2.07,2.03,2.6,2.6,2.6,2.6,
                   1.9,1.93,1.91,3.2,3.2,3.2,3.2),
                 ncol=2,nrow=7)

# keep sds the same...
c_sds <- matrix(c(rep(.06,7), rep(.06,7)),
                ncol=2,nrow=7)
#################################################################################

# reorganize diet for analysis
s_names <- c("UGO","KGO","Penguin","Cricket","Diddle-dee",
             "Mtn Berry","Limpet")
s_means <- matrix(c(ugo$d13C,kgo$d13C,penguin$d13C,cricket$d13C,
                   diddle_dee$d13C,mtn_berry$d13C,limpet$d13C,
                   ugo$d15N,kgo$d15N,penguin$d15N,cricket$d15N,
                   diddle_dee$d15N,mtn_berry$d15N,limpet$d15N),
                 ncol=2,nrow=7)
s_sds <- matrix(c(ugo$sdC,kgo$sdC,penguin$sdC,cricket$sdC,
                  diddle_dee$sdC,mtn_berry$sdC,limpet$sdC,
                  ugo$sdN,kgo$sdN,penguin$sdN,cricket$sdN,
                  diddle_dee$sdN,mtn_berry$sdN,limpet$sdN),
                ncol=2,nrow=7)


# load data for SAGF and plot
simmr_in <- simmr_load(mixtures=mix,
                       source_names = s_names,
                       source_means = s_means,
                       source_sds = s_sds,
                       correction_means = c_means,
                       correction_sds = c_sds)

plot(simmr_in,
        xlab=expression(paste(delta^13, "C (\u2030)",sep="")),
        ylab=expression(paste(delta^15, "N (\u2030)",sep="")),
        title="Isospace plot of data")

#################################################################################

#### Running the Model
# run without priors, initially
simmr_out1 <- simmr_mcmc(simmr_in)

# check stuff
summary(simmr_out1,type="diagnostics") # looks converged
posterior_predictive(simmr_out1) # looks OK...

# check priors
windows();prior_viz(simmr_out1)

# check %
summary(simmr_out1, type="quantiles")
# looks like lots of inland berries, big mix of stuff...
# good amount of uncertainty, which is to be expected

# plot everything

windows();plot(simmr_out1, type="density")
# looks like limpets, KGO, penguins

windows();plot(simmr_out1, type="matrix")
# nothing is clear to combine - as the authors say:
# "are an unavoidable part of stable isotope mixing models"

windows();plot(simmr_out1,type="boxplot")

##################################################################################

#### Running the Model with priors from Bugge's data
# based on % volume found in scat, January 2019
#s_names <- c("UGO","KGO","Penguin","Cricket","Diddle-dee",
#             "Mtn Berry","Limpet")

proportion_means <- c(.0709,.0709,.060,.622,.008,.165,.002)
proportion_sds <- rep(0.1,7)

prior <- simmr_elicit(7, proportion_means, proportion_sds)

# run model with new priors
simmr_out2 <- simmr_mcmc(simmr_in,
                         prior_control=list(means=prior$mean,
                                            sd=prior$sd))

# new quantiles
summary(simmr_out2, type="quantiles")

# plot
windows();plot(simmr_out2, type="boxplot")
windows();plot(simmr_out2, type="density")

##################################################################################

#### COMBINE STUFF - add in SAMPLES
# added samples need to look like
# mtn_berry <- samp %>% filter(sp == "mountain_berry") %>%
#   select(individual, d13C, d15N) %>%
#   summarize(sdC = sd(d13C), sdN = sd(d15N),
#             d13C = mean(d13C), d15N = mean(d15N))

# add lamb, sea lion, sea birds
# reorganize diet for analysis

# lamb = -23.8 C, 0.57SE; 5.3 N, 1.03SE // Piasentier et al. (2003) Meat Science

# sea lion1 = -12.48 C, 0.68SE; 20.97 N, 0.77SE // Huckstadt et al. (2007) Journal of Experimental Marine Biology and Ecology
# sea lion2 = -14 C, 0.65SE; 16.4 N, 0.5SE // Baylis et al. (2015) Oecologia - female sea lion vibrissae

# Wilson's storm petrel = -19.1 C, 2.24SE; 13.68 N, 1.44SE // Quillfeldt et al. (2005) MEPS
# Seabirds = -17.2 C, 2.28SE; 14.63 N, 2.5SE // Quillfeldt et al. (2009) Polar Bio

s_names <- c("UGO","KGO","Penguin","Terrestrial Invertebrates","Diddle-dee",
             "Mtn Berry","Limpet","Sea Lion","Other Seabirds",
             "Sheep")

s_means <- matrix(c(ugo$d13C,kgo$d13C,penguin$d13C,cricket$d13C,
                    diddle_dee$d13C,mtn_berry$d13C,limpet$d13C,
                    -13.06,-17.2,-23.8,
                    ugo$d15N,kgo$d15N,penguin$d15N,cricket$d15N,
                    diddle_dee$d15N,mtn_berry$d15N,limpet$d15N,
                    18.75,14.63,5.3),
                  ncol=2,nrow=10)

s_sds <- matrix(c(ugo$sdC,kgo$sdC,penguin$sdC,cricket$sdC,
                  diddle_dee$sdC,mtn_berry$sdC,limpet$sdC,
                  0.73,2.28,0.57,
                  ugo$sdN,kgo$sdN,penguin$sdN,cricket$sdN,
                  diddle_dee$sdN,mtn_berry$sdN,limpet$sdN,
                  0.9,2.5,1.03),
                ncol=2,nrow=10)

#0.55 C; 1.28 N
c_means <-matrix(c(2.05,2.07,2.03,2.6,2.6,2.6,2.6,1.8,2.04,2.6, #C
                   1.9,1.93,1.91,3.2,3.2,3.2,3.2,3.16,1.93,3.2), #N
                 ncol=2,nrow=10)

# keep sds the same...
c_sds <- matrix(c(rep(.06,10), rep(.06,10)),
                ncol=2,nrow=10)

# load data for SAGF and plot
simmr_in <- simmr_load(mixtures=mix,
                       source_names = s_names,
                       source_means = s_means,
                       source_sds = s_sds,
                       correction_means = c_means,
                       correction_sds = c_sds)

windows();
plot(simmr_in,
     xlab=expression(paste(delta^13, "C (\u2030)",sep="")),
     ylab=expression(paste(delta^15, "N (\u2030)",sep="")),
     title="Isospace plot of data")

# without priors
simmr_out3 <- simmr_mcmc(simmr_in)
windows();plot(simmr_out3, type="boxplot")

# combine sources - without priors s_names <- c("UGO","KGO","Penguin","Cricket","Diddle-dee","Mtn Berry","Limpet","Sea Lion","Sea Birds","Sheep")

combined_1 <- combine_sources(simmr_out3,
                              to_combine = c("Diddle-dee","Mtn Berry"),
                              new_source_name = "Plants")
windows();plot(combined_1,type="boxplot")
windows();plot(combined_1$input)

combined_2 <- combine_sources(combined_1,
                              to_combine = c("UGO","Sheep"),
                              new_source_name = "Terrestrial Grazers")
windows();plot(combined_2,type="boxplot")
windows();plot(combined_2$input)

combined_3 <- combine_sources(combined_2,
                              to_combine = c("KGO","Limpet"),
                              new_source_name = "Marine Grazers")
windows();plot(combined_3,type="boxplot")
windows();plot(combined_3$input)

combined_4 <- combine_sources(combined_3,
                              to_combine = c("Penguin","Sea Lion"),
                              new_source_name = "Marine Piscivores")
windows();plot(combined_4,type="boxplot")
windows();plot(combined_4$input)

combined_5 <- combine_sources(combined_4,
                              to_combine = c("Marine Piscivores","Other Seabirds"),
                              new_source_name = "Marine Piscivores")
windows();plot(combined_5,type="boxplot")
windows();plot(combined_5$input)

#################################################################################

# okay, now add priors from Bugge's....

#### Running the Model with priors from Bugge's data
# names = Terrestrial Grazers, Marine Grazers, Marine Piscivores, Terrestrial Invertebrates, Plants
# assume 1/4 of bird species are seabirds (besides penguins and geese = 1%)
# TG = .035; MG = 035; MP = .04; TI = .625; P = .2 = .935
proportion_means <- c(.037,.037,.043,.6715,.2115)
proportion_sds <- rep(.1,5)

prior <- simmr_elicit(5, proportion_means, proportion_sds)

# run model with new priors
simmr_out_fin <- simmr_mcmc(combined_5$input,
                         prior_control=list(means=prior$mean,
                                            sd=prior$sd))

windows();plot(simmr_out_fin, type="boxplot")
summary(simmr_out_fin, type="quantiles")

#################################################################################

#### Add in the FIW

fiw <- samp %>% filter(sp=="FIW" & HairPart!="Underfur") %>%
  select(individual,d15N,d13C) %>%
  group_by(individual) %>%
  summarize(d15N = mean(d15N), d13C = mean(d13C)) # average over individual

# adjust for SUESS effect
fiw$d13C <- fiw$d13C - 3.198

# format
# add in ancient warrah from Kit Hamley's thesis, -3.198 on carbon for Suess Effect
# also presume that TEFs are the same for warrahs
mix <- matrix(c(c(sagf$d13C,fiw$d13C,-9.6-3.198-1.4,-12.3-3.198-1.4,-12.2-3.198-1.4), # subtracting 1.4
                c(sagf$d15N,fiw$d15N,22.1,17.3,17.7)), ncol=2, nrow=43)
colnames(mix) <- c("d13C","d15N")

# add groups
grp <- as.factor(c(rep(1,32),rep(2,8),rep(3,3)))
#grp <- as.factor(c(rep(1,32),rep(2,8)))

# load data
simmr_in_fiw <- simmr_load(mixtures=mix,
                       source_names = s_names,
                       source_means = s_means,
                       source_sds = s_sds,
                       correction_means = c_means,
                       correction_sds = c_sds,
                       group = grp)

windows();
plot(simmr_in_fiw, group=1:3,
     xlab=expression(paste(delta^13, "C (\u2030)",sep="")),
     ylab=expression(paste(delta^15, "N (\u2030)",sep="")),
     title="Isospace plot of data")

simmr_out_fiw <- simmr_mcmc(simmr_in_fiw)


combined_1 <- combine_sources(simmr_out_fiw,
                              to_combine = c("Diddle-dee","Mtn Berry"),
                              new_source_name = "Plants")
windows();plot(combined_1,type="boxplot")
windows();plot(combined_1$input)

combined_2 <- combine_sources(combined_1,
                              to_combine = c("UGO","Sheep"),
                              new_source_name = "Terrestrial Grazers")
windows();plot(combined_2,type="boxplot")
windows();plot(combined_2$input)

combined_3 <- combine_sources(combined_2,
                              to_combine = c("KGO","Limpet"),
                              new_source_name = "Marine Grazers")
windows();plot(combined_3,type="boxplot")
windows();plot(combined_3$input)

combined_4 <- combine_sources(combined_3,
                              to_combine = c("Penguin","Sea Lion"),
                              new_source_name = "Marine Piscivores")
windows();plot(combined_4,type="boxplot")
windows();plot(combined_4$input)

combined_5 <- combine_sources(combined_4,
                              to_combine = c("Marine Piscivores","Other Seabirds"),
                              new_source_name = "Marine Piscivores")
windows();plot(combined_5,type="boxplot")
windows();plot(combined_5$input,group=1:3,
               xlab=expression(paste(delta^13, "C (\u2030)",sep="")),
               ylab=expression(paste(delta^15, "N (\u2030)",sep="")),
               title="Isospace plot of data")

####################################################################################
