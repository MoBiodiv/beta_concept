# investigation of Bateman birds diversity: 
# standardise effort for a comparison of sites in riparian habitats from 
# either engineered or natural landscapes
library(tidyverse)

# downloaded 3rd Oct 2022   

dat <- read_csv('./data/46_core_birds_ee23527b9fad8b2ead1a6f0b4471ab1e.csv')


# separate date in year, month, day
dat <- dat %>% 
  separate(survey_date, into = c('year', 'month', 'day'), sep = '-', remove = FALSE) %>% 
  mutate(year = as.numeric(year),
         month = as.numeric(month),
         day = as.numeric(day))

# want to look at spring and winter
dat %>% 
  distinct(month)

dat <- dat %>% 
  mutate(season = case_when(month %in% c(12,1,2) ~ 'winter',
                            month %in% 3:5 ~ 'spring',
                            month %in% 6:8 ~ 'summer',
                            month %in% 9:11 ~ 'autumn',
                            TRUE ~ as.character(month)))

# to start, look at similar data to that in Banville et al. 2017
# discard year==2000, and retain only spring samples from perennial riparian habitat
riparian <- dat %>% 
  filter(location_type=='Riparian' & year > 2000) %>% 
  filter(season %in% c('spring')) %>% 
  separate(site_code, into = c('sub_type', 'site'), remove=FALSE) %>% 
  filter(sub_type %in% c('PE', 'PN')) %>% 
  # there was no sample in spring of 2003, drop whole year
  filter(year != 2003)

# Riparian habitat sub-types include: 
# (3) perennial-engineered (PE, n=3), and 
# (4) perennial-natural (PN, n=3).

# unbalanced
riparian %>% 
  distinct(site_code)

# TODO: examine sensitivity to this choice
# select 3 of each type 
riparian2 <- riparian %>% 
  filter(site_code != 'PE-13A') %>% 
  mutate(habitat = case_when(sub_type == 'EE' ~ 'ephemeral-engineered',
                             sub_type == 'EN' ~ 'ephemeral-natural',
                             sub_type == 'PE' ~ 'perennial-engineered',
                             sub_type == 'PN' ~ 'perennial-natural')) %>% 
  select(-c(sub_type, site))

# check we've got three sites for each treatment (engineered and natural)
riparian2 %>% 
  distinct(site_code)

# how many sites per year? 2009 is missing one site
riparian2 %>% 
  group_by(habitat, year) %>% 
  summarise(n_sites = n_distinct(site_code)) %>% 
  ungroup() %>% 
  filter(n_sites < 3)

# how many samples per site per season per year?
riparian2 %>% 
  group_by(site_code, year, season) %>% 
  summarise(n_samps = n_distinct(survey_date)) %>% 
  ungroup() %>% 
  # we're most interested in the minimum
  filter(n_samps == min(n_samps))


# get two samples from each site/year/season
rip2 <- riparian2 %>% 
  # some birds have NA counts?
  filter(!is.na(bird_count)) %>% 
  group_by(site_code, habitat, year, season, survey_date) %>% 
  # use the code as species ID; rename for clarity
  rename(species = code) %>% 
  nest() %>% 
  group_by(site_code, habitat, year, season) %>% 
  slice(1:2) %>% 
  # now want to combine these two samples as our alpha-grain
  unnest(data) %>% 
  group_by(site_code, habitat, year, season, species) %>% 
  summarise(N = sum(bird_count)) %>% 
  # nest these data again
  nest() %>% 
  ungroup() %>% 
  # drop 2009 (one treatment has one site less than the other)
  filter(year != 2009)

# some balance checks
rip2 %>% 
  group_by(habitat, season) %>% 
  summarise(n_sites = n_distinct(site_code))
# 2009 
rip2 %>% 
  group_by(habitat, season, year) %>% 
  summarise(n_sites = n_distinct(site_code)) %>% 
  ungroup() %>% 
  filter(n_sites != 3)
