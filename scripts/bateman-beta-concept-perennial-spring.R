# investigation of Bateman birds diversity: 
library(tidyverse)

# downloaded 3rd Oct 2022   

dat <- read_csv('~/Dropbox/MoBD (Measurements of Beta diversity)/Beta_paper/data/knb-lter-cap.46.15/46_core_birds_ee23527b9fad8b2ead1a6f0b4471ab1e.csv')


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
  ungroup()
  filter(n_sites < 3)

# how many samples per site per season per year?
riparian2 %>% 
  group_by(site_code, year, season) %>% 
  summarise(n_samps = n_distinct(survey_date)) %>% 
  ungroup() %>% 
  # we're most interested in the minimum
  filter(n_samps == min(n_samps))

# TODO: sensitivity to the samples selected
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

# first total spatiotemporal variation in community composition
# gamma-scale is all years and sites combined; alpha is a single site-year
# combination

# want to do leave-one-out (jackknife) resampling to visualise variation
counter <- rip2 %>% 
  distinct(site_code, year)

total_beta_jk <- tibble()

for(i in 1:nrow(counter)){
  print(paste(i, ' of ', nrow(counter), ' resamples'))
  temp <- rip2 %>% 
    # drop one block & year combination from each treatment
    group_by(habitat) %>% 
    arrange(-desc(year)) %>% 
    slice(-i) %>% 
    mutate(resamp = i) %>% 
    unnest(data) %>% 
    nest(data = c(site_code, year, species, N)) %>% 
    mutate(wide_data = map(data, ~pivot_wider(data = .,
                                              names_from = species,
                                              values_from = N,
                                              values_fill = 0))) %>% 
    ungroup()
  
  total_beta_jk <- bind_rows(total_beta_jk, 
                             temp)
}

# the maximum number of distinct communities in the jackknife data is: 41
# it is equal to the number of site-year combinations
total_beta_jk$wide_data[[1]] %>% nrow()

# calculate target coverage for standardisation
total_targetC_resamps <- total_beta_jk %>% 
  mutate(target_C = map(wide_data, ~mobr::C_target(x = .[,-c(1,2)], factor = 2))) %>% 
  unnest(target_C) %>% 
  ungroup() %>% 
  summarise(C_target = min(target_C))

# calculate number of individuals for beta_S_n
total_targetN_resamps <- total_beta_jk %>% 
  unnest(data) %>% 
  group_by(habitat, season, site_code, year, resamp) %>% 
  summarise(J = sum(N)) %>% 
  ungroup() %>% 
  filter(J == min(J)) %>% 
  pull(J) %>% 
  unique()

total_beta_resamps_calcs <- total_beta_jk %>% 
  mutate(beta_C = map(wide_data, possibly(~mobr::beta_C(x = .[,-c(1,2)], 
                                                        C = total_targetC_resamps$C_target, 
                                                        extrapolation = TRUE), 
                                          otherwise = NULL)),
         beta_S = map(wide_data, ~mobr::calc_comm_div(abund_mat = .[,-c(1,2)], 
                                                      index = 'S', 
                                                      scales = 'beta', 
                                                      coverage = FALSE)),
         beta_S_PIE = map(wide_data, ~mobr::calc_comm_div(abund_mat = .[,-c(1,2)], 
                                                          index = 'S_PIE', 
                                                          scales = 'beta', 
                                                          coverage = FALSE)),
         beta_S_n = map(wide_data, ~mobr::calc_comm_div(abund_mat = .[,-c(1,2)], 
                                                          index = 'S_n', 
                                                          effort = total_targetN_resamps,
                                                          scales = 'beta', 
                                                          coverage = FALSE)))

total_beta_resamps <- total_beta_resamps_calcs %>% 
  unnest(cols = beta_S) %>% 
  rename(beta_S = value) %>% 
  dplyr::select(-c(scale, index, sample_size, effort, coverage)) %>% 
  unnest(cols = beta_S_PIE) %>% 
  rename(beta_S_PIE = value) %>% 
  dplyr::select(-c(scale, index, sample_size, effort, coverage)) %>% 
  unnest(cols = beta_S_n) %>% 
  rename(beta_S_n = value) %>% 
  dplyr::select(-c(scale, index, sample_size, effort, coverage)) %>% 
  unnest(cols = beta_C) %>% 
  dplyr::select(-c(data, wide_data)) %>% 
  pivot_longer(cols = beta_C:beta_S_n, names_to = 'metric', 
               values_to = 'value') %>% 
  mutate(metric = factor(metric, levels = c('beta_S', 'beta_S_PIE', 
                                          'beta_S_n', 'beta_C'))) 

# summarise for figure
total_beta_summary <- total_beta_resamps %>% 
  group_by(habitat, metric) %>% 
  summarise(location = mean(value),
            Q5 = quantile(value, probs = c(0.05)),
            Q95 = quantile(value, probs = c(0.95))) %>% 
  ungroup()

ggplot() +
  geom_point(data = total_beta_summary,
             aes(x = metric, y = location, colour = habitat, shape = habitat,
                 group = habitat),
             size = 1.5, 
             position = position_dodge(width = 0.1)) +
  geom_linerange(data = total_beta_summary,
                 aes(x = metric, ymin = Q5, ymax = Q95,
                     colour = habitat, group = habitat),
                 position = position_dodge(width = 0.1)) +
  scale_x_discrete(limits = c('beta_S', 'beta_S_PIE', 'beta_S_n', 'beta_C'),
                   labels = c(expression(beta[S]), expression(beta[S[PIE]]), 
                              expression(beta[S[n]]), expression(beta[C])),
                   name = 'Metric') +
  scale_colour_manual(name = 'Habitat',
                      values = c('perennial-engineered' = '#a6611a',
                                 'perennial-natural' = '#80cdc1'),
                      label = c('engineered',
                                'natural'),
                      # guide = 'none'
                      )  +
  scale_shape_manual(name = 'Habitat',
                     values = c('perennial-engineered' = 19,
                                'perennial-natural' = 17),
                     label = c('engineered',
                               'natural'),
                     # guide = 'none'
                     )  +
  labs(y = 'Metric value') +
  theme_bw() +
  theme(legend.position = c(1,1),
        legend.direction = 'vertical',
        legend.justification = c(1,1),
        legend.background = element_blank())

ggsave('~/Dropbox/MoBD (Measurements of Beta diversity)/Beta_paper/figs/case-study-total-spatiotemporal-results.pdf',
       width = 100, height = 100, units = 'mm')

# and the rarefaction visualisation of total spatiotemporal
# ibr: gamma-space-time
gamma_spat_temp <- rip2 %>% 
  unnest(data) %>% 
  group_by(habitat, season, species) %>% 
  summarise(N = sum(N)) %>% 
  nest() %>% 
  mutate(ibr = map(data, ~mobr::rarefaction(.x$N, method = 'IBR'))) %>% 
  unnest(ibr) %>% 
  # still grouped, so
  mutate(N = 1:n())

# all the alphas
alpha_ibr_dat <- rip2 %>% 
  mutate(ibr = map(data, ~mobr::rarefaction(.x$N, method = 'IBR'))) %>% 
  unnest(ibr) %>% 
  group_by(site_code, habitat, year, season) %>% 
  mutate(N = 1:n())

# need a target N for the comparison
# focus on perennial habitat in spring
targetN <- alpha_ibr_dat %>% 
  select(-data) %>% 
  filter(N == max(N)) %>% 
  ungroup() %>% 
  group_by(habitat, season) %>%
  summarise(targetN = min(N))

targetN_perennial_spring <- targetN %>% 
  ungroup() %>% 
  filter(targetN == min(targetN)) %>% 
  pull(targetN)

average_alpha_ibr_dat <- alpha_ibr_dat %>% 
  select(-data) %>% 
  slice(1:targetN_perennial_spring) %>% 
  group_by(habitat, season, N) %>% 
  summarise(ibr_hat = mean(ibr),
            ibr_Q95 = quantile(ibr, probs = 0.95),
            ibr_Q05 = quantile(ibr, probs = 0.05))

rarefaction_zoom_inset <- ggplot() +
  geom_ribbon(data = average_alpha_ibr_dat,
              aes(x = N, ymin = ibr_Q05, ymax = ibr_Q95, fill = habitat),
              alpha = 0.33) +
  geom_line(data = average_alpha_ibr_dat,
            aes(x = N, y = ibr_hat, linetype = 'alpha', colour = habitat),
            size = 1.5) +
  geom_line(data = gamma_spat_temp %>% 
              filter(N < 50),
            aes(x = N, y = ibr, linetype = 'gamma', colour = habitat),
            size = 1.5) +
  scale_colour_manual(name = 'Habitat',
                      values = c('perennial-engineered' = '#a6611a',
                                 'perennial-natural' = '#80cdc1'),
                      guide = 'none')  +
  scale_fill_manual(name = 'Habitat',
                    values = c('perennial-engineered' = '#a6611a',
                               'perennial-natural' = '#80cdc1'),
                    guide = 'none')  +
  scale_linetype_manual(name = 'Scale',
                        values = c('alpha' = 2,
                                   'gamma' = 1),
                        label = c(expression(alpha),
                                  expression(gamma)), ) +
  labs(x = '',
       y = '') +
  theme_classic() +
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_blank(),
        legend.direction = 'horizontal') +
  guides(linetype = guide_legend(title.position = 'top',reverse=TRUE))


spatiotemporal_gamma_ibr <-
ggplot() +
  geom_line(data = gamma_spat_temp ,
            aes(x = N, y = ibr, colour = habitat),
            size = 1.5) +
  scale_colour_manual(name = 'Habitat',
                      values = c('perennial-engineered' = '#a6611a',
                                 'perennial-natural' = '#80cdc1'),
                      labels = c('engineered', 'natural'))  +
  labs(x = 'Number of individuals',
       y = 'Expected number of species') +
  theme_bw() +
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_blank()) +
  guides(colour = guide_legend(reverse=TRUE))

spatiotemporal_gamma_ibr + 
  annotation_custom(ggplotGrob(rarefaction_zoom_inset), 
                    xmin = 3000, xmax = 7750, ymin = 1, ymax = 100)
ggsave('~/Dropbox/MoBD (Measurements of Beta diversity)/Beta_paper/figs/case-study-total-spatiotemporal-rarefaction.pdf',
       width = 200, height = 150, units = 'mm')

# spatial variation in total temporal beta-diversity
# leave-one-out (jackknife) resampling to visualise variation
temporal_beta_calcs <- rip2 %>% 
  unnest(data) %>% 
  group_by(site_code, habitat, season) %>% 
  nest(data = c(year, code, N)) %>% 
  mutate(wide_data = map(data, ~pivot_wider(data = ., 
                                            names_from = code, 
                                            values_from = N,
                                            values_fill = 0))) %>% 
  ungroup()

temporal_targetC <- temporal_beta_calcs %>% 
  mutate(target_C = map(wide_data, ~mobr::C_target(x = .[,-1], 
                                                   factor = 2))) %>%
  unnest(target_C) %>% 
  ungroup() %>% 
  summarise(C_target = min(target_C)) %>% 
  ungroup()

temporal_targetN <- rip2 %>% 
  filter(season=='spring' & habitat %in% c('perennial-engineered', 
                                           'perennial-natural')) %>%
  unnest(data) %>% 
  group_by(site_code, habitat, season, year) %>% 
  summarise(N = sum(N)) %>% 
  ungroup() %>% 
  filter(N == min(N)) %>% 
  pull(N)
  

temporal_beta_calcs <- temporal_beta_calcs %>% 
  mutate(beta_C = map(wide_data, ~mobr::beta_C(x = .[,-1], 
                                               C = temporal_targetC$C_target, 
                                               extrapolation = TRUE)),
         beta_S = map(wide_data, ~mobr::calc_comm_div(abund_mat = .[,-1], 
                                                      index = 'S', 
                                                      coverage = FALSE, 
                                                      scales = 'beta')),
         beta_S_PIE = map(wide_data, ~mobr::calc_comm_div(abund_mat = .[,-1], 
                                                          index = 'S_PIE', 
                                                          coverage = FALSE, 
                                                          scales = 'beta')),
         beta_S_n = map(wide_data, ~mobr::calc_comm_div(abund_mat = .[,-1], 
                                                          index = 'S_n', 
                                                          coverage = FALSE, 
                                                          effort = temporal_targetN,
                                                          scales = 'beta')))

temp_beta_dat <- temporal_beta_calcs %>% 
  unnest(cols = beta_S) %>% 
  rename(beta_S = value) %>% 
  dplyr::select(-c(scale, index, sample_size, effort, coverage)) %>% 
  unnest(cols = beta_S_PIE) %>% 
  rename(beta_S_PIE = value) %>% 
  dplyr::select(-c(scale, index, sample_size, effort, coverage)) %>% 
  unnest(cols = beta_S_n) %>% 
  rename(beta_S_n = value) %>% 
  dplyr::select(-c(scale, index, sample_size, effort, coverage)) %>% 
  unnest(cols = beta_C) %>% 
  dplyr::select(-c(data, wide_data)) %>% 
  pivot_longer(cols = beta_C:beta_S_n, names_to = 'index', 
               values_to = 'value') %>% 
  mutate(index = factor(index, levels = c('beta_S', 'beta_S_PIE', 
                                          'beta_S_n', 'beta_C'))) 

temp_beta_dat %>% 
  ggplot() +
  facet_grid(habitat ~ season) +
  geom_point(aes(x = index, y = value, colour = site_code)) +
  geom_hline(yintercept = 1, lty = 2)


spatial_beta_dat %>% 
  ggplot() +
  facet_grid(habitat~season) +
  stat_smooth(aes(x = year, y = value, colour = index),
              method = 'gam') +
  geom_point(aes(x = year, y = value, colour = index))

# contrast beta-diversity
contrast_perennial <- rip2 %>% 
  ungroup() %>% 
  filter(habitat %in% c('perennial-engineered', 'perennial-natural') &
           season == 'spring') %>% 
  unnest(data) %>% 
  # first, pool the two sites so as the maximum contrast diversity == 2
  group_by(year, habitat, season, code) %>% 
  summarise(N = sum(N)) %>% 
  nest(data = c(habitat, code, N)) %>% 
  mutate(wide_data = map(data, ~pivot_wider(data = .,
                                            names_from = code,
                                            values_from = N,
                                            values_fill = 0))) %>% 
  ungroup()

targetC_contr_perennial <- contrast_perennial %>% 
  mutate(target_C = map(wide_data, ~mobr::C_target(x = .[,-1], factor = 2))) %>% 
  unnest(target_C) %>% 
  ungroup() %>% 
  summarise(C_target = min(target_C)) %>% 
  pull(C_target)

targetN_contr_perennial <- rip2 %>% 
  ungroup() %>% 
  filter(habitat %in% c('perennial-engineered', 'perennial-natural') &
           season=='spring') %>% 
  unnest(data) %>% 
  # first, pool the two sites so as the maximum contrast diversity == 2
  group_by(year, habitat, season) %>% 
  summarise(N = sum(N)) %>% 
  ungroup() %>% 
  filter(N == min(N)) %>% 
  pull(N)
  
contr_beta_calcs_perennial <- contrast_perennial %>% 
  mutate(beta_C = map(wide_data, possibly(~mobr::beta_C(x = .[,-1], 
                                                        C = targetC_contr_perennial, 
                                                        extrapolation = TRUE), 
                                          otherwise = NULL)),
         beta_S = map(wide_data, ~mobr::calc_comm_div(abund_mat = .[,-1], 
                                                      index = 'S', 
                                                      scales = 'beta', 
                                                      coverage = FALSE)),
         beta_S_PIE = map(wide_data, ~mobr::calc_comm_div(abund_mat = .[,-1], 
                                                          index = 'S_PIE', 
                                                          scales = 'beta', 
                                                          coverage = FALSE)),
         beta_S_n = map(wide_data, ~mobr::calc_comm_div(abund_mat = .[,-1], 
                                                          index = 'S_n',
                                                          effort = targetN_contr_perennial,
                                                          scales = 'beta', 
                                                          coverage = FALSE)))


contr_beta_dat <- contr_beta_calcs_perennial %>%
                              unnest(cols = beta_S) %>%
                              rename(beta_S = value) %>%
                              dplyr::select(-c(scale, index, sample_size, effort, coverage)) %>%
                              unnest(cols = beta_S_PIE) %>%
                              rename(beta_S_PIE = value) %>%
                              dplyr::select(-c(scale, index, sample_size, effort, coverage)) %>%
                              unnest(cols = beta_S_n) %>%
                              rename(beta_S_n = value) %>%
                              dplyr::select(-c(scale, index, sample_size, effort, coverage)) %>%
                              unnest(cols = beta_C) %>%
                              dplyr::select(-c(data, wide_data)) %>%
                              pivot_longer(cols = beta_C:beta_S_n, names_to = 'index', values_to = 'value') %>%
                              mutate(index = factor(index, levels = c('beta_S', 'beta_S_PIE', 
                                                                      'beta_S_n', 'beta_C')),
                                     contrast = 'perennial')

contr_beta_dat %>% 
  ggplot() +
  facet_grid(contrast~season) +
  stat_smooth(aes(x = year, y = value, colour = index),
              method = 'gam') +
  geom_point(aes(x = year, y = value, colour = index))


# present total spatial-temporal beta-diversity using the rarefaction curves?
# doesn't work, the gamma-scale swamps the alpha


# calculations for total beta-diversity
total_beta <- rip2 %>% 
  ungroup() %>% 
  filter(season=='spring' & habitat %in% c('perennial-engineered', 
                                           'perennial-natural')) %>% 
  unnest(data) %>% 
  dplyr::select(-season) %>% 
  group_by(habitat) %>% 
  nest(data = c(site_code, year, code, N)) %>% 
  mutate(wide_data = map(data, ~pivot_wider(data = .,
                                            names_from = code,
                                            values_from = N,
                                            values_fill = 0))) %>% 
  ungroup()

total_targetC <- total_beta %>% 
  mutate(target_C = map(wide_data, ~mobr::C_target(x = .[,-c(1,2)], factor = 2))) %>% 
  unnest(target_C) %>% 
  ungroup() %>% 
  summarise(C_target = min(target_C))

total_beta_calcs <- total_beta %>% 
  mutate(beta_C = map(wide_data, possibly(~mobr::beta_C(x = .[,-c(1,2)], 
                                                        C = total_targetC$C_target, 
                                                        extrapolation = TRUE), 
                                          otherwise = NULL)),
         beta_S = map(wide_data, ~mobr::calc_comm_div(abund_mat = .[,-c(1,2)], 
                                                      index = 'S', 
                                                      scales = 'beta', 
                                                      coverage = FALSE)),
         beta_S_PIE = map(wide_data, ~mobr::calc_comm_div(abund_mat = .[,-c(1,2)], 
                                                          index = 'S_PIE', 
                                                          scales = 'beta', 
                                                          coverage = FALSE)),
         beta_S_n = map(wide_data, ~mobr::calc_comm_div(abund_mat = .[,-c(1,2)], 
                                                          index = 'S_n', 
                                                          scales = 'beta',
                                                        effort = targetN_perennial_spring,
                                                          coverage = FALSE)))

total_beta_results <- total_beta_calcs %>% 
  unnest(cols = beta_S) %>% 
  rename(beta_S = value) %>% 
  dplyr::select(-c(scale, index, sample_size, effort, coverage)) %>% 
  unnest(cols = beta_S_PIE) %>% 
  rename(beta_S_PIE = value) %>% 
  dplyr::select(-c(scale, index, sample_size, effort, coverage)) %>% 
  unnest(cols = beta_S_n) %>% 
  rename(beta_S_n = value) %>% 
  dplyr::select(-c(scale, index, sample_size, effort, coverage)) %>% 
  unnest(cols = beta_C) %>% 
  dplyr::select(-c(data, wide_data)) %>% 
  pivot_longer(cols = beta_C:beta_S_n, names_to = 'index', 
               values_to = 'value') %>% 
  mutate(index = factor(index, levels = c('beta_S', 'beta_S_PIE', 
                                          'beta_S_n', 'beta_C'))) 

# shorten labels for facets
hab = c('perennial-engineered' = 'engineered',
        'perennial-natural' = 'natural')



panelA <- temp_beta_dat %>% 
  filter(season=='spring' & habitat %in% c('perennial-engineered', 
                                           'perennial-natural')) %>%
  ggplot() +
  facet_grid(habitat ~ season) +
  geom_point(aes(x = index, y = value, group = site_code, colour = index)) +
  geom_hline(yintercept = 1, lty = 2) +
  scale_x_discrete(limits = c('beta_S', 'beta_S_PIE', 'beta_S_n', 'beta_C'),
                   labels = c(expression(beta[S]), expression(beta[S[PIE]]), 
                              expression(beta[S[n]]), expression(beta[C])),
                   name = 'Metric') +
  scale_colour_manual(values = c('beta_S' = '#d01c8b',
                                 'beta_S_PIE' = '#f1b6da',
                                 'beta_S_n' = '#b8e186',
                                 'beta_C' = '#4dac26'),
                      labels = c(expression(beta[S]), expression(beta[S[PIE]]), 
                                 expression(beta[S[n]]), expression(beta[C])),
                      name = 'Metric') +
  labs(tag = 'A',
      subtitle = expression(paste('Spatial variation in temporal ', beta,
                                  '-diversity'))) +
  theme_bw() +
  theme(legend.position = 'none')

# panelB <- 
spatial_beta_dat %>% 
  filter(index != 'beta_S_n') %>% 
  ggplot() +
  facet_grid(~index) +
  stat_smooth(aes(x = year, y = value, colour = index, linetype = habitat),
              method = 'gam') +
  geom_point(aes(x = year, y = value, colour = index)) +
  scale_colour_manual(values = c('beta_S' = '#d01c8b',
                                 'beta_S_PIE' = '#f1b6da',
                                 'beta_S_n' = '#b8e186',
                                 'beta_C' = '#4dac26'),
                      labels = c(expression(beta[S]), expression(beta[S[PIE]]), 
                                 expression(beta[S[n]]), expression(beta[C])),
                      name = '') +
  labs(tag = 'B',
       subtitle = expression(paste('Temporal variation in spatial ', 
                                   beta, '-diversity'))) +
  theme_bw() +
  theme(legend.position = c(0,0),
        legend.justification = c(0,0),
        legend.background = element_blank()) +
  guides(colour = guide_legend(nrow = 1))

panelC <- contr_beta_dat %>% 
  ggplot() +
  facet_grid(contrast~season) +
  stat_smooth(aes(x = year, y = value, colour = index),
              method = 'gam') +
  geom_point(aes(x = year, y = value, colour = index)) +
  scale_colour_manual(values = c('beta_S' = '#d01c8b',
                                 'beta_S_PIE' = '#f1b6da',
                                 'beta_S_n' = '#b8e186',
                                 'beta_C' = '#4dac26'),
                      labels = c(expression(beta[S]), expression(beta[S[PIE]]), 
                                 expression(beta[S[n]]), expression(beta[C])),
                      name = '') +
  labs(tag = 'C',
       subtitle = expression(paste('Temporal variation in natural vs engineered ', 
                                   beta, '-diversity'))) +
  theme_bw() +
  theme(legend.position = c(1,0),
        legend.justification = c(1,0),
        legend.background = element_blank()) +
  guides(colour = guide_legend(nrow = 1))

panelD <- cowplot::plot_grid(total_rarefaction_fig +
                               theme(plot.margin = unit(c(2,0,2,2), 'mm')),
                             cowplot::plot_grid(NULL,
                                                total_beta_metrics +
                                                  theme(plot.margin = unit(c(0,0,0,0), 'mm'))                   ,
                                                NULL,
                                                nrow = 3,
                                                rel_heights = c(0.2,1,0.2)),
                             ncol = 2,
                             rel_widths = c(1, 0.6),
                             align = 'hv', axis = 'lrtb' )

cowplot::plot_grid(panelA,
                   panelB,
                   panelC,
                   panelD)

ggsave('~/Dropbox/MoBD (Measurements of Beta diversity)/Beta_paper/figs/bateman-draft.png',
       width = 290, height = 250, units = 'mm')

# temporal beta-diversity  rarefaction figure
gamma_time <- rip2 %>% 
  unnest(data) %>% 
  group_by(site_code, habitat, season, code) %>% 
  summarise(N = sum(N)) %>% 
  nest() %>% 
  mutate(ibr = map(data, ~mobr::rarefaction(.x$N, method = 'IBR'))) %>% 
  unnest(ibr) %>% 
  # still grouped, so
  mutate(N = 1:n())

# need a target N for the comparison
# focus on perennial habitat in spring
targetN_alphaTime <- alpha_ibr_dat %>% 
  select(-data) %>% 
  filter(N == max(N)) %>% 
  ungroup() %>% 
  group_by(site_code, habitat, season) %>%
  summarise(targetN = min(N))

targetN_perennial_spring <- targetN_alphaTime %>% 
  ungroup() %>% 
  filter(season=='spring' & habitat %in% c('perennial-engineered', 
                                           'perennial-natural')) %>% 
  filter(targetN == min(targetN)) %>% 
  pull(targetN)


average_alpha_time <- alpha_ibr_dat %>% 
  select(-data) %>% 
  filter(season=='spring' & habitat %in% c('perennial-engineered', 
                                           'perennial-natural')) %>%
  slice(1:targetN_perennial_spring) %>% 
  group_by(site_code, habitat, season, N) %>% 
  summarise(ibr_hat = mean(ibr),
            ibr_Q95 = quantile(ibr, probs = 0.95),
            ibr_Q05 = quantile(ibr, probs = 0.05))


temp_beta_rarefaction <- ggplot() +
  facet_grid(habitat ~ season) + 
  geom_ribbon(data = average_alpha_time,
              aes(x = N, ymin = ibr_Q05, ymax = ibr_Q95, 
                  group = interaction(site_code, habitat, season),
                  fill = site_code, colour = NULL),
              alpha = 0.33) +
  geom_line(data = average_alpha_time,
            aes(x = N, y = ibr_hat, linetype = 'alpha', 
                colour = site_code),
            size = 1.5) +
  geom_line(data = gamma_time %>% 
              filter(N < 50) %>% 
              filter(season=='spring' & habitat %in% c('perennial-engineered', 
                                                       'perennial-natural')),
            aes(x = N, y = ibr, linetype = 'gamma', colour = site_code),
            size = 1.5) +
  scale_fill_manual(name = 'Site',
                      values = c('PE-11A' = '#a6611a',
                                 'PE-1D' = '#dfc27d',
                                 'PN-1B' = '#80cdc1',
                                 'PN-7A' = '#018571'),
                    guide = 'none') +
  scale_colour_manual(name = 'Site',
                    values = c('PE-11A' = '#a6611a',
                               'PE-1D' = '#dfc27d',
                               'PN-1B' = '#80cdc1',
                               'PN-7A' = '#018571'),
                    guide = 'none') +
  scale_linetype_manual(name = 'Scale',
                        values = c('alpha' = 2,
                                   'gamma' = 1),
                        label = c(expression(alpha),
                                  expression(gamma)), ) +
  labs(x = 'Number of individuals',
       y = 'Expected number of species',
       subtitle = expression(paste('Total temporal ', beta,
                                   ' - diversity')),
       tag = 'A') +
  theme_bw() +
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.background = element_blank(),
        legend.direction = 'horizontal',
        plot.margin = unit(c(2,0,2,2), units = 'mm')) +
  guides(linetype = guide_legend(title.position = 'top',
                                 reverse=TRUE))#,
         # fill = guide_legend(title.position = 'top'),
         # colour = guide_legend(title.position = 'top'))

temp_beta_metrics <- ggplot() +
  # facet_wrap(~habitat, ncol = 1, labeller = labeller(habitat = hab)) +
  geom_point(data = temp_beta_dat,
             aes(x = index, y = value, colour = site_code, shape = habitat, 
                 group = habitat),
             position = position_dodge(width = 0.1), size = 1.5) +
  scale_x_discrete(limits = c('beta_S', 'beta_S_PIE', 'beta_S_n', 'beta_C'),
                   labels = c(expression(beta[S]), expression(beta[S[PIE]]), 
                              expression(beta[S[n]]), expression(beta[C])),
                   name = 'Metric') +
  scale_colour_manual(name = 'Site',
                      values = c('PE-11A' = '#a6611a',
                                 'PE-1D' = '#dfc27d',
                                 'PN-1B' = '#80cdc1',
                                 'PN-7A' = '#018571'),
                      guide = 'none')  +
  scale_shape_manual(name = '',
                     values = c('perennial-engineered' = 19,
                                 'perennial-natural' = 17),
                     label = c('engineered',
                               'natural'),
                     guide = 'none')  +
  theme_bw() +
  theme(legend.position = 'top',
        legend.direction = 'vertical',
        legend.justification = c(1,1),
        plot.margin = unit(x = c(0,0,0,0), units = 'mm'))


panelA_alt <- cowplot::plot_grid(temp_beta_rarefaction,
                             cowplot::plot_grid(NULL,
                                                temp_beta_metrics,
                                                NULL,
                                                nrow = 3, 
                                                rel_heights = c(0.2,1,0.2)),
                             ncol = 2,
                             rel_widths = c(1, 0.6),
                             align = 'hv', axis = 'lrtb' )

cowplot::plot_grid(panelA_alt,
                   panelB,
                   panelD,
                   panelC)

ggsave('~/Dropbox/MoBD (Measurements of Beta diversity)/Beta_paper/figs/bateman-draft2.png',
       width = 290, height = 200, units = 'mm')
