## temporal variation in spatial beta-diversity
library(mobr)
# get standardised data
source('./scripts/bateman-beta-concept-perennial-spring-standardise-effort.R')

spatial_beta_calcs <- rip2 %>% 
  unnest(data) %>% 
  group_by(habitat, year, season) %>% 
  nest(data = c(site_code, species, N)) %>% 
  mutate(wide_data = map(data, ~pivot_wider(data = ., 
                                            names_from = species, 
                                            values_from = N,
                                            values_fill = 0))) %>% 
  ungroup()

spatial_targetC <- spatial_beta_calcs %>% 
  mutate(target_C = map(wide_data, ~mobr::calc_C_target(x = .[,-c(1)], 
                                                   factor = 2))) %>%
  unnest(target_C) %>% 
  ungroup() %>% 
  summarise(C_target = min(target_C)) %>% 
  ungroup()

spatial_targetN <- rip2 %>% 
  filter(season=='spring' & habitat %in% c('perennial-engineered', 
                                           'perennial-natural')) %>%
  unnest(data) %>% 
  group_by(site_code, habitat, season, year) %>% 
  summarise(N = sum(N)) %>% 
  ungroup() %>% 
  filter(N == min(N)) %>% 
  pull(N)


spatial_beta_calcs <- spatial_beta_calcs %>% 
  mutate(beta_S_C = map(wide_data, ~mobr::calc_beta_div(abund_mat = .[ , -1], index = 'S_C', 
                                                      C_target_gamma = spatial_targetC$C_target)),
         beta_S = map(wide_data, ~mobr::calc_beta_div(abund_mat = .[ , -1], 
                                                      index = 'S')),
         beta_S_PIE = map(wide_data, ~mobr::calc_beta_div(abund_mat = .[ , -1], 
                                                          index = 'S_PIE')),
         beta_S_n = map(wide_data, ~mobr::calc_beta_div(abund_mat = .[ , -1], 
                                                        index = 'S_n', 
                                                        effort = spatial_targetN)))


spatial_beta_dat <- spatial_beta_calcs %>% 
  unnest(cols = beta_S) %>% 
  rename(beta_S = value) %>% 
  dplyr::select(-c(scale, index, sample_size, effort, gamma_coverage)) %>% 
  unnest(cols = beta_S_PIE) %>% 
  rename(beta_S_PIE = value) %>% 
  dplyr::select(-c(scale, index, sample_size, effort, gamma_coverage)) %>% 
  unnest(cols = beta_S_n) %>% 
  rename(beta_S_n = value) %>% 
  dplyr::select(-c(scale, index, sample_size, effort, gamma_coverage)) %>% 
  unnest(cols = beta_S_C) %>%
  rename(beta_S_C = value) %>% 
  dplyr::select(-c(scale, index, sample_size, effort, gamma_coverage)) %>% 
  pivot_longer(cols = beta_S_C:beta_S_n, names_to = 'index', 
               values_to = 'value') %>% 
  mutate(index = factor(index, levels = c('beta_S',  
                                          'beta_S_n', 
                                          'beta_S_C', 
                                          'beta_S_PIE'))) 

spatial_beta_dat$index <- factor(spatial_beta_dat$index,
                                 labels = c(expression(beta[S]),
                                            expression(beta[S[n]]),
                                            expression(beta[C]),
                                            expression(beta[S[PIE]])))
panel_B <- spatial_beta_dat %>% 
  filter(index %in% c('beta[S]', 'beta[C]')) %>% 
  ggplot() +
  facet_grid(~index, labeller = label_parsed) +
  stat_smooth(aes(x = year, y = value, colour = habitat, linetype = habitat,
                  fill = habitat),
              method = 'gam') +
  geom_point(aes(x = year, y = value, colour = habitat)) +
  scale_colour_manual(name = 'Habitat',
                      values = c('perennial-engineered' = '#a6611a',
                                 'perennial-natural' = '#80cdc1'),
                      label = c('engineered',
                                'natural'),
                      # guide = 'none'
  )  +
  scale_linetype_manual(name = 'Habitat',
                        values = c('perennial-engineered' = 1,
                                   'perennial-natural' = 1),
                        label = c('engineered',
                                  'natural')) +
  scale_fill_manual(name = 'Habitat',
                    values = c('perennial-engineered' = '#a6611a',
                               'perennial-natural' = '#80cdc1'),
                        label = c('engineered',
                                  'natural')) +
  labs(x = 'Year',
       y = 'Metric value') +
  scale_x_continuous(breaks = c(2004, 2008, 2012, 2016),
                     labels = c(2004, '', 2012, '')) +
  theme_bw() +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_blank()) +
  guides(colour = guide_legend(nrow = 2, title = NULL, reverse = TRUE),
         linetype = guide_legend(nrow = 2, title = NULL, reverse = TRUE),
         fill = guide_legend(nrow = 2, title = NULL, reverse = TRUE))

panel_B
ggsave('./figs/case-study-total-spatial-beta-results.pdf',
       width = 150, height = 90, units = 'mm')
