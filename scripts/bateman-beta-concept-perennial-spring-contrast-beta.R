


# contrast beta-diversity
contrast_perennial <- rip2 %>% 
  ungroup() %>% 
  unnest(data) %>% 
  # first, pool the two sites so as the maximum contrast diversity == 2
  group_by(year, habitat, season, species) %>% 
  summarise(N = sum(N)) %>% 
  nest(data = c(habitat, species, N)) %>% 
  mutate(wide_data = map(data, ~pivot_wider(data = .,
                                            names_from = species,
                                            values_from = N,
                                            values_fill = 0))) %>% 
  ungroup()

targetC_contr_perennial <- contrast_perennial %>% 
  mutate(target_C = map(wide_data, ~mobr::calc_C_target(x = .[,-1], factor = 2))) %>% 
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
  mutate(beta_S_C = map(wide_data, possibly(~mobr::calc_beta_div(abund_mat = .[,-1], 
                                                                 index = 'S_C',
                                                        C_target_gamma = targetC_contr_perennial), 
                                          otherwise = NULL)),
         beta_S = map(wide_data, ~mobr::calc_beta_div(abund_mat = .[,-1], 
                                                      index = 'S')), 
         beta_S_PIE = map(wide_data, ~mobr::calc_beta_div(abund_mat = .[,-1], 
                                                          index = 'S_PIE')), 
         beta_S_n = map(wide_data, ~mobr::calc_beta_div(abund_mat = .[,-1], 
                                                        index = 'S_n',
                                                        effort = targetN_contr_perennial)))


contr_beta_dat <- contr_beta_calcs_perennial %>%
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
  pivot_longer(cols = beta_S_C:beta_S_n, names_to = 'index', values_to = 'value') %>%
  mutate(index = factor(index, levels = c('beta_S', 'beta_S_PIE', 
                                          'beta_S_n', 'beta_S_C')),
         contrast = 'perennial')


contr_beta_dat$index <- factor(contr_beta_dat$index,
                                 labels = c(expression(beta[S]),
                                            expression(beta[S[PIE]]),
                                            expression(beta[S[n]]),
                                            expression(beta[C])))

panel_D <-
contr_beta_dat %>% 
  filter(index %in% c('beta[S]', 'beta[C]')) %>% 
  ggplot() +
  # facet_grid(contrast~season) +
  stat_smooth(aes(x = year, y = value, colour = index, fill = index),
              method = 'gam') +
  geom_point(aes(x = year, y = value, colour = index)) +
  scale_colour_manual(values = c('beta[S]' = '#d01c8b',
                                 # 'beta[S[n]]' = '#b8e186',
                                 'beta[C]' = '#4dac26'
                                 # 'beta[S[PIE]]' = '#f1b6da'
                                 ),
                      labels = c(expression(beta[S]),
                                 # expression(beta[S[n]]), 
                                 expression(beta[C])),
                                 # expression(beta[S[PIE]])),
                      name = '') +
  scale_fill_manual(values = c('beta[S]' = '#d01c8b',
                               # 'beta[S[n]]' = '#b8e186',
                               'beta[C]' = '#4dac26'
                               # 'beta[S[PIE]]' = '#f1b6da'
                               ),
                    labels = c(expression(beta[S]),
                               # expression(beta[S[n]]), 
                               expression(beta[C])),
                               # expression(beta[S[PIE]])),
                    name = '') +
  labs(x = 'Year', 
       y = 'Metric value') +
  theme_bw() +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.background = element_blank()) +
  guides(colour = guide_legend(nrow = 1),
         fill = guide_legend(nrow = 1))


panel_D
ggsave('./figs/case-study-composition_diff_time_results.pdf',
       width = 120, height = 90, units = 'mm')

       

top <- cowplot::plot_grid(panel_A, panel_B,
                   nrow = 1, labels = c("A", "B"))
                   
bottom <- cowplot::plot_grid(panel_C, panel_D,
                             rel_widths = c(1, 0.5),
                             nrow = 1,
                             labels = c("C", "D"))

cowplot::plot_grid(top, bottom,
                   nrow = 2)

ggsave('./figs/case-study-results.pdf', width = 270, height = 180, units = 'mm')


