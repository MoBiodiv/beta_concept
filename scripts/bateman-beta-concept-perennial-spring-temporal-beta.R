# leave-one-out (jackknife) for temporal beta-diversity

# get standardised data
source('~/Dropbox/MoBD (Measurements of Beta diversity)/Beta_paper/bateman-beta-concept-perennial-spring-standardise-effort.R')

counter <- rip2 %>% 
  distinct(site_code, year)

temporal_beta_jk <- tibble()

for(i in 1:nrow(counter)){
  print(paste(i, ' of ', nrow(counter), ' resamples'))
  temp <- rip2 %>% 
    # drop one year from each site
    group_by(site_code) %>% 
    arrange(-desc(year)) %>% 
    slice(-i) %>% 
    mutate(resamp = i) %>% 
    ungroup() %>% 
    unnest(data) %>% 
    nest(data = c(year, species, N)) %>% 
    mutate(wide_data = map(data, ~pivot_wider(data = .,
                                              names_from = species,
                                              values_from = N,
                                              values_fill = 0))) %>% 
    ungroup()
  
  temporal_beta_jk <- bind_rows(temporal_beta_jk, 
                             temp)
}

# the maximum number of distinct communities in the jackknife data is: 13
# it is equal to the number of years sampled at site
temporal_beta_jk$wide_data[[1]] %>% nrow()

# calculate target coverage for standardisation
temporal_targetC_resamps <- temporal_beta_jk %>% 
  mutate(target_C = map(wide_data, ~mobr::C_target(x = .[,-c(1)], factor = 2))) %>% 
  unnest(target_C) %>% 
  ungroup() %>% 
  summarise(C_target = min(target_C))

# calculate number of individuals for beta_S_n
temporal_targetN_resamps <- temporal_beta_jk %>% 
  unnest(data) %>% 
  group_by(habitat, season, site_code, year, resamp) %>% 
  summarise(J = sum(N)) %>% 
  ungroup() %>% 
  filter(J == min(J)) %>% 
  pull(J) %>% 
  unique()

temporal_beta_resamps_calcs <- temporal_beta_jk %>% 
  mutate(beta_C = map(wide_data, possibly(~mobr::beta_C(x = .[,-c(1)], 
                                                        C = temporal_targetC_resamps$C_target, 
                                                        extrapolation = TRUE), 
                                          otherwise = NULL)),
         beta_S = map(wide_data, ~mobr::calc_comm_div(abund_mat = .[,-c(1)], 
                                                      index = 'S', 
                                                      scales = 'beta', 
                                                      coverage = FALSE)),
         beta_S_PIE = map(wide_data, ~mobr::calc_comm_div(abund_mat = .[,-c(1)], 
                                                          index = 'S_PIE', 
                                                          scales = 'beta', 
                                                          coverage = FALSE)),
         beta_S_n = map(wide_data, ~mobr::calc_comm_div(abund_mat = .[,-c(1)], 
                                                        index = 'S_n', 
                                                        effort = temporal_targetN_resamps,
                                                        scales = 'beta', 
                                                        coverage = FALSE)))

temporal_beta_resamps <- temporal_beta_resamps_calcs %>% 
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
temporal_beta_summary <- temporal_beta_resamps %>% 
  group_by(habitat, site_code, metric) %>% 
  summarise(location = mean(value),
            Q5 = quantile(value, probs = c(0.05)),
            Q95 = quantile(value, probs = c(0.95))) %>% 
  ungroup()

temporal_beta_summary2 <- temporal_beta_resamps %>% 
  group_by(habitat, metric) %>% 
  summarise(location = mean(value),
            Q5 = quantile(value, probs = c(0.05)),
            Q95 = quantile(value, probs = c(0.95))) %>% 
  ungroup() %>% 
  separate(habitat, c('water_avail', 'habitat'))

temporal_beta_results <-
ggplot() +
  geom_point(data = temporal_beta_summary2 %>% 
               filter(metric %in% c('beta_S', 'beta_C')),
             aes(x = metric, y = location, colour = habitat, shape = habitat,
                 group = habitat),
             size = 1.5, alpha = 1,
             position = position_dodge(width = 0.25)) +
  geom_text(data = temporal_beta_summary2 %>% 
              filter(metric %in% c('beta_S')),
            aes(x = metric, y = location, colour = habitat, 
                label = habitat),
            nudge_x = 0.4, hjust = 0,
  ) +
  geom_linerange(data = temporal_beta_summary2 %>% 
                   filter(metric %in% c('beta_S', 'beta_C')),
                 aes(x = metric, ymin = Q5, ymax = Q95,
                     colour = habitat, group = habitat),
                 position = position_dodge(width = 0.25)) +
  scale_x_discrete(limits = c('beta_S','beta_C'),#, 'beta_S_PIE'
                   labels = c(expression(beta[S]), 
                              expression(beta[C])
                              # expression(beta[S[PIE]])
                              ),
                   name = 'Metric') +
  scale_colour_manual(name = 'Habitat',
                      values = c('engineered' = '#a6611a',
                                 'natural' = '#80cdc1'),
                      label = c('engineered',
                                'natural'),
                      # guide = 'none'
  )  +
  scale_shape_manual(name = 'Habitat',
                     values = c('engineered' = 19,
                                'natural' = 17),
                     label = c('engineered',
                               'natural'),
                     # guide = 'none'
  )  +
  labs(y = 'Metric value',
       tag = 'ii.') +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'vertical',
        legend.justification = c(1,1),
        legend.background = element_blank(),
        plot.tag = element_text(face = 'bold', hjust = 0, vjust = 0),
        plot.tag.position = c(0.2, 0.95)) +
  guides(colour = guide_legend(reverse = TRUE, title = NULL),
         shape = guide_legend(reverse = TRUE, title = NULL))

# ggsave('~/Dropbox/MoBD (Measurements of Beta diversity)/Beta_paper/figs/case-study-total-temporal-results.pdf',
#        width = 100, height = 100, units = 'mm')

# and the rarefaction visualisation
gamma_time <- rip2 %>% 
  unnest(data) %>% 
  group_by(site_code, habitat, season, species) %>% 
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
  summarise(targetN = min(N)) %>% 
  ungroup() %>% 
  filter(targetN == min(targetN)) %>% 
  pull(targetN)


average_alpha_time <- alpha_ibr_dat %>% 
  select(-data) %>% 
  slice(1:targetN_alphaTime) %>% 
  group_by(site_code, habitat, season, N) %>% 
  summarise(ibr_hat = median(ibr),
            ibr_Q95 = quantile(ibr, probs = 0.95),
            ibr_Q05 = quantile(ibr, probs = 0.05))

average_alpha_time$habitat <- factor(average_alpha_time$habitat,
                                     levels = c('perennial-natural',
                                                'perennial-engineered'))
gamma_time$habitat <- factor(gamma_time$habitat,
                                     levels = c('perennial-natural',
                                                'perennial-engineered'))

temp_beta_rarefaction_inset <-
ggplot() +
  facet_wrap(habitat ~ site_code, scales = 'free') +
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
                    values = c('PE-10B' = '#8c510a',
                               'PE-11A' = '#d8b365',
                               'PE-1D' = '#f6e8c3',
                               'PN-1B' = '#c7eae5',
                               'PN-2A' = '#5ab4ac',
                               'PN-7A' = '#01665e'),
                    guide = 'none') +
  scale_colour_manual(name = 'Site',
                      values = c('PE-10B' = '#8c510a',
                                 'PE-11A' = '#d8b365',
                                 'PE-1D' = '#f6e8c3',
                                 'PN-1B' = '#c7eae5',
                                 'PN-2A' = '#5ab4ac',
                                 'PN-7A' = '#01665e'),
                      guide = 'none') +
  scale_linetype_manual(name = 'Scale',
                        values = c('alpha' = 2,
                                   'gamma' = 1),
                        label = c(expression(alpha),
                                  expression(gamma))) +
  scale_x_continuous(limits = c(0, 50), breaks = c(0, 25, 50)) +
  scale_y_continuous(limits = c(0, 25), breaks = c(0,5,10,15,20,25), labels = c(0,'',10,'', 20, '')) +
  labs(x = 'Number of individuals',
       y = 'Expected number of species') +
  theme_classic() +
  theme(legend.position = 'top',
        legend.key.width = unit(10, units = 'mm'),
        legend.justification = c(0,1),
        legend.background = element_blank(),
        legend.direction = 'horizontal',
        strip.text = element_blank(),
        plot.margin = unit(x = c(0,0,0,0), units = 'mm'),
        strip.background = element_blank(),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7)) +
  guides(linetype = guide_legend(title = NULL, reverse=TRUE, 
                                 label.position = 'left', label.hjust = 0))


temporal_gamma_ibr <-
  ggplot() +
  geom_line(data = gamma_time ,
            aes(x = N, y = ibr, colour = site_code),
            size = 1.5) +
  scale_colour_manual(name = 'Sites',   
                      values = c('PE-10B' = '#8c510a',
                                 'PE-11A' = '#d8b365',
                                 'PE-1D' = '#f6e8c3',
                                 'PN-1B' = '#c7eae5',
                                 'PN-2A' = '#5ab4ac',
                                 'PN-7A' = '#01665e'),
                      labels = c('engineered','engineered', 'engineered',
                                 'natural', 'natural', 'natural'),
                      guide = guide_legend(reverse = TRUE,  direction = "horizontal",
                                           title.position = "top",
                                           label.position = "bottom",
                                           label.hjust = 0,
                                           label.vjust = 1,
                                           byrow = TRUE,
                                           label.theme = element_text(angle = 90, size = 7)))  +
  scale_x_continuous(breaks = c(0, 1000, 2000, 3000, 4000),
                     labels = c(0, 1000, '', 3000, '')) +
  labs(x = 'Number of individuals',
       y = 'Expected number of species',
       tag = 'i.') +
  theme_bw() +
  theme(legend.position = c(0.15,0.55),
        legend.justification = c(0,1),
        legend.background = element_blank(),
        legend.key.width = unit(4, units = 'mm'),
        # legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        plot.tag = element_text(face = 'bold', hjust = 0, vjust = 0),
        plot.tag.position = c(0.125, 0.95)) 

temporal_beta_rarefaction <- temporal_gamma_ibr + 
  annotation_custom(ggplotGrob(temp_beta_rarefaction_inset), 
                    xmin = 1375, xmax = 3900, ymin = -2, ymax = 70)

panel_C <- cowplot::plot_grid(temporal_beta_rarefaction, 
                              temporal_beta_results,
                              nrow = 1, rel_widths = c(1, 0.5))
