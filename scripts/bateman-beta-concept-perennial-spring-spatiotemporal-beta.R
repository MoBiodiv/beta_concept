library(mobr)

source('./scripts/bateman-beta-concept-perennial-spring-standardise-effort.R')

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
  mutate(target_C = map(wide_data, ~mobr::calc_C_target(x = .[ , -c(1,2)], factor = 2))) %>% 
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
  mutate(beta_S_C = map(wide_data, possibly(~mobr::calc_beta_div(abund_mat = .[ , -c(1,2)],
                                                              index = 'S_C', 
                                                              C_target_gamma = total_targetC_resamps$C_target), 
                                          otherwise = NULL)),
         beta_S = map(wide_data, ~mobr::calc_beta_div(abund_mat = .[ , -c(1,2)], 
                                                      index = 'S')), 
         beta_S_PIE = map(wide_data, ~mobr::calc_beta_div(abund_mat = .[,-c(1,2)], 
                                                          index = 'S_PIE')), 
         beta_S_n = map(wide_data, ~mobr::calc_beta_div(abund_mat = .[,-c(1,2)], 
                                                        index = 'S_n', 
                                                        effort = total_targetN_resamps)))

total_beta_resamps <- total_beta_resamps_calcs %>% 
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
  pivot_longer(cols = beta_S_C:beta_S_n, names_to = 'metric', 
               values_to = 'value') %>% 
  mutate(metric = factor(metric, levels = c('beta_S', 'beta_S_PIE', 
                                          'beta_S_n', 'beta_S_C'))) 

# summarise for figure
total_beta_summary <- total_beta_resamps %>% 
  group_by(habitat, metric) %>% 
  summarise(location = mean(value),
            Q5 = quantile(value, probs = c(0.05)),
            Q95 = quantile(value, probs = c(0.95))) %>% 
  ungroup() %>% 
  separate(habitat, c('water_availability', 'habitat'))

results_plot <-
ggplot() +
  geom_point(data = total_beta_summary %>% 
               filter(metric %in% c('beta_S', 'beta_S_C')),
             aes(x = metric, y = location, colour = habitat, shape = habitat,
                 group = habitat),
             size = 2, 
             position = position_dodge(width = 0.1)) +
  geom_text(data = total_beta_summary %>% 
              filter(metric %in% c('beta_S')),
             aes(x = metric, y = location, colour = habitat, 
                 label = habitat),
             nudge_x = 0.2, hjust = 0, size = 3
              ) +
  geom_linerange(data = total_beta_summary %>% 
                   filter(metric %in% c('beta_S', 'beta_S_C')),
                 aes(x = metric, ymin = Q5, ymax = Q95,
                     colour = habitat, group = habitat),
                 position = position_dodge(width = 0.1)) +
  scale_x_discrete(limits = c('beta_S', 'beta_S_C'),#, 'beta_S_PIE'
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
                      guide = 'none'
                      )  +
  scale_shape_manual(name = 'Habitat',
                     values = c('engineered' = 19,
                                'natural' = 17),
                     label = c('engineered',
                               'natural'),
                     guide = 'none'
                     )  +
  labs(y = 'Metric value',
       title = 'b)') +
  theme_bw() +
  theme(legend.position = 'none',
        legend.direction = 'vertical',
        legend.justification = c(1,1),
        legend.background = element_blank(),
        legend.text = element_text(size = 8),
        plot.tag = element_text(face = 'bold', hjust = 0, vjust = 0),
        plot.tag.position = c(0.3, 0.95)) +
  guides(colour = guide_legend(reverse = TRUE, title = NULL),
         shape = guide_legend(reverse = TRUE, title = NULL))


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
  theme(legend.position = 'top',
        legend.key.width = unit(10, units = 'mm'),
        legend.justification = c(0,1),
        legend.background = element_blank(),
        legend.direction = 'horizontal') +
  guides(linetype = guide_legend(title = NULL, reverse=TRUE, 
                                 label.position = 'left', label.hjust = 0))


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
       y = 'Expected number of species',
       title = 'a)') +
  theme_bw() +
  theme(legend.position = c(0.15,0.5),
        legend.justification = c(0,1),
        legend.background = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        plot.tag = element_text(face = 'bold', hjust = 0, vjust = 0),
        plot.tag.position = c(0.15, 0.95)) +
  guides(colour = guide_legend(reverse=TRUE))

rarefaction_curves <- spatiotemporal_gamma_ibr + 
  annotation_custom(ggplotGrob(rarefaction_zoom_inset), 
                    xmin = 3000, xmax = 7750, ymin = 1, ymax = 100)

panel_A <- cowplot::plot_grid(rarefaction_curves, results_plot,
                   nrow = 1,
                   rel_widths = c(1, 0.5), 
                   align = 'hv')
panel_A
ggsave('./figs/case-study-total-spatiotemporal-results.pdf',
       width = 150, height = 90, units = 'mm')

