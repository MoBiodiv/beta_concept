# demonstration that beta-C captures spatial structure
# simulate three regional communities: all with the same number of species, 
# and variation in log abundance of species; two communities with have no spatial structure (i.e.,
# individuals are distributed randomly in space [poisson process]), one with more individuals than the other;
# the third community has within species aggregation (thomas process for the distribution of individuals of each species),
# and same number of individuals as the smaller (spatially random) community.

# devtools::install_github("MoBiodiv/mobsim")
rm(list=ls())
library(mobsim)
library(mobr)
library(tidyverse)
library(cowplot)
# devtools::install_github("thomasp85/patchwork")
library(patchwork)
source('./scripts/gg_legend.R')


sad1 <- sim_sad(s_pool = 50, n_sim = N_pool_big, sad_type = 'lnorm',
                sad_coef = list(meanlog = log(N_pool_big/S_pool), sdlog = sd_log),
                fix_s_sim = TRUE)

sad2 <- sim_sad(s_pool = NULL, n_sim = N_pool_big, sad_type = 'ls',
                sad_coef = list(N = N_pool_big, alpha = 11.5))

length(sad2)

plot(sad1, method = "rank", ylim = c(1,224))
lines(1:length(sad2), sort(sad2, decreasing = TRUE), type = "b", col = 2)


set.seed(199)
S_pool = 50
N_pool_small = 400
N_pool_big = 1000
sd_log = 0.01
sd_log2 = 0.5

panel1 = sim_poisson_community(s_pool = S_pool, n_sim = N_pool_big, sad_type = 'lnorm',
                               sad_coef = list(meanlog = log(N_pool_big/S_pool), sdlog = sd_log),
                               fix_s_sim = TRUE)

panel2 = sim_poisson_community(s_pool = S_pool, n_sim = N_pool_small, sad_type = 'lnorm',
                               sad_coef = list(meanlog = log(N_pool_small/S_pool), sdlog = sd_log),
                               fix_s_sim = TRUE)

panel3 = sim_thomas_community(s_pool = S_pool, n_sim = N_pool_big, sad_type = 'lnorm',
                               sad_coef = list(meanlog = log(N_pool_big/S_pool), sdlog = sd_log),
                              fix_s_sim = TRUE)

panel4 = sim_poisson_community(s_pool = NULL, n_sim = N_pool_big, sad_type = 'ls',
                               sad_coef = list(N = N_pool_big, alpha = 12))

length(table(panel4$census$species))


# sample 4 quadrats (same as Dan's conceptual figure of perfect sample)
sample1 <- sample_quadrats(panel1, n_quadrats = 4, plot = FALSE,
                           quadrat_area = 0.25, method = 'grid', 
                           x0 = 0, y0 = 0, delta_x = 0.5, delta_y = 0.5)

sample2 <- sample_quadrats(panel2, n_quadrats = 4, plot = FALSE,
                           quadrat_area = 0.25, method = 'grid', 
                           x0 = 0, y0 = 0, delta_x = 0.5, delta_y = 0.5)

sample3 <- sample_quadrats(panel3, n_quadrats = 4, plot = FALSE,
                           quadrat_area = 0.25, method = 'grid', 
                           x0 = 0, y0 = 0, delta_x = 0.5, delta_y = 0.5)

sample4 <- sample_quadrats(panel4, n_quadrats = 4, plot = FALSE,
                           quadrat_area = 0.25, method = 'grid', 
                           x0 = 0, y0 = 0, delta_x = 0.5, delta_y = 0.5)

# calculate target coverage (use same target )
ct1 = calc_C_target(sample1$spec_dat, factor = 2)
ct2 = calc_C_target(sample2$spec_dat, factor = 2)
ct3 = calc_C_target(sample3$spec_dat, factor = 2)
ct4 = calc_C_target(sample4$spec_dat, factor = 2)

ctarget = min(ct1, ct2, ct3, ct4)

bC1 = calc_beta_div(sample1$spec_dat, 'S_C', C_target_gamma = ctarget)
bC2 = calc_beta_div(sample2$spec_dat, 'S_C', C_target_gamma = ctarget)
bC3 = calc_beta_div(sample3$spec_dat, 'S_C', C_target_gamma = ctarget)
bC4 = calc_beta_div(sample4$spec_dat, 'S_C', C_target_gamma = ctarget)

# also need target for beta_Sn
Ntarget = min(rowSums(sample1$spec_dat), rowSums(sample2$spec_dat), 
              rowSums(sample3$spec_dat), rowSums(sample4$spec_dat))
# calculate betaS
bS1 <- calc_comm_div(sample1$spec_dat,
              index = c('S', 'S_PIE', 'S_n'),
              scales = 'beta', effort = Ntarget)
bS2 <- calc_comm_div(sample2$spec_dat,
                     index = c('S', 'S_PIE', 'S_n'),
                     scales = 'beta', effort = Ntarget)
bS3 <- calc_comm_div(sample3$spec_dat,
                     index = c('S', 'S_PIE', 'S_n'),
                     scales = 'beta', effort = Ntarget)

bS4 <- calc_comm_div(sample4$spec_dat,
                     index = c('S', 'S_PIE', 'S_n'),
                     scales = 'beta', effort = Ntarget)

panel1_map <- ggplot(panel1$census) +
  geom_point(aes(x = x , y = y, colour = species)) +
  labs(subtitle = paste0('S =  ', S_pool, ', N = ', N_pool_big, '\nrandom')) +
  scale_x_continuous(breaks = c(0, 0.5, 1), expand = c(0,0)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0,0)) +
  scale_color_viridis_d() +
  labs(tag = 'A') +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size=20),
        plot.subtitle = element_text(size=16),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = 'light grey', linetype = 2, size = 1),
        panel.border = element_rect(linetype = 'solid', colour = 'black', size = 1.5)) 

panel2_map <- ggplot(panel2$census) +
  geom_point(aes(x = x , y = y, colour = species)) +
  labs(subtitle = paste0('S =  ', S_pool, ', N = ', N_pool_small, '\nrandom; lower N')) +
  scale_x_continuous(breaks = c(0, 0.5, 1), expand = c(0,0)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0,0)) +
  scale_color_viridis_d() +
  labs(tag = 'C') +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size=20),
        plot.subtitle = element_text(size=16),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = 'light grey', linetype = 2, size = 1),
        panel.border = element_rect(linetype = 'solid', colour = 'black', size = 1.5))

panel3_map <- ggplot(panel3$census) +
  geom_point(aes(x = x , y = y, colour = species)) +
  labs(subtitle = paste0('S =  ', S_pool, ', N = ', N_pool_big, '\nclumped within species')) +
  scale_x_continuous(breaks = c(0, 0.5, 1), expand = c(0,0)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0,0)) +
  scale_color_viridis_d() +
  labs(tag = 'D') +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size=20),
        plot.subtitle = element_text(size=16),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = 'light grey', linetype = 2, size = 1),
        panel.border = element_rect(linetype = 'solid', colour = 'black', size = 1.5))

panel4_map <- ggplot(panel4$census) +
  geom_point(aes(x = x , y = y, colour = species)) +
  labs(subtitle = paste0('S =  ', S_pool, ', N = ', N_pool_big, '\nrandom; less even SAD')) +
  scale_x_continuous(breaks = c(0, 0.5, 1), expand = c(0,0)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0,0)) +
  scale_color_viridis_d() +
  labs(tag = 'B') +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size=20),
        plot.subtitle = element_text(size=16),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = 'light grey', linetype = 2, size = 1),
        panel.border = element_rect(linetype = 'solid', colour = 'black', size = 1.5))

# # not used (use Dan's simple cartoon)
# cowplot::plot_grid(panel1_map,
#                   panel2_map,
#                   panel3_map, nrow = 1)

# rarefaction curves for results
minN1 = min(rowSums(sample1$spec_dat))
alphabar_curve1 <- bind_rows(rarefaction(sample1$spec_dat[1,], method = 'IBR', effort = 1:minN1),
          rarefaction(sample1$spec_dat[2,], method = 'IBR', effort = 1:minN1),
          rarefaction(sample1$spec_dat[3,], method = 'IBR', effort = 1:minN1),
          rarefaction(sample1$spec_dat[4,], method = 'IBR', effort = 1:minN1)) %>% 
  colMeans() %>% 
  as.numeric() %>% 
  as_tibble() %>% 
  rename(S = value) %>% 
  mutate(N = 1:n())

gamma_curve1 = rarefaction(sample1$spec_dat, method = 'IBR') %>% 
  as.numeric() %>% 
  as_tibble() %>% 
  rename(S = value) %>% 
  mutate(N = 1:n())

minN2 = min(rowSums(sample2$spec_dat))
alphabar_curve2 <- bind_rows(rarefaction(sample2$spec_dat[1,], method = 'IBR', effort = 1:minN2),
                             rarefaction(sample2$spec_dat[2,], method = 'IBR', effort = 1:minN2),
                             rarefaction(sample2$spec_dat[3,], method = 'IBR', effort = 1:minN2),
                             rarefaction(sample2$spec_dat[4,], method = 'IBR', effort = 1:minN2)) %>% 
  colMeans() %>% 
  as.numeric() %>% 
  as_tibble() %>% 
  rename(S = value) %>% 
  mutate(N = 1:n())

gamma_curve2 = rarefaction(sample2$spec_dat, method = 'IBR') %>% 
  as.numeric() %>% 
  as_tibble() %>% 
  rename(S = value) %>% 
  mutate(N = 1:n())


minN3 = min(rowSums(sample3$spec_dat))
alphabar_curve3 <- bind_rows(rarefaction(sample3$spec_dat[1,], method = 'IBR', effort = 1:minN3),
                             rarefaction(sample3$spec_dat[2,], method = 'IBR', effort = 1:minN3),
                             rarefaction(sample3$spec_dat[3,], method = 'IBR', effort = 1:minN3),
                             rarefaction(sample3$spec_dat[4,], method = 'IBR', effort = 1:minN3)) %>% 
  colMeans() %>% 
  as.numeric() %>% 
  as_tibble() %>% 
  rename(S = value) %>% 
  mutate(N = 1:n())

gamma_curve3 = rarefaction(sample3$spec_dat, method = 'IBR') %>% 
  as.numeric() %>% 
  as_tibble() %>% 
  rename(S = value) %>% 
  mutate(N = 1:n())


minN4 = min(rowSums(sample4$spec_dat))
alphabar_curve4 <- bind_rows(rarefaction(sample4$spec_dat[1,], method = 'IBR', effort = 1:minN4),
                             rarefaction(sample4$spec_dat[2,], method = 'IBR', effort = 1:minN4),
                             rarefaction(sample4$spec_dat[4,], method = 'IBR', effort = 1:minN4),
                             rarefaction(sample4$spec_dat[4,], method = 'IBR', effort = 1:minN4)) %>% 
  colMeans() %>% 
  as.numeric() %>% 
  as_tibble() %>% 
  rename(S = value) %>% 
  mutate(N = 1:n())

gamma_curve4 = rarefaction(sample4$spec_dat, method = 'IBR') %>% 
  as.numeric() %>% 
  as_tibble() %>% 
  rename(S = value) %>% 
  mutate(N = 1:n())

rare1 <- ggplot() +
  geom_line(data = gamma_curve1,
            aes(x = N, y = S), size = 1.5) +
  geom_line(data = alphabar_curve1,
            aes(x = N, y = S), size = 1.5, lty = 2, colour = 'light grey') +
  # annotate('text', x = 700, y = 35,
  #          label = paste('beta[W]==', round(bS1$value[[1]], digits = 1)),
  #          parse = T, size = 6) +
  # annotate('text', x = 700, y = 25,
  #          label = paste('beta[C]==', round(bC1[[1]], digits = 1)),
  #          parse = T, size = 6) +
  labs(tag = 'E') +
  theme_classic() +
  theme(text=element_text(size=20)) +
  theme(#axis.text = element_blank(),
        axis.title = element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_text(size=13))

rare2 <- ggplot() +
  geom_line(data = gamma_curve2,
            aes(x = N, y = S), size = 1.5) +
  geom_line(data = alphabar_curve2,
            aes(x = N, y = S), size = 1.5, lty = 2, colour = 'light grey') +
  # annotate('text', x = 300, y = 32,
  #          label = paste('beta[W]==', round(bS2$value[[1]], digits = 1)),
  #          parse = T, size = 6) +
  # annotate('text', x = 300, y = 22,
  #          label = paste('beta[C]==', round(bC2[[1]], digits = 1)),
  #          parse = T, size = 6) +
  labs(tag = 'G') +
  theme_classic() +
  theme(text=element_text(size=20)) +
  theme(#axis.text = element_blank(),
        axis.title = element_blank(), 
        axis.ticks = element_blank(),
        axis.text = element_text(size=13))


rare3 <- ggplot() +
  geom_line(data = gamma_curve3,
            aes(x = N, y = S), size = 1.5) +
  geom_line(data = alphabar_curve3,
            aes(x = N, y = S), size = 1.5, lty = 2, colour = 'light grey') +
  # annotate('text', x = 700, y = 35,
  #          label = paste('beta[W]==', round(bS3$value[[1]], digits = 1)),
  #          parse = T, size = 6) +
  # annotate('text', x = 700, y = 25,
  #          label = paste('beta[C]==', round(bC3[[1]], digits = 1)),
  #          parse = T, size = 6) +
  
  labs(tag = 'H') +
  theme_classic() +
  theme(text=element_text(size=20)) +
  theme(#axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(size=13))


rare4 <- ggplot() +
  geom_line(data = gamma_curve4,
            aes(x = N, y = S), size = 1.5) +
  geom_line(data = alphabar_curve4,
            aes(x = N, y = S), size = 1.5, lty = 2, colour = 'light grey') +
  # annotate('text', x = 700, y = 35,
  #          label = paste('beta[W]==', round(bS3$value[[1]], digits = 1)),
  #          parse = T, size = 6) +
  # annotate('text', x = 700, y = 25,
  #          label = paste('beta[C]==', round(bC3[[1]], digits = 1)),
  #          parse = T, size = 6) +
  
  labs(tag = 'F') +
  theme_classic() +
  theme(text=element_text(size=20)) +
  theme(#axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_text(size=13))

# legend for linetype and colour 
curve_legend = ggplot() +
  geom_line(data = gamma_curve3,
            aes(x = N, y = S, 
                colour = 'gamma', linetype = 'gamma'), size = 1.5) +
  geom_line(data = alphabar_curve3,
            aes(x = N, y = S, colour = 'alpha', linetype = 'alpha'), size = 1.5) +
  scale_colour_manual(name = 'Scale',
                      values = c('alpha' = 'light grey',
                                 'gamma' = 'black'),
                      label = c(expression(alpha),
                                expression(gamma))) +
  scale_linetype_manual(name = 'Scale',
                        values = c('alpha' = 2,
                                   'gamma' = 1),
                        label = c(expression(alpha),
                                  expression(gamma))) +
  annotate('text', x = 700, y = 35,
           label = paste('beta[W]==', round(bS3$value[[1]], digits = 1)),
           parse = T, size = 6) +
  annotate('text', x = 700, y = 25,
           label = paste('beta[C]==', round(bC3$value[[1]], digits = 1)),
           parse = T, size = 6) +
  
  theme_classic() +
  theme(text=element_text(size=18)) +
  theme(#axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.position = 'top')

curve_leg = gg_legend(curve_legend)

# insets plots (to replace reporting values with text)
inset1 <- ggplot() +
  geom_point(data = bind_rows(bS1 %>% 
                                filter(index=='beta_S'), 
                              bC1),
             aes(x = index, y = value)) +
  scale_x_discrete(limits = c('beta_S', 'beta_S_C'),
                   labels = c(expression(beta[S]), expression(beta[C])), #expression(beta[S[PIE]]), expression(beta[S[n]]),
                   name = '') +
  scale_y_continuous(breaks= c(1, 1.1, 1.2, 1.3, 1.4),
                     labels = c(1, '', 1.2, '', 1.4),
                     limits = c(0.95, 1.4), name = '') +
  theme_minimal() +
  theme(text=element_text(size=20)) +
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=16))
        


inset2 <- ggplot() +
  geom_point(data = bind_rows(bS2 %>% 
                                filter(index=='beta_S'), 
                              bC2),
             aes(x = index, y = value)) +
  scale_x_discrete(limits = c('beta_S', 'beta_S_C'),
                   labels = c(expression(beta[S]), expression(beta[C])),
                   name = '') +
  scale_y_continuous(breaks= c(1, 1.1, 1.2, 1.3, 1.4),
                     labels = c(1, '', 1.2, '', 1.4),
                     limits = c(0.95, 1.4), name = '') +
  theme_minimal() +
  theme(text=element_text(size=20)) +
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=16))

inset3 <- ggplot() +
  geom_point(data = bind_rows(bS3 %>% 
                                filter(index=='beta_S'), 
                              bC3),
             aes(x = index, y = value)) +
  scale_x_discrete(limits = c('beta_S', 'beta_S_C'),
                   labels = c(expression(beta[S]), expression(beta[C])),
                   name = '') +
  scale_y_continuous(breaks= c(1, 1.1, 1.2, 1.3, 1.4),
                     labels = c(1, '', 1.2, '', 1.4),
                     limits = c(0.95, 1.4), name = '') +
  theme_minimal() +
  theme(text=element_text(size=20)) +
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=16))

inset4 <- ggplot() +
  geom_point(data = bind_rows(bS4 %>% 
                                filter(index=='beta_S'), 
                              bC4),
             aes(x = index, y = value)) +
  scale_x_discrete(limits = c('beta_S', 'beta_S_C'),
                   labels = c(expression(beta[S]), expression(beta[C])),
                   name = '') +
  scale_y_continuous(breaks= c(1, 1.1, 1.2, 1.3, 1.4),
                     labels = c(1, '', 1.2, '', 1.4),
                     limits = c(0.95, 1.4), name = '') +
  theme_minimal() +
  theme(text=element_text(size=20)) +
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=16))

all_curves <- plot_grid(curve_leg,
                        plot_grid(NULL, 
                                  rare1 + annotation_custom(ggplotGrob(inset1), xmin = 250, xmax = 1000, ymin = -5, ymax = 40),
                                  rare4 + annotation_custom(ggplotGrob(inset4), xmin = 250, xmax = 1000, ymin = -5, ymax = 40),
                                  rare2 + annotation_custom(ggplotGrob(inset2), xmin = 100, xmax = 390, ymin = -5, ymax = 40),
                                  rare3 + annotation_custom(ggplotGrob(inset3), xmin = 250, xmax = 1000, ymin = -5, ymax = 40),
                                  rel_widths = c(0.05, 1, 1, 1, 1),
                                  nrow = 1),
                        NULL,
                        rel_heights = c(0.08, 1, 0.05), ncol = 1) +
  draw_label(label = 'Number of individuals', y = 0.03, size = 18) +
  draw_label(label = 'Species richness', x = 0.01, angle = 90, size = 18) 


plot_grid(cowplot::plot_grid(panel1_map,
                             panel4_map,
                             panel2_map,
                             panel3_map, 
                             nrow = 1), 
          all_curves,
          rel_heights = c(1, 0.9),
          nrow = 2) 

ggsave('./figs/non-spatial-vs-spatial-beta-diversity.pdf',
        width = 293, height = 160, units = 'mm')

# for making additional edits in powerpoint 
ggsave('./figs/non-spatial-vs-spatial-beta-diversity.svg',
       width = 293, height = 160, units = 'mm')
