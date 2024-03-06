
### method shape table 
###	 prop. 0   soft 1   hard 2 adaptive 3 POET 4
### glasso 5 banded 6 BayesG 7 BayesKBC 8

### method line-type table 
###	prop. "solid"   soft "longdash" hard "dotted" adaptive "4C88C488" 
### POET "12345678"
### glasso "twodash" banded "dotdash" BayesG "dashed" BayesKBC "1F"

### method color table 
### prop. "#F8766D" soft "#D39200" hard "#93AA00" adaptive "#00BA38" 
### POET  "#00C19F"
### glasso "#00B9E3" banded "#619CFF" BayesG "#DB72FB" BayesKBC "#FF61C3" 



require(ggplot2)
require(readxl)
require(devEMF)
require(svglite)

PLOT_DATA_A <- read_excel(path = "Results For Figure 2 B1.xlsx", sheet = "Scenario 1")
PLOT_DATA_A_TIME <- read_excel(path = "Results For Figure 2 B1.xlsx", sheet = "Scenario 1 Time")
PLOT_DATA_B <- read_excel(path = "Results For Figure 2 B1.xlsx", sheet = "Scenario 2")
PLOT_DATA_C <- read_excel(path = "Results For Figure 2 B1.xlsx", sheet = "Scenario 3")

PLOT_DATA_A <- data.frame(PLOT_DATA_A)
PLOT_DATA_A_TIME <- data.frame(PLOT_DATA_A_TIME)
PLOT_DATA_B <- data.frame(PLOT_DATA_B)
PLOT_DATA_C <- data.frame(PLOT_DATA_C)

PLOT_DATA_A_SIGMA <- PLOT_DATA_A[1 : 90, ]
PLOT_DATA_A_OMEGA <- PLOT_DATA_A[91: 180, ]
PLOT_DATA_B_SIGMA <- PLOT_DATA_B[1 : 30, ]
PLOT_DATA_B_OMEGA <- PLOT_DATA_B[31: 60, ]

PLOT_DATA_A_SIGMA$method <- factor(PLOT_DATA_A_SIGMA$method, levels = c("prop.", "soft", "hard", "adaptive", "POET"))
PLOT_DATA_A_OMEGA$method <- factor(PLOT_DATA_A_OMEGA$method, levels = c("prop.", "glasso", "banded", "BayesG", "BayesKBC"))
PLOT_DATA_A_TIME$method <- factor(PLOT_DATA_A_TIME$method, levels = c("prop.", "soft", "hard", "adaptive", "POET", "glasso", "banded", "BayesG", "BayesKBC"))
PLOT_DATA_A_TIME$matrix <- factor(PLOT_DATA_A_TIME$matrix, levels = c("Sigma", "Omega"))
PLOT_DATA_B_SIGMA$method <- factor(PLOT_DATA_B_SIGMA$method, levels = c("prop.", "soft", "hard", "adaptive", "POET"))
PLOT_DATA_B_OMEGA$method <- factor(PLOT_DATA_B_OMEGA$method, levels = c("prop.", "glasso", "banded", "BayesG", "BayesKBC"))
PLOT_DATA_C$method <- factor(PLOT_DATA_C$method, levels = c("prop.", "soft", "hard", "adaptive", "POET"))

##################
### Scenario 1 ###
##################

sample_size <- 50 ### 50, 100, 150
line_size <- 2
shape_size <- 5
font_size <- 38

ggplot(data = subset(PLOT_DATA_A_SIGMA, n == sample_size), 
	mapping = aes(x = p, 
				  y = log(loss),
		  	  group = method)) +
  geom_line(mapping = aes(color = method, linetype = method), 
			   size = line_size,
			  alpha = 0.9) +
  geom_point(mapping = aes(color = method, shape = method), 
				size = shape_size) +
  theme_classic() +
  scale_x_continuous(breaks = c(150, 225, 300)) +
  labs(x = expression(paste('dimension ', italic(p))),
       y = expression(paste('log loss of ', Sigma, ' estimator'))) +
  theme(legend.position = c(0.90, 0.2),
				   text = element_text(size = font_size),
	   strip.background = element_blank()) +
  scale_shape_manual(values = seq(0, 4)) +
  scale_linetype_manual(values = c("solid", "longdash", "dotted", "4C88C488", "12345678")) + 
  scale_color_manual(values = c("#F8766D", "#D39200", "#93AA00", "#00BA38", "#00C19F")) + 
  facet_wrap(~ norm, scales = "free_y")
  
ggsave("Scenario_1_Sigma_SampleSize_50_138_by_915.svg", width = 13.8, height = 9.15, units = "in") # 13.8 by 9.15

ggplot(data = subset(PLOT_DATA_A_OMEGA, n == sample_size), 
	mapping = aes(x = p, 
				  y = log(loss),
		  	  group = method)) +
  geom_line(mapping = aes(color = method, linetype = method), 
			   size = line_size,
			  alpha = 0.9) +
  geom_point(mapping = aes(color = method, shape = method), 
				size = shape_size) +
  theme_classic() +
  scale_x_continuous(breaks = c(150, 225, 300)) +
  labs(x = expression(paste('dimension ', italic(p))),
       y = expression(paste('log loss of ', Omega, ' estimator'))) +
  theme(legend.position = c(0.90, 0.3),
				   text = element_text(size = font_size),
	   strip.background = element_blank()) +
  scale_shape_manual(values = c(0, seq(5, 8))) +
  scale_linetype_manual(values = c("solid", "twodash", "dotdash", "dashed", "1F")) + 
  scale_color_manual(values = c("#F8766D", "#00B9E3", "#619CFF", "#DB72FB", "#FF61C3")) + 
  facet_wrap(~ norm, scales = "free_y")

ggsave("Scenario_1_Omega_SampleSize_50_138_by_915.svg", width = 13.8, height = 9.15, units = "in")

#######################
### Scenario 1 Time ###
#######################

ggplot(data = subset(PLOT_DATA_A_TIME, n == sample_size), 
	mapping = aes(x = p, 
				  y = log(time),
		  	  group = method)) +
  geom_line(mapping = aes(color = method, linetype = method), 
			   size = line_size,
			  alpha = 0.9) +
  geom_point(mapping = aes(color = method, shape = method), 
				size = shape_size) +
  theme_classic() +
  scale_x_continuous(breaks = c(150, 225, 300)) +
  labs(x = expression(paste('dimension ', italic(p))),
       y = "log of the execution time (in second)") +
  theme(legend.position = c(0.90, 0.7),
				   text = element_text(size = font_size),
	   strip.background = element_blank()) +
  scale_shape_manual(values = seq(0, 8)) +
  scale_linetype_manual(values = c("solid", "longdash", "dotted", "4C88C488", "12345678", "twodash", "dotdash", "dashed", "1F")) + 
  scale_color_manual(values = c("#F8766D", "#D39200", "#93AA00", "#00BA38", "#00C19F", "#00B9E3", "#619CFF", "#DB72FB", "#FF61C3")) + 
  facet_wrap(~ matrix, scales = "free_y")

ggsave("Scenario_1_Time_SampleSize_50_138_by_915.svg", width = 13.8, height = 9.15, units = "in")

##################
### Scenario 2 ###
##################

ggplot(data = PLOT_DATA_B_SIGMA, 
	mapping = aes(x = sigma, 
				  y = log(loss),
		  	  group = method)) +
  geom_line(mapping = aes(color = method, linetype = method), 
			   size = line_size,
			  alpha = 0.9) +
  geom_point(mapping = aes(color = method, shape = method), 
				size = shape_size) +
  theme_classic() +
  scale_x_continuous(breaks = c(0.1, 0.5, 0.8)) +
  labs(x = expression(paste('coefficient ', sigma)),
       y = expression(paste('log loss of ', Sigma, ' estimator'))) +
  theme(legend.position = c(0.90, 0.2),
				   text = element_text(size = font_size),
	   strip.background = element_blank()) +
  scale_shape_manual(values = seq(0, 4)) +
  scale_linetype_manual(values = c("solid", "longdash", "dotted", "4C88C488", "12345678")) + 
  scale_color_manual(values = c("#F8766D", "#D39200", "#93AA00", "#00BA38", "#00C19F")) + 
  facet_wrap(~ norm, scales = "free_y")

ggsave("Scenario_2_Sigma_138_by_915.svg", width = 13.8, height = 9.15, units = "in")


ggplot(data = PLOT_DATA_B_OMEGA, 
	mapping = aes(x = sigma, 
				  y = log(loss),
		  	  group = method)) +
  geom_line(mapping = aes(color = method, linetype = method), 
			   size = line_size,
			  alpha = 0.9) +
  geom_point(mapping = aes(color = method, shape = method), 
				size = shape_size) +
  theme_classic() +
  scale_x_continuous(breaks = c(0.1, 0.5, 0.8)) +
  labs(x = expression(paste('coefficient ', sigma)),
       y = expression(paste('log loss of ', Omega, ' estimator'))) +
  theme(legend.position = c(0.90, 0.6),
				   text = element_text(size = font_size),
       strip.background = element_blank()) +
  scale_shape_manual(values = c(0, seq(5, 8))) +
  scale_linetype_manual(values = c("solid", "twodash", "dotdash", "dashed", "1F")) +
  scale_color_manual(values = c("#F8766D", "#00B9E3", "#619CFF", "#DB72FB", "#FF61C3")) + 
  facet_wrap(~ norm, scales = "free_y")

ggsave("Scenario_2_Omega_138_by_915.svg", width = 13.8, height = 9.15, units = "in")

##################
### Scenario 3 ###
##################

ggplot(data = PLOT_DATA_C, 
	mapping = aes(x = K, 
				  y = log(loss),
		  	  group = method)) +
  geom_line(mapping = aes(color = method, linetype = method), 
			   size = line_size,
			  alpha = 0.9) +
  geom_point(mapping = aes(color = method, shape = method), 
				size = shape_size) +
  theme_classic() +
  scale_x_continuous(breaks = c(30, 40, 50)) +
  labs(x = expression(paste('number of blocks ', italic(K))),
       y = expression(paste('log loss of ', Sigma, ' estimator'))) +
  theme(legend.position = c(0.90, 0.4),
				   text = element_text(size = font_size),
	   strip.background = element_blank()) +
  scale_shape_manual(values = seq(0, 4)) +
  scale_linetype_manual(values = c("solid", "longdash", "dotted", "4C88C488", "12345678")) + 
  scale_color_manual(values = c("#F8766D", "#D39200", "#93AA00", "#00BA38", "#00C19F")) + 
  facet_wrap(~ norm, scales = "free_y")

ggsave("Scenario_3_138_by_915.svg", width = 13.8, height = 9.15, units = "in")







