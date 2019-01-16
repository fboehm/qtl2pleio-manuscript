# from page 197-198 on google drive research notebook

trace <- rep(1:4, each = 5)
distance <- rep(c(0, 0.5, 1, 2, 3), times = 4)
power <- c(0.0425, 0.61, 0.8575, 0.935, 0.985,
           0.0575, 0.96, 1, 1, 1,
           0.05, 0.4175, 0.6375, 0.7725, 0.8775,
           0.0425, 0.8125, 0.98, 1, 1
           )
library(dplyr)
tibble(trace, distance, power) -> mytib
library(ggplot2)
mytib %>%
  group_by(trace) %>%
  ggplot() + geom_line(aes(x = distance, y = power, colour = as.factor(trace))) +
  broman::theme_karl(legend.position = "none") +
  xlab("Distance (cM)") + ylab("Power")
ggsave("R/power-curves.eps", height=4, width=7.5)
ggsave("R/power-curves.png")
