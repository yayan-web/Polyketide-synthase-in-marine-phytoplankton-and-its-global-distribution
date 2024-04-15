setwd("D:/冈比亚藻综述/浮游植物PKS基因的多样性/全球分布")
library(ggplot2)
library(reshape2)
library(ggpubr)
data = read.csv("domin_env_reads_2.csv", check.names=FALSE, encoding = "UTF-8")
colnames(data)
unique(data$Depth)

length(data_long)


data_long = melt(data, measure.vars = c("PF14765 (DH)",
                                        "PF02801 (KS-C)",
                                        "PF00109 (KS-N)",
                                        "PF00975 (TE)",
                                        "PF08659 (KR)",
                                        "PF00698 (AT)",
                                        "PF00550 (PP-ACP)",
                                        "PF16197 (CE)"))

X11()
length(subset(data_long, Size %in% c("0.8 - 5", 
                                     "5 - 20", 
                                     "20 - 180")
              &Depth == "SRF"))
data_long$Size = factor(data_long$Size, levels = unique(data_long$Size))

p1 = ggplot(subset(data_long, Size %in% c("0.8 - 5", 
                                         "5 - 20", 
                                         "20 - 180")
                   &variable %in% c("PF14765 (DH)",
                                    "PF02801 (KS-C)",
                                    "PF00109 (KS-N)",
                                    "PF00975 (TE)")),
           aes(x = `Temperature (°C)`,  y = log2(value), col = Depth)) +
  geom_point(size = 0.5) +
  facet_wrap(variable~Size, 
             scales = "free", 
             ncol = 3) +
  geom_smooth(se = F, method = "lm", formula = y~x+I(x^2)) +
  theme_bw() +
  xlim(10, 30) +
  labs(x = "Temperature (°C)",
       y = "log2(Relative abundance)") +
  theme(strip.text = element_text(size = 8),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))

p2 = ggplot(subset(data_long, Size %in% c("0.8 - 5", 
                                          "5 - 20", 
                                          "20 - 180")
                   &variable %in% c("PF08659 (KR)",
                                    "PF00698 (AT)",
                                    "PF00550 (PP-ACP)",
                                    "PF16197 (CE)")),
            aes(x = `Temperature (°C)`, y = log2(value), col = Depth)) +
  geom_point(size = 0.5) +
  facet_wrap(variable~Size, 
             scales = "free", 
             ncol = 3) +
  geom_smooth(se = F, method = "lm", formula = y~x+I(x^2)) +
  theme_bw() +
  labs(x = "Temperature (°C)",
       y = "log2(Relative abundance") +
  xlim(10, 30)+
  theme(strip.text = element_text(size = 8),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))


ggarrange(p1, p2, ncol = 2, common.legend = T, legend = "bottom")

ggsave("domain_T.tiff", dpi = 300)






























































