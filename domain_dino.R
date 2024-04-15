library(reshape2)
library(ggplot2)
setwd("D:/论文/PKS/浮游植物PKS基因的多样性/全球分布")
phy_ann = read.csv("phy_ann.csv", header = T,
                   encoding = "UTF-8",
                   check.names = F)
colnames(phy_ann)
unique(phy_ann$Phylum)
dino_ann = subset(phy_ann, Class == "Dinophyceae")
dino_ann = dino_ann[c("Phylum", "Order")]
dino_ann = dino_ann[!duplicated(dino_ann),]
phy_data = read.csv("phy_order_data.csv", header = T, 
                    encoding = "UTF-8", check.names = F)
colnames(phy_data)[2] = "Order"

phy_data2 = as.data.frame(t(phy_data[3:length(phy_data)]))
phy_data2_sum = rowSums(phy_data2)
phy_data3 = phy_data2/phy_data2_sum
phy_data4 = as.data.frame(t(phy_data3))
phy_data5 = cbind(phy_data[2], phy_data4)


dino_data = merge(dino_ann, phy_data5, by = "Order")
head(dino_data)
rownames(dino_data) =  dino_data$Order
dino_data$Order = NULL
dino_data$Phylum = NULL
dino_data$TARA_N000000077
dino_data = as.data.frame(t(dino_data))









prd_data = read.csv("domin_env_reads_2.csv", header = T, check.names = F,
                    encoding = "UTF-8")
prd_data[1] = NULL
prd_data[12] = NULL
head(prd_data)

id_barcode = read.csv("Carrabec_et_al_Suppl_data_6.txt", se = "\t", check.names = F)

head(id_barcode)
colnames(prd_data)[1] = "Sample Code"
prd_data = merge(prd_data, id_barcode[c("Tara Oceans ID", "Sample Code")], by = "Sample Code")

colnames(prd_data)[35] = "Barcode" 




colnames(dino_data)
dino_data = merge(sample_ids[c("Barcode", "Stations")], dino_data, by = "Barcode")

head(prd_data)
head(dino_data)


length(dino_data$Barcode)
length(unique(dino_data$Barcode))

length(prd_data$Barcode)
length(unique(prd_data$Barcode))

prd_env_dino = merge(prd_data, dino_data, by = "Barcode")


colnames(prd_env_dino)

tail(prd_env_dino)
colnames(prd_env_dino[40:53])

prd_env_dino["Dinophyceae"] = rowSums(prd_env_dino[40:53])

prd_env_dino_long = melt(prd_env_dino, measure.vars = c("PF14765 (DH)", "PF02801 (KS-C)", 
                                                        "PF00109 (KS-N)", "PF00975 (TE)", 
                                                        "PF08659 (KR)", "PF00698 (AT)", 
                                                        "PF16197 (CE)"))
x11()


prd_env_dino_long$value = as.numeric(prd_env_dino_long$value)
unique(prd_env_dino_long$Depth)

library(GGally)
colnames(prd_env_dino)

head(prd_env_dino)

write.csv(prd_env_dino, "prd_env_dino.csv")


# [1] "Blastodiniales"    "Dinophysiales"     "Gonyaulacales"     "Gymnodiniales"     "Lophodiniales"    
# [6] "Noctilucales"      "Oxyrrhinales"      "Peridiniales"      "Phytodiniales"     "Prorocentrales"   
# [11] "Pyrocystales"      "Suessiales"        "Thoracosphaerales" "Uncultured"   


head(prd_env_dino_long)

p = ggplot(data = subset(prd_env_dino_long,
                         # variable %in% c("PF14765 (DH)", "PF00975 (TE)")
                         ), 
           aes(x = log10(Gymnodiniales*100), y = log10(value*100), col = Depth)) +
  scale_color_brewer(palette = "Dark2") +
  geom_point(size = 1, alpha = 0.5) +
  theme_bw() +
  labs(x = "log10(Relative abundance of Gymnodiniales (%))",
       y = "log10(Relative abundance of Domains (%))") +
  geom_smooth(method = "lm", alpha = 0.1, linewidth = 1) +
  facet_wrap(~variable, scale = "free_y", ncol = 4) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.position = c(0.85, 0.2),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15),
        strip.text = element_text(size = 13, face = "bold"))
p

ggsave("Gym_domins.tiff", dpi = 300)

# 正向的的：Gymnodiniales，

for (depth in c("DCM", "SRF")) {
  for (domain in c("PF14765 (DH)", "PF02801 (KS-C)", "PF00109 (KS-N)", "PF00975 (TE)", 
                   "PF08659 (KR)", "PF00698 (AT)", "PF16197 (CE)")){
    df <- subset(prd_env_dino_long, Depth == depth & variable == domain)
    result <- cor.test(log10(df$Gymnodiniales*100), log10(df$value*100), method = "pearson")
    print(paste(depth, domain, result$estimate, result$p.value, sep = "  "))
  }
}
# [1] "DCM  PF14765 (DH)  0.449213131826998  3.70041773445925e-05"
# [1] "DCM  PF02801 (KS-C)  0.482760395428485  7.59656698658141e-06"
# [1] "DCM  PF00109 (KS-N)  0.236543147126754  0.0370645033972327"
# [1] "DCM  PF00975 (TE)  0.476418851021601  1.03799532951903e-05"
# [1] "DCM  PF08659 (KR)  0.111283817660217  0.332055857926491"
# [1] "DCM  PF00698 (AT)  0.486870316542147  6.1844052475931e-06"
# [1] "DCM  PF16197 (CE)  0.165870457495458  0.146681985171587"
# [1] "SRF  PF14765 (DH)  0.340102323961185  5.44003998802522e-05"
# [1] "SRF  PF02801 (KS-C)  0.370941845519861  9.46976830305529e-06"
# [1] "SRF  PF00109 (KS-N)  0.358002325072597  2.01605761817796e-05"
# [1] "SRF  PF00975 (TE)  0.358675375419002  1.93990627021343e-05"
# [1] "SRF  PF08659 (KR)  0.20563188775279  0.0167277660440958"
# [1] "SRF  PF00698 (AT)  0.428886351825851  2.10369006604336e-07"
# [1] "SRF  PF16197 (CE)  0.0774506962525274  0.371926887205708"



# 新的环境因子

head(prd_data)
prd_data = prd_data[c("PF14765 (DH)", "PF02801 (KS-C)", "PF00109 (KS-N)", "PF00975 (TE)", 
                      "PF08659 (KR)", "PF00698 (AT)", "PF16197 (CE)", "Barcode")]
env_data = read.csv("enviro_18SV9v1.csv", header = T, check.names = F, encoding = "utf-8")

env_data = env_data[c("Barcode", "Temperature", "Chlorophyll_A", "Shannon_Darwin*", 
                      "Fraction_Lower", "Fraction_Upper", "Depth_Nominal")]

prd_env = merge(prd_data, env_data, by = "Barcode")
prd_env["Fraction"] = paste(prd_env$Fraction_Lower, prd_env$Fraction_Upper, sep = "-")
prd_env$Fraction_Lower = NULL
prd_env$Fraction_Upper = NULL

colnames(prd_env)

library(GGally)

head(prd_env)
tail(prd_env)

x11()
prd_env = subset(prd_env, Fraction %in% c("0.8-5", "5-20", "20-180"))
prd_env$Fraction = factor(prd_env$Fraction, levels = c("0.8-5", "5-20", "20-180"))

ggpairs(prd_env, columns = c("PF14765 (DH)", "PF00975 (TE)", "PF00698 (AT)", 
                             "PF16197 (CE)", 
                             "Chlorophyll_A", "Shannon_Darwin*"), 
        alpha = 0.5,
        ggplot2::aes(colour = interaction(Depth_Nominal , Fraction))) +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw()

prd_env_long = melt(prd_env, measure.vars = c("PF14765 (DH)", 
                                              "PF02801 (KS-C)", 
                                              "PF00109 (KS-N)",
                                              "PF00975 (TE)", 
                                              "PF08659 (KR)", 
                                              "PF00698 (AT)", 
                                              "PF16197 (CE)"))

prd_env_long = subset(prd_env_long, Fraction %in% c("0.8-5", "5-20", "20-180"))


prd_env_long$Fraction = factor(prd_env_long$Fraction, levels = c("0.8-5",
                                                                       "5-20",
                                                                       "20-180"))

p = ggplot(subset(prd_env_long, 
                  variable %in% c("PF14765 (DH)", "PF00698 (AT)", "PF00975 (TE)")
                  ), 
           aes(x = Temperature, y = log10(value*100), col = Depth_Nominal)) + 
  scale_color_brewer(palette = "Dark2") +
  geom_point(size = 3, alpha = 0.2) +
  theme_bw() +
  geom_smooth(data = subset(prd_env_long, 
                            variable %in% c("PF14765 (DH)", "PF00698 (AT)", "PF00975 (TE)") &
                              Depth_Nominal == "SRF"),
    method = "glm", size = 2,
    formula =  y ~ x + I((x - 26.0) * (x > 26.0)), se = F) +
  geom_vline(xintercept = 26, linewidth = 1, linetype = "dashed") +
  geom_smooth(data = subset(prd_env_long, 
                            variable %in% c("PF14765 (DH)", "PF00698 (AT)", "PF00975 (TE)") &
                              Depth_Nominal == "DCM"),
              method = "glm", size = 2,
              formula =  y ~ x , se = F) +
  facet_wrap(variable~Fraction , scale = "free_y", ncol = 3) +
  labs(x = "Temperature [°C]",
       y = "log10(Relative abundance of Domains (%))") +
  theme(strip.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.position = "top",
        legend.title = element_text(size = 20))
p

ggsave("domain_t_2.tiff", dpi = 300)



for (depth in c("DCM", "SRF")) {
  for (domain in c("PF14765 (DH)", "PF00975 (TE)", 
                   "PF00698 (AT)")){
    for (size in c("0.8-5", "5-20", "20-180")) {
      df = subset(prd_env_long, Depth_Nominal == depth & variable == domain & Fraction == size)
      df2 = subset(df, Temperature < 26)
      result = cor.test(df2$Temperature, log10(df2$value*100), method = "pearson")
      print(paste(depth, domain, size, result$estimate, result$p.value, sep = "  "))
    }
  }
}


dino_env = merge(dino_data, env_data, by = "Barcode")

colnames(dino_env)
dino_env[2:5] = NULL


dino_env["Fraction"] = paste(dino_env$Fraction_Lower, dino_env$Fraction_Upper, sep = "-")
dino_env$Fraction_Lower = NULL
dino_env$Fraction_Upper = NULL

dino_env = subset(dino_env, Fraction %in% c("0.8-5", "5-20", "20-180"))


dino_env$Fraction = factor(dino_env$Fraction, levels = c("0.8-5",
                                                                       "5-20",
                                                                       "20-180"))


p = ggplot(dino_env, aes(x = Temperature, y = Gymnodiniales*100)) +
  geom_point() +
  facet_wrap(~Fraction, scale = "free_y") +
  theme_bw() +
  geom_smooth(method = "glm", size = 2,
              formula =  y ~ x + I((x - 25.0) * (x > 25.0)), se = F) +
  labs(x = "Temperature [°C]", y = "Gymnodiniales (%)") +
  theme(strip.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.position = "top",
        legend.title = element_text(size = 20))
p

head(prd_env)

colnames(prd_data)
p = ggplot(prd_data, aes(x = `PF02801 (KS-C)`*100, y = `PF00109 (KS-N)`*100)) + 
  geom_point(alpha = 0.2) +
  theme_bw() +
  scale_x_continuous(trans = "log10") +
  scale_y_continuous(trans = "log10") +
  geom_smooth(method = "lm", linewidth = 2) +
  annotate("text", x = 0.000002, y = 0.0002, 
           label = "p-value < 2.2e-16",
           size = 8) +
  labs(x = "Relative abudance of KS-C (%)", y = "Relative abudance of KS-N (%)") +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15),)
p

ggsave("ks_c_n_corr.tiff", dpi = 300)
cor.test(prd_data$`PF02801 (KS-C)`, prd_data$`PF00109 (KS-N)`, method = "spearman")


colnames(prd_env_long)

p = ggplot(prd_env_long, aes(x = Fraction, 
                             y = value,
                             group = Fraction,
                             col = Fraction)) +
  scale_color_brewer(palette = "Dark2") +
  geom_jitter(alpha = 0.4) +
  geom_boxplot(alpha = 0.2) +
  theme_bw() +
  facet_wrap(~variable, scale = "free_y") +
  scale_y_continuous(trans = "log10") +
  labs(y = "Relative abudance of domains (%)", x = "Fraction [μm]") +
  theme(axis.title = element_text(size = 20),
        strip.text = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        legend.position = "none",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15))
p
ggsave("domains percent in diff fraction.tiff", dpi = 300)


# 与海岸线距离的关系
prd_data = prd_data[c("PF14765 (DH)", "PF02801 (KS-C)", "PF00109 (KS-N)", "PF00975 (TE)", 
                      "PF08659 (KR)", "PF00698 (AT)", "PF16197 (CE)", "Barcode")]
env_data = read.csv("enviro_18SV9v1.csv", header = T, check.names = F, encoding = "utf-8")
colnames(env_data)

env_data = env_data[c("Barcode", "Distance_coast", 
                      "Fraction_Lower", "Fraction_Upper", "Depth_Nominal")]

prd_env = merge(prd_data, env_data, by = "Barcode")
prd_env["Fraction"] = paste(prd_env$Fraction_Lower, prd_env$Fraction_Upper, sep = "-")
prd_env$Fraction_Lower = NULL
prd_env$Fraction_Upper = NULL

colnames(prd_env)

prd_env = subset(prd_env, Fraction %in% c("0.8-5", "5-20", "20-180"))
prd_env$Fraction = factor(prd_env$Fraction, levels = c("0.8-5", "5-20", "20-180"))

x11()
prd_env_long = melt(prd_env, measure.vars = c("PF14765 (DH)", 
                                              "PF02801 (KS-C)", 
                                              "PF00109 (KS-N)",
                                              "PF00975 (TE)", 
                                              "PF08659 (KR)", 
                                              "PF00698 (AT)", 
                                              "PF16197 (CE)"))

prd_env_long = subset(prd_env_long, Fraction %in% c("0.8-5", "5-20", "20-180"))


prd_env_long$Fraction = factor(prd_env_long$Fraction, levels = c("0.8-5",
                                                                 "5-20",
                                                                 "20-180"))
p = ggplot(prd_env_long, aes(x = log10(Distance_coast), 
                             y = value,
                             group = Fraction,
                             col = Fraction)) +
  scale_color_brewer(palette = "Dark2") +
  geom_jitter(alpha = 0.4) +
  theme_bw() +
  geom_smooth(method = "lm") +
  facet_wrap(Fraction~variable, scale = "free_y") +
  scale_y_continuous(trans = "log10") +
  labs(y = "Relative abudance of domains (%)") +
  theme(axis.title = element_text(size = 20),
        strip.text = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15),
        legend.position = "none",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15))
p



