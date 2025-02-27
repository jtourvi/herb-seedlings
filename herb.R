###Code for Tourville & Dovciak (2025): Patterns of herbivory and resource utilization 
#of tree seedlings along an altitudinal gradient

#Set-up working directory here

#required packages
library(plyr)
library(ggplot2)
library(cowplot)
library(reshape)
library(patchwork)
library(afex)
library(nlme)
library(lme4)
library(ggstance)
library(sjPlot)

#Foliar nutrients
tissues = read.csv("tissues_new.csv", header = TRUE)
head(tissues)

#Foliar Data management
tissues_a = subset(tissues, spp == "ACSA")
tissues_f = subset(tissues, spp == "FAGR")

tissue = ddply(tissues, ~spp + site, summarise, 
               P = mean(leaf_P), SEP = sd(leaf_P)/sqrt((length(leaf_P))), 
               K = mean(leaf_K), SEK = sd(leaf_K)/sqrt((length(leaf_K))),
               CA = mean(leaf_Ca), SEC = sd(leaf_Ca)/sqrt((length(leaf_Ca))), 
               MG = mean(leaf_Mg), SEMg = sd(leaf_Mg)/sqrt((length(leaf_Mg))),
               AL = mean(leaf_Al), SEA = sd(leaf_Al)/sqrt((length(leaf_Al))),
               S = mean(leaf_S), SES = sd(leaf_S)/sqrt((length(leaf_S))),
               MN = mean(leaf_Mn), SEMn = sd(leaf_Mn)/sqrt((length(leaf_Mn))),
               FE = mean(leaf_Fe), SEF = sd(leaf_Fe)/sqrt((length(leaf_Fe))),
               CU = mean(leaf_Cu), SECu = sd(leaf_Cu)/sqrt((length(leaf_Cu))),
               B = mean(leaf_B), SEB = sd(leaf_B)/sqrt((length(leaf_B))),
               ZN = mean(leafZn), SEZ = sd(leafZn)/sqrt((length(leafZn))),
               Na = mean(leaf_Na), SEN = sd(leaf_Na)/sqrt((length(leaf_Na))))

tissue

tissue1 = ddply(tissues, ~spp + treat, summarise, N = mean(leaf_N), SEN = sd(leaf_N)/sqrt((length(leaf_N))),
                P = mean(leaf_P), SEP = sd(leaf_P)/sqrt((length(leaf_P))), 
                K = mean(leaf_K), SEK = sd(leaf_K)/sqrt((length(leaf_K))),
                CA = mean(leaf_Ca), SEC = sd(leaf_Ca)/sqrt((length(leaf_Ca))), 
                MG = mean(leaf_Mg), SEM = sd(leaf_Mg)/sqrt((length(leaf_Mg))),
                AL = mean(leaf_Al), SEA = sd(leaf_Al)/sqrt((length(leaf_Al))))

tissue_a = subset(tissue, spp == "ACSA")
tissue_f = subset(tissue, spp == "FAGR")

tissue1_a = subset(tissue1, spp == "ACSA")
tissue1_f = subset(tissue1, spp == "FAGR")

tissue1_a$order = factor(tissue1_a$treat, levels=c("N", "D", "C"))
tissue1_f$order = factor(tissue1_f$treat, levels=c("N", "D", "C"))

tissue_a$order = factor(tissue_a$site, levels=c("L", "M", "H"))
tissue_f$order = factor(tissue_f$site, levels=c("L", "M", "H"))

tissue$order = factor(tissue$site, levels=c("L", "M", "H"))
tissue$spp = factor(tissue$spp, levels=c("ACSA", "FAGR"))

###maple leaf - site differences LMM

N_null = lme(leaf_N ~ site + treat, random = 1~mountain, data = tissues_a)
summary(N_null)

N_null = lmer(leaf_N ~ site + treat + (1|mountain), data = tissues_a)
summary(N_null)
anova(N_null)
Anova(N_null)

P_null = lmer(leaf_P ~ site + treat + (1|mountain), data = tissues_a)
summary(P_null)
anova(P_null)
Anova(P_null)

K_null = lmer(leaf_K ~ site + treat + (1|mountain), data = tissues_a)
summary(K_null)
anova(K_null)
Anova(K_null)

CA_null = lmer(leaf_Ca ~ site + treat + (1|mountain), data = tissues_a)
summary(CA_null)
anova(CA_null)
Anova(CA_null)

MG_null = lmer(leaf_Mg ~ site + treat + (1|mountain), data = tissues_a)
summary(MG_null)
anova(MG_null)
Anova(MG_null)

AL_null = lmer(leaf_Al ~ site + treat + (1|mountain), data = tissues_a)
summary(AL_null)
anova(AL_null)
Anova(AL_null)

###Beech leaf - site differences LMM

N_null = lmer(leaf_N ~ site + treat + (1|mountain), data = tissues_f)
summary(N_null)
anova(N_null)
Anova(N_null)

P_null = lmer(leaf_P ~ site + treat + (1|mountain), data = tissues_f)
summary(P_null)
anova(P_null)
Anova(P_null)

K_null = lmer(leaf_K ~ site + treat + (1|mountain), data = tissues_f)
summary(K_null)
anova(K_null)
Anova(K_null)

CA_null = lmer(leaf_Ca ~ site + treat + (1|mountain), data = tissues_f)
summary(CA_null)
anova(CA_null)
Anova(CA_null)

MG_null = lmer(leaf_Mg ~ site + treat + (1|mountain), data = tissues_f)
summary(MG_null)
anova(MG_null)
Anova(MG_null)

AL_null = lmer(leaf_Al ~ site + treat + (1|mountain), data = tissues_f)
summary(AL_null)
anova(AL_null)
Anova(AL_null)

N_null = glmer(leaf_N ~ site + (1|tissues$mountain), data = tissues, family="gaussian", na.action=na.omit)
summary(N_null)

P_null = glmer(leaf_P ~ site + treat + (1|tissues$mountain), data = tissues, family="gaussian", na.action=na.omit)
summary(P_null)

K_null = glmer(leaf_K ~ site + (1|tissues$mountain), data = tissues, family="gaussian", na.action=na.omit)
summary(K_null)

CA_null = glmer(leaf_Ca ~ site + (1|tissue_a$mountain), data = tissue_a, family="gaussian", na.action=na.omit)
summary(CA_null)

MG_null = glmer(leaf_Mg ~ site + (1|tissues$mountain), data = tissues, family="gaussian", na.action=na.omit)
summary(MG_null)

AL_null = glmer(leaf_Al ~ site + (1|tissue_a$mountain), data = tissue_a, family="gaussian", na.action=na.omit)
summary(AL_null)

###Soil data
soils = read.csv("VT_soil_analysis.csv", header = TRUE)

soil = ddply(soils, ~mountain + site, summarise, soil_pH = mean(pH), som = mean(OM_Pct), soil_P = mean(Avail_P), soil_K = mean(K), soil_Ca = mean(Ca),
             soil_Mg = mean(Mg), soil_Zn = mean(Zn), soil_B = mean(B), soil_Mn = mean(Mn), soil_Cu = mean(Cu), soil_Fe = mean(Fe), soil_Al = mean(Al), soil_Na = mean(Na),
             soil_S = mean(S), Ex_acid = mean(Exch_Acid), CEC = mean(ECEC))

lm_ph = lm(soil_pH ~ site, data = soil)
aovp = aov(lm_ph)
TukeyHSD(aovp)

lm_som = lm(som ~ site, data = soil)
aovs = aov(lm_som)
TukeyHSD(aovs)

lm_cec = lm(CEC ~ site, data = soil)
aovc = aov(lm_cec)
TukeyHSD(aovc)

lm_al = lm(soil_Al ~ site, data = soil)
aova = aov(lm_al)
TukeyHSD(aova)

lm_p = lm(soil_P ~ site, data = soil)
aovpp = aov(lm_p)
TukeyHSD(aovpp)

lm_k = lm(soil_K ~ site, data = soil)
aovk = aov(lm_k)
TukeyHSD(aovk)

lm_ca = lm(soil_Ca ~ site, data = soil)
aovcc = aov(lm_ca)
TukeyHSD(aovcc)

lm_mg = lm(soil_Mg ~ site, data = soil)
aovm = aov(lm_mg)
TukeyHSD(aovm)

###Resource utilization
#Read data
#Note correlations and resource utilization calculations done in excel - can also be done
#with data from above
cor_a = read.csv("corr_a.csv", header = TRUE)
cor_a_s = read.csv("corr_a_s.csv", header = TRUE)
cor_f = read.csv("corr_f.csv", header = TRUE)
cor_f_s = read.csv("corr_f_s.csv", header = TRUE)

dfa <- melt(cor_a, id = "Foliar_Nutrient")
colnames(dfa) <- c("x", "y", "Value")
dfa$order = factor(dfa$x, levels=c("N", "P", "K", "Ca", "Mg", "Al"))

dfas <- melt(cor_a_s, id = "Foliar_Nutrient")
colnames(dfas) <- c("x", "y", "Value")
dfas$order = factor(dfas$x, levels=c("N", "P", "K", "Ca", "Mg", "Al"))
dfas$order1 = factor(dfas$y, levels=c("Hardwoods", "Ecotone", "Spruce.fir"),
                     labels = c("Hardwoods", "Ecotone", "Spruce-fir"))

dff <- melt(cor_f, id = "Foliar_Nutrient")
colnames(dff) <- c("x", "y", "Value")
dff$order = factor(dff$x, levels=c("N", "P", "K", "Ca", "Mg", "Al"))

dffs <- melt(cor_f_s, id = "Foliar_Nutrient")
colnames(dffs) <- c("x", "y", "Value")
dffs$order = factor(dffs$x, levels=c("N", "P", "K", "Ca", "Mg", "Al"))
dffs$order1 = factor(dffs$y, levels=c("Hardwoods", "Ecotone", "Spruce.fir"),
                     labels = c("Hardwoods", "Ecotone", "Spruce-fir"))
#Plotting
tile_a = ggplot(dfa, aes(x = order, y = y, fill = Value)) +
  geom_tile(color = "black") +
  geom_text(aes(label = Value), color = "white", size = 3) +
  scale_fill_viridis_b() +
  xlab("Maple Foliar Nutrient") +
  ylab("Soil Variable") +
  ggtitle("A") +
  coord_fixed() +
  theme_bw() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))
tile_a

tile_as  = ggplot(dfas, aes(x = order, y = order1, fill = Value)) +
  geom_tile(color = "black") +
  geom_text(aes(label = Value), color = "white", size = 3) +
  scale_fill_viridis_b() +
  xlab("Maple Foliar Nutrient") +
  ylab("Planting Site NRUS") +
  ggtitle("A") +
  coord_fixed() +
  theme_bw() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))
tile_as

tile_f = ggplot(dff, aes(x = order, y = y, fill = Value)) +
  geom_tile(color = "black") +
  geom_text(aes(label = Value), color = "white", size = 3) +
  scale_fill_viridis_b() +
  xlab("Beech Foliar Nutrient") +
  ylab("Soil Variable") +
  ggtitle("B") +
  coord_fixed() +
  theme_bw() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))
tile_f

tile_fs  = ggplot(dffs, aes(x = order, y = order1, fill = Value)) +
  geom_tile(color = "black") +
  geom_text(aes(label = Value), color = "white", size = 3) +
  scale_fill_viridis_b() +
  xlab("Beech Foliar Nutrient") +
  ylab("") +
  ggtitle("B") +
  coord_fixed() +
  theme_bw() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"))
tile_fs

n = tile_a + tile_as + tile_f + tile_fs + 
  plot_layout(ncol = 2)
n

ni = tile_a + tile_f + 
  plot_layout(ncol = 1)
ni

ns = tile_as + tile_fs + 
  plot_layout(ncol = 2)
ns

tiff(file = "ru_plots.tif", width = 7.5, height = 2, units = 'in', res = 600, pointsize = 11)

ns

dev.off()

tiff(file = "ru_corr_plots.tif", width = 4.5, height = 6, units = 'in', res = 600, pointsize = 11)

ni

dev.off()

###Herb data (foliar damage and browse)

herb_h = read.csv("herb_h.csv", header = TRUE)#herbivory formatted from above for plotting
herb_b = read.csv("herb_b.csv", header = TRUE)#browse formatted from above for plotting

herb_h$site <- factor(herb_h$site, levels = c("Hardwoods", "Ecotone", "Spruce-fir"))
herb_h$species <- factor(herb_h$species, levels = c("Sugar maple", "American beech"))

herb_b$site <- factor(herb_b$site, levels = c("Hardwoods", "Ecotone", "Spruce-fir"))
herb_b$species <- factor(herb_b$species, levels = c("Sugar maple", "American beech"))

Number = c("a", "ab", "b", "c", "b", "b")
Number1 = c("a", "b", "b", "a", "a", "b")

#Plotting
h = ggplot(herb_h, aes(x = site, y = perc, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  ylab("Foliar Herbivory (% of Pots)") +
  scale_fill_manual(values = c("blue", "darkorange")) +
  xlab("") +
  theme_bw() +
  ylim(0,90) +
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(size = 10)) +
  geom_text(aes(label=Number), position=position_dodge(width=0.9), vjust=-0.3, size = 5)
h

b = ggplot(herb_b, aes(x = site, y = perc, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  ylab("Browse (% of Pots)") +
  scale_fill_manual(name = "Species", values = c("blue", "darkorange")) +
  xlab("") +
  ylim(0,40) +
  theme_bw() +
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  theme(axis.text.x = element_text(size = 10)) +
  theme(legend.position='bottom') +
  geom_text(aes(label=Number1), position=position_dodge(width=0.9), vjust=-0.3, size = 5)
b

herb_plot = h / b
herb_plot

tiff(file = "herb_plots.tif", width = 4, height = 6, units = 'in', res = 600, pointsize = 11)

herb_plot

dev.off()

#Foliar nuttrient Plotting
foliar = read.csv("foliar_n.csv", header = TRUE)#data formatted from above for plotting

foliar$site <- factor(foliar$site, levels = c("Hardwoods", "Ecotone", "Spruce-fir"))
foliar$species <- factor(foliar$species, levels = c("Sugar maple", "American beech"))
foliar$nutrient <- factor(foliar$nutrient, levels = c("N (%)", "P (%)", "K (%)", "Ca (%)", "Mg (%)", "Al (mg/kg)"))
Number = c("a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a",
           "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a", "a")

f = ggplot(foliar, aes(x = site, y = mean_n, fill = species, group = species)) +
  geom_line(position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = mean_n - mean_nse, ymax = mean_n + mean_nse), width = 0.2, position = position_dodge(0.5)) +
  geom_point(pch = 21, size = 3, position = position_dodge(0.5)) +
  #geom_errorbar(aes(ymin = scaled_s, ymax = scaled_s + scaled_sse), width = 0.2, color = "grey50", alpha = 0.5) +
  #geom_point(aes(y = scaled_s), alpha = 0.3, color = "grey50") +
  ylab("Foliar Nutrient Concentration") +
  xlab("") +
  scale_fill_manual(name = "Species", values = c("blue", "darkorange")) +
  facet_wrap(~ nutrient, ncol = 3, scales = "free") +
  theme_bw() +
  theme(legend.position='bottom')
  #geom_text(aes(label=Number), position=position_dodge(width=0.5), vjust=-2, size = 5)
f

tiff(file = "f_plots.tif", width = 7, height = 5, units = 'in', res = 600, pointsize = 11)

f

dev.off()

soliar = read.csv("soil_fig.csv", header = TRUE)#soil data formatted from above for plotting

soliar$site <- factor(soliar$site, levels = c("Hardwoods", "Ecotone", "Spruce-fir"))
#soliar$species <- factor(soliar$species, levels = c("Sugar maple", "American beech"))
soliar$soil_name <- factor(soliar$soil_name, levels = c("pH", "P (mg/kg)", "K (mg/kg)", "Ca (mg/kg)", "Mg (mg/kg)", "Al (mg/kg)"))

s = ggplot(soliar, aes(x = site, y = mean_n)) +
  geom_errorbar(aes(ymin = mean_n - mean_nse, ymax = mean_n + mean_nse), width = 0.2) +
  geom_point(pch = 16, size = 3) +
  #geom_errorbar(aes(ymin = scaled_s, ymax = scaled_s + scaled_sse), width = 0.2, color = "grey50", alpha = 0.5) +
  #geom_point(aes(y = scaled_s), alpha = 0.3, color = "grey50") +
  ylab("Soil Nutrient Concentration") +
  xlab("") +
  #scale_fill_manual(name = "Species", values = c("blue", "darkorange")) +
  facet_wrap(~ soil_name, ncol = 3, scales = "free") +
  theme_bw()
  #theme(legend.position='bottom')
s

tiff(file = "s_plots.tif", width = 7, height = 5, units = 'in', res = 600, pointsize = 11)

s

dev.off()

fig3 = (s / f) + plot_annotation(tag_levels = 'A')
fig3

tiff(file = "fig3.tif", width = 7, height = 9, units = 'in', res = 600, pointsize = 11)

fig3

dev.off()

#Survival Modeling with all Fixed Effects

trans_all = read.csv("trans_all.csv", header = TRUE)

trans_ac = subset(trans_all, spp == "ACSA")
trans_fa = subset(trans_all, spp == "FAGR")

glmm1 = glmer(survival ~ herb * elev + (1|mountain), family = "binomial", data = trans_ac)
summary(glmm1)
anova(glmm1)
Anova(glmm1)

glmm2 = glmer(survival ~ herb * elev + (1|mountain), family = "binomial", data = trans_fa)
summary(glmm2)
anova(glmm2)
Anova(glmm2)

glmm3 = glmer(survival ~ browse * elev + (1|mountain), family = "binomial", data = trans_ac)
summary(glmm3)
anova(glmm3)
Anova(glmm3)

glmm4 = glmer(survival ~ browse * elev + (1|mountain), family = "binomial", data = trans_fa)
summary(glmm4)
anova(glmm4)
Anova(glmm4)

glmm5 = glmer(survival ~ N * elev + (1|mountain), family = "binomial", data = trans_ac)
summary(glmm5)
anova(glmm5)
Anova(glmm5)

glmm6 = glmer(survival ~ P * elev + (1|mountain), family = "binomial", data = trans_ac)
summary(glmm6)
anova(glmm6)
Anova(glmm6)

glmm7 = glmer(survival ~ K * elev + (1|mountain), family = "binomial", data = trans_ac)
summary(glmm7)
anova(glmm7)
Anova(glmm7)

glmm8 = glmer(survival ~ Ca * elev + (1|mountain), family = "binomial", data = trans_ac)
summary(glmm8)
anova(glmm8)
Anova(glmm8)

glmm9 = glmer(survival ~ Mg * elev + (1|mountain), family = "binomial", data = trans_ac)
summary(glmm9)
anova(glmm9)
Anova(glmm9)

glmm10 = glmer(survival ~ Al * elev + (1|mountain), family = "binomial", data = trans_ac)
summary(glmm10)
anova(glmm10)
Anova(glmm10)

glmm11 = glmer(survival ~ N * elev + (1|mountain), family = "binomial", data = trans_fa)
summary(glmm11)
anova(glmm11)
Anova(glmm11)

glmm12 = glmer(survival ~ P * elev + (1|mountain), family = "binomial", data = trans_fa)
summary(glmm12)
anova(glmm12)
Anova(glmm12)

glmm13 = glmer(survival ~ K * elev + (1|mountain), family = "binomial", data = trans_fa)
summary(glmm13)
anova(glmm13)
Anova(glmm13)

glmm14 = glmer(survival ~ Ca * elev + (1|mountain), family = "binomial", data = trans_fa)
summary(glmm14)
anova(glmm14)
Anova(glmm14)

glmm15 = glmer(survival ~ Mg * elev + (1|mountain), family = "binomial", data = trans_fa)
summary(glmm15)
anova(glmm15)
Anova(glmm15)

glmm16 = glmer(survival ~ Al * elev + (1|mountain), family = "binomial", data = trans_fa)
summary(glmm16)
anova(glmm16)
Anova(glmm16)

#Survival Plotting

coef_ac = read.csv("coef_ac.csv", header = TRUE)#maple coefficients from modeling in new dataframe
coef_fa = read.csv("coef_fa.csv", header = TRUE)#beech coefficients from modeling in new dataframe

coef_ac$species <- factor(coef_ac$species, levels = c("Sugar maple", "American beech"))
coef_ac$variable <- factor(coef_ac$variable, levels = c("Al", "Mg", "Ca", "K", "P", "N", "Browse x Elevation", "Browse", "Foliar Herbivory x Elevation", "Foliar Herbivory", "Planting Site (Elevation)"))

coef_fa$species <- factor(coef_fa$species, levels = c("Sugar maple", "American beech"))
coef_fa$variable <- factor(coef_fa$variable, levels = c("Al", "Mg", "Ca", "K", "P", "N", "Browse x Elevation", "Browse", "Foliar Herbivory x Elevation", "Foliar Herbivory", "Planting Site (Elevation)"))
                                   
cof1 = ggplot(coef_ac, aes(x = estimate, y = variable, color = sign, shape = group)) +
  geom_linerangeh(aes(xmin = estimate - se, xmax = estimate + se)) +
  geom_point(size = 3) +
  ylab("") +
  xlab("Coefficient Estimate") +
  ggtitle("Sugar maple") +
  geom_vline(xintercept = 0, lty = 2) +
  scale_shape_manual(name = "", values=c(1,19)) +
  scale_color_manual(name = "", values = c("blue", "darkorange")) +
  theme_bw() +
  theme(legend.position='none')
cof1

cof2 = ggplot(coef_fa, aes(x = estimate, y = variable, color = sign, shape = group)) +
  geom_linerangeh(aes(xmin = estimate - se, xmax = estimate + se)) +
  geom_point(size = 3) +
  ylab("") +
  xlab("Coefficient Estimate") +
  ggtitle("American beech") +
  geom_vline(xintercept = 0, lty = 2) +
  scale_shape_manual(name = "", values=c(1,19)) +
  scale_color_manual(name = "", values = c("blue", "darkorange")) +
  theme_bw() +
  theme(legend.position='none')
cof2

cof = cof1 + cof2
cof

tiff(file = "cof_plots.tif", width = 7, height = 4, units = 'in', res = 600, pointsize = 11)

cof

dev.off()

#Biomass Modeling and Plotting - in SI only

bglmm1 = glmer(total_abv_avg ~ herb * elev + (1|mountain), family = "gaussian", data = trans_ac)
summary(bglmm1)
anova(bglmm1)
Anova(bglmm1)

bglmm2 = glmer(total_abv_avg ~ herb * elev + (1|mountain), family = "gaussian", data = trans_fa)
summary(bglmm2)
anova(bglmm2)
Anova(bglmm2)

bglmm3 = glmer(total_abv_avg ~ browse * elev + (1|mountain), family = "gaussian", data = trans_ac)
summary(bglmm3)
anova(bglmm3)
Anova(bglmm3)

bglmm4 = glmer(total_abv_avg ~ browse * elev + (1|mountain), family = "gaussian", data = trans_fa)
summary(bglmm4)
anova(bglmm4)
Anova(bglmm4)

bglmm5 = glmer(total_abv_avg ~ N * elev + (1|mountain), family = "gaussian", data = trans_ac)
summary(bglmm5)
anova(bglmm5)
Anova(bglmm5)

bglmm6 = glmer(total_abv_avg ~ P * elev + (1|mountain), family = "gaussian", data = trans_ac)
summary(bglmm6)
anova(bglmm6)
Anova(bglmm6)

bglmm7 = glmer(total_abv_avg ~ K * elev + (1|mountain), family = "gaussian", data = trans_ac)
summary(bglmm7)
anova(bglmm7)
Anova(bglmm7)

bglmm8 = glmer(total_abv_avg ~ Ca * elev + (1|mountain), family = "gaussian", data = trans_ac)
summary(bglmm8)
anova(bglmm8)
Anova(bglmm8)

bglmm9 = glmer(total_abv_avg ~ Mg * elev + (1|mountain), family = "gaussian", data = trans_ac)
summary(bglmm9)
anova(bglmm9)
Anova(bglmm9)

bglmm10 = glmer(total_abv_avg ~ Al * elev + (1|mountain), family = "gaussian", data = trans_ac)
summary(bglmm10)
anova(bglmm10)
Anova(bglmm10)

bglmm11 = glmer(total_abv_avg ~ N * elev + (1|mountain), family = "gaussian", data = trans_fa)
summary(bglmm11)
anova(bglmm11)
Anova(bglmm11)

bglmm12 = glmer(total_abv_avg ~ P * elev + (1|mountain), family = "gaussian", data = trans_fa)
summary(bglmm12)
anova(bglmm12)
Anova(bglmm12)

bglmm13 = glmer(total_abv_avg ~ K * elev + (1|mountain), family = "gaussian", data = trans_fa)
summary(bglmm13)
anova(bglmm13)
Anova(bglmm13)

bglmm14 = glmer(total_abv_avg ~ Ca * elev + (1|mountain), family = "gaussian", data = trans_fa)
summary(bglmm14)
anova(bglmm14)
Anova(bglmm14)

bglmm15 = glmer(total_abv_avg ~ Mg * elev + (1|mountain), family = "gaussian", data = trans_fa)
summary(bglmm15)
anova(bglmm15)
Anova(bglmm15)

bglmm16 = glmer(total_abv_avg ~ Al * elev + (1|mountain), family = "gaussian", data = trans_fa)
summary(bglmm16)
anova(bglmm16)
Anova(bglmm16)

coef_ac_b = read.csv("coef_ac_bio.csv", header = TRUE)#maple coefficients from biomass modeling in new dataframe
coef_fa_b = read.csv("coef_fa_bio.csv", header = TRUE)#beech coefficients from biomass modeling in new dataframe

coef_ac_b$species <- factor(coef_ac_b$species, levels = c("Sugar maple", "American beech"))
coef_ac_b$variable <- factor(coef_ac_b$variable, levels = c("Al", "Mg", "Ca", "K", "P", "N", "Browse x Elevation", "Browse", "Foliar Herbivory x Elevation", "Foliar Herbivory", "Planting Site (Elevation)"))

coef_fa_b$species <- factor(coef_fa_b$species, levels = c("Sugar maple", "American beech"))
coef_fa_b$variable <- factor(coef_fa_b$variable, levels = c("Al", "Mg", "Ca", "K", "P", "N", "Browse x Elevation", "Browse", "Foliar Herbivory x Elevation", "Foliar Herbivory", "Planting Site (Elevation)"))

cof1b = ggplot(coef_ac_b, aes(x = estimate, y = variable, color = sign, shape = group)) +
  geom_linerangeh(aes(xmin = estimate - se, xmax = estimate + se)) +
  geom_point(size = 3) +
  ylab("") +
  xlab("Coefficient Estimate") +
  ggtitle("Sugar maple - Biomass") +
  geom_vline(xintercept = 0, lty = 2) +
  scale_shape_manual(name = "", values=c(1,19)) +
  scale_color_manual(name = "", values = c("blue", "darkorange")) +
  theme_bw() +
  theme(legend.position='none')
cof1b

cof2b = ggplot(coef_fa_b, aes(x = estimate, y = variable, color = sign, shape = group)) +
  geom_linerangeh(aes(xmin = estimate - se, xmax = estimate + se)) +
  geom_point(size = 3) +
  ylab("") +
  xlab("Coefficient Estimate") +
  ggtitle("American beech - Biomass") +
  geom_vline(xintercept = 0, lty = 2) +
  scale_shape_manual(name = "", values=c(1,19)) +
  scale_color_manual(name = "", values = c("blue", "darkorange")) +
  theme_bw() +
  theme(legend.position='none')
cof2b

cofb = cof1b + cof2b
cofb

tiff(file = "cofb_plots.tif", width = 8.5, height = 4, units = 'in', res = 600, pointsize = 11)

cofb

dev.off()