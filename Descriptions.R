
load("MEMENTO_DWS2026.Rdata")

library(ggbeeswarm)
library(patchwork)
library(tidyverse)
library(lcmm)

vols_fs2 <- vols_fs %>% 
            select(id:eTIV, TotalGrayVol, Hippocampus_L, Hippocampus_R, CortexVol, WMH)


g1 <- ggplot(data=vols_fs2 %>% 
         pivot_longer(eTIV:WMH) %>% 
         filter(visit=="M0"),
       aes(x=value/1000, y=pipeline, color=pipeline)) + 
  geom_quasirandom() + 
  geom_boxplot(fill=NA, color="black", width=0.2, outlier.shape = NA) +
  scale_color_manual(values=c("#81b29a", "#ffafcc"))+
  labs(y="FS pipeline", x="cm3 (eq. mL)") +
  facet_grid("M0"~name, scale="free") + 
  theme_bw() + 
  theme(legend.position = "none")

g2 <- ggplot(data=vols_fs2 %>% 
               pivot_longer(eTIV:WMH) %>% 
               filter(visit=="M24"),
             aes(x=value/1000, y=pipeline, color=pipeline)) + 
  geom_quasirandom() + 
  geom_boxplot(fill=NA, color="black", width=0.2, outlier.shape = NA) +
  scale_color_manual(values=c("#81b29a", "#ffafcc"))+
  labs(y="FS pipeline", x="cm3 (eq. mL)") +
  facet_grid("M24"~name, scale="free") + 
  theme_bw() + 
  theme(legend.position = "none")

g1/g2 + plot_layout(guides = "collect", axes = "collect")


z1 <- vols_fs2 %>% 
      filter(pipeline=="cross") %>% 
      pivot_longer(eTIV:WMH, values_to = "cross")

z2 <- vols_fs2 %>% 
  filter(pipeline=="longit") %>% 
  pivot_longer(eTIV:WMH, values_to = "longit")


temp <- full_join(z1, 
                  z2 %>% select(id, visit, name, longit), 
                  by=c("id", "visit", "name"), relationship = "many-to-many") %>% 
        mutate(cross = cross/1000,
               longit = longit/1000)

ggplot(data=temp, 
       aes(x=cross, y=longit, color=visit)) + 
  geom_point() + 
  geom_smooth(method="loess", color="blue") +
  geom_abline(slope=1, intercept=0) +
  scale_color_manual(values=c("gray", "#9a8c98"))+
  labs(x="Transveral pipeline",
       y="Longitudinal pipeline") +
  facet_wrap(~name, scales="free", nrow=1) + 
  theme_bw()


z1 <- vols_fs2 %>% 
  filter(visit=="M0") %>% 
  pivot_longer(eTIV:WMH, values_to = "M0")

z2 <- vols_fs2 %>% 
  filter(visit=="M24") %>% 
  pivot_longer(eTIV:WMH, values_to = "M24")

temp <- full_join(z1, 
                  z2 %>% select(id, pipeline, name, M24, fu), by=c("id", "pipeline", "name")) %>% 
        mutate(M0 = M0/1000,
               M24 = M24/1000,
               delta = (M24 - M0) / (fu.y - fu.x) )


ggplot(data=temp,
       aes(x=pipeline, y=delta, color=pipeline) ) +
  geom_quasirandom() +
  geom_boxplot(fill=NA, color="black", width=0.2, outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype=2) +
  labs(y="Annual slope") +
  facet_wrap(~name, scales="free", nrow=1) + 
  scale_color_manual(values=c("#81b29a", "#ffafcc"))+
  theme_bw() + 
  theme(legend.position = "none")


z <- temp %>% 
     filter(pipeline=="cross") %>% 
     select(id, name, delta) %>% 
     rename(cross = delta) %>% 
     left_join(temp %>% 
                 filter(pipeline=="longit") %>% 
                 select(id, name, delta) %>% 
                 rename(longit = delta), by=c("id", "name"))

ggplot(data=z, 
       aes(x=cross, y=longit)) + 
  geom_point() + 
  geom_smooth(method="lm", color="blue") +
  geom_hline(yintercept = 0, linetype=2) + 
  geom_vline(xintercept = 0, linetype=2) + 
  labs(x="Transveral pipeline",
       y="Longitudinal pipeline") +
  facet_wrap(~name, scales="free", nrow=1) + 
  theme_bw()



temp <- vols_fs2 %>% 
        filter(visit=="M0") %>% 
  pivot_longer(cols=c(eTIV:WMH)) %>% 
  mutate(Template = str_replace(Template, "avgtemplate_", ""))

ggplot(data=temp,
       aes(x=Template, y=value/1000, color=Template) ) +
  geom_quasirandom() +
  geom_boxplot(fill=NA, color="black", width=0.2, outlier.shape = NA) +
  facet_wrap(~name, scales="free", nrow=1) + 
  scale_color_manual(values=c("#335c67", "#81b29a", "#ffafcc"))+
  theme_bw() + 
  theme(legend.position = "none")



z1 <- vols_fs2 %>% 
         select(id, age, pipeline, eTIV:WMH) %>% 
         pivot_longer(eTIV:WMH) %>% 
         mutate(age_c = (age - 70)/10,
                age2= age_c * age_c,
                vol = value / 1000)

tomodel <- z1 %>% 
           filter(pipeline=="cross") %>% 
           rename(cross = vol) %>% 
           select(-pipeline, -value) %>% 
           left_join(z1 %>% 
                       filter(pipeline=="longit") %>% 
                       rename(longit = vol) %>% 
                       select(id, age, name, longit),
                     by=c("id", "age", "name")) %>% 
           arrange(id, name, age) %>% 
           na.omit()
  

topred <- data.frame("age" = seq(50,90,by=1)) %>% 
          mutate(age_c = (age - 70)/10,
                 age2= age_c * age_c)
  

fun_preds <- function(ROI) {  
  
tt <- subset(tomodel, name==ROI)
  
m_cross <- hlme(data=tt, cross ~ age_c + age2, random=~age_c + age2, subject="id")
m_longit <- hlme(data=tt,longit ~ age_c + age2, random=~age_c + age2, subject="id")


p1 <- predictY(m_cross, topred, draws = T)$pred %>% 
      data.frame() 
p1$age <- topred$age

p2 <- predictY(m_longit, topred, draws = T)$pred %>% 
      data.frame()
p2$age <- topred$age

preds <- bind_rows(list("cross" = p1,
                        "longit" = p2), .id="pipeline") %>% 
         mutate(name=ROI)

return(preds)
}

p_cortex <- fun_preds("CortexVol")
p_eTIV <- fun_preds("eTIV")
p_hipp <- fun_preds("Hippocampus_L")
p_gray <- fun_preds("TotalGrayVol")
p_WMH <- fun_preds("WMH")


preds <- bind_rows(mget(ls(pattern = "^p_")))

ggplot(data=preds, aes(x=age, y=Ypred, ymin=lower.Ypred, ymax=upper.Ypred,
                       fill=pipeline, color=pipeline)) + 
  geom_ribbon(alpha=0.2) +
  geom_path() +
  scale_color_manual(values=c("#81b29a", "#ffafcc")) +
  scale_fill_manual(values=c("#81b29a", "#ffafcc")) +
  facet_wrap(~name, scale="free_y", nrow=1) +
  theme_bw() + 
  theme(legend.position = "bottom")


tt <- subset(tomodel, name=="Hippocampus_L")

m_cross <- hlme(data=tt, cross ~ age_c + age2, random=~age_c + age2, subject="id")
m_longit <- hlme(data=tt,longit ~ age_c + age2, random=~age_c + age2, subject="id")

summary(m_cross)


z1 <- vols_fs2 %>% 
  select(id, fu, pipeline, eTIV:WMH) %>% 
  pivot_longer(eTIV:WMH) %>% 
  left_join(patsum %>% select(id, AGE_CONS), by="id") %>% 
  mutate(vol = value / 1000)

tomodel <- z1 %>% 
  filter(pipeline=="cross") %>% 
  rename(cross = vol) %>% 
  select(-pipeline, -value) %>% 
  left_join(z1 %>% 
              filter(pipeline=="longit") %>% 
              rename(longit = vol) %>% 
              select(id, fu, name, longit),
            by=c("id", "fu", "name")) %>% 
  arrange(id, name, fu) %>% 
  na.omit()



fun_mods <- function(ROI) {  
  
  tt <- subset(tomodel, name==ROI)
  
  m_cross <- hlme(data=tt, cross ~ AGE_CONS + fu, random=~1, subject="id")
  m_longit <- hlme(data=tt,longit ~ AGE_CONS + fu, random=~1, subject="id")

  z1 <- summary(m_cross) %>% data.frame() %>% rownames_to_column("term")
  z2 <- summary(m_longit) %>% data.frame() %>% rownames_to_column("term")
  
  parms <- bind_rows(list("cross" = z1,
                          "longit" = z2), .id="pipeline") %>% 
           mutate(name = ROI)

  return(parms)
}

parms_cortex <- fun_mods("CortexVol")
parms_eTIV <- fun_mods("eTIV")
parms_hipp <- fun_mods("Hippocampus_L")
parms_gray <- fun_mods("TotalGrayVol")
parms_WMH <- fun_mods("WMH")


temp <- bind_rows(parms_cortex, parms_eTIV, parms_hipp, parms_gray, parms_WMH) %>% 
        filter(term != "intercept")

ggplot(data=temp, aes(x=term, y=coef, ymin = coef - 1.96*Se, ymax = coef + 1.96*Se, color=pipeline)) + 
  geom_pointrange(position = position_dodge(0.2)) +
  scale_color_manual(values=c("#81b29a", "#ffafcc")) +
  facet_wrap(~name, scale = "free_y", nrow=1) + 
  theme_bw() + 
  theme(legend.position = "bottom")









z1 <- vols_fs %>% 
  select(id, fu, pipeline, eTIV:WMH) %>% 
  pivot_longer(eTIV:WMH) %>% 
  left_join(patsum %>% select(id, AGE_CONS), by="id") %>% 
  mutate(vol = ifelse(name=="BrainSegVol_to_eTIV", value, value / 1000))


tomodel <- z1 %>% 
  filter(pipeline=="cross") %>% 
  rename(cross = vol) %>% 
  select(-pipeline, -value) %>% 
  left_join(z1 %>% 
              filter(pipeline=="longit") %>% 
              rename(longit = vol) %>% 
              select(id, fu, name, longit),
            by=c("id", "fu", "name")) %>% 
  arrange(id, name, fu) %>% 
  na.omit()

tt <- subset(tomodel, name=="Hippocampus_R") %>% 
      # mutate(cross = log(cross),
      #        longit = log(longit)) %>% 
      na.omit()

m_cross <- hlme(data=tt, cross ~ AGE_CONS + fu, random = ~1 + fu, subject="id", idiag = T)
m_longit <- hlme(data=tt,longit ~ AGE_CONS + fu, random = ~1 + fu, subject="id", idiag = T)

#summary(m_cross)

z1 <-  data.frame(id=m_cross$predRE$id,
                  intX=m_cross$predRE$intercept,
                  intL=m_longit$predRE$intercept,
                  slopeX=m_cross$predRE$fu,
                  slopeL=m_longit$predRE$fu)


ggplot(z1, aes(x=slopeX,y=slopeL))+
  geom_point() +
  geom_smooth(method="loess") +
  geom_abline(slope = 1,intercept = 0,color="red")



m_cross <- nlme::lme(data=tt, fixed = cross ~ AGE_CONS + fu, random = ~1 + fu | id)
m_longit <- nlme::lme(data=tt, fixed = longit ~ AGE_CONS + fu, random = ~1 + fu | id)

m_longit <- hlme(data=tt,longit ~ AGE_CONS + fu, random = ~1 + fu, subject="id")



