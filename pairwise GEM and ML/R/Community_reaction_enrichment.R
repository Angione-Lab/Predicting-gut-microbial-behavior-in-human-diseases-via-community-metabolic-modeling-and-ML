library(tidyverse)
library(viridis)
library(stringr)
#library(fgsea)
library(plotly)
library(heatmaply)
library(R.matlab)
library(xlsx)
library(tdr) #tdStats
library(Metrics)
library(ggthemes)
library(extrafont)
library(likert)

library(RColorBrewer)
library(BBmisc)

library(hrbrthemes)
library(reshape2)
library(plyr)
library(ggrepel)

library(dplyr)    # alternatively, this also loads %>%
library(tidyverse)

library(data.table)
library(stringr)


library(ggpubr)
library(rstatix)
library(ggplot2)
library(introdataviz)


ex_rxns <- read.csv('../community_extracted_exRxns/community_samples_all.csv')
rownames(ex_rxns) <- ex_rxns$Model
ex_rxns <- ex_rxns %>% select(-c("Model"))

ex_rxn_meta <- ex_rxns %>%
  select(c(Patient, interaction_type))

ex_rxn_meta$ID <- rownames(ex_rxns)

variances <- apply(ex_rxns %>% select(-c(Patient, interaction_type)), 2, var)
zero_variance_cols <- names(variances[variances <= 1e-4])
ex_rxns <- ex_rxns[, !names(ex_rxns) %in% zero_variance_cols]


formulas <- paste(names(ex_rxns)[1:(ncol(ex_rxns)-2)], "~ Patient")

output <- t(sapply(formulas, function(f) {
  pairwise_t_test(ex_rxns, as.formula(f))
}))

feature_name = row.names(as.data.frame(output))
df_out =  data.frame()
res = data.frame()
for (i in 1:dim(output)[1]) {
  feature = gsub( " ", "", str_split(feature_name[i], "~")[[1]][1])
  res = cbind(feature, as.data.frame(output[i,]))
  df_out = rbind(df_out, res)
  
}

df_out <- df_out[df_out$p.adj < 0.005, ]
df_out <- df_out[with(df_out, order(p.adj)), ]

df_out <- df_out[1:10,]

ex_rxns <- ex_rxns %>% select(c(df_out$feature, Patient, interaction_type))


labels = ex_rxn_meta$Patient
fluxes <- ex_rxns # %>% select(-c("X")) #,"label"

sep_into_cols <- colnames(fluxes %>% select(-c(Patient, interaction_type)))

pf <- fluxes %>%
  pivot_longer(
    sep_into_cols,
    names_to = "Reactions",
    values_to = "flux rate",
    values_drop_na = TRUE
  )


pf$`flux rate` <- as.numeric(pf$`flux rate`)

pwc <- pf %>%
  group_by(Reactions) %>%
  pairwise_t_test(
    `flux rate` ~ Patient, paired = FALSE,
    p.adjust.method = "bonferroni"
  )
pwc
#file_name = paste('Results/metabolic_analysis/',sub('/','_',trimws(supermodules$Super.Module.class[i])), "_significance.csv", sep = "")
#write.csv(pwc, file_name)

#pwc <- filter(pwc, p.adj < 0.005)
#pwc <- pwc %>% filter(p.adj.signif != "ns")
pwc <- pwc %>% add_xy_position(x = "Reactions")


p <- ggboxplot(pf, x = "Reactions", y = "flux rate", color = "Patient", alpha= 1,
          palette = c("#F79B9B", "#9BB9F4"),
          add = "mean_se",   short.panel.labs = TRUE,
          ggtheme = theme_bw(base_size = 20),
          title = "Community exchange reactions" ) +
    rotate_x_text(angle = 0)+
    geom_hline(yintercept = mean(pf$`flux rate`), linetype = 2)+
    stat_pvalue_manual(pwc, hide.ns = TRUE)

p <- p + theme(axis.text.x = element_text(angle = 0, hjust=1), 
                 plot.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank(),
                 axis.text.y = element_text(size = 12), 
                 plot.title = element_text(size=18,  hjust = 0.5)) #+ stat_pvalue_manual(pwc, hide.ns = TRUE)

p <-p + coord_flip()
print(p)


ggsave(
  "images/Community_Exrxns_sig_boxplot.pdf",
  plot = p,
  device = NULL,
  #path = 'Results/metabolic_analysis/',
  scale = 1,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  bg = NULL
)

fluxes <- ex_rxns %>% select(c("modelB_EX_k_e","modelB_EX_mg2_e", "modelB_EX_cl_e","modelA_EX_ca2_e", "modelB_EX_ca2_e",  "Patient", "interaction_type" )) #,"label"

sep_into_cols <- colnames(fluxes %>% select(-c(Patient, interaction_type)))

pf <- fluxes %>%
  pivot_longer(
    sep_into_cols,
    names_to = "Reactions",
    values_to = "flux rate",
    values_drop_na = TRUE
  )


pf$`flux rate` <- as.numeric(pf$`flux rate`)

pwc <- pf %>%
  group_by(Reactions) %>%
  pairwise_t_test(
    `flux rate` ~ Patient, paired = FALSE,
    p.adjust.method = "bonferroni"
  )
pwc

pwc <- pwc %>% add_xy_position(x = "Reactions")

p <- ggboxplot(pf, x = "Reactions", y = "flux rate", color = "Patient", alpha= 1,
               palette = c("#F79B9B", "#9BB9F4"),
               add = "mean_se",   short.panel.labs = TRUE,
               ggtheme = theme_bw(base_size = 20),
               title = "Community exchange reactions" ) +
  rotate_x_text(angle = 0)+
  geom_hline(yintercept = mean(pf$`flux rate`), linetype = 2)+
  stat_pvalue_manual(pwc, hide.ns = TRUE)

p <- p + theme(axis.text.x = element_text(angle = 0, hjust=1), 
               plot.background = element_blank(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_blank(),
               panel.background = element_blank(),
               axis.text.y = element_text(size = 12), 
               plot.title = element_text(size=18,  hjust = 0.5)) #+ stat_pvalue_manual(pwc, hide.ns = TRUE)

p <-p + coord_flip()
print(p)

ggsave(
  "images/community_exrxns_small_range_1_potassium.pdf",
  plot = p,
  device = NULL,
  #path = 'Results/metabolic_analysis/',
  scale = 1,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  bg = NULL
)

