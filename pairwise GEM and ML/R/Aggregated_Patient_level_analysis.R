library(tidyverse)
library(umap)
library("factoextra")
library(rstatix)
library(ggplot2)
library(ggpubr)
theme_set(theme_bw(18))

ex_rxns <- read.csv('..\\aggregated_patient_level_features.csv')
rownames(ex_rxns) <- ex_rxns$key_0
ex_rxns <- ex_rxns %>% select(-c("key_0"))

ex_rxn_meta <- ex_rxns %>%
  select(label)

ex_rxn_meta$ID <- rownames(ex_rxns)


ex_rxns <- ex_rxns %>% 
  select(-label)%>%
  mutate(ID=row_number()) 

custom.config <- umap.defaults
custom.config$n_neighbors <- 5
custom.config$min_dist <- 0.9

set.seed(142)
umap_fit <- ex_rxns %>%
  select(where(is.numeric)) %>%
  scale() %>% 
  umap(config=custom.config)


umap_df <- umap_fit$layout %>%
  as.data.frame()%>%
  # rename(UMAP1="V1",
  #        UMAP2="V2") %>%
  mutate(ID=rownames(umap_fit$layout)) %>%
  inner_join(ex_rxn_meta, by="ID")


umap_df %>%
  ggplot(aes(x = V1, 
             y = V2, 
             color = label
             ))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 0, hjust=0.5), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.text.y = element_text(size = 18), 
        plot.title = element_text(size=20,  hjust = 0.5))+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "Exchange reactions UMAP") +
  scale_color_manual(values = c("C" = "#F79B9B",
                                "H" = "#9BB9F4"))




ggsave("Patient_level_umap-all_features.pdf")


#%%

ex_rxns <- read.csv('..\\aggregated_patient_level_features.csv')
rownames(ex_rxns) <- ex_rxns$key_0
ex_rxns <- ex_rxns %>% select(-c("key_0"))

ex_rxn_meta <- ex_rxns %>%
  select(c(label))

ex_rxn_meta$ID <- rownames(ex_rxns)

ex_rxns <- ex_rxns %>% select(c('EX_h_e', 'EX_co2_e', 'EX_fe3_e', 'EX_pyovd_kt_e', 'EX_fe3pyovd_kt_e',
                                'EX_o2_e', 'label'))


formulas <- paste(names(ex_rxns)[1:(ncol(ex_rxns)-1)], "~ label")

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


labels = ex_rxn_meta$label
fluxes <- ex_rxns # %>% select(-c("X")) #,"label"

sep_into_cols <- colnames(fluxes %>% select(-c(label)))

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
    `flux rate` ~ label, paired = FALSE,
    p.adjust.method = "bonferroni"
  )
pwc
pwc <- pwc %>% add_xy_position(x = "Reactions")



p <- ggboxplot(pf, x = "Reactions", y = "flux rate", color = "label", alpha= 1,
               palette = c("#F79B9B", "#9BB9F4"),
               add = "mean_se",   short.panel.labs = TRUE,
               ggtheme = theme_bw(base_size = 20),
               title = "Patient level exchange reactions" ) +
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

ggsave("images/Patient_level_top_shap_imp_reactions.pdf")




#%%

interaction_type <- read.csv('..\\covid_healthy.csv')
rownames(interaction_type) <- interaction_type$X
interaction_type <- interaction_type %>% select(-c("X"))

interaction_meta <- interaction_type %>%
  select(label)

interaction_meta$ID <- rownames(interaction_type)


interaction_type <- interaction_type %>% 
  select(-label)%>%
  mutate(ID=row_number()) 

custom.config <- umap.defaults
custom.config$n_neighbors <- 10
custom.config$min_dist <- 0.5

set.seed(142)
umap_fit <- interaction_type %>%
  select(where(is.numeric)) %>%
  scale() %>% 
  umap(config=custom.config)


umap_df <- umap_fit$layout %>%
  as.data.frame()%>%
  # rename(UMAP1="V1",
  #        UMAP2="V2") %>%
  mutate(ID=rownames(umap_fit$layout)) %>%
  inner_join(interaction_meta, by="ID")


umap_df %>%
  ggplot(aes(x = V1, 
             y = V2, 
             color = label
  ))+
  geom_point()+
  theme(axis.text.x = element_text(angle = 0, hjust=0.5), 
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #axis.text.y = element_text(size = 18), 
        plot.title = element_text(size=20,  hjust = 0.5))+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle = "Community interaction") +
  scale_color_manual(values = c("C" = "#F79B9B",
                                "H" = "#9BB9F4"))

ggsave("Patient_level_interaction_type.pdf")



ex_ftr <- ex_rxns %>% select(-c(ID))
pca <- prcomp(ex_ftr,
              scale = TRUE)
fviz_pca_biplot(pca,
                label="var",
                #addEllipses=TRUE,
                #ellipse.level=0.60,
                habillage = ex_rxn_meta$label)






