library(tidyverse)
library(umap)
library("factoextra")
library(rstatix)
theme_set(theme_bw(18))

ex_rxns <- read.csv('../community_extracted_exRxns/community_samples_all.csv')
rownames(ex_rxns) <- ex_rxns$Model
ex_rxns <- ex_rxns %>% select(-c("Model"))


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

df_out <- df_out[df_out$p.adj.signif < 0.005, ]
df_out <- df_out[with(df_out, order(p.adj)), ]
df_out <- df_out[1:10,]

ex_rxn_meta <- ex_rxns %>%
  select(c(Patient, interaction_type))


ex_rxn_meta$ID <- rownames(ex_rxns)


ex_rxns <- ex_rxns %>% 
  select(-c(Patient, interaction_type))#%>%
  # mutate(ID=row_number()) 

#df_out <- df_out[df_out$feature %in% imp_feat$X,]

ex_rxns <- ex_rxns %>% select(df_out$feature) %>% mutate(ID=row_number()) 

custom.config <- umap.defaults
custom.config$n_neighbors <- 15
custom.config$min_dist <- 0.2

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
             #color = interaction_type,
             color = Patient,
             # shape = interaction_type
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
       subtitle = "Community exchange reactions UMAP") +
  # scale_color_manual(values = c(" Neutralism" = "red",
  #                               " Amensalism" = "blue",
  #                               " Competition" = "green"),
  #                    )
  # 
  scale_color_manual(values = c("COVID" = "#F79B9B",
                               "Healthy" = "#9BB9F4"))
  

ggsave("images/Community_exrxn_UMAP.pdf")


ex_ftr <- ex_rxns %>% select(-c(ID))
pca <- prcomp(ex_ftr, 
              scale = TRUE)
fviz_pca_biplot(pca,
                label="var",
                addEllipses=TRUE, 
                ellipse.level=0.60,
                habillage = ex_rxn_meta$interaction_type)



# imp_feat <- read.csv(('../MetabolicEnrichmentAnalysis/shap_imp.csv'))
# 
# imp_feat <- imp_feat[1:6,]
# 
# imp_feat <- imp_feat %>% 
#   filter(str_detect(X, "^model"))
# 
